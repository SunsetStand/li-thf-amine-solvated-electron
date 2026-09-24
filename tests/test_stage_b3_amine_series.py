from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from solvelec.rendering import render_stage_b3_amine_series_cp2k

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "stage_b3_amine_series.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("stage_b3_amine_series", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


class StageB3AmineSeriesTests(unittest.TestCase):
    def test_renderer_keeps_pair_fixed_nuclei_li_free_and_unconstrained(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coordinates = root / "coordinates.xyz"
            cell = root / "cell.inc"
            coordinates.write_text("2\ntest\nH 1 1 1\nGh 5 5 5\n", encoding="utf-8")
            cell.write_text(
                "&CELL\n  A 10 0 0\n  B 0 10 0\n  C 0 0 10\n  PERIODIC XYZ\n&END CELL\n",
                encoding="utf-8",
            )
            methods = json.loads((ROOT / "configs" / "methods.yaml").read_text())
            settings = methods["stage_b3_amine_series"]
            states = {record["id"]: record for record in settings["states"]}
            for state_id in ("neutral", "anion"):
                output = root / state_id / "cp2k.inp"
                render_stage_b3_amine_series_cp2k(
                    ROOT
                    / "workflow"
                    / "templates"
                    / "cp2k"
                    / "stage_b2c_preferential_smoke.inp.tpl",
                    output,
                    project=f"test_{state_id}_stage_b3",
                    coordinates_path=coordinates,
                    cell_path=cell,
                    method=settings,
                    state=states[state_id],
                )
                text = output.read_text(encoding="utf-8")
                self.assertIn("RUN_TYPE ENERGY", text)
                self.assertNotIn("&KIND LI", text.upper())
                self.assertNotIn("&CDFT", text.upper())
                self.assertNotIn("&CONSTRAINT", text.upper())
                self.assertEqual("&E_DENSITY_CUBE" in text, state_id == "anion")

    def test_pilot_summary_treats_tmeda_as_nh_free_control(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths: list[str] = []
            systems = ("deta_3m", "tmeda_1p5m", "tmeda_3m")
            for system in systems:
                for role in ("nh_facing", "nh_control"):
                    if system == "deta_3m":
                        energy = 1.10 if role == "nh_facing" else 1.00
                        degenerate = False
                    else:
                        energy = 0.95
                        degenerate = True
                    path = root / f"{system}-{role}.json"
                    _write_json(
                        path,
                        {
                            "ready": True,
                            "system_id": system,
                            "amine": system.split("_")[0],
                            "replica": 1,
                            "seed_role": role,
                            "nh_degenerate_reference": degenerate,
                            "energies": {"vertical_attachment_proxy_ev": energy},
                        },
                    )
                    paths.append(str(path))
            output = root / "summary.json"
            self.assertEqual(
                module.run_pilot_summary(
                    SimpleNamespace(
                        records=paths,
                        expected_systems=list(systems),
                        expected_replicas=[1],
                        energy_tolerance_ev=0.05,
                        campaign="pilot",
                        output=str(output),
                    )
                ),
                0,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(summary["ready"])
            comparisons = {
                record["system_id"]: record for record in summary["comparisons"]
            }
            self.assertEqual(
                comparisons["deta_3m"]["pilot_result"],
                "nh_facing_has_larger_attachment_proxy",
            )
            self.assertEqual(
                comparisons["tmeda_1p5m"]["pilot_result"],
                "nh_free_or_geometrically_degenerate_reference",
            )

    def test_reanalysis_summary_rejects_missing_records(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            record = root / "record.json"
            _write_json(
                record,
                {
                    "ready": True,
                    "system_id": "eda_1p5m",
                    "amine": "eda",
                    "replica": 1,
                    "seed_role": "eda_rich",
                    "nh_donor_count_total": 36,
                    "positive_spin_fraction_on_nh_bearing_amine_molecules": 0.1,
                },
            )
            output = root / "summary.json"
            module.run_reanalysis_summary(
                SimpleNamespace(
                    records=[str(record)],
                    expected_count=2,
                    campaign="pilot",
                    output=str(output),
                )
            )
            summary = json.loads(output.read_text(encoding="utf-8"))
            self.assertFalse(summary["ready"])
            self.assertIn("expected 2 records", summary["problems"][0])


if __name__ == "__main__":
    unittest.main()
