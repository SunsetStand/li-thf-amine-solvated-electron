from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

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
    def test_source_handoff_rejects_superseded_failed_or_mismatched_sources(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            spec_path, validation_path = root / "spec.json", root / "validation.json"
            tpr, trajectory = root / "production.tpr", root / "production.xtc"
            tpr.write_bytes(b"immutable-tpr")
            trajectory.write_bytes(b"immutable-trajectory")
            args = SimpleNamespace(campaign="pilot", source_campaign="amine_series_pilot",
                                   spec=str(spec_path), validation=str(validation_path),
                                   tpr=str(tpr), trajectory=str(trajectory),
                                   methods=str(ROOT / "configs/methods.yaml"), output_dir=str(root / "handoff"))
            spec = {"system_id": "12pda_1p5m", "replica": 1}
            validation = {**spec, "metrics": {"ready": True},
                          "inputs": {key: str(getattr(args, key)) for key in ("spec", "tpr", "trajectory")}}
            _write_json(spec_path, spec)
            _write_json(validation_path, validation)
            self.assertEqual(module.validate_source_handoff(args), spec)
            for bad_spec, bad_validation in (
                ({**spec, "system_id": "12pda_3m"}, validation),
                (spec, {**validation, "metrics": {"ready": False}}),
                (spec, {**validation, "replica": 2}),
                (spec, {**validation, "inputs": {**validation["inputs"], "trajectory": str(tpr)}}),
            ):
                with self.subTest(spec=bad_spec, validation=bad_validation):
                    _write_json(spec_path, bad_spec)
                    _write_json(validation_path, bad_validation)
                    with self.assertRaises(ValueError), patch.object(module.subprocess, "run") as execute:
                        module.run_source_handoff(args)
                    execute.assert_not_called()
            self.assertEqual(tpr.read_bytes(), b"immutable-tpr")
            self.assertEqual(trajectory.read_bytes(), b"immutable-trajectory")

    def test_source_handoff_runs_only_analysis_snapshot_and_candidate_tools(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            validation = root / "validation.json"
            _write_json(validation, {"metrics": {"ready": True}})
            args = SimpleNamespace(campaign="pilot", source_campaign="amine_series_pilot",
                                   spec=str(root / "spec.json"), validation=str(validation),
                                   tpr="immutable.tpr", trajectory="immutable.xtc", methods="methods.json",
                                   output_dir=str(root / "handoff"))
            commands = []
            def execute(command, *, check):
                self.assertTrue(check)
                commands.append(command)
                output = Path(command[command.index("--output") + 1])
                output.parent.mkdir(parents=True, exist_ok=True)
                if command[2] == "gate":
                    summary = Path(command[command.index("--summary") + 1])
                    output.write_text(f"sha256 {module.sha256_file(summary)} {summary.name}\n")
                else:
                    _write_json(output, {"ready": True})
            with patch.object(module, "validate_source_handoff"), patch.object(module.subprocess, "run", side_effect=execute):
                self.assertEqual(module.run_source_handoff(args), 0)
            self.assertEqual([command[2] for command in commands], ["analyze", "select", "prepare", "summary", "gate"])
            self.assertTrue(all(Path(command[1]).name in ("analyze_classical_ensemble.py", "prepare_stage_b.py") for command in commands))
            summary = json.loads((root / "handoff/stage_b_candidates.summary.json").read_text())
            self.assertEqual(summary["source_validation"]["sha256"], module.sha256_file(validation))

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
