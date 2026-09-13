from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from solvelec.localization_analysis import MolecularTopology
from solvelec.preferential_solvation import (
    local_composition,
    molecular_heavy_atom_centers,
    select_preferential_sites,
)
from solvelec.rendering import render_stage_b2c_preferential_cp2k

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "prepare_stage_b2c_preferential.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("prepare_stage_b2c_preferential", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


class PreferentialSolvationTests(unittest.TestCase):
    def test_periodic_molecule_centers_and_rich_poor_selection(self) -> None:
        cell = np.diag([20.0, 20.0, 20.0])
        elements = ["C", "O", "H", "N", "C", "H"]
        positions = np.asarray(
            [
                [2.0, 2.0, 2.0],
                [2.8, 2.0, 2.0],
                [2.0, 2.8, 2.0],
                [15.0, 15.0, 15.0],
                [15.8, 15.0, 15.0],
                [15.0, 15.8, 15.0],
            ],
            dtype=float,
        )
        topology = MolecularTopology(
            np.asarray([0, 0, 0, 1, 1, 1], dtype=np.int64),
            (
                {
                    "molecule_id": 1,
                    "component": "thf",
                    "atom_indices_1based": [1, 2, 3],
                },
                {
                    "molecule_id": 2,
                    "component": "eda",
                    "atom_indices_1based": [4, 5, 6],
                },
            ),
            (),
        )
        centers = molecular_heavy_atom_centers(elements, positions, cell, topology)
        around_eda = local_composition(
            [15.4, 15.0, 15.0],
            centers,
            cell,
            [2.0, 4.0, 6.0],
            {"thf": 1, "eda": 1},
            amine_component="eda",
        )
        self.assertEqual(around_eda["shells"][1]["amine_mole_fraction"], 1.0)
        sites = [
            {"rank": 1, "cartesian_angstrom": [15.4, 15.0, 15.0]},
            {"rank": 2, "cartesian_angstrom": [2.4, 2.0, 2.0]},
            {"rank": 3, "cartesian_angstrom": [10.0, 10.0, 10.0]},
        ]
        rich, poor = select_preferential_sites(
            sites,
            centers,
            cell,
            {"thf": 1, "eda": 1},
            amine_component="eda",
            shell_radii_angstrom=[2.0, 4.0, 6.0],
            selection_radius_angstrom=4.0,
            minimum_separation_angstrom=2.0,
        )
        self.assertEqual((rich["seed_role"], rich["site"]["rank"]), ("eda_rich", 1))
        self.assertEqual((poor["seed_role"], poor["site"]["rank"]), ("eda_poor", 2))
        self.assertGreater(
            rich["selection_amine_mole_fraction"], poor["selection_amine_mole_fraction"]
        )

    def test_matched_pair_render_has_no_li_or_localization_constraint(self) -> None:
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
            settings = methods["stage_b2c_preferential_smoke"]
            states = {record["id"]: record for record in settings["states"]}
            for state_id in ("neutral", "anion"):
                output = root / state_id / "cp2k.inp"
                render_stage_b2c_preferential_cp2k(
                    ROOT
                    / "workflow"
                    / "templates"
                    / "cp2k"
                    / "stage_b2c_preferential_smoke.inp.tpl",
                    output,
                    project=f"test_{state_id}",
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

    def test_summary_interprets_larger_attachment_proxy_as_rich_preference(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records: list[str] = []
            for system, fractions, energies in (
                ("pure_thf", (0.0, 0.0), (1.0, 1.0)),
                ("eda_1p5m", (0.5, 0.1), (1.2, 1.0)),
                ("eda_3m", (0.7, 0.2), (1.4, 1.0)),
            ):
                for role, fraction, energy in zip(
                    ("eda_rich", "eda_poor"), fractions, energies, strict=True
                ):
                    path = root / f"{system}-{role}.json"
                    _write_json(
                        path,
                        {
                            "ready": True,
                            "system_id": system,
                            "amine": None if system == "pure_thf" else "eda",
                            "replica": 1,
                            "seed_role": role,
                            "composition_degenerate_reference": system == "pure_thf",
                            "local_composition_at_seed": {
                                "shells": [
                                    {
                                        "radius_angstrom": radius,
                                        "amine_mole_fraction": fraction,
                                    }
                                    for radius in (4.0, 6.0, 8.0)
                                ]
                            },
                            "energies": {"vertical_attachment_proxy_ev": energy},
                            "positive_spin_fraction_by_component": {},
                        },
                    )
                    records.append(str(path))
            output = root / "summary.json"
            self.assertEqual(
                module.run_summary(
                    SimpleNamespace(
                        records=records,
                        expected_systems=["pure_thf", "eda_1p5m", "eda_3m"],
                        expected_replicas=[1],
                        expected_seed_roles=["eda_rich", "eda_poor"],
                        methods=str(ROOT / "configs" / "methods.yaml"),
                        campaign="pilot",
                        output=str(output),
                    )
                ),
                0,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(summary["ready"])
            by_system = {record["system_id"]: record for record in summary["comparisons"]}
            self.assertEqual(
                by_system["pure_thf"]["preference_at_smoke_tolerance"],
                "composition_degenerate_reference",
            )
            self.assertEqual(
                by_system["eda_1p5m"]["preference_at_smoke_tolerance"], "eda_rich"
            )
            self.assertGreater(
                by_system["eda_3m"]["vertical_attachment_proxy_ev"]["rich_minus_poor"],
                0,
            )


if __name__ == "__main__":
    unittest.main()
