from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pypdf import PdfReader

from solvelec.cube import CubeData, write_cube

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "reports" / "stage_b" / "build_report.py"
PUBLISHED = ROOT / "reports" / "stage_b"
SYSTEMS = ("pure_thf", "eda_1p5m")
STATES = ("li0_diabatic", "li_plus_e_diabatic")


def _load_script():
    spec = importlib.util.spec_from_file_location("build_stage_b_report", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _write_gate(summary: Path, gate: Path) -> None:
    gate.write_text(f"sha256 {_sha256(summary)}  {summary.name}\n", encoding="utf-8")


def _source(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _write_cube(path: Path, center: tuple[float, float, float], sign_pair: bool = False) -> None:
    shape = (12, 12, 12)
    grid = np.indices(shape, dtype=float).transpose(1, 2, 3, 0)
    center_array = np.asarray(center)
    values = np.exp(-np.sum((grid - center_array) ** 2, axis=-1) / 4.0)
    if sign_pair:
        values -= np.exp(-np.sum((grid - np.asarray((2.0, 2.0, 2.0))) ** 2, axis=-1) / 3.0)
    write_cube(
        path,
        CubeData(
            comments=("synthetic Stage-B report fixture", "angstrom grid"),
            origin=np.zeros(3),
            axes=np.eye(3),
            atoms=np.asarray(
                [
                    [3.0, 0.0, 2.0, 2.0, 2.0],
                    [0.0, 0.0, 8.0, 8.0, 8.0],
                ]
            ),
            values=values,
            coordinate_unit="angstrom",
        ),
    )


def create_fixture(root: Path) -> tuple[Path, Path]:
    run_root = root / "runs"
    base = run_root / "pilot"
    base.mkdir(parents=True)
    targets = {"compact": 3.5, "separated": 6.0, "distant": 9.0}
    systems_summary: dict[str, dict] = {}
    for system_index, system in enumerate(SYSTEMS):
        systems_summary[system] = {"replicas": [1, 2, 3]}
        for replica in (1, 2, 3):
            candidate_dir = base / "stage_b" / system / f"r{replica}" / "candidates"
            records = []
            for candidate_index, (candidate, target) in enumerate(targets.items()):
                leaf = candidate_dir / candidate
                leaf.mkdir(parents=True, exist_ok=True)
                xyz = leaf / "coordinates.xyz"
                cell = leaf / "cell.inc"
                xyz.write_text(
                    "2\nsynthetic\nLi 2 2 2\nGh 8 8 8\n", encoding="utf-8"
                )
                cell.write_text(
                    "&CELL\n  ABC 12 12 12\n  PERIODIC XYZ\n&END CELL\n", encoding="utf-8"
                )
                achieved = target + 0.08 * (replica - 2) + 0.03 * system_index
                records.append(
                    {
                        "candidate_id": candidate,
                        "ready": True,
                        "target_li_cavity_distance_angstrom": target,
                        "achieved_li_cavity_distance_angstrom": achieved,
                        "li_site": {"radius_angstrom": 2.1 + 0.05 * candidate_index},
                        "cavity_basis_site": {
                            "radius_angstrom": 2.5 + 0.05 * candidate_index
                        },
                        "structure": {"xyz": _source(xyz), "cell": _source(cell)},
                    }
                )
            _write_json(
                candidate_dir / "manifest.json",
                {
                    "ready": True,
                    "system_id": system,
                    "replica": replica,
                    "candidate_count": 3,
                    "candidates": records,
                },
            )

    candidate_summary = base / "stage_b_candidates.summary.json"
    _write_json(
        candidate_summary,
        {"ready": True, "record_count": 6, "systems": systems_summary},
    )
    _write_gate(candidate_summary, base / "stage_b_candidates.done")

    legacy_records = []
    mechanism_records = []
    gaps = {"pure_thf": 1.9497304227957517, "eda_1p5m": 0.5081838643066795}
    state_values = {
        ("pure_thf", "li0_diabatic"): (-2793.041822276325, 3.0, 2.992804421057),
        ("pure_thf", "li_plus_e_diabatic"): (-2792.970171004862, 2.0, 2.035839353097),
        ("eda_1p5m", "li0_diabatic"): (-3117.3402362892894, 3.0, 2.969287908158),
        ("eda_1p5m", "li_plus_e_diabatic"): (-3117.3215608767355, 2.0, 2.045021367901),
    }
    for system in SYSTEMS:
        state_records = []
        for state in STATES:
            directory = base / "stage_b" / system / "r1" / "mechanism_smoke" / "separated" / state
            directory.mkdir(parents=True, exist_ok=True)
            cp2k_input = directory / "cp2k.inp"
            cp2k_output = directory / "cp2k.out"
            cp2k_input.write_text(f"# synthetic {system} {state}\n", encoding="utf-8")
            cp2k_output.write_text("PROGRAM ENDED AT synthetic\n", encoding="utf-8")
            energy, target, current = state_values[(system, state)]
            state_records.append(
                {
                    "system_id": system,
                    "replica": 1,
                    "candidate_id": "separated",
                    "state_id": state,
                    "converged": True,
                    "normal_termination": True,
                    "energy_hartree": energy,
                    "cdft_constraint_gate": {
                        "target_electrons": target,
                        "current_electrons": current,
                        "deviation_electrons": current - target,
                        "tolerance_electrons": 0.05,
                        "iteration": 3 if system == "pure_thf" else 4,
                    },
                    "input": _source(cp2k_input),
                    "output": _source(cp2k_output),
                }
            )

        mechanism_records.append(
            {
                "system_id": system,
                "replica": 1,
                "candidate_id": "separated",
                "ready": True,
                "complete_pair": True,
                "states": state_records,
                "fixed_geometry_diabatic_gap_ev": gaps[system],
            }
        )
        legacy_state = state_records[1]
        legacy_records.append(
            {
                key: legacy_state[key]
                for key in (
                    "system_id",
                    "replica",
                    "candidate_id",
                    "converged",
                    "normal_termination",
                    "energy_hartree",
                    "cdft_constraint_gate",
                    "input",
                    "output",
                )
            }
        )

    legacy_summary = base / "stage_b_cp2k_smoke.summary.json"
    _write_json(
        legacy_summary,
        {
            "ready": True,
            "scientific_status": "NUMERICAL_SMOKE_ONLY_NOT_A_LOCALIZATION_RESULT",
            "records": legacy_records,
        },
    )
    _write_gate(legacy_summary, base / "stage_b.done")

    mechanism_summary = base / "stage_b_mechanism_smoke.summary.json"
    _write_json(
        mechanism_summary,
        {
            "campaign": "pilot",
            "ready": True,
            "energy_definition": "E(li_plus_e_diabatic) - E(li0_diabatic)",
            "records": mechanism_records,
        },
    )
    _write_gate(mechanism_summary, base / "stage_b_mechanism_smoke.done")

    localization_values = {
        ("eda_1p5m", "li0_diabatic"): (3.3870, 0.4682, 0.1547, 0.0974, 0.0052),
        ("eda_1p5m", "li_plus_e_diabatic"): (4.3920, 0.0180, 0.0785, 0.4531, 0.0506),
        ("pure_thf", "li0_diabatic"): (4.2410, 0.4489, 0.1014, 0.1374, 0.0003),
        ("pure_thf", "li_plus_e_diabatic"): (4.5542, 0.0076, 0.1028, 0.4474, 0.0285),
    }
    states = []
    pairs = []
    for system_index, system in enumerate(SYSTEMS):
        for state_index, state in enumerate(STATES):
            directory = base / "stage_b" / system / "r1" / "mechanism_smoke" / "separated" / state
            spin_cube = directory / "spin.cube"
            electron_cube = directory / "electron.cube"
            center = (2.0, 2.0, 2.0) if state_index == 0 else (7.0, 7.0 - system_index, 7.0)
            _write_cube(spin_cube, center)
            _write_cube(electron_cube, center)
            radius, li, maximum_molecule, interstitial, cavity = localization_values[
                (system, state)
            ]
            states.append(
                {
                    "ready": True,
                    "system_id": system,
                    "replica": 1,
                    "state_id": state,
                    "spin_density": {
                        "signed_integral": 1.0,
                        "radius_angstrom": radius,
                        "inverse_participation_ratio": 0.01,
                    },
                    "geometric_partition": {
                        "li_positive_spin_fraction": li,
                        "maximum_solvent_molecule_positive_spin_fraction": maximum_molecule,
                        "interstitial_positive_spin_fraction": interstitial,
                        "ghost_cavity_probe_positive_spin_fraction": cavity,
                    },
                    "localization_proxy_flags": {
                        "li_centered": False,
                        "molecule_centered": False,
                        "interstitial": False,
                        "ghost_cavity_centered": False,
                    },
                    "inputs": {
                        "spin_cube": _source(spin_cube),
                        "electron_cube": _source(electron_cube),
                    },
                }
            )
        pair_dir = base / "stage_b" / system / "r1" / "mechanism_smoke" / "separated"
        difference_cube = pair_dir / "electron_density_difference.cube"
        _write_cube(difference_cube, (7.0, 7.0 - system_index, 7.0), sign_pair=True)
        pair_values = (
            (-0.0001, 1.8112, -1.0154, 0.2374, 0.0163)
            if system == "pure_thf"
            else (0.0005, 1.8385, -1.1366, 0.2348, 0.0251)
        )
        signed, rearranged, li_change, interstitial_acc, cavity_acc = pair_values
        pairs.append(
            {
                "ready": True,
                "system_id": system,
                "replica": 1,
                "inputs": {"density_difference_cube": _source(difference_cube)},
                "density_difference": {
                    "ready": True,
                    "signed_integral_electrons": signed,
                    "rearranged_electrons_half_l1": rearranged,
                    "li_region": {"density_change_electrons": li_change},
                    "interstitial_region": {
                        "fraction_of_all_accumulation": interstitial_acc
                    },
                    "ghost_cavity_probe": {
                        "fraction_of_all_accumulation": cavity_acc
                    },
                },
            }
        )

    localization_summary = base / "stage_b_localization.summary.json"
    _write_json(
        localization_summary,
        {
            "campaign": "pilot",
            "ready": True,
            "state_count": 4,
            "pair_count": 2,
            "source_mechanism_summary": _source(mechanism_summary),
            "states": states,
            "pairs": pairs,
        },
    )
    _write_gate(localization_summary, base / "stage_b_localization.done")
    for name in (
        "stage_b_localization.states.csv",
        "stage_b_localization.molecules.csv",
        "stage_b_localization.pairs.csv",
    ):
        (base / name).write_text("synthetic,ready\nfixture,true\n", encoding="utf-8")
    return run_root, root / "report"


class StageBReportTests(unittest.TestCase):
    def test_published_report_matches_server_provenance(self) -> None:
        metrics = json.loads((PUBLISHED / "stage_b_metrics.json").read_text(encoding="utf-8"))
        provenance = json.loads(
            (PUBLISHED / "report_provenance.json").read_text(encoding="utf-8")
        )
        self.assertTrue(metrics["ready"])
        self.assertTrue(provenance["ready"])
        self.assertEqual(metrics["source_file_count"], 93)
        self.assertEqual(
            metrics["scientific_status"],
            "COMPLETED_NUMERICAL_STAGE_B_PILOT_NOT_PRODUCTION_MECHANISM",
        )
        self.assertAlmostEqual(
            metrics["mechanism"]["gaps_ev"]["pure_thf"], 1.9497304227957517
        )
        self.assertAlmostEqual(
            metrics["mechanism"]["gaps_ev"]["eda_1p5m"], 0.5081838643066795
        )

        published_files = {path.name: path for path in PUBLISHED.rglob("*") if path.is_file()}
        self.assertEqual(len(provenance["products"]), 7)
        for product in provenance["products"]:
            with self.subTest(product=product["path"]):
                local = published_files[Path(product["path"]).name]
                self.assertEqual(local.stat().st_size, product["size_bytes"])
                self.assertEqual(_sha256(local), product["sha256"])
        self.assertEqual(len(PdfReader(str(PUBLISHED / "stage_b_report_zh.pdf")).pages), 10)

    def test_complete_report_builds_with_hash_validated_fixture(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            run_root, output = create_fixture(Path(directory))
            metrics = module.build(run_root, "pilot", output)

            self.assertTrue(metrics["ready"])
            self.assertEqual(metrics["source_file_count"], 68)
            self.assertAlmostEqual(metrics["mechanism"]["gaps_ev"]["pure_thf"], 1.9497304)
            self.assertEqual(len(metrics["localization"]["states"]), 4)
            self.assertEqual(len(metrics["localization"]["pairs"]), 2)
            self.assertEqual(
                metrics["candidate_bank"]["pure_thf"]["compact"][
                    "achieved_distances_angstrom"
                ],
                [3.42, 3.5, 3.58],
            )

            pdf = output / "stage_b_report_zh.pdf"
            self.assertGreaterEqual(len(PdfReader(str(pdf)).pages), 9)
            provenance = json.loads(
                (output / "report_provenance.json").read_text(encoding="utf-8")
            )
            self.assertTrue(provenance["ready"])
            self.assertEqual(len(provenance["products"]), 7)
            self.assertTrue((output / "figures" / "cube_maps.png").is_file())
            self.assertTrue(
                (output / "data" / "cp2k_inputs" / "pure_thf" / "li0_diabatic.inp").is_file()
            )

    def test_report_workflow_is_analysis_only_and_immutable(self) -> None:
        rules = (ROOT / "workflow" / "rules" / "58_stage_b_report.smk").read_text(
            encoding="utf-8"
        )
        build_rule, gate_rule = rules.split("rule stage_b_report:", 1)
        inputs, parameters = build_rule.split("    params:", 1)
        self.assertNotIn("stage_b_localization.done", inputs)
        self.assertIn("accepted_stage_b=stage_b_report_handoff", parameters)
        self.assertIn("trajectory_analysis", build_rule)
        self.assertNotIn("cp2k.psmp", rules.lower())
        self.assertNotIn("gmx", rules.lower())
        self.assertIn("--summary {input.provenance:q}", gate_rule)


if __name__ == "__main__":
    unittest.main()
