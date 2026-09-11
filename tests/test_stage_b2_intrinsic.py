from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from solvelec.candidates import read_xyz, write_xyz
from solvelec.parsers import HARTREE_TO_EV
from solvelec.provenance import sha256_file

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "prepare_stage_b2_intrinsic.py"

THF_ELEMENTS = ["C", "C", "C", "O", "C", "H", "H", "H", "H", "H", "H", "H", "H"]
THF_POSITIONS = np.asarray(
    [
        [3.08, 3.12, 5.13],
        [3.74, 1.91, 4.42],
        [5.13, 2.45, 3.97],
        [5.26, 3.76, 4.53],
        [4.42, 3.80, 5.67],
        [2.49, 2.67, 5.94],
        [2.53, 3.67, 4.37],
        [3.17, 1.39, 3.65],
        [3.99, 1.11, 5.14],
        [5.21, 2.48, 2.87],
        [6.03, 1.98, 4.37],
        [4.09, 4.84, 5.79],
        [4.87, 3.21, 6.48],
    ],
    dtype=float,
)


def _load_script():
    spec = importlib.util.spec_from_file_location("prepare_stage_b2_intrinsic", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _write_gate(summary: Path, gate: Path) -> None:
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    gate.write_text(f"sha256 {digest}  {summary.name}\n", encoding="utf-8")


class StageB2IntrinsicTests(unittest.TestCase):
    def test_prepare_removes_li_and_builds_two_independent_void_seeds(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_xyz = root / "candidate.xyz"
            source_cell = root / "cell.inc"
            write_xyz(
                source_xyz,
                ["LI", *THF_ELEMENTS, "Gh"],
                np.vstack([[1.0, 1.0, 1.0], THF_POSITIONS, [9.0, 9.0, 9.0]]),
                "synthetic Stage-B candidate",
            )
            source_cell.write_text(
                "&CELL\n  A 12 0 0\n  B 0 12 0\n  C 0 0 12\n  PERIODIC XYZ\n&END CELL\n",
                encoding="utf-8",
            )
            candidate_manifest = root / "candidate-manifest.json"
            _write_json(
                candidate_manifest,
                {
                    "ready": True,
                    "system_id": "pure_thf",
                    "amine": None,
                    "replica": 1,
                    "snapshot_id": "pure_thf-r1-test",
                    "void_search": {
                        "sites": [
                            {
                                "rank": 1,
                                "radius_angstrom": 2.1,
                                "fractional": [0.75, 0.75, 0.75],
                                "cartesian_angstrom": [9.0, 9.0, 9.0],
                            },
                            {
                                "rank": 2,
                                "radius_angstrom": 1.9,
                                "fractional": [0.75, 0.25, 0.75],
                                "cartesian_angstrom": [9.0, 3.0, 9.0],
                            },
                        ]
                    },
                    "candidates": [
                        {
                            "ready": True,
                            "candidate_id": "separated",
                            "structure": {
                                "xyz": {
                                    "path": str(source_xyz),
                                    "sha256": sha256_file(source_xyz),
                                },
                                "cell": {
                                    "path": str(source_cell),
                                    "sha256": sha256_file(source_cell),
                                },
                            },
                        }
                    ],
                },
            )
            candidate_summary = root / "candidate-summary.json"
            _write_json(candidate_summary, {"ready": True})
            candidate_gate = root / "stage_b_candidates.done"
            _write_gate(candidate_summary, candidate_gate)
            spec = root / "spec.json"
            _write_json(
                spec,
                {
                    "system_id": "pure_thf",
                    "amine": None,
                    "replica": 1,
                    "component_counts": {"thf": 1},
                },
            )
            output_dir = root / "stage_b2"
            manifest = output_dir / "manifest.json"
            args = SimpleNamespace(
                candidate_manifest=str(candidate_manifest),
                candidate_summary=str(candidate_summary),
                candidate_gate=str(candidate_gate),
                spec=str(spec),
                methods=str(ROOT / "configs" / "methods.yaml"),
                output_dir=str(output_dir),
                output=str(manifest),
            )

            self.assertEqual(module.run_prepare(args), 0)
            result = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertTrue(result["ready"])
            self.assertEqual(result["seed_ids"], ["void_01", "void_02"])
            self.assertFalse(result["experiment_context"]["pfas_present"])
            for record in result["seeds"]:
                elements, _positions, comment = read_xyz(record["structure"]["xyz"]["path"])
                self.assertNotIn("LI", elements)
                self.assertEqual(elements.count("GH"), 1)
                self.assertIn("Li_absent", comment)
                self.assertEqual(record["electronic_state"]["charge"], -1)

            cp2k_input = root / "cp2k.inp"
            self.assertEqual(
                module.run_render(
                    SimpleNamespace(
                        manifest=str(manifest),
                        methods=str(ROOT / "configs" / "methods.yaml"),
                        template=str(
                            ROOT
                            / "workflow"
                            / "templates"
                            / "cp2k"
                            / "stage_b2_intrinsic_smoke.inp.tpl"
                        ),
                        seed="void_01",
                        project="pure_thf_r1_void_01_stage_b2_intrinsic_smoke",
                        output=str(cp2k_input),
                    )
                ),
                0,
            )
            text = cp2k_input.read_text(encoding="utf-8")
            self.assertIn("CHARGE -1", text)
            self.assertIn("MULTIPLICITY 2", text)
            self.assertIn("&KIND Gh", text)
            self.assertNotIn("&KIND LI", text.upper())
            self.assertNotIn("&CDFT", text.upper())
            self.assertNotIn("&CONSTRAINT", text.upper())

    def test_summary_compares_seed_outcomes_and_writes_checksum_gate(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records = []
            for system_index, system in enumerate(("pure_thf", "eda_1p5m")):
                for seed_index, seed in enumerate(("void_01", "void_02")):
                    record = root / f"{system}-{seed}.json"
                    _write_json(
                        record,
                        {
                            "ready": True,
                            "system_id": system,
                            "amine": None if system == "pure_thf" else "eda",
                            "replica": 1,
                            "seed_id": seed,
                            "energy_hartree": -100.0 - system_index - 0.01 * seed_index,
                            "localization_proxy_flags": {
                                "li_centered": False,
                                "interstitial": True,
                            },
                        },
                    )
                    records.append(str(record))
            summary = root / "summary.json"
            self.assertEqual(
                module.run_summary(
                    SimpleNamespace(
                        records=records,
                        expected_systems=["pure_thf", "eda_1p5m"],
                        expected_replicas=[1],
                        expected_seeds=["void_01", "void_02"],
                        campaign="pilot",
                        methods=str(ROOT / "configs" / "methods.yaml"),
                        output=str(summary),
                    )
                ),
                0,
            )
            result = json.loads(summary.read_text(encoding="utf-8"))
            self.assertTrue(result["ready"])
            self.assertEqual(result["record_count"], 4)
            self.assertAlmostEqual(result["systems"][0]["energy_spread_ev"], 0.01 * HARTREE_TO_EV)
            self.assertIn("PFAS is absent", result["method_limitations"][-1])

            gate = root / "stage_b2_intrinsic_smoke.done"
            self.assertEqual(
                module.run_gate(SimpleNamespace(summary=str(summary), output=str(gate))), 0
            )
            self.assertEqual(
                gate.read_text(encoding="utf-8").split(),
                ["sha256", hashlib.sha256(summary.read_bytes()).hexdigest(), summary.name],
            )


if __name__ == "__main__":
    unittest.main()
