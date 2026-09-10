from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from solvelec.cube import CubeData, read_cube, write_cube
from solvelec.localization_analysis import (
    cube_compatibility_problems,
    density_difference_record,
    infer_molecular_topology,
    state_localization_record,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "analyze_stage_b_localization.py"

THF_ELEMENTS = ["C", "C", "C", "O", "C", "H", "H", "H", "H", "H", "H", "H", "H"]
THF_POSITIONS = np.asarray(
    [
        [13.0800008774, 13.1200008392, 5.1300001144],
        [13.7400007248, 11.9100008011, 4.4200000763],
        [15.1300010681, 12.4499998093, 3.9700002670],
        [15.2600002289, 13.7600002289, 4.5300002098],
        [14.4200000763, 13.8000011444, 5.6700000763],
        [12.4900007248, 12.6700010300, 5.9400005341],
        [12.5300006866, 13.6700010300, 4.3699998856],
        [13.1700000763, 11.3900003433, 3.6500000954],
        [13.9900007248, 11.1100006104, 5.1399998665],
        [15.2100000381, 12.4800004959, 2.8699998856],
        [16.0300006866, 11.9800004959, 4.3699998856],
        [14.0900001526, 14.8400011063, 5.7900004387],
        [14.8700008392, 13.2100009918, 6.4800000191],
    ],
    dtype=float,
)

SETTINGS = {
    "signed_spin_integral_min": 0.75,
    "signed_spin_integral_max": 1.25,
    "cube_grid_tolerance_angstrom": 0.001,
    "density_difference_charge_tolerance_electrons": 0.10,
    "covalent_bond_scale": 1.25,
    "vdw_region_scale": 1.0,
    "cavity_probe_radius_angstrom": 2.5,
    "li_positive_spin_fraction_threshold": 0.50,
    "molecular_positive_spin_fraction_threshold": 0.50,
    "cavity_positive_spin_fraction_threshold": 0.55,
    "interstitial_positive_spin_fraction_threshold": 0.55,
}


def _cube(values: np.ndarray) -> CubeData:
    return CubeData(
        comments=("synthetic", "test"),
        origin=np.zeros(3),
        axes=np.eye(3),
        atoms=np.empty((0, 5)),
        values=np.asarray(values, dtype=float),
        coordinate_unit="angstrom",
    )


def _load_script():
    spec = importlib.util.spec_from_file_location("analyze_stage_b_localization", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_cube(path: Path, values: np.ndarray, elements: list[str], positions: np.ndarray) -> None:
    atomic_numbers = {"H": 1, "LI": 3, "C": 6, "N": 7, "O": 8}
    lines = ["synthetic", "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"]
    lines.append(f"{len(elements):5d} 0.0 0.0 0.0")
    for count, vector in zip(values.shape, np.eye(3), strict=True):
        lines.append(f"{-count:5d} {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")
    for element, position in zip(elements, positions, strict=True):
        number = atomic_numbers[element.upper()]
        lines.append(
            f"{number:5d} 0.0 {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}"
        )
    flat = values.reshape(-1)
    for start in range(0, len(flat), 6):
        lines.append(" ".join(f"{value:.8E}" for value in flat[start : start + 6]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class StageBLocalizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.cell = np.diag([22.0, 22.0, 22.0])
        self.elements = ["LI", *THF_ELEMENTS, "GH"]
        self.positions = np.vstack([[1.0, 1.0, 1.0], THF_POSITIONS, [10.0, 10.0, 10.0]])
        self.topology = infer_molecular_topology(
            self.elements,
            self.positions,
            self.cell,
            {"thf": 1},
            {"thf": {"element_counts": {"C": 4, "H": 8, "O": 1}}},
            bond_scale=1.25,
        )

    def test_infers_solvent_molecule_and_excludes_li_and_ghost(self) -> None:
        self.assertEqual(self.topology.problems, ())
        self.assertEqual(len(self.topology.molecules), 1)
        self.assertEqual(self.topology.molecules[0]["component"], "thf")
        self.assertEqual(self.topology.molecule_ids[0], -1)
        self.assertEqual(self.topology.molecule_ids[-1], -1)

    def test_state_partition_distinguishes_li_and_ghost_regions(self) -> None:
        values = np.zeros((22, 22, 22))
        values[1, 1, 1] = 1.0
        result = state_localization_record(
            _cube(values),
            self.elements,
            self.positions,
            self.cell,
            self.topology,
            self.positions[-1],
            SETTINGS,
        )
        self.assertAlmostEqual(result["spin_density"]["signed_integral"], 1.0)
        self.assertAlmostEqual(
            result["geometric_partition"]["li_positive_spin_fraction"], 1.0
        )
        self.assertTrue(result["localization_proxy_flags"]["li_centered"])
        self.assertFalse(result["localization_proxy_flags"]["ghost_cavity_centered"])

        values.fill(0.0)
        values[10, 10, 10] = 1.0
        cavity = state_localization_record(
            _cube(values),
            self.elements,
            self.positions,
            self.cell,
            self.topology,
            self.positions[-1],
            SETTINGS,
        )
        self.assertTrue(cavity["localization_proxy_flags"]["ghost_cavity_centered"])
        self.assertTrue(cavity["localization_proxy_flags"]["interstitial"])

    def test_density_difference_is_same_grid_charge_redistribution(self) -> None:
        li0 = np.zeros((22, 22, 22))
        separated = np.zeros((22, 22, 22))
        li0[1, 1, 1] = 1.0
        separated[10, 10, 10] = 1.0
        result = density_difference_record(
            _cube(li0),
            _cube(separated),
            self.elements,
            self.positions,
            self.topology,
            self.positions[-1],
            SETTINGS,
        )
        self.assertTrue(result["ready"])
        self.assertAlmostEqual(result["signed_integral_electrons"], 0.0)
        self.assertAlmostEqual(result["rearranged_electrons_half_l1"], 1.0)
        self.assertGreater(
            result["ghost_cavity_probe"]["density_accumulation_electrons"], 0.9
        )
        self.assertLess(result["li_region"]["density_change_electrons"], -0.9)

    def test_cube_grid_mismatch_is_rejected(self) -> None:
        first = _cube(np.zeros((2, 2, 2)))
        second = CubeData(
            comments=("synthetic", "test"),
            origin=np.zeros(3),
            axes=np.diag([1.0, 1.0, 1.1]),
            atoms=np.empty((0, 5)),
            values=np.zeros((2, 2, 2)),
            coordinate_unit="angstrom",
        )
        self.assertIn(
            "cube grid vectors differ",
            cube_compatibility_problems(first, second, tolerance_angstrom=0.001),
        )

    def test_cube_writer_round_trip(self) -> None:
        values = np.arange(8, dtype=float).reshape((2, 2, 2)) / 10.0
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "roundtrip.cube"
            write_cube(path, _cube(values))
            restored = read_cube(path)
        self.assertEqual(restored.coordinate_unit, "angstrom")
        np.testing.assert_allclose(restored.values, values)

    def test_script_writes_state_pair_summary_and_csv_audits(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coordinates = root / "coordinates.xyz"
            coordinate_lines = [str(len(self.elements)), "synthetic candidate"]
            coordinate_lines.extend(
                f"{element} {position[0]} {position[1]} {position[2]}"
                for element, position in zip(self.elements, self.positions, strict=True)
            )
            coordinates.write_text("\n".join(coordinate_lines) + "\n", encoding="utf-8")
            cell = root / "cell.inc"
            cell.write_text(
                "&CELL\n  A 22 0 0\n  B 0 22 0\n  C 0 0 22\n  PERIODIC XYZ\n&END CELL\n",
                encoding="utf-8",
            )
            from solvelec.provenance import sha256_file

            metadata = root / "metadata.json"
            metadata.write_text(
                json.dumps(
                    {
                        "ready": True,
                        "system_id": "pure_thf",
                        "replica": 1,
                        "candidate_id": "separated",
                        "cavity_basis_site": {"cartesian_angstrom": [10.0, 10.0, 10.0]},
                        "structure": {
                            "xyz": {"sha256": sha256_file(coordinates)},
                            "cell": {"sha256": sha256_file(cell)},
                        },
                    }
                ),
                encoding="utf-8",
            )
            spec = root / "spec.json"
            spec.write_text(
                json.dumps(
                    {
                        "system_id": "pure_thf",
                        "replica": 1,
                        "amine": None,
                        "component_counts": {"thf": 1},
                    }
                ),
                encoding="utf-8",
            )
            methods = root / "methods.json"
            methods.write_text(json.dumps({"stage_b_localization": SETTINGS}), encoding="utf-8")
            systems = root / "systems.json"
            systems.write_text(
                json.dumps(
                    {
                        "thf": {"element_counts": {"C": 4, "H": 8, "O": 1}},
                        "amines": {},
                    }
                ),
                encoding="utf-8",
            )
            real_elements = self.elements[:-1]
            real_positions = self.positions[:-1]
            li0_spin_values = np.zeros((22, 22, 22))
            separated_spin_values = np.zeros((22, 22, 22))
            li0_density_values = np.zeros((22, 22, 22))
            separated_density_values = np.zeros((22, 22, 22))
            li0_spin_values[1, 1, 1] = 1.0
            separated_spin_values[10, 10, 10] = 1.0
            li0_density_values[1, 1, 1] = 1.0
            separated_density_values[10, 10, 10] = 1.0
            cubes = {}
            for name, values in (
                ("li0_spin", li0_spin_values),
                ("separated_spin", separated_spin_values),
                ("li0_density", li0_density_values),
                ("separated_density", separated_density_values),
            ):
                path = root / f"{name}.cube"
                _write_cube(path, values, real_elements, real_positions)
                cubes[name] = path

            records = []
            for state, spin, electron in (
                ("li0_diabatic", cubes["li0_spin"], cubes["li0_density"]),
                (
                    "li_plus_e_diabatic",
                    cubes["separated_spin"],
                    cubes["separated_density"],
                ),
            ):
                output = root / f"{state}.json"
                args = SimpleNamespace(
                    campaign="pilot",
                    system="pure_thf",
                    replica=1,
                    candidate="separated",
                    candidate_metadata=str(metadata),
                    coordinates=str(coordinates),
                    cell=str(cell),
                    spec=str(spec),
                    methods=str(methods),
                    systems=str(systems),
                    output=str(output),
                    state=state,
                    spin_cube=str(spin),
                    electron_cube=str(electron),
                )
                self.assertEqual(module.run_state(args), 0)
                self.assertTrue(json.loads(output.read_text(encoding="utf-8"))["ready"])
                records.append(output)

            pair = root / "pair.json"
            pair_args = SimpleNamespace(
                campaign="pilot",
                system="pure_thf",
                replica=1,
                candidate="separated",
                candidate_metadata=str(metadata),
                coordinates=str(coordinates),
                cell=str(cell),
                spec=str(spec),
                methods=str(methods),
                systems=str(systems),
                output=str(pair),
                state_records=[str(path) for path in records],
                li0_electron_cube=str(cubes["li0_density"]),
                li_plus_e_electron_cube=str(cubes["separated_density"]),
                difference_cube=str(root / "difference.cube"),
            )
            self.assertEqual(module.run_pair(pair_args), 0)
            self.assertTrue(json.loads(pair.read_text(encoding="utf-8"))["ready"])
            self.assertTrue((root / "difference.cube").is_file())

            mechanism = root / "mechanism.json"
            mechanism.write_text(json.dumps({"ready": True}), encoding="utf-8")
            summary = root / "summary.json"
            summary_args = SimpleNamespace(
                campaign="pilot",
                mechanism_summary=str(mechanism),
                states=[str(path) for path in records],
                pairs=[str(pair)],
                states_csv=str(root / "states.csv"),
                molecules_csv=str(root / "molecules.csv"),
                pairs_csv=str(root / "pairs.csv"),
                output=str(summary),
            )
            self.assertEqual(module.run_summary(summary_args), 0)
            result = json.loads(summary.read_text(encoding="utf-8"))
            self.assertTrue(result["ready"])
            self.assertEqual(result["state_count"], 2)
            self.assertTrue((root / "molecules.csv").is_file())


if __name__ == "__main__":
    unittest.main()
