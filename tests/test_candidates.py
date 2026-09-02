from __future__ import annotations

import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

from solvelec.candidates import ranked_void_sites, read_xyz, select_candidate_pairs, write_xyz
from solvelec.config import load_repository_configs
from solvelec.rendering import render_stage_b_cp2k

ROOT = Path(__file__).resolve().parents[1]


class CandidateTests(unittest.TestCase):
    def test_ranked_void_sites_are_deterministic_and_separated(self) -> None:
        cell = np.diag([12.0, 12.0, 12.0])
        atoms = np.asarray([[0.0, 0.0, 0.0], [6.0, 6.0, 6.0]])
        radii = np.asarray([1.5, 1.5])
        arguments = dict(
            count=6,
            points_per_axis=6,
            refinement_levels=1,
            minimum_separation_angstrom=2.0,
            minimum_clearance_angstrom=0.1,
        )
        first = ranked_void_sites(atoms, radii, cell, **arguments)
        second = ranked_void_sites(atoms, radii, cell, **arguments)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 6)
        self.assertTrue(all(site["radius_angstrom"] > 0 for site in first))

    def test_candidate_pair_targets_are_selected_without_pair_reuse(self) -> None:
        sites = [
            {
                "rank": index + 1,
                "radius_angstrom": 2.0 - index * 0.05,
                "fractional": [index / 10.0, 0.0, 0.0],
                "cartesian_angstrom": [float(index), 0.0, 0.0],
            }
            for index in range(10)
        ]
        definitions = [
            {
                "id": "compact",
                "target_li_cavity_distance_angstrom": 3.0,
                "tolerance_angstrom": 0.1,
            },
            {
                "id": "separated",
                "target_li_cavity_distance_angstrom": 6.0,
                "tolerance_angstrom": 0.1,
            },
        ]
        selected = select_candidate_pairs(sites, np.diag([20.0, 20.0, 20.0]), definitions)
        self.assertEqual([record["id"] for record in selected], ["compact", "separated"])
        self.assertAlmostEqual(selected[0]["achieved_li_cavity_distance_angstrom"], 3.0)
        self.assertAlmostEqual(selected[1]["achieved_li_cavity_distance_angstrom"], 6.0)

    def test_stage_b_xyz_and_cp2k_input_preserve_explicit_roles(self) -> None:
        _, _, methods = load_repository_configs(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            xyz = directory / "candidate.xyz"
            cell = directory / "cell.inc"
            output = directory / "cp2k.inp"
            write_xyz(
                xyz,
                ["LI", "O", "H", "Gh"],
                np.asarray(
                    [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [2.5, 2.0, 2.0], [7.0, 7.0, 7.0]]
                ),
                "test",
            )
            cell.write_text(
                "&CELL\n  A 10 0 0\n  B 0 10 0\n  C 0 0 10\n  PERIODIC XYZ\n&END CELL\n",
                encoding="utf-8",
            )
            elements, positions, _ = read_xyz(xyz)
            self.assertEqual(elements, ["LI", "O", "H", "GH"])
            self.assertEqual(positions.shape, (4, 3))
            render_stage_b_cp2k(
                ROOT / "workflow" / "templates" / "cp2k" / "stage_b_smoke.inp.tpl",
                output,
                project="candidate_smoke",
                coordinates_path=xyz,
                cell_path=cell,
                method=methods["stage_b_smoke"],
                li_atom_index=1,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn("TARGET 2.0", text)
            self.assertIn("ATOMS 1", text)
            self.assertIn("POTENTIAL GTH-PBE-q3", text)
            self.assertIn("&KIND Gh", text)
            self.assertIn("GHOST TRUE", text)
            self.assertIn("&E_DENSITY_CUBE", text)
            self.assertNotIn("&SPIN_DENSITY_CUBE", text)
            self.assertNotIn("&HF", text)

            inconsistent = deepcopy(methods["stage_b_smoke"])
            inconsistent["li_target_valence_electrons"] = 0.0
            with self.assertRaisesRegex(ValueError, "pseudopotential valence - 1"):
                render_stage_b_cp2k(
                    ROOT / "workflow" / "templates" / "cp2k" / "stage_b_smoke.inp.tpl",
                    output,
                    project="invalid_li_target",
                    coordinates_path=xyz,
                    cell_path=cell,
                    method=inconsistent,
                    li_atom_index=1,
                )


if __name__ == "__main__":
    unittest.main()
