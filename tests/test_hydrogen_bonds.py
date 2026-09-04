from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from solvelec.cavity_visualization import cavity_centered_positions, render_cavity_hbond_svg
from solvelec.hydrogen_bonds import (
    find_hydrogen_bonds,
    hydrogen_bond_counts,
    point_segment_distance,
)


class HydrogenBondTests(unittest.TestCase):
    def test_point_segment_distance(self) -> None:
        distance = point_segment_distance(
            np.asarray([0.0, 1.0, 0.0]),
            np.asarray([-1.0, 0.0, 0.0]),
            np.asarray([1.0, 0.0, 0.0]),
        )
        self.assertAlmostEqual(distance, 1.0)

    def test_eda_eda_bond_excludes_same_residue_and_labels_cavity(self) -> None:
        positions = np.asarray(
            [
                [5.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [3.2, 0.0, 0.0],
            ]
        )
        bonds = find_hydrogen_bonds(
            positions,
            {0: [1]},
            {"eda_n": np.asarray([0, 2, 3]), "thf_o": np.asarray([], dtype=int)},
            np.asarray([0, 0, 1, 0]),
            np.diag([20.0, 20.0, 20.0]),
            np.asarray([4.0, 1.0, 0.0]),
            0.75,
            distance_cutoff_angstrom=3.5,
            angle_cutoff_degree=150.0,
            cavity_shell_thickness_angstrom=0.5,
            cavity_bridge_margin_angstrom=0.25,
        )
        self.assertEqual(len(bonds), 1)
        self.assertEqual(bonds[0].bond_type, "eda_nh_eda_n")
        self.assertTrue(bonds[0].cavity_associated)
        self.assertTrue(bonds[0].cavity_bridging)
        counts = hydrogen_bond_counts(bonds)
        self.assertEqual(counts["eda_eda_hydrogen_bonds"], 1)
        self.assertEqual(counts["eda_thf_hydrogen_bonds"], 0)

    def test_svg_contains_cavity_and_bond_legend(self) -> None:
        atoms = [
            {
                "index": 0,
                "element": "N",
                "position_angstrom": [-1.0, 0.0, 0.0],
            },
            {
                "index": 1,
                "element": "H",
                "position_angstrom": [0.0, 0.0, 0.0],
            },
            {
                "index": 2,
                "element": "O",
                "position_angstrom": [1.0, 0.0, 0.0],
            },
        ]
        bonds = [
            {
                "hydrogen_index": 1,
                "acceptor_index": 2,
                "bond_type": "eda_nh_thf_o",
                "cavity_associated": True,
                "cavity_bridging": True,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "view.svg"
            render_cavity_hbond_svg(
                output,
                title="test cavity",
                atoms=atoms,
                covalent_bonds=[(0, 1)],
                hydrogen_bonds=bonds,
                cavity_radius_angstrom=1.0,
                view_radius_angstrom=5.0,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn("test cavity", text)
            self.assertIn("EDA N-H...O(THF)", text)
            self.assertEqual(text.count("projection"), 3)

    def test_cavity_centering_reassembles_boundary_crossing_residue(self) -> None:
        positions = np.asarray(
            [
                [9.8, 5.0, 5.0],
                [0.2, 5.0, 5.0],
                [0.4, 5.0, 5.0],
            ]
        )
        centered = cavity_centered_positions(
            positions,
            np.asarray([0, 0, 0]),
            np.asarray([0, 1]),
            np.asarray([0.0, 5.0, 5.0]),
            np.diag([10.0, 10.0, 10.0]),
        )
        self.assertLess(float(np.ptp(centered[:, 0])), 1.0)
        self.assertLess(float(np.max(np.abs(centered[:, 0]))), 1.0)


if __name__ == "__main__":
    unittest.main()
