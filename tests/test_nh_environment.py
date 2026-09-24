from __future__ import annotations

import unittest

import numpy as np

from solvelec.localization_analysis import MolecularTopology
from solvelec.nh_environment import (
    infer_nh_donors,
    nh_environment,
    select_nh_contrast_sites,
)


class NhEnvironmentTests(unittest.TestCase):
    def test_infers_real_nh_bonds_but_not_nearby_carbon_hydrogens(self) -> None:
        cell = np.diag([20.0, 20.0, 20.0])
        elements = ["N", "H", "H", "C", "N", "C", "H"]
        positions = np.asarray(
            [
                [1.0, 1.0, 1.0],
                [1.95, 1.0, 1.0],
                [1.0, 1.95, 1.0],
                [2.4, 1.0, 1.0],
                [10.0, 10.0, 10.0],
                [11.4, 10.0, 10.0],
                [12.45, 10.0, 10.0],
            ],
            dtype=float,
        )
        topology = MolecularTopology(
            np.asarray([0, 0, 0, 0, 1, 1, 1], dtype=np.int64),
            (
                {
                    "molecule_id": 1,
                    "component": "eda",
                    "atom_indices_1based": [1, 2, 3, 4],
                },
                {
                    "molecule_id": 2,
                    "component": "tmeda",
                    "atom_indices_1based": [5, 6, 7],
                },
            ),
            (),
        )
        donors = infer_nh_donors(
            elements,
            positions,
            cell,
            topology,
            maximum_nh_bond_angstrom=1.25,
        )
        self.assertEqual(len(donors), 2)
        self.assertEqual({record["component"] for record in donors}, {"eda"})
        self.assertEqual(
            {record["hydrogen_atom_index_1based"] for record in donors}, {2, 3}
        )

    def test_periodic_nh_angle_and_shell_counts(self) -> None:
        cell = np.diag([10.0, 10.0, 10.0])
        elements = ["N", "H"]
        positions = np.asarray([[9.8, 5.0, 5.0], [0.7, 5.0, 5.0]], dtype=float)
        topology = MolecularTopology(
            np.asarray([0, 0], dtype=np.int64),
            (
                {
                    "molecule_id": 1,
                    "component": "eda",
                    "atom_indices_1based": [1, 2],
                },
            ),
            (),
        )
        donors = infer_nh_donors(
            elements,
            positions,
            cell,
            topology,
            maximum_nh_bond_angstrom=1.25,
        )
        environment = nh_environment(
            [1.5, 5.0, 5.0],
            donors,
            positions,
            cell,
            shell_radii_angstrom=[1.0, 3.0],
            inward_angle_minimum_degrees=150.0,
            distance_decay_angstrom=2.5,
        )
        self.assertEqual(environment["nh_donor_count_total"], 1)
        self.assertEqual(environment["shells"][0]["inward_nh_count"], 1)
        nearest = environment["nearest_donors"][0]
        self.assertAlmostEqual(nearest["nh_bond_length_angstrom"], 0.9)
        self.assertAlmostEqual(nearest["n_h_center_angle_degrees"], 180.0)

    def test_selects_facing_and_control_sites_and_marks_nh_free_degeneracy(self) -> None:
        cell = np.diag([20.0, 20.0, 20.0])
        elements = ["N", "H"]
        positions = np.asarray([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], dtype=float)
        topology = MolecularTopology(
            np.asarray([0, 0], dtype=np.int64),
            (
                {
                    "molecule_id": 1,
                    "component": "eda",
                    "atom_indices_1based": [1, 2],
                },
            ),
            (),
        )
        donors = infer_nh_donors(
            elements,
            positions,
            cell,
            topology,
            maximum_nh_bond_angstrom=1.25,
        )
        sites = [
            {"rank": 1, "cartesian_angstrom": [3.0, 1.0, 1.0]},
            {"rank": 2, "cartesian_angstrom": [10.0, 10.0, 10.0]},
            {"rank": 3, "cartesian_angstrom": [1.0, 10.0, 10.0]},
        ]
        facing, control = select_nh_contrast_sites(
            sites,
            donors,
            positions,
            cell,
            shell_radii_angstrom=[3.0, 6.0],
            selection_radius_angstrom=3.0,
            inward_angle_minimum_degrees=120.0,
            distance_decay_angstrom=2.5,
            minimum_separation_angstrom=2.0,
        )
        self.assertEqual((facing["seed_role"], facing["site"]["rank"]), ("nh_facing", 1))
        self.assertEqual(control["seed_role"], "nh_control")
        self.assertFalse(facing["nh_degenerate_reference"])

        no_nh_facing, no_nh_control = select_nh_contrast_sites(
            sites,
            [],
            positions,
            cell,
            shell_radii_angstrom=[3.0, 6.0],
            selection_radius_angstrom=3.0,
            inward_angle_minimum_degrees=120.0,
            distance_decay_angstrom=2.5,
            minimum_separation_angstrom=2.0,
        )
        self.assertEqual(no_nh_facing["site"]["rank"], 1)
        self.assertEqual(no_nh_control["site"]["rank"], 2)
        self.assertTrue(no_nh_facing["nh_degenerate_reference"])


if __name__ == "__main__":
    unittest.main()
