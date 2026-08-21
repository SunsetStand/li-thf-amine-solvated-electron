from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from solvelec.cube import analyze_spin_density, read_cube
from solvelec.pbc import minimum_image_displacement, periodic_weighted_centroid

ROOT = Path(__file__).resolve().parents[1]


class PbcCubeTests(unittest.TestCase):
    def test_minimum_image(self) -> None:
        delta = minimum_image_displacement([9.8, 0, 0], [0.2, 0, 0], [10, 10, 10])
        np.testing.assert_allclose(delta, [-0.4, 0, 0], atol=1.0e-12)

    def test_centroid_crosses_boundary(self) -> None:
        centroid = periodic_weighted_centroid([[0.1, 1, 1], [9.9, 1, 1]], [1, 1], [10, 10, 10])
        self.assertLess(min(centroid[0], 10.0 - centroid[0]), 0.2)

    def test_cube_metrics(self) -> None:
        cube = read_cube(ROOT / "tests" / "fixtures" / "boundary_spin.cube")
        metrics = analyze_spin_density(cube)
        self.assertAlmostEqual(metrics.electron_count, 1.0, places=12)
        self.assertAlmostEqual(metrics.signed_integral, 1.0, places=12)
        self.assertEqual(metrics.positive_voxels, 4)
        self.assertGreater(metrics.radius, 0)
        self.assertTrue(np.isfinite(metrics.inverse_participation_ratio))


if __name__ == "__main__":
    unittest.main()
