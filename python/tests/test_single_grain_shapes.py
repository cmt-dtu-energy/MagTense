from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np

SINGLE_GRAIN_DIR = Path(__file__).resolve().parents[1] / "experiments" / "single_grain"
sys.path.insert(0, str(SINGLE_GRAIN_DIR))

from single_grain_coercivity import (  # noqa: E402
    build_shape_mesh,
    resolve_shape_spec,
)
from utils.simulation import create_single_grain_problem  # noqa: E402


class SingleGrainShapeTests(unittest.TestCase):
    def test_cube_keeps_uniform_grid_interpretation(self) -> None:
        spec = resolve_shape_spec("Cube", None, 20e-9)

        self.assertEqual(spec.shape, "cube")
        self.assertEqual(spec.variant, "cube")
        self.assertIsNone(build_shape_mesh(spec, 4))

        problem = create_single_grain_problem(20e-9, 4, cuda=False, cvode=False)
        self.assertEqual(problem.ntot, 64)
        self.assertEqual(problem.grid_type, 1)
        np.testing.assert_allclose(problem.grid_L, [20e-9, 20e-9, 20e-9])

    def test_sphere_uses_outer_size_as_diameter(self) -> None:
        spec = resolve_shape_spec("sphere", None, 20e-9)
        mesh = build_shape_mesh(spec, 3)

        self.assertEqual(spec.parameters["radius"], 10e-9)
        self.assertIsNotNone(mesh)
        self.assertEqual(mesh.shape, "sphere")
        self.assertEqual(mesh.target_tiles, 27)
        self.assertGreater(mesh.achieved_tiles, 0)
        np.testing.assert_allclose(mesh.root_bounds, [[-10e-9] * 3, [10e-9] * 3])

    def test_ellipsoid_variants_use_mild_axis_elongation(self) -> None:
        x_spec = resolve_shape_spec("elipsoid", "x", 20e-9)
        z_spec = resolve_shape_spec("ellipsoid_z", None, 20e-9)

        self.assertEqual(x_spec.variant, "ellipsoid_x")
        self.assertEqual(z_spec.variant, "ellipsoid_z")
        np.testing.assert_allclose(x_spec.parameters["semi_axes"], [10e-9, 7e-9, 7e-9])
        np.testing.assert_allclose(z_spec.parameters["semi_axes"], [7e-9, 7e-9, 10e-9])

        x_mesh = build_shape_mesh(x_spec, 3)
        z_mesh = build_shape_mesh(z_spec, 3)
        self.assertEqual(x_mesh.shape, "ellipsoid")
        self.assertEqual(z_mesh.shape, "ellipsoid")
        self.assertGreater(x_mesh.achieved_tiles, 0)
        self.assertGreater(z_mesh.achieved_tiles, 0)

    def test_other_supported_shapes_resolve_outer_size(self) -> None:
        rectangle = resolve_shape_spec("recantgle", None, 10e-9)
        cylinder = resolve_shape_spec("cylinder", None, 10e-9)
        hexagon = resolve_shape_spec("hexagon", None, 10e-9)

        np.testing.assert_allclose(rectangle.parameters["dimensions"], [10e-9, 7e-9, 7e-9])
        self.assertEqual(cylinder.parameters["radius"], 5e-9)
        self.assertEqual(cylinder.parameters["length"], 10e-9)
        self.assertEqual(hexagon.parameters["side_length"], 5e-9)
        self.assertEqual(hexagon.parameters["height"], 10e-9)

    def test_invalid_shape_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown shape"):
            resolve_shape_spec("pyramid", None, 10e-9)


if __name__ == "__main__":
    unittest.main()
