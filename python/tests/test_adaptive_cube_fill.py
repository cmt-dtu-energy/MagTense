"""Tests for conservative adaptive cube filling."""

from __future__ import annotations

from math import sqrt
from pathlib import Path
import sys
import unittest

import numpy as np

SINGLE_GRAIN_DIR = Path(__file__).resolve().parents[1] / "experiments" / "single_grain"
sys.path.insert(0, str(SINGLE_GRAIN_DIR))

from utils.geometry import (  # noqa: E402
    AxisAlignedBox,
    Cube,
    Cylinder,
    Ellipsoid,
    H0_REFERENCE_TABLE,
    HexagonalPrism,
    Sphere,
    adaptive_cube_fill,
    assert_no_cube_overlaps,
    calibrate_h0_sweep,
    conservative_cube_condition,
    cubes_overlap,
    estimate_h0_from_volume,
    fill_domain_mesh,
    generate_cylinder,
    generate_hexagonal_prism,
    recommended_h0,
    select_best_h0_records,
    _off_axis_priority,
)


class AdaptiveCubeFillTests(unittest.TestCase):
    def _assert_domain_mesh_is_inside_and_non_overlapping(self, domain, mesh) -> None:
        cubes = [
            Cube(
                center=center,
                h=float(dimensions[0]),
                level=int(level),
                i=index,
                j=0,
                k=0,
            )
            for index, (center, dimensions, level) in enumerate(
                zip(mesh.centers, mesh.dimensions, mesh.levels)
            )
        ]
        self.assertGreater(len(cubes), 0)
        for cube in cubes:
            self.assertEqual(domain.classify_cube(cube, 1e-12).name, "INSIDE")
        assert_no_cube_overlaps(cubes, tol=1e-14)

    def _assert_center_mirror_orbits_are_complete(self, cubes, center) -> None:
        center = np.asarray(center, dtype=float)
        counts: dict[tuple[int, float, tuple[int, int, int]], int] = {}
        for cube in cubes:
            relative = np.abs(cube.center - center)
            key = (
                cube.level,
                round(cube.h, 12),
                tuple(int(round(value / cube.h * 1_000_000_000)) for value in relative),
            )
            counts[key] = counts.get(key, 0) + 1

        for (level, h, scaled), count in counts.items():
            expected = 2 ** sum(value != 0 for value in scaled)
            self.assertEqual(
                count,
                expected,
                f"level={level}, h={h}, orbit={scaled} is not mirror-complete",
            )

    def test_sphere_cubes_are_conservatively_inside_and_non_overlapping(self) -> None:
        eps = 1e-12
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=140,
            h0=0.5,
            h_min=0.03125,
            max_depth=5,
            eps=eps,
            overshoot_policy="soft",
        )

        self.assertGreater(len(cubes), 0)
        self.assertGreater(len(cubes), 0)
        self.assertLessEqual(len(cubes), int(1.1 * 140))
        self.assertEqual(stats.accepted_count, len(cubes))
        for cube in cubes:
            required = sqrt(3.0) * cube.h / 2.0 + eps
            self.assertGreaterEqual(domain.sdf(cube.center), required)
            self.assertTrue(conservative_cube_condition(domain, cube, eps))
        assert_no_cube_overlaps(cubes, tol=1e-14)

    def test_box_aligned_fill_has_exact_volume(self) -> None:
        domain = AxisAlignedBox(center=[0.0, 0.0, 0.0], dimensions=[1.0, 1.0, 1.0])
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            bbox=domain.bounding_box,
            N_target=1000,
            h0=0.25,
            h_min=0.25,
            max_depth=0,
            eps=1e-12,
            overshoot_policy="soft",
        )

        self.assertEqual(len(cubes), 64)
        self.assertAlmostEqual(stats.accepted_volume, 1.0)
        self.assertAlmostEqual(stats.fill_fraction, 1.0)
        self.assertEqual(stats.boundary_cells_left_unresolved, 0)
        assert_no_cube_overlaps(cubes, tol=1e-14)

    def test_refinement_produces_dyadic_levels_near_boundary(self) -> None:
        domain = Ellipsoid(center=[0.0, 0.0, 0.0], semi_axes=[1.0, 0.75, 0.5])
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=220,
            h0=0.5,
            h_min=0.03125,
            max_depth=5,
            eps=1e-12,
            overshoot_policy="soft",
        )

        levels = sorted({cube.level for cube in cubes})
        sizes = sorted({round(cube.h, 12) for cube in cubes}, reverse=True)
        self.assertGreater(stats.max_level_reached, 0)
        self.assertGreaterEqual(len(levels), 2)
        for level, size in zip(levels, sizes):
            self.assertAlmostEqual(size, 0.5 / (2**level))

        finest = [cube for cube in cubes if cube.level == max(levels)]
        coarsest = [cube for cube in cubes if cube.level == min(levels)]
        finest_boundary_distance = min(abs(domain.sdf(cube.center)) for cube in finest)
        coarsest_boundary_distance = min(abs(domain.sdf(cube.center)) for cube in coarsest)
        self.assertLess(finest_boundary_distance, coarsest_boundary_distance)

    def test_termination_reports_max_depth_and_unresolved_boundary(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            bbox=[[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]],
            N_target=10,
            h0=1.0,
            h_min=0.0,
            max_depth=0,
            eps=1e-12,
            overshoot_policy="soft",
        )

        self.assertEqual(cubes, [])
        self.assertEqual(stats.termination_reason, "max_depth reached")
        self.assertGreater(stats.boundary_cells_left_unresolved, 0)
        self.assertEqual(stats.accepted_count, 0)

    def test_boundary_crossing_root_cube_is_not_accepted(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            bbox=[[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]],
            N_target=1,
            h0=2.0,
            h_min=2.0,
            max_depth=0,
            eps=1e-12,
            overshoot_policy="soft",
        )

        self.assertEqual(cubes, [])
        self.assertEqual(stats.accepted_count, 0)
        self.assertEqual(stats.boundary_cells_left_unresolved, 1)

    def test_overlap_helper_allows_touching_faces_only(self) -> None:
        first = Cube(center=[0.0, 0.0, 0.0], h=1.0, level=0, i=0, j=0, k=0)
        touching = Cube(center=[1.0, 0.0, 0.0], h=1.0, level=0, i=1, j=0, k=0)
        overlapping = Cube(center=[0.99, 0.0, 0.0], h=1.0, level=0, i=1, j=0, k=0)

        self.assertFalse(cubes_overlap(first, touching, tol=1e-12))
        self.assertTrue(cubes_overlap(first, overlapping, tol=1e-12))

    def test_legacy_mesh_wrapper_returns_prism_mesh_metadata(self) -> None:
        domain = AxisAlignedBox(center=[0.0, 0.0, 0.0], dimensions=[1.0, 1.0, 1.0])
        mesh = fill_domain_mesh(
            domain,
            N_target=100,
            h0=0.5,
            h_min=0.5,
            max_depth=0,
            overshoot_policy="soft",
        )

        self.assertEqual(mesh.achieved_tiles, 8)
        self.assertEqual(mesh.dimensions.shape, mesh.centers.shape)
        self.assertIn("cube_fill_stats", mesh.refinement_metadata)
        self.assertEqual(mesh.to_micromag_kwargs()["res"], [8, 1, 1])

    def test_cylinder_fill_is_conservative_and_mixed_level(self) -> None:
        domain = Cylinder(center=[0.0, 0.0, 0.0], radius=1.0, length=2.0, axis="z")
        mesh = fill_domain_mesh(
            domain,
            N_target=180,
            h0=0.5,
            h_min=0.03125,
            max_depth=5,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_priority",
        )

        self._assert_domain_mesh_is_inside_and_non_overlapping(domain, mesh)
        self.assertGreater(len(set(mesh.levels.tolist())), 1)
        self.assertLess(np.min(mesh.levels), np.max(mesh.levels))

    def test_hexagonal_prism_fill_is_conservative(self) -> None:
        for rotation in (0.0, 15.0):
            with self.subTest(rotation=rotation):
                domain = HexagonalPrism(
                    center=[0.0, 0.0, 0.0],
                    side_length=1.0,
                    height=0.75,
                    axis="z",
                    rotation_degrees=rotation,
                )
                mesh = fill_domain_mesh(
                    domain,
                    N_target=160,
                    h0=0.25,
                    h_min=0.03125,
                    max_depth=4,
                    eps=1e-12,
                    overshoot_policy="soft",
                    queue_policy="symmetric_priority",
                )

                self._assert_domain_mesh_is_inside_and_non_overlapping(domain, mesh)
                self.assertLessEqual(mesh.represented_volume, domain.volume + 1e-12)

    def test_symmetric_breadth_first_preserves_sphere_mirror_orbits(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=90,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_breadth_first",
        )

        self.assertGreater(len(cubes), 0)
        self.assertEqual(stats.queue_policy, "symmetric_breadth_first")
        self._assert_center_mirror_orbits_are_complete(cubes, domain.center)
        assert_no_cube_overlaps(cubes, tol=1e-14)

    def test_target_limited_symmetric_refinement_stops_on_complete_orbits(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=24,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_breadth_first",
            complete_layers=False,
        )

        self.assertGreaterEqual(len(cubes), 24)
        self.assertIn(stats.termination_reason, {"target reached", "boundary queue empty"})
        self._assert_center_mirror_orbits_are_complete(cubes, domain.center)

    def test_complete_layer_fill_commits_a_full_level_within_tolerance(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=80,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_breadth_first",
            complete_layers=True,
        )

        self.assertEqual(len(cubes), 80)
        self.assertTrue(stats.complete_layers)
        self.assertTrue(stats.target_accepted)
        self.assertEqual(stats.target_status, "at_target")
        self.assertEqual(stats.termination_reason, "target tolerance reached")
        self.assertEqual(stats.layers_used, 1)
        self._assert_center_mirror_orbits_are_complete(cubes, domain.center)

    def test_complete_layer_fill_discards_a_full_level_above_tolerance(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=24,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            complete_layers=True,
            queue_policy="breadth_first",
        )

        self.assertEqual(len(cubes), 8)
        self.assertEqual({cube.level for cube in cubes}, {0})
        self.assertFalse(stats.target_accepted)
        self.assertEqual(stats.target_status, "below_tolerance")
        self.assertEqual(
            stats.termination_reason,
            "next complete layer exceeds target tolerance",
        )

    def test_complete_layers_do_not_depend_on_symmetry_queue_policy(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        common = dict(
            domain=domain,
            N_target=80,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            complete_layers=True,
        )
        breadth_first, _ = adaptive_cube_fill(
            **common, queue_policy="breadth_first"
        )
        symmetric, _ = adaptive_cube_fill(
            **common, queue_policy="symmetric_priority"
        )

        breadth_geometry = [(cube.level, *cube.center) for cube in breadth_first]
        symmetric_geometry = [(cube.level, *cube.center) for cube in symmetric]
        self.assertEqual(breadth_geometry, symmetric_geometry)

    def test_complete_layer_accepts_below_target_within_tolerance(self) -> None:
        domain = AxisAlignedBox(center=[0.0, 0.0, 0.0], dimensions=[1.5, 1.0, 1.0])
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=13,
            h0=0.5,
            h_min=0.5,
            max_depth=0,
            eps=1e-12,
            complete_layers=True,
        )

        self.assertEqual(len(cubes), 12)
        self.assertTrue(stats.target_accepted)
        self.assertEqual(stats.target_status, "within_lower_tolerance")

    def test_strict_min_inside_fraction_matches_default_conservative_fill(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        kwargs = dict(
            domain=domain,
            N_target=90,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            complete_layers=True,
        )

        strict_default, default_stats = adaptive_cube_fill(**kwargs)
        strict_explicit, explicit_stats = adaptive_cube_fill(
            **kwargs,
            min_inside_fraction=1.0,
        )

        self.assertEqual(len(strict_default), len(strict_explicit))
        self.assertAlmostEqual(default_stats.accepted_volume, explicit_stats.accepted_volume)

    def test_min_inside_fraction_validation(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        with self.assertRaises(ValueError):
            adaptive_cube_fill(domain=domain, min_inside_fraction=0.49)
        with self.assertRaises(ValueError):
            adaptive_cube_fill(domain=domain, min_inside_fraction=1.01)

    def test_partial_inside_fraction_can_accept_more_volume(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        strict, strict_stats = adaptive_cube_fill(
            domain=domain,
            N_target=90,
            h0=0.5,
            h_min=0.125,
            max_depth=2,
            eps=1e-12,
            complete_layers=True,
            min_inside_fraction=1.0,
        )
        partial, partial_stats = adaptive_cube_fill(
            domain=domain,
            N_target=90,
            h0=0.5,
            h_min=0.125,
            max_depth=2,
            eps=1e-12,
            complete_layers=True,
            min_inside_fraction=0.5,
            inside_fraction_samples=3,
        )

        self.assertGreaterEqual(len(partial), len(strict))
        self.assertGreaterEqual(partial_stats.accepted_volume, strict_stats.accepted_volume)

    def test_grid_shift_candidates_are_ranked_by_target_then_fill(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=90,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            complete_layers=True,
            grid_shifts="half_step",
        )

        self.assertGreater(len(cubes), 0)
        self.assertEqual(len(stats.grid_shift_candidates), 27)
        self.assertIn(
            stats.selected_grid_shift,
            [row["grid_shift"] for row in stats.grid_shift_candidates],
        )

    def test_legacy_overshoot_policy_available_when_complete_layers_disabled(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        _cubes, stats = adaptive_cube_fill(
            domain=domain,
            N_target=24,
            h0=0.5,
            h_min=0.0625,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
            complete_layers=False,
        )

        self.assertFalse(stats.complete_layers)
        self.assertIn(stats.termination_reason, {"target reached", "boundary queue empty"})

    def test_cylinder_priority_favors_radial_boundary_before_axial_caps(self) -> None:
        domain = Cylinder(center=[0.0, 0.0, 0.0], radius=1.0, length=2.0, axis="z")
        radial_boundary = Cube(center=[0.75, 0.0, 0.0], h=0.5, level=0, i=0, j=0, k=0)
        axial_cap = Cube(center=[0.0, 0.0, 0.75], h=0.5, level=0, i=0, j=0, k=1)

        self.assertLess(
            _off_axis_priority(radial_boundary, domain),
            _off_axis_priority(axial_cap, domain),
        )

    def test_h0_volume_estimate_matches_formula(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        target = 125
        eta = 0.8

        expected = (eta * domain.volume / target) ** (1.0 / 3.0)
        self.assertAlmostEqual(estimate_h0_from_volume(domain, target, eta=eta), expected)

    def test_h0_sweep_returns_one_record_per_target_and_multiplier(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        records = calibrate_h0_sweep(
            domain,
            target_values=[16, 32],
            h0_multipliers=[0.75, 1.0],
            h_min_ratio=0.25,
            max_depth=2,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_breadth_first",
        )

        self.assertEqual(len(records), 4)
        for record in records:
            self.assertEqual(record["shape"], "sphere")
            self.assertIn(record["target"], {16, 32})
            self.assertIn(record["h0_multiplier"], {0.75, 1.0})
            self.assertGreater(record["h0"], 0.0)
            self.assertGreaterEqual(record["accepted_count"], 0)
            self.assertIn("count_ratio", record)
            self.assertIn("within_target_band", record)
            self.assertIn("short_circuited", record)
            self.assertEqual(
                record["within_target_band"],
                record["count_ratio"] >= 0.9,
            )

    def test_h0_sweep_short_circuits_over_resolved_level0_grid(self) -> None:
        domain = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        records = calibrate_h0_sweep(
            domain,
            target_values=[16],
            h0_multipliers=[0.05],
            h_min_ratio=0.25,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
            queue_policy="symmetric_breadth_first",
            target_count_max_ratio=1.25,
        )

        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertTrue(record["short_circuited"])
        self.assertIsNone(record["fill_fraction"])
        self.assertEqual(record["termination_reason"], "level-0 count exceeds target band")
        self.assertEqual(record["selection_warning"], "increase h0 or increase target cells")
        self.assertGreater(record["count_ratio"], record["target_count_max_ratio"])

    def test_best_h0_records_ignore_high_fill_massive_overshoot(self) -> None:
        records = [
            {
                "shape": "sphere",
                "target": 100,
                "h0_multiplier": 1.0,
                "fill_fraction": 0.80,
                "count_error": 1,
                "count_ratio": 1.01,
                "within_target_band": True,
                "short_circuited": False,
                "unresolved_boundary_cells": 0,
                "classified_cells": 100,
            },
            {
                "shape": "sphere",
                "target": 100,
                "h0_multiplier": 1.5,
                "fill_fraction": 0.99,
                "count_error": 100000,
                "count_ratio": 1001.0,
                "within_target_band": False,
                "short_circuited": True,
                "unresolved_boundary_cells": 10,
                "classified_cells": 200,
                "selection_warning": "increase h0 or increase target cells",
            },
        ]

        best = select_best_h0_records(records)
        self.assertEqual(best[0]["h0_multiplier"], 1.0)
        self.assertEqual(best[0]["selection_warning"], "")

    def test_best_h0_records_rank_by_count_error_before_fill_fraction(self) -> None:
        records = [
            {
                "shape": "sphere",
                "target": 200,
                "h0_multiplier": 0.75,
                "fill_fraction": 0.75,
                "count_error": 5,
                "count_ratio": 1.025,
                "within_target_band": True,
                "target_accepted": True,
                "short_circuited": False,
                "unresolved_boundary_cells": 2,
                "classified_cells": 100,
            },
            {
                "shape": "sphere",
                "target": 200,
                "h0_multiplier": 1.25,
                "fill_fraction": 0.78,
                "count_error": 10,
                "count_ratio": 0.95,
                "within_target_band": True,
                "target_accepted": True,
                "short_circuited": False,
                "unresolved_boundary_cells": 8,
                "classified_cells": 100,
            },
        ]

        best = select_best_h0_records(records)
        by_target = {record["target"]: record for record in best}
        self.assertEqual(by_target[200]["h0_multiplier"], 0.75)

    def test_best_h0_records_fallback_when_no_row_is_inside_target_band(self) -> None:
        records = [
            {
                "shape": "sphere",
                "target": 100,
                "h0_multiplier": 0.25,
                "fill_fraction": 0.99,
                "count_error": 1000,
                "count_ratio": 11.0,
                "within_target_band": False,
                "short_circuited": True,
                "unresolved_boundary_cells": 0,
                "classified_cells": 1000,
                "selection_warning": "increase h0 or increase target cells",
            },
            {
                "shape": "sphere",
                "target": 100,
                "h0_multiplier": 2.0,
                "fill_fraction": 0.40,
                "count_error": 10,
                "count_ratio": 0.90,
                "within_target_band": False,
                "short_circuited": False,
                "unresolved_boundary_cells": 20,
                "classified_cells": 100,
            },
        ]

        best = select_best_h0_records(records)
        self.assertEqual(best[0]["h0_multiplier"], 2.0)
        self.assertEqual(best[0]["selection_warning"], "no rows within target-count band")

    def test_recommended_h0_uses_table_and_falls_back_for_unknown_shape(self) -> None:
        sphere = Sphere(center=[0.0, 0.0, 0.0], radius=1.0)
        analytic = estimate_h0_from_volume(sphere, 100)
        table_multiplier = H0_REFERENCE_TABLE["sphere"][0].h0_multiplier

        self.assertAlmostEqual(recommended_h0(sphere, 100), analytic * table_multiplier)
        self.assertGreater(recommended_h0(sphere, 100), 0.0)
        self.assertAlmostEqual(recommended_h0(sphere, 100, shape_key="unknown"), analytic)

    def test_cylinder_and_hexagon_wrappers_return_metadata(self) -> None:
        cylinder = generate_cylinder(
            1.0,
            2.0,
            "z",
            120,
            h0=0.5,
            h_min=0.03125,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
        )
        hexagon = generate_hexagonal_prism(
            1.0,
            0.75,
            "z",
            120,
            rotation_degrees=15.0,
            h0=0.25,
            h_min=0.03125,
            max_depth=4,
            eps=1e-12,
            overshoot_policy="soft",
        )

        self.assertEqual(cylinder.shape, "cylinder")
        self.assertAlmostEqual(cylinder.shape_metadata["radius"], 1.0)
        self.assertAlmostEqual(cylinder.shape_metadata["length"], 2.0)
        self.assertEqual(cylinder.shape_metadata["axis"], "z")
        self.assertIn("volume", cylinder.shape_metadata)
        self.assertIn("cube_fill_stats", cylinder.refinement_metadata)

        self.assertEqual(hexagon.shape, "hexagonal_prism")
        self.assertAlmostEqual(hexagon.shape_metadata["side_length"], 1.0)
        self.assertAlmostEqual(hexagon.shape_metadata["height"], 0.75)
        self.assertEqual(hexagon.shape_metadata["axis"], "z")
        self.assertAlmostEqual(hexagon.shape_metadata["rotation_degrees"], 15.0)
        self.assertIn("vertices", hexagon.shape_metadata)
        self.assertIn("volume", hexagon.shape_metadata)


if __name__ == "__main__":
    unittest.main()
