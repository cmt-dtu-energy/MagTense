"""Tests for reusable single-grain result analysis helpers."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

SINGLE_GRAIN_DIR = Path(__file__).resolve().parents[1] / "experiments" / "single_grain"
sys.path.insert(0, str(SINGLE_GRAIN_DIR))

from utils.magnetization import (  # noqa: E402
    angle_between_vectors_deg,
    extract_tile_magnetization,
    extract_tile_positions,
)
from utils.metrics import MU0, calculate_hysteresis_metrics  # noqa: E402
from utils.results import (  # noqa: E402
    discover_result_files,
    write_summary_csv,
    load_result,
    load_results,
)
from utils.reporting import material_summary  # noqa: E402
import single_grain_coercivity as experiment  # noqa: E402


class MetricTests(unittest.TestCase):
    def test_known_linear_demagnetization_curve(self) -> None:
        """Interpolation and segment optimization reproduce analytic values."""
        H = np.array([1.0e5, 0.0, -5.0e5, -1.0e6])
        M = np.array([1.0e6, 1.0e6, 5.0e5, 0.0])
        metrics = calculate_hysteresis_metrics(H, M, 1.0e6)

        self.assertAlmostEqual(metrics.Hc_A_per_m, -1.0e6)
        self.assertAlmostEqual(metrics.Mr_A_per_m, 1.0e6)
        self.assertAlmostEqual(metrics.Mr_over_Ms, 1.0)
        expected_bh = MU0 * 1.25e11
        self.assertAlmostEqual(metrics.BH_max_J_per_m3, expected_bh)
        self.assertEqual(metrics.coercivity_status, "ok")
        self.assertEqual(metrics.remanence_status, "ok")
        self.assertEqual(metrics.energy_product_status, "ok")

    def test_missing_crossings_have_diagnostics(self) -> None:
        H = np.array([2.0, 1.0])
        M = np.array([3.0, 2.0])
        metrics = calculate_hysteresis_metrics(H, M, 4.0)
        self.assertTrue(np.isnan(metrics.Hc_A_per_m))
        self.assertTrue(np.isnan(metrics.Mr_A_per_m))
        self.assertTrue(np.isnan(metrics.BH_max_J_per_m3))
        self.assertEqual(metrics.coercivity_status, "no_crossing")
        self.assertEqual(metrics.remanence_status, "no_crossing")
        self.assertEqual(metrics.energy_product_status, "no_second_quadrant")

    def test_insufficient_points_are_reported(self) -> None:
        metrics = calculate_hysteresis_metrics(np.array([0.0]), np.array([1.0]), 1.0)
        self.assertEqual(metrics.coercivity_status, "insufficient_points")
        self.assertEqual(metrics.remanence_status, "insufficient_points")
        self.assertEqual(metrics.energy_product_status, "insufficient_points")


class ResultLoadingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _write_result(
        self,
        path: Path,
        *,
        L: float = 40e-9,
        Ms: float = 1.0e6,
        **extra: object,
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        H = np.array([1.0e5, 0.0, -5.0e5, -1.0e6])
        M = np.array([Ms, Ms, 0.5 * Ms, 0.0])
        np.savez(
            path,
            H_array=MU0 * H,
            H_array_A_per_m=H,
            M_array=M,
            Mx_array=np.zeros_like(M),
            My_array=np.zeros_like(M),
            Mz_array=M,
            L=L,
            size_factor=4.0,
            n=12,
            ntot=12**3,
            Ms=Ms,
            A0=7e-12,
            K0=8.9e6,
            runtime=2.5,
            **extra,
        )

    def test_batch_manifest_and_physical_size_are_loaded(self) -> None:
        batch = self.root / "res_sm" / "batch_1"
        result = batch / "results" / "single_grain_size40nm_n12_P_A_FS1.0e-03.npz"
        self._write_result(result)
        (batch / "manifest.json").write_text(
            json.dumps(
                {
                    "preset": "sm2fe17n3",
                    "material": {
                        "mu0_Ms_T": 1.54,
                        "A0_J_per_m": 7e-12,
                        "K0_J_per_m3": 8.9e6,
                    },
                    "batch": {"label": "batch_1"},
                }
            ),
            encoding="utf-8",
        )

        record = load_result(result)
        self.assertEqual(record.material, "sm2fe17n3")
        self.assertEqual(record.batch, "batch_1")
        self.assertEqual(record.shape, "cube")
        self.assertEqual(record.shape_variant, "cube")
        self.assertAlmostEqual(record.size_nm, 40.0)
        self.assertTrue(record.periodic)
        self.assertTrue(record.adaptive)
        self.assertAlmostEqual(record.adaptive_dh_min_t, 0.001)
        self.assertTrue(np.isfinite(record.metrics.BH_max_kJ_per_m3))

    def test_legacy_filename_and_material_fallback(self) -> None:
        result = self.root / "results" / "single_grain_sf4_n12_A_FS1.0e-02.npz"
        self._write_result(result)
        record = load_result(result)
        self.assertFalse(record.periodic)
        self.assertEqual(record.preset, None)
        self.assertIn("Ms=", record.material)
        self.assertEqual(record.shape, "cube")
        self.assertEqual(record.shape_variant, "cube")
        self.assertAlmostEqual(record.adaptive_dh_min_t, 0.01)

    def test_shape_fields_are_loaded_from_result_file(self) -> None:
        result = self.root / "results" / "single_grain_sphere_size40nm_n12_A_FS1.0e-03.npz"
        self._write_result(
            result,
            shape="sphere",
            shape_variant="sphere",
            shape_metadata=np.array({"radius": 20e-9}, dtype=object),
        )

        record = load_result(result)
        self.assertEqual(record.shape, "sphere")
        self.assertEqual(record.shape_variant, "sphere")
        self.assertEqual(record.shape_metadata, {"radius": 20e-9})
        self.assertIn(" / sphere", record.comparison_label)

    def test_shape_falls_back_to_manifest_and_shape_directory(self) -> None:
        manifest_batch = self.root / "res_fe" / "batch_1"
        manifest_result = manifest_batch / "results" / "single_grain_size40nm_n12_A_FS1.0e-03.npz"
        self._write_result(manifest_result)
        (manifest_batch / "manifest.json").write_text(
            json.dumps({"shape": "ellipsoid_x", "batch": {"label": "batch_1"}}),
            encoding="utf-8",
        )

        path_result = (
            self.root
            / "res_fe"
            / "shape_ellipsoid_z"
            / "batch_2"
            / "results"
            / "single_grain_size40nm_n12_A_FS1.0e-03.npz"
        )
        self._write_result(path_result)

        manifest_record = load_result(manifest_result)
        path_record = load_result(path_result)
        self.assertEqual(manifest_record.shape_variant, "ellipsoid_x")
        self.assertEqual(path_record.shape_variant, "ellipsoid_z")

    def test_summary_and_export_include_shape_fields(self) -> None:
        result = self.root / "results" / "single_grain_sphere_size40nm_n12_A_FS1.0e-03.npz"
        self._write_result(result, shape="sphere", shape_variant="sphere")
        record = load_result(result)

        rows = material_summary([record])
        self.assertEqual(rows[0]["shape"], "sphere")
        self.assertEqual(rows[0]["shape_variant"], "sphere")
        self.assertEqual(rows[0]["comparison_label"], record.comparison_label)

        csv_path = write_summary_csv([record], self.root / "summary.csv")
        csv_text = csv_path.read_text(encoding="utf-8")
        self.assertIn("shape,shape_variant,comparison_label", csv_text)
        self.assertIn(",sphere,sphere,", csv_text)

    def test_discovery_deduplicates_overlapping_roots_and_reports_bad_files(self) -> None:
        valid = self.root / "results" / "single_grain_valid.npz"
        invalid = self.root / "results" / "single_grain_invalid.npz"
        self._write_result(valid)
        invalid.write_text("not an npz", encoding="utf-8")
        paths = discover_result_files([self.root, self.root / "results"])
        self.assertEqual(paths, [invalid, valid])
        records, errors = load_results(self.root)
        self.assertEqual(len(records), 1)
        self.assertEqual(errors[0][0], invalid)


class MagnetizationExtractionTests(unittest.TestCase):
    def test_extracts_standard_four_dimensional_layout(self) -> None:
        magnetization = np.arange(2 * 3 * 4 * 3, dtype=float).reshape(2, 3, 4, 3)
        positions = np.arange(1 * 3 * 1 * 3, dtype=float).reshape(1, 3, 1, 3)
        res = np.empty(3, dtype=object)
        res[0] = None
        res[1] = magnetization
        res[2] = positions
        np.testing.assert_array_equal(extract_tile_magnetization(res, 2), magnetization[0, :, 2, :])
        np.testing.assert_array_equal(extract_tile_positions(res), positions[0, :, 0, :])

    def test_vector_angles(self) -> None:
        first = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        second = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]])
        np.testing.assert_allclose(angle_between_vectors_deg(first, second), [90.0, 180.0])


class SimulationOutputTests(unittest.TestCase):
    def test_run_preserves_old_keys_and_adds_derived_metrics(self) -> None:
        """The CLI implementation writes the extended schema without a solver run."""

        class FakeProblem:
            def run_hysteresis(self, H_ext):
                return [H_ext]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mean_m = (
                np.zeros(3),
                np.zeros(3),
                np.array([1.0e6, 5.0e5, 0.0]),
            )
            with (
                patch.object(experiment, "create_single_grain_problem", return_value=FakeProblem()),
                patch.object(experiment, "extract_mean_magnetisation", return_value=mean_m),
            ):
                experiment.run_single_grain_coercivity(
                    L=40e-9,
                    size_factor=4.0,
                    n=2,
                    field_max_t=0.0,
                    field_min_t=-1.0,
                    field_step_t=-0.5,
                    output_dir=root / "results",
                    timer_log_dir=root / "timers",
                    plotting=False,
                    Ms=1.0e6,
                    output_stem_override="schema_test",
                )

            with np.load(root / "results" / "schema_test.npz", allow_pickle=True) as data:
                for old_key in (
                    "res", "H_array", "H_array_A_per_m", "M_array",
                    "Mx_array", "My_array", "Mz_array", "Hc", "Hc_T",
                    "runtime", "L", "n", "ntot", "size_factor", "Ms", "K0", "A0",
                ):
                    self.assertIn(old_key, data.files)
                for new_key in (
                    "Mr_A_per_m", "mu0_Mr_T", "Mr_over_Ms",
                    "BH_max_J_per_m3", "BH_max_kJ_per_m3", "BH_max_MGOe",
                    "coercivity_status", "remanence_status", "energy_product_status",
                    "adaptive", "periodic", "adaptive_dh_min_t",
                ):
                    self.assertIn(new_key, data.files)


if __name__ == "__main__":
    unittest.main()
