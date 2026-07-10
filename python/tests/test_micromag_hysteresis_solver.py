"""Regression tests for the unified micromagnetic Fortran entry point.

These tests focus on the Python/Fortran interface rather than the numerical
solution. The compiled routine is replaced with a small fake so the tests can
inspect every argument sent by ``MicromagProblem`` and emulate the exact f2py
return order without running a full micromagnetic simulation.

The suite verifies four public contracts:

* ``hysteresis_solver`` maps ``static`` and ``adaptive`` to the Fortran flags.
* Both Python hysteresis methods call the same ``runmicromagsimulation`` symbol.
* Static callers retain the historical 13-item result layout.
* Adaptive callers receive only accepted field steps plus the accepted count.
"""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import magtense.micromag as micromag_module
from magtense.micromag import MicromagProblem


def _problem(hysteresis_solver: str = "static") -> MicromagProblem:
    """Create the smallest problem needed to exercise wrapper logic."""
    problem = MicromagProblem(
        res=[1, 1, 1],
        solver="explicit",
        hysteresis_solver=hysteresis_solver,
    )
    problem.nt = 2
    problem.t = np.linspace(0.0, 1e-9, problem.nt)
    problem.nt_conv = 2
    problem.t_conv = problem.t.copy()
    return problem


def _fortran_result(n_fields: int, n_accepted: int) -> list:
    """Build a recognizable result in the unified f2py output order.

    ``RunMicroMagSimulation`` now always returns 14 values. Index 7 is the new
    accepted-field count, while indices 8-13 contain the exchange-matrix data
    that historically started at index 7 for static simulations.
    """
    field_shape = (2, 1, n_fields, 3)
    return [
        np.full(2, 10.0),                         # t_out
        np.full(field_shape, 11.0),               # M_out
        np.full((1, 3), 12.0),                    # points
        np.full(field_shape, 13.0),               # exchange field
        np.full(field_shape, 14.0),               # external field
        np.full(field_shape, 15.0),               # demagnetizing field
        np.full(field_shape, 16.0),               # anisotropy field
        n_accepted,
        2,                                        # exchange nonzero count
        np.arange(12) + 20,                       # exchange rows
        np.arange(12) + 40,                       # exchange columns
        np.arange(12, dtype=float) + 60.0,        # exchange values
        70,                                       # exchange row dimension
        80,                                       # exchange column dimension
    ]


class HysteresisSolverTests(unittest.TestCase):
    def test_defaults_to_static_and_validates(self) -> None:
        """Static is backward compatible and invalid mode names fail early."""
        self.assertEqual(MicromagProblem(res=[1, 1, 1]).hysteresis_solver, 1)
        self.assertEqual(_problem("static").hysteresis_solver, 1)
        self.assertEqual(_problem("adaptive").hysteresis_solver, 2)
        with self.assertRaisesRegex(ValueError, "static.*adaptive"):
            MicromagProblem(res=[1, 1, 1], hysteresis_solver="invalid")

    def test_static_hysteresis_uses_unified_fortran_entry_point(self) -> None:
        """Static mode passes inert adaptive inputs and restores old ordering."""
        captured = {}

        def fake_run(**kwargs):
            captured.update(kwargs)
            return _fortran_result(n_fields=3, n_accepted=0)

        fake_source = SimpleNamespace(
            fortrantopythonio=SimpleNamespace(runmicromagsimulation=fake_run)
        )
        h_ext = np.arange(12, dtype=np.float64).reshape(3, 4)
        with patch.object(micromag_module, "magtensesource", fake_source):
            result = _problem("static").run_hysteresis(h_ext)

        self.assertEqual(captured["hysteresis_solver"], 1)
        np.testing.assert_array_equal(captured["hext"], h_ext)
        self.assertEqual(captured["nt_hext_out"], 3)
        self.assertEqual(captured["maxhextsteps"], 0)
        self.assertEqual(captured["use_switch_refine"], 0)
        np.testing.assert_array_equal(captured["h_start"], np.zeros(3))
        np.testing.assert_array_equal(captured["h_end"], np.zeros(3))

        # The internal accepted-count value is removed for static callers. All
        # following exchange values shift back to their historical positions.
        self.assertEqual(len(result), 13)
        self.assertEqual(result[7], 2)
        np.testing.assert_array_equal(result[8], [20, 21])
        np.testing.assert_array_equal(result[9], [40, 41])
        np.testing.assert_array_equal(result[10], [60.0, 61.0])
        self.assertEqual(result[11], 70)
        self.assertEqual(result[12], 80)

    def test_adaptive_hysteresis_uses_unified_fortran_entry_point(self) -> None:
        """Adaptive mode forwards controls and slices all field-dependent data."""
        captured = {}

        def fake_run(**kwargs):
            captured.update(kwargs)
            return _fortran_result(n_fields=4, n_accepted=2)

        fake_source = SimpleNamespace(
            fortrantopythonio=SimpleNamespace(runmicromagsimulation=fake_run)
        )
        with patch.object(micromag_module, "magtensesource", fake_source):
            result = _problem("adaptive").run_hysteresis_adaptive(
                H_start=np.array([0.0, 0.0, 1.0]),
                H_end=np.array([0.0, 0.0, -1.0]),
                dH_initial=0.5,
                dH_min=0.1,
                dH_max=1.0,
                max_steps=4,
                dM_min=0.002,
                dM_target=0.02,
                dM_reject=0.08,
                dH_grow=1.8,
                dH_shrink=0.6,
                switch_refine_dH=0.15,
            )

        self.assertEqual(captured["hysteresis_solver"], 2)
        self.assertEqual(captured["maxhextsteps"], 4)
        self.assertEqual(captured["nt_hext_out"], 4)
        np.testing.assert_array_equal(captured["h_start"], [0.0, 0.0, 1.0])
        np.testing.assert_array_equal(captured["h_end"], [0.0, 0.0, -1.0])
        self.assertEqual(captured["dh_initial"], 0.5)
        self.assertEqual(captured["dh_min"], 0.1)
        self.assertEqual(captured["dh_max"], 1.0)
        self.assertEqual(captured["dm_min"], 0.002)
        self.assertEqual(captured["dm_target"], 0.02)
        self.assertEqual(captured["dm_reject"], 0.08)
        self.assertEqual(captured["dh_grow"], 1.8)
        self.assertEqual(captured["dh_shrink"], 0.6)
        self.assertEqual(captured["switch_refine_dh"], 0.15)
        self.assertEqual(captured["use_switch_refine"], 1)

        # Fortran allocates four field slots, but reports that only two were
        # accepted. Every field-dependent output must therefore be sliced to 2.
        self.assertEqual(len(result), 14)
        self.assertEqual(result[-1], 2)
        for index in (1, 3, 4, 5, 6):
            self.assertEqual(result[index].shape[2], 2)
        np.testing.assert_array_equal(result[8], [20, 21])
        np.testing.assert_array_equal(result[9], [40, 41])
        np.testing.assert_array_equal(result[10], [60.0, 61.0])

    def test_methods_reject_the_wrong_mode(self) -> None:
        """Each public method rejects a problem configured for the other mode."""
        with self.assertRaisesRegex(ValueError, "hysteresis_solver='static'"):
            _problem("adaptive").run_hysteresis(np.zeros((2, 4)))
        with self.assertRaisesRegex(ValueError, "hysteresis_solver='adaptive'"):
            _problem("static").run_hysteresis_adaptive(
                H_start=np.array([0.0, 0.0, 1.0]),
                H_end=np.array([0.0, 0.0, -1.0]),
                dH_initial=0.5,
                dH_min=0.1,
                dH_max=1.0,
                max_steps=4,
            )

    def test_adaptive_requires_explicit_time_solver(self) -> None:
        """Adaptive field stepping remains restricted to equilibrium solves."""
        problem = MicromagProblem(
            res=[1, 1, 1],
            solver="dynamic",
            hysteresis_solver="adaptive",
        )
        with self.assertRaisesRegex(ValueError, "explicit solver"):
            problem.run_hysteresis_adaptive(
                H_start=np.array([0.0, 0.0, 1.0]),
                H_end=np.array([0.0, 0.0, -1.0]),
                dH_initial=0.5,
                dH_min=0.1,
                dH_max=1.0,
                max_steps=4,
            )

    def test_adaptive_input_validation_happens_before_fortran(self) -> None:
        """Malformed adaptive controls fail in Python before compiled code runs."""
        problem = _problem("adaptive")
        valid = {
            "H_start": np.array([0.0, 0.0, 1.0]),
            "H_end": np.array([0.0, 0.0, -1.0]),
            "dH_initial": 0.5,
            "dH_min": 0.1,
            "dH_max": 1.0,
            "max_steps": 4,
        }
        invalid_cases = (
            ({"H_start": np.zeros(2)}, "3-vectors"),
            ({"max_steps": 0}, "max_steps must be positive"),
            ({"dH_initial": 0.0}, "must be positive"),
            ({"dH_min": 2.0}, "smaller than or equal"),
            ({"H_end": valid["H_start"]}, "must differ"),
            ({"dH_grow": 1.0}, "larger than 1"),
            ({"dH_shrink": 1.0}, "between 0 and 1"),
            ({"switch_refine_dH": 0.0}, "must be positive"),
        )
        for overrides, message in invalid_cases:
            with self.subTest(overrides=overrides):
                arguments = {**valid, **overrides}
                with self.assertRaisesRegex(ValueError, message):
                    problem.run_hysteresis_adaptive(**arguments)


if __name__ == "__main__":
    unittest.main()
