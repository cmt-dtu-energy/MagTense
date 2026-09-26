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
* The energies and relaxation diagnostics that Fortran returns after the historical
  outputs land on the problem object, sliced to the accepted steps for adaptive runs.
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

    ``RunMicroMagSimulation`` returns 19 values. Index 7 is the accepted-field
    count, indices 8-13 contain the exchange-matrix data that historically started
    at index 7 for static simulations, and indices 14-18 are the energies and the
    relaxation diagnostics, which the wrapper moves onto the problem object.
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
        np.full((2, n_fields, 4), 90.0),          # energies
        np.arange(n_fields) + 100,                # field evaluations
        np.arange(n_fields) + 200,                # minimizer iterations
        np.full(n_fields, 1e-6),                  # final torque
        np.full(n_fields, -1),                    # status
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
        problem = _problem("static")
        with patch.object(micromag_module, "magtensesource", fake_source):
            result = problem.run_hysteresis(h_ext)

        self.assertEqual(captured["hysteresis_solver"], 1)
        # The minimizer settings travel with every call, at their defaults here
        self.assertEqual(captured["min_tol"], 1e-5)
        self.assertEqual(captured["min_maxiter"], 10000)
        self.assertEqual(captured["min_maxrot"], 0.3)
        self.assertEqual(captured["min_fallback"], 1)
        self.assertEqual(captured["min_saddle_check"], 1)
        self.assertEqual(captured["min_predictor"], 0)
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

        # The diagnostics are not in the list but on the problem, unsliced for a static run
        self.assertEqual(problem.E_out.shape, (2, 3, 4))
        np.testing.assert_array_equal(problem.n_feval, [100, 101, 102])
        np.testing.assert_array_equal(problem.min_iter, [200, 201, 202])
        np.testing.assert_array_equal(problem.min_status, [-1, -1, -1])

    def test_adaptive_hysteresis_uses_unified_fortran_entry_point(self) -> None:
        """Adaptive mode forwards controls and slices all field-dependent data."""
        captured = {}

        def fake_run(**kwargs):
            captured.update(kwargs)
            return _fortran_result(n_fields=4, n_accepted=2)

        fake_source = SimpleNamespace(
            fortrantopythonio=SimpleNamespace(runmicromagsimulation=fake_run)
        )
        problem = _problem("adaptive")
        with patch.object(micromag_module, "magtensesource", fake_source):
            result = problem.run_hysteresis_adaptive(
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
        # The Fortran loop stores the starting field in slot 1 and then one slot per accepted
        # step, so it needs max_steps + 1 slots. Sizing these to max_steps overran
        # problem%Hext by one row on the last step and made the copy-out of M/H a
        # non-conforming array assignment.
        self.assertEqual(captured["nt_hext"], 5)
        self.assertEqual(captured["nt_hext_out"], 5)
        self.assertEqual(captured["hext"].shape, (5, 4))
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

        # The diagnostics are sliced to the accepted steps as well
        self.assertEqual(problem.E_out.shape, (2, 2, 4))
        np.testing.assert_array_equal(problem.n_feval, [100, 101])
        np.testing.assert_array_equal(problem.min_torque, [1e-6, 1e-6])

    def test_minimizer_settings_are_stored(self) -> None:
        """The minimizer settings land on the problem as given."""
        problem = MicromagProblem(
            res=[1, 1, 1], solver="explicit", min_tol=2e-6, min_maxiter=50,
            min_maxrot=0.1, min_fallback=False, min_saddle_check=False, min_predictor=True,
        )
        self.assertEqual(problem.min_tol, 2e-6)
        self.assertEqual(problem.min_maxiter, 50)
        self.assertEqual(problem.min_maxrot, 0.1)
        self.assertEqual(problem.min_fallback, 0)
        self.assertEqual(problem.min_saddle_check, 0)
        self.assertEqual(problem.min_predictor, 1)

    def test_explicit_uses_the_minimizer(self) -> None:
        """'explicit' relaxes with the minimizer, 'explicit_ll' with the LL integration."""
        self.assertEqual(MicromagProblem(res=[1, 1, 1]).solver, 2)
        self.assertEqual(MicromagProblem(res=[1, 1, 1], solver="explicit").solver, 3)
        self.assertEqual(MicromagProblem(res=[1, 1, 1], solver="explicit_ll").solver, 1)
        for removed in ("minimizer", "implicit"):
            with self.assertRaisesRegex(ValueError, "explicit_ll"):
                MicromagProblem(res=[1, 1, 1], solver=removed)

        captured = {}

        def fake_run(**kwargs):
            captured.update(kwargs)
            return _fortran_result(n_fields=2, n_accepted=0)

        fake_source = SimpleNamespace(
            fortrantopythonio=SimpleNamespace(runmicromagsimulation=fake_run)
        )
        h_ext = np.zeros((2, 4))
        with patch.object(micromag_module, "magtensesource", fake_source):
            _problem().run_hysteresis(h_ext)
            self.assertEqual(captured["solver"], 3)

            # The minimizer cannot include the thermal field, so a finite temperature, even one
            # set after the solver, stops 'explicit' before Fortran is called and points the
            # user to 'explicit_ll' ...
            captured.clear()
            problem = _problem()
            problem.T = 300.0
            with self.assertRaisesRegex(ValueError, "solver='explicit_ll'"):
                problem.run_hysteresis(h_ext)
            self.assertEqual(captured, {})

            # ... which accepts it.
            problem.solver = "explicit_ll"
            problem.run_hysteresis(h_ext)
            self.assertEqual(captured["solver"], 1)

    def test_adaptive_accepts_the_minimizer(self) -> None:
        """The adaptive field stepping works with either equilibrium solver."""
        captured = {}

        def fake_run(**kwargs):
            captured.update(kwargs)
            return _fortran_result(n_fields=3, n_accepted=1)

        fake_source = SimpleNamespace(
            fortrantopythonio=SimpleNamespace(runmicromagsimulation=fake_run)
        )
        problem = MicromagProblem(
            res=[1, 1, 1], solver="explicit", hysteresis_solver="adaptive", min_tol=3e-6
        )
        problem.nt = 2
        problem.t = np.linspace(0.0, 1e-9, problem.nt)
        with patch.object(micromag_module, "magtensesource", fake_source):
            problem.run_hysteresis_adaptive(
                H_start=np.array([0.0, 0.0, 1.0]),
                H_end=np.array([0.0, 0.0, -1.0]),
                dH_initial=0.5,
                dH_min=0.1,
                dH_max=1.0,
                max_steps=2,
            )
        self.assertEqual(captured["solver"], 3)
        self.assertEqual(captured["min_tol"], 3e-6)

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
        with self.assertRaisesRegex(ValueError, "solver='explicit' or 'explicit_ll'"):
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
