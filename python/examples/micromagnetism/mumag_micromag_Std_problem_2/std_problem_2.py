"""muMag standard problem 2: the quasi-static hysteresis loop of a rectangular bar.

The bar has the aspect ratio 5 : 1 : 0.1 and is swept over its width d in units of the exchange
length. For every d the remanent magnetization components and the coercive field are read off the
loop and compared with the published solutions. This is the Python counterpart of
matlab/examples/Micromagnetism/mumag_micromag_Std_problem_2/Standard_problem_2.m and uses the same
grid, field schedule and damping ramp, and takes the same options.

The equilibrium at each field is found either by integrating the Landau-Lifshitz equation in time
(solver 'explicit') or by the energy minimizer (solver 'minimizer'), selected with use_minimizer.
Either way the number of effective-field evaluations spent is printed, which is the cost to compare
between the two methods.

The field sweep is either the fixed table of 40 fields of the Matlab example or, with use_adaptive,
the adaptive field stepping of MagTense, which takes large steps where the loop is smooth and
refines the step to dH_min across the coercive field, where the sign of the mean magnetization
along the field changes.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from magtense.micromag import MicromagProblem

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
MUMAG_DIR = (
    REPOSITORY_ROOT / "documentation" / "examples_mumag_validation" / "Validation_standard_problem_2"
)


def std_prob_2(
    res: tuple[int, int, int] = (100, 20, 1),
    d_loop: np.ndarray | None = None,
    cuda: bool = True,
    cvode: bool = False,
    use_minimizer: bool = True,
    use_adaptive: bool = True,
    plotting: bool = True,
    figpath: Path | None = None,
) -> dict:
    """Run standard problem 2 for a range of bar widths.

    Args:
        res: cells along x, y and z.
        d_loop: scale factors of the bar; the width is 1e-6 m times the factor. The default is
            the ten values of the Matlab example.
        cuda: use CUDA for the calculations.
        cvode: use CVODE for the numerical time evolution.
        use_minimizer: relax with the energy minimizer instead of the time integration. The
            minimizer ignores the time window and the damping ramp and stops when the largest
            torque is below problem.min_tol.
        use_adaptive: sweep the field with the adaptive stepping instead of the fixed table of
            40 fields. The step starts at the fixed-table spacing of 0.005 T, grows to at most
            0.02 T where the magnetization hardly changes, and is refined to 0.0005 T across the
            coercive field, which locates Hc ten times more precisely for a few times the cost.
        plotting: show (or save, with figpath) the remanence and the coercive field against d/l_ex
            together with the published solutions.
        figpath: directory to save the figure in. None shows it interactively.

    Returns:
        A dict with d/l_ex, the remanence components Mxr and Myr, the coercive field Hc (all in
        units of Ms), the number of effective-field evaluations per d, and the last loop.
    """
    mu0 = 4 * np.pi * 1e-7
    if d_loop is None:
        d_loop = np.linspace(0.05, 0.5, 10)
    d_loop = np.atleast_1d(np.asarray(d_loop, dtype=np.float64))

    Ms = 1000e3
    A0 = 1.74532925199e-10
    alpha = 1e3
    max_H = 0.1  # T
    n_fields = 40
    hyst_dir = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    lex = np.sqrt(A0 / (0.5 * mu0 * Ms**2))

    # The damping is ramped up from alpha to 1e5 alpha over the first 2 ns of each relaxation, which
    # kills the precession quickly and lets the time integration settle faster. alpha=0 tells the
    # solver to use the table. The minimizer does not use either.
    t_alpha = np.linspace(0, 10e-9, 100)

    problem = MicromagProblem(
        res=res,
        solver="minimizer" if use_minimizer else "explicit",
        hysteresis_solver="adaptive" if use_adaptive else "static",
        A0=A0,
        Ms=Ms,
        K0=0.0,
        alpha=0.0,
        t_alpha=t_alpha,
        alpha_fct=lambda t: alpha * 10.0 ** (5 * np.minimum(t, 2e-9) / 2e-9),
        m0=np.full((int(np.prod(res)), 3), 1 / np.sqrt(3.0)),
        cuda=cuda,
        cvode=cvode,
        # The adaptive sweep reports the fields it accepted through the returned H_ext array
        usereturnhall=use_adaptive,
    )

    # Two output times per field and a convergence check at the second one, as in the Matlab
    # example. The explicit solver integrates the 40 ns window; the check compares the end state
    # with the start state, so it only stops the integration early when nothing moved at all.
    problem.t = np.linspace(0, 40e-9, 2)
    problem.nt = 2
    problem.t_conv = problem.t.copy()
    problem.nt_conv = 2
    problem.conv_tol = np.repeat(1e-6, 2)

    # The field table: one row per constant field, walking from +max_H to -max_H along hyst_dir.
    # Column 0 carries the signed field magnitude in T, which the solver ignores but which is
    # convenient for reading the loop back.
    H_T = np.linspace(max_H, -max_H, n_fields)
    H_ext = np.zeros((n_fields, 4))
    H_ext[:, 0] = H_T
    H_ext[:, 1:4] = (H_T / mu0)[:, None] * hyst_dir[None, :]
    H_norm = H_T / (mu0 * Ms)  # signed applied field in units of Ms

    results = {
        "dlex": np.zeros(len(d_loop)), "Mxr": np.zeros(len(d_loop)), "Myr": np.zeros(len(d_loop)),
        "Hc": np.zeros(len(d_loop)), "n_feval": np.zeros(len(d_loop), dtype=int),
        "n_fields": np.zeros(len(d_loop), dtype=int), "H_norm": H_norm, "M": None,
    }

    for i, d in enumerate(d_loop):
        problem.grid_L = np.array([5e-6, 1e-6, 1e-7]) * d
        results["dlex"][i] = problem.grid_L[1] / lex
        print(f"Running d/l_ex = {results['dlex'][i]:.2f}, i.e. {i + 1}/{len(d_loop)}")

        if use_adaptive:
            result = problem.run_hysteresis_adaptive(
                H_start=(max_H / mu0) * hyst_dir,
                H_end=(-max_H / mu0) * hyst_dir,
                dH_initial=0.005 / mu0,
                dH_min=0.0005 / mu0,
                dH_max=0.02 / mu0,
                max_steps=400,
                # The step-control thresholds on the change of the mean magnetization. The defaults
                # suit a hard grain that barely moves between switching events; this soft bar
                # changes by a few percent per step everywhere, so looser thresholds keep the step
                # at the table spacing in the smooth parts and leave the refinement to the switch
                # test. At 40 x 8 x 1 cells this gives about 50 fields against 300 with the defaults,
                # with the same coercive field to four digits.
                dM_min=0.01,
                dM_target=0.05,
                dM_reject=0.15,
                switch_refine_dH=0.0005 / mu0,
            )
            n_acc = int(result[-1])
            M_out = result[1]
            # The accepted fields, signed along the sweep direction, in T and in units of Ms
            H_acc = np.asarray(result[4][0, 0, :n_acc, :]) @ hyst_dir
            H_T = mu0 * H_acc
            H_norm = H_acc / Ms
        else:
            M_out = problem.run_hysteresis(H_ext)[1]
            n_acc = n_fields
        results["n_fields"][i] = n_acc

        # The reduced magnetization at the end of each relaxation, averaged over the cells
        m = M_out[-1]  # (ntot, n_fields, 3)
        m = m / np.linalg.norm(m, axis=2, keepdims=True)
        m_avg = m.mean(axis=0)  # (n_fields, 3)
        M_par = m_avg @ hyst_dir

        method = "minimizer" if use_minimizer else "LL relaxation"
        results["n_feval"][i] = int(problem.n_feval.sum())
        extra = (
            f", {int(problem.min_iter.sum())} iterations, "
            f"{int(np.sum(problem.min_status == 2))} fields not converged"
            if use_minimizer else ""
        )
        sweep = f"{n_acc} adaptive fields" if use_adaptive else f"{n_acc} fixed fields"
        print(f"   {method}, {sweep}: {results['n_feval'][i]} field evaluations{extra}")

        # Remanence at H = 0 and coercive field at M = 0. np.interp needs increasing abscissae and
        # the sweep runs downwards, so the arrays are reversed.
        results["Mxr"][i] = np.interp(0.0, H_norm[::-1], m_avg[::-1, 0])
        results["Myr"][i] = np.interp(0.0, H_norm[::-1], m_avg[::-1, 1])
        results["Hc"][i] = np.interp(0.0, M_par[::-1], H_norm[::-1])
        results["M"] = M_par
        results["H_norm"] = H_norm

    if plotting:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        axes[0].plot(results["dlex"], results["Mxr"], "k.", markersize=12, label="MagTense")
        axes[1].plot(results["dlex"], results["Myr"], "k.", markersize=12, label="MagTense")
        axes[2].plot(results["dlex"], np.abs(results["Hc"]), "k.", markersize=12, label="MagTense")
        for name in ("Streibl", "McMichael", "Lopez-Diaz", "Donahue"):
            data = np.loadtxt(MUMAG_DIR / f"{name}.txt", comments="%")
            axes[0].plot(data[:, 0], data[:, 1], "d", label=name)
            axes[1].plot(data[:, 0], data[:, 2], "d", label=name)
            axes[2].plot(data[:, 0], data[:, 3], "d", label=name)
        for ax, label in zip(axes, ("$M_{xr}/M_s$", "$M_{yr}/M_s$", "$|H_c|/M_s$"), strict=True):
            ax.set_xlabel("$d/l_{ex}$")
            ax.set_ylabel(label)
            ax.grid(True)
        axes[2].legend()
        fig.suptitle(
            "Standard problem 2, " + ("energy minimizer" if use_minimizer else "LL relaxation")
            + (", adaptive field steps" if use_adaptive else "")
        )
        fig.tight_layout()
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            fig.savefig(figpath / "2_remanence_coercivity.png")

    return results


if __name__ == "__main__":
    std_prob_2(
        cuda=True,
        cvode=False,
        use_minimizer=True,   # set False to relax by integrating the LL equation in time instead
        use_adaptive=True,    # set False to sweep the fixed table of 40 fields instead
        plotting=True,
        figpath=None,
    )
