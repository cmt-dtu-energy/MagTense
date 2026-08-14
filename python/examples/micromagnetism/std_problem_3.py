from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from magtense.micromag import MicromagProblem
from magtense.utils import plot_M_avg_seq, plot_M_thin_film
from magtense.computeMagneticEnergy import compute_magnetic_energy

def std_prob_3(
    res: tuple[int, int, int] = (10, 10, 10),
    L_loop: np.ndarray | None = None,
    cuda: bool = False,
    cvode: bool = False,
    plotting: bool = True,
    figpath: Path | None = None,
    plot_details: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Run the muMag standard problem 3 for a range of cube sizes.

    Returns:
        L_loop: the cube edge lengths in units of the exchange length that were simulated.
        E_arr: the equilibrium energy terms with shape (4, len(L_loop), 2), where the first axis
            is (exchange, external, demagnetisation, anisotropy) and the last axis selects the
            flower state (0) or the vortex state (1). The single domain limit is the L where the
            two total energies cross.
    """
    mu0 = 4 * np.pi * 1e-7
    A0 = 1.74532925199e-10
    Ms = 1e6
    if L_loop is None:
        L_loop = np.linspace(8, 9, 10)

    problem = MicromagProblem(
        res=res,
        A0=A0,
        Ms=Ms,
        K0=0.1 * 0.5 * mu0 * Ms**2,
        alpha=1e3,
        cuda=cuda,
        cvode=cvode,
        usereturnhall=True ,
    )

    #--------- disable fmm -----
    problem.use_fmm = 0
    #----------------------------
    
    #--------------- set trace/timing options -------------
    problem.use_fmm = 0
    problem.window_enabled = 0
    problem.window_interval = 30.0
    problem.trace_enabled = 0
    problem.flush_each = 1
    problem.trace_verbose = 2
    problem.timer_log_file =  "std_3_timer.log"
    problem.trace_log_file =  "std_3_trace.log"
    #------------------------------------------------------

    problem.u_ea[:, 2] = 1
    lex = np.sqrt(A0 / (0.5 * mu0 * Ms**2))

    def h_ext_fct(t) -> np.ndarray:
        return np.atleast_2d(t).T * np.array([0, 0, 0])

    E_arr = np.zeros(shape=(4, len(L_loop), 2))

    for i in range(len(L_loop)):
        print(f"ITERATION: {i} / {len(L_loop)}")
        for j in range(2):
            if j == 0:
                print("Flower state")
                problem.m0[:, 0:2] = 0
                problem.m0[:, 2] = 1
                t_end = 10e-9

            else:
                print("Vortex state")
                xv = np.linspace(-1, 1, res[0])
                yv = np.linspace(-1, 1, res[1])
                zv = np.linspace(-1, 1, res[2])
                [x, y, z] = np.meshgrid(xv, yv, zv, indexing="ij")
                xvec = np.sin(np.arctan2(z, x))
                yvec = -np.cos(np.arctan2(z, x))
                problem.m0[:, 0] = xvec.swapaxes(0, 2).reshape(-1)
                problem.m0[:, 2] = yvec.swapaxes(0, 2).reshape(-1)
                problem.m0 = problem.m0 / np.tile(
                    np.expand_dims(np.sqrt(np.sum(problem.m0**2, axis=1)), axis=1),
                    (1, 3),
                )
                t_end = 200e-9

            problem.grid_L = np.array([lex, lex, lex]) * L_loop[i]
            tile_volumes = np.full(int(np.prod(res)), np.prod(problem.grid_L) / np.prod(res))

            t, M_out, _, H_exc, H_ext, H_dem, H_ani = (
                problem.run_simulation(
                    t_end=t_end, nt=50, fct_h_ext=h_ext_fct, nt_h_ext=2
                )
            )[:7]

            if plot_details:
                M_sq = np.squeeze(M_out, axis=2)
                plot_M_avg_seq(t, M_sq, figpath=figpath)
                plot_M_thin_film(M_sq[0], res, title="3_start", figpath=figpath)
                plot_M_thin_film(M_sq[-1], res, title="3_end", figpath=figpath)

            # Calculate the energy terms
            E_exc, E_ext, E_dem, E_ani = compute_magnetic_energy(
                M=M_out, H_exc=H_exc, H_ext=H_ext, H_dem=H_dem, H_ani=H_ani, Vols=tile_volumes, problem=problem, Ms=Ms
            )

            E_arr[:, i, j] = np.array([E_exc[-1].item(), E_ext[-1].item(), E_dem[-1].item(), E_ani[-1].item()])

            if plot_details:
                plt.clf()
                for E_x in [E_exc, E_ext, E_dem, E_ani]:
                    plt.plot(t, mu0 * E_x - mu0 * E_x[0], ".")
                plt.xlabel("Time [s]")
                plt.ylabel("Energy [-]")
                plt.legend([r"$E_{exc}$", r"$E_{ext}$", r"$E_{dem}$", r"$E_{ani}$"])
                if figpath is None:
                    plt.show()
                else:
                    figpath.mkdir(parents=True, exist_ok=True)
                    plt.savefig(figpath / "3_details.png")

    if plotting:
        plt.clf()
        plt.plot(L_loop, np.sum(E_arr[:, :, 0], axis=0), ".", color="blue", markersize=8)
        plt.plot(L_loop, np.sum(E_arr[:, :, 1], axis=0), ".", color="orange", markersize=8)
        plt.xlabel("L [l_ex]")
        plt.ylabel("Energy [-]")
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / "3_solution.png")

    return L_loop, E_arr


if __name__ == "__main__":
    std_prob_3(
        cuda=True,
        cvode=False,
        plotting=True,
        plot_details=False,
        #figpath=Path(__file__).parent.absolute().joinpath("..", "figs"),
        figpath=None,
    )
