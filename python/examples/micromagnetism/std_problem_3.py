# %%
import numpy as np
import matplotlib.pyplot as plt

from typing import List

from magtense.micromag import MicromagProblem
from magtense.utils import plot_M_avg_seq, plot_M_thin_film


def std_prob_3(
    res: List[int] = [10, 10, 10],
    L_loop: np.ndarray = np.linspace(8, 9, 10),
    cuda: bool = False,
    show: bool = True,
    show_details: bool = False,
) -> None:
    mu0 = 4 * np.pi * 1e-7
    Ms = 1e6

    problem = MicromagProblem(
        res,
        A0=1.74532925199e-10,
        Ms=Ms,
        K0=0.1 * 0.5 * mu0 * Ms**2,
        alpha=1e3,
        cuda=cuda,
    )

    problem.u_ea[:, 2] = 1
    lex = np.sqrt(problem.A0 / (0.5 * mu0 * Ms**2))

    def Hext_fct(t):
        return np.atleast_2d(t).T * [0, 0, 0]

    E_arr = np.zeros(shape=(4, len(L_loop), 2))

    for i in range(len(L_loop)):
        print(f"ITERATION: {i} / {len(L_loop)}")
        for j in range(2):
            if j == 0:
                print("Flower state")
                problem.m0[:, 0:2] = 0
                problem.m0[:, 2] = 1
                t_end = 10e-9

            elif j == 1:
                print("Vortex state")
                xv = np.linspace(-1, 1, res[0])
                yv = np.linspace(-1, 1, res[1])
                zv = np.linspace(-1, 1, res[2])
                [x, y, z] = np.meshgrid(xv, yv, zv, indexing="ij")
                xvec = np.sin(np.arctan2(z, x))
                yvec = -np.cos(np.arctan2(z, x))
                problem.m0[:, 0] = xvec.swapaxes(0, 2).reshape((-1))
                problem.m0[:, 2] = yvec.swapaxes(0, 2).reshape((-1))
                problem.m0 = problem.m0 / np.tile(
                    np.expand_dims(np.sqrt(np.sum(problem.m0**2, axis=1)), axis=1),
                    (1, 3),
                )
                t_end = 200e-9

            problem.grid_L = np.array([lex, lex, lex]) * L_loop[i]
            t, M_out, _, H_exc, H_ext, H_dem, H_ani = problem.run_simulation(
                t_end, 50, Hext_fct, 2
            )

            if show_details:
                M_sq = np.squeeze(M_out, axis=2)
                plot_M_avg_seq(t, M_sq)
                plot_M_thin_film(M_sq[0], res)
                plot_M_thin_film(M_sq[-1], res)

            # Calculate the energy terms
            E_exc = np.sum(
                (1 / 2)
                * (
                    M_out[:, :, 0, 0] * H_exc[:, :, 0, 0]
                    + M_out[:, :, 0, 1] * H_exc[:, :, 0, 1]
                    + M_out[:, :, 0, 2] * H_exc[:, :, 0, 2]
                ),
                axis=1,
            )

            E_ext = np.sum(
                M_out[:, :, 0, 0] * H_ext[:, :, 0, 0]
                + M_out[:, :, 0, 1] * H_ext[:, :, 0, 1]
                + M_out[:, :, 0, 2] * H_ext[:, :, 0, 2],
                axis=1,
            )

            E_dem = np.sum(
                (1 / 2)
                * (
                    M_out[:, :, 0, 0] * H_dem[:, :, 0, 0]
                    + M_out[:, :, 0, 1] * H_dem[:, :, 0, 1]
                    + M_out[:, :, 0, 2] * H_dem[:, :, 0, 2]
                ),
                axis=1,
            )

            E_ani = np.sum(
                (1 / 2)
                * (
                    M_out[:, :, 0, 0] * H_ani[:, :, 0, 0]
                    + M_out[:, :, 0, 1] * H_ani[:, :, 0, 1]
                    + M_out[:, :, 0, 2] * H_ani[:, :, 0, 2]
                ),
                axis=1,
            )

            E_arr[:, i, j] = mu0 * np.array(
                [E_exc[-1], E_ext[-1], E_dem[-1], E_ani[-1]]
            )

            if show_details:
                plt.clf()
                for E_x in [E_exc, E_ext, E_dem, E_ani]:
                    plt.plot(t, mu0 * E_x - mu0 * E_x[0], ".")
                plt.xlabel("Time [s]")
                plt.ylabel("Energy [-]")
                plt.legend([r"$E_{exc}$", r"$E_{ext}$", r"$E_{dem}$", r"$E_{ani}$"])
                plt.show()

    if show:
        plt.clf()
        plt.plot(L_loop, np.sum(E_arr[:, :, 0], axis=0), ".")
        plt.plot(L_loop, np.sum(E_arr[:, :, 1], axis=0), "o")
        plt.xlabel("L [l_ex]")
        plt.ylabel("Energy [-]")
        plt.show()


# %%

if __name__ == "__main__":
    std_prob_3(show=False, cuda=False, show_details=True)
