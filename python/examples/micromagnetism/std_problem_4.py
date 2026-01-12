from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from magtense.micromag import MicromagProblem
from magtense.utils import plot_M_thin_film


def std_prob_4(
    NIST_field: int = 1,
    res: tuple[int, int, int] = (36, 9, 1),
    cuda: bool = False,
    cvode: bool = False,
    unstructured: bool = False,
    mesh_file: str = "unstructured_grains_6_res_80_20_ref_2",
    plotting: bool = True,
    figpath: Path | None = None,
) -> list[float]:
    mu0 = 4 * np.pi * 1e-7
    grid_L = [500e-9, 125e-9, 3e-9]

    if unstructured:
        grid_type = "unstructuredPrisms"
        # Load unstructured mesh points and grain data
        mesh_data = np.loadtxt(
            Path(__file__).parent.absolute() / ".." / "meshes" / f"{mesh_file}.txt"
        )
        grid_pts = np.array(mesh_data[:, :3])
        grid_abc = np.array(mesh_data[:, 3:])
        res = (len(grid_pts), 1, 1)
    else:
        grid_type = "uniform"
        grid_pts = None
        grid_abc = None

    ### Magnetization to s-state
    problem_ini = MicromagProblem(
        res=res,
        grid_L=grid_L,
        m0=1 / np.sqrt(3),
        alpha=4.42e3,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        grid_type=grid_type,
        cuda=cuda,
        cvode=cvode,
    )
    h_ext = np.array([1, 1, 1]) / mu0

    def h_ext_fct_init(t) -> np.ndarray:
        return np.expand_dims(np.where(t < 1e-09, 1e-09 - t, 0), axis=1) * h_ext

    result = problem_ini.run_simulation(
        t_end=100e-9,
        nt=200,
        fct_h_ext=h_ext_fct_init,
        nt_h_ext=2000,
    )
    M_sq_ini = np.squeeze(result[1], axis=2)
    if unstructured:
        exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = result[7:]
        passexch = 1
    else:
        exch_nval = 1
        exch_nrow = 1
        exch_ncols = 1
        ExchMat_r = None
        ExchMat_c = None
        ExchMat_v = None
        passexch = 0

    ### Time-dependent solver
    problem_dym = MicromagProblem(
        res=res,
        grid_L=grid_L,
        m0=M_sq_ini[-1],
        alpha=4.42e3,
        gamma=2.21e5,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        grid_type=grid_type,
        exch_rows=ExchMat_r,
        exch_col=ExchMat_c,
        exch_val=ExchMat_v,
        exch_nval=exch_nval,
        exch_nrow=exch_nrow,
        exch_ncols=exch_ncols,
        passexch=passexch,
        cuda=cuda,
        cvode=cvode,
    )

    # Two applied external fields of std problem 4
    if NIST_field == 1:
        h_ext_nist = np.array([-24.6, 4.3, 0])
    elif NIST_field == 2:
        h_ext_nist = np.array([-35.5, -6.3, 0])
    else:
        raise NotImplementedError

    def h_ext_fct(t) -> np.ndarray:
        return np.expand_dims(t > -1, axis=1) * (h_ext_nist / 1000 / mu0)

    t_dym, M_out = problem_dym.run_simulation(
        t_end=1e-9,
        nt=200,
        fct_h_ext=h_ext_fct,
        nt_h_ext=2000,
    )[:2]

    # Copy as otherwise later a NumPy view is created into the Fortran-owned memory
    M_sq_dym = np.squeeze(M_out.copy(), axis=2)
    Mx = np.mean(M_sq_dym[:, :, 0], axis=1)
    My = np.mean(M_sq_dym[:, :, 1], axis=1)
    Mz = np.mean(M_sq_dym[:, :, 2], axis=1)

    ## Compare with published solutions available from mumag webpage
    mumag_eval_path = (
        Path(__file__).parent.absolute()
        / ".."
        / ".."
        / ".."
        / "documentation"
        / "examples_mumag_validation"
        / "Validation_standard_problem_4"
    )

    with Path.open(
        Path(mumag_eval_path, f"Field_{NIST_field}_mumag_mean_solution.txt"), "r"
    ) as file:
        T = file.readlines()[1:]

    M_mumag = np.asarray([line.split() for line in T], dtype=np.float64)

    # Interpolate the MagTense solution to the NIST-published solutions and
    # calculate the difference between the results as an integral.
    t = np.linspace(0, 1e-9, 1000)
    Magtense_Mx_interpolated = np.interp(t, t_dym, Mx)
    Magtense_My_interpolated = np.interp(t, t_dym, My)
    Magtense_Mz_interpolated = np.interp(t, t_dym, Mz)
    int_error = [
        np.trapezoid(np.abs(M_mumag[:, 0] - Magtense_Mx_interpolated), t),
        np.trapezoid(np.abs(M_mumag[:, 2] - Magtense_My_interpolated), t),
        np.trapezoid(np.abs(M_mumag[:, 4] - Magtense_Mz_interpolated), t),
    ]

    if plotting:
        _, ax1 = plt.subplots()

        ax1.plot(t_dym, Mx, "rx")
        ax1.plot(t_dym, My, "gx")
        ax1.plot(t_dym, Mz, "bx")

        ax1.plot(t, M_mumag[:, 0], "r-")
        ax1.plot(t, M_mumag[:, 2], "g-")
        ax1.plot(t, M_mumag[:, 4], "b-")

        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="x",
                color="r",
                label=r"MagTense $M_x$",
                linestyle="None",
            ),
            Line2D(
                [0],
                [0],
                marker="x",
                color="g",
                label=r"MagTense $M_y$",
                linestyle="None",
            ),
            Line2D(
                [0],
                [0],
                marker="x",
                color="b",
                label=r"MagTense $M_z$",
                linestyle="None",
            ),
            Line2D([0], [0], marker="none", color="r", label=r"$\mu{}mag$ $<M_x>$"),
            Line2D([0], [0], marker="none", color="g", label=r"$\mu{}mag$ $<M_y>$"),
            Line2D([0], [0], marker="none", color="b", label=r"$\mu{}mag$ $<M_z>$"),
        ]
        ax1.legend(handles=legend_elements)
        plt.setp(plt.gca().get_legend().get_texts(), fontsize="14")
        plt.xlabel("Time [s]", fontsize="14")
        plt.ylabel(r"$M_i$" + " [-]", fontsize="14")
        plt.title(f"Standard problem 4, Field {NIST_field}")
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / f"4_field_{NIST_field}_mt_vs_mumag.png")

        if not unstructured:
            plot_M_thin_film(M_sq_dym[0], res, "Start_state", figpath=figpath)
            plot_M_thin_film(M_sq_dym[-1], res, "Final_state", figpath=figpath)

    return int_error


if __name__ == "__main__":
    int_error = std_prob_4(
        NIST_field=1,
        cuda=True,
        cvode=False,
        unstructured=True,
        plotting=True,
        #figpath=Path(__file__).parent.absolute().joinpath("..", "figs"),
        figpath=None,
    )
