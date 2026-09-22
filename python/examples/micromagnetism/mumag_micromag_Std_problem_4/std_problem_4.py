"""muMag standard problem 4: the switching dynamics of a thin film under a reversed field.

It is a two-stage calculation: first an s-state is relaxed in a decaying field, then that state
is the initial condition for the dynamic run in one of the two fields of the problem, whose mean
magnetization is compared with the published mean solutions. This is the Python counterpart of
matlab/examples/Micromagnetism/mumag_micromag_Std_problem_4/Standard_problem_4.m and takes the
same options.

The film can be run on two kinds of mesh, selected with mesh_type:

  'uniform'            a regular grid of res[0] x res[1] x res[2] cells. This is the default and
                       the only mesh type that uses the res argument.
  'unstructuredPrisms' a grid of unstructured Cartesian prisms read from the text file mesh_file
                       (one row per prism: centre [x,y,z] and side lengths [a,b,c] in metres).
                       The exchange operator that MagTense builds from the mesh in the first stage
                       is handed to the second stage, so the mesh is analysed only once.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from magtense.micromag import MicromagProblem
from magtense.utils import plot_M_thin_film

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
MUMAG_DIR = (
    REPOSITORY_ROOT / "documentation" / "examples_mumag_validation" / "Validation_standard_problem_4"
)
DEFAULT_MESH_FILE = MUMAG_DIR / "Std_prob_4_unstructured_mesh_grains_6_res_80_20_ref_2.txt"
# MUMAG_DIR / "Std_prob_4_unstructured_mesh_grains_6_res_100_25_ref_3.txt" is the finer mesh
MESH_TYPES = ("uniform", "unstructuredPrisms")

# Entries of the exchange matrix that MagTense reserves per cell for returning it. The stencil of
# the unstructured operator reaches every cell sharing a vertex with a face, so it is far wider
# than the seven point stencil of the uniform grid.
EXCH_PRESIZE = {"uniform": 12, "unstructuredPrisms": 64}


def std_prob_4(
    mumag_field: int = 1,
    res: tuple[int, int, int] = (36, 9, 1),
    mesh_type: str = "uniform",
    mesh_file: Path | str = DEFAULT_MESH_FILE,
    cuda: bool = True,
    cvode: bool = False,
    use_avgn: bool = True,
    plotting: bool = True,
    figpath: Path | None = None,
) -> tuple[list[float], list[float]]:
    """Run the muMag standard problem 4 and compare with the published mean solutions.

    Args:
        mumag_field: 1 or 2, which of the two applied fields of the problem description to use.
        res: cells along x, y and z of the uniform grid. Ignored for the unstructured mesh.
        mesh_type: 'uniform' or 'unstructuredPrisms', see the module docstring.
        mesh_file: the text file with the unstructured Cartesian mesh.
        cuda: use CUDA for the calculations.
        cvode: use CVODE for the numerical time evolution.
        use_avgn: use the averaged prism tensor for the demag field.
        plotting: show (or save, with figpath) the starting state and the comparison with mumag.
        figpath: directory to save the figures in. None shows them interactively.

    Returns:
        int_error: the integral of |M_MagTense - M_mumag| over the simulated second, for the
            x, y and z component. Has units of seconds.
        rel_int_error: the same integral divided by the integral of |M_mumag| and expressed in
            percent. This is the same measure as calculate_relative_integral_error.m uses in the
            MATLAB test suite, so the two can be compared directly.
    """
    if mesh_type not in MESH_TYPES:
        raise ValueError(f"mesh_type must be one of {MESH_TYPES}, got {mesh_type!r}")
    if mumag_field not in (1, 2):
        raise ValueError(f"mumag_field must be 1 or 2, got {mumag_field!r}")

    mu0 = 4 * np.pi * 1e-7
    grid_L = [500e-9, 125e-9, 3e-9]

    # Setup the mesh
    if mesh_type == "unstructuredPrisms":
        mesh_data = np.loadtxt(mesh_file, comments="%")
        grid_pts = mesh_data[:, 0:3]
        grid_abc = mesh_data[:, 3:6]
        res = (len(grid_pts), 1, 1)
        volumes = np.prod(grid_abc, axis=1)
        print(f"Prisms N_grid = {len(grid_pts)}")
    else:
        grid_pts = None
        grid_abc = None
        volumes = np.ones(int(np.prod(res)))

    ### Magnetization to s-state
    problem_ini = MicromagProblem(
        res=res,
        grid_L=grid_L,
        grid_type=mesh_type,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        exch_presize=EXCH_PRESIZE[mesh_type],
        m0=1 / np.sqrt(3),
        alpha=4.42e3,
        gamma=0,
        Ms=8e5,
        K0=0,
        A0=1.3e-11,
        cuda=cuda,
        cvode=cvode,
        useavgn=use_avgn,
    )
    h_ext = np.array([1, 1, 1]) / mu0

    def h_ext_fct_init(t) -> np.ndarray:
        return np.expand_dims(np.where(t < 1e-09, 1e-09 - t, 0), axis=1) * h_ext

    problem_ini.setTimeDis = 100
    problem_ini.timer_log_file = "std_4_ini_timer.log"
    problem_ini.trace_log_file = "std_4_ini_trace.log"

    result = problem_ini.run_simulation(
        t_end=100e-9,
        nt=200,
        fct_h_ext=h_ext_fct_init,
        nt_h_ext=2000,
    )
    M_sq_ini = np.squeeze(result[1], axis=2)
    pts = result[2]

    ### Time-dependent solver
    # The exchange operator MagTense built from the unstructured mesh in the first stage is
    # passed on, as the Matlab example does with setExchangeMatrixCOO, so the mesh is only
    # analysed once. On the uniform grid the operator is trivially rebuilt.
    exch_kwargs = {}
    if mesh_type == "unstructuredPrisms":
        exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = result[7:]
        if exch_nval > EXCH_PRESIZE[mesh_type] * np.prod(res):
            raise RuntimeError(
                f"the exchange matrix has {exch_nval} entries, more than the "
                f"{EXCH_PRESIZE[mesh_type] * np.prod(res)} reserved for returning it; raise EXCH_PRESIZE"
            )
        exch_kwargs = {
            "exch_nval": exch_nval, "exch_nrow": exch_nrow, "exch_ncols": exch_ncols,
            "exch_rows": ExchMat_r, "exch_cols": ExchMat_c, "exch_val": ExchMat_v, "passexch": 1,
        }

    problem_dym = MicromagProblem(
        res=res,
        grid_L=grid_L,
        grid_type=mesh_type,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        m0=M_sq_ini[-1],
        alpha=4.42e3,
        gamma=2.21e5,
        Ms=8e5,
        K0=0,
        A0=1.3e-11,
        cuda=cuda,
        cvode=cvode,
        useavgn=use_avgn,
        **exch_kwargs,
    )

    # Two applied external fields of std problem 4
    if mumag_field == 1:
        h_ext_mumag = np.array([-24.6, 4.3, 0])
    else:
        h_ext_mumag = np.array([-35.5, -6.3, 0])

    def h_ext_fct(t) -> np.ndarray:
        return np.expand_dims(t > -1, axis=1) * (h_ext_mumag / 1000 / mu0)

    problem_dym.setTimeDis = 10
    problem_dym.timer_log_file = "std_4_dym_timer.log"
    problem_dym.trace_log_file = "std_4_dym_trace.log"

    t_dym, M_out = problem_dym.run_simulation(
        t_end=1e-9,
        nt=200,
        fct_h_ext=h_ext_fct,
        nt_h_ext=2000,
    )[:2]

    # Copy as otherwise later a NumPy view is created into the Fortran-owned memory. The mean
    # magnetization is weighted by the cell volumes, which matters for the unstructured mesh
    M_sq_dym = np.squeeze(M_out.copy(), axis=2)
    weights = volumes / volumes.sum()
    Mx, My, Mz = np.einsum("tnc,n->tc", M_sq_dym, weights).T

    ## Compare with published solutions available from mumag webpage
    M_mumag = np.loadtxt(MUMAG_DIR / f"Field_{mumag_field}_mumag_mean_solution.txt", skiprows=1)

    # Interpolate the MagTense solution to the mumag-published solutions and
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
    # Normalised the way the MATLAB test suite does it, so the acceptance limits carry over
    reference = [
        np.trapezoid(np.abs(M_mumag[:, 0]), t),
        np.trapezoid(np.abs(M_mumag[:, 2]), t),
        np.trapezoid(np.abs(M_mumag[:, 4]), t),
    ]
    rel_int_error = [
        float(err / ref * 100) for err, ref in zip(int_error, reference, strict=True)
    ]

    if plotting:
        # The starting state of the dynamical simulation
        plt.figure()
        M_end = M_sq_ini[-1]
        plt.quiver(pts[:, 0], pts[:, 1], M_end[:, 0], M_end[:, 1], pivot="mid")
        plt.axis("equal")
        plt.title("Starting state of dynamical simulation")
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / f"4_start_state_{mesh_type}.png")

        _, ax1 = plt.subplots()

        ax1.plot(t_dym, Mx, "rx")
        ax1.plot(t_dym, My, "gx")
        ax1.plot(t_dym, Mz, "bx")

        ax1.plot(t, M_mumag[:, 0], "r-")
        ax1.plot(t, M_mumag[:, 2], "g-")
        ax1.plot(t, M_mumag[:, 4], "b-")

        legend_elements = [
            Line2D([0], [0], marker="x", color="r", label=r"MagTense $M_x$", linestyle="None"),
            Line2D([0], [0], marker="x", color="g", label=r"MagTense $M_y$", linestyle="None"),
            Line2D([0], [0], marker="x", color="b", label=r"MagTense $M_z$", linestyle="None"),
            Line2D([0], [0], marker="none", color="r", label=r"$\mu{}mag$ $<M_x>$"),
            Line2D([0], [0], marker="none", color="g", label=r"$\mu{}mag$ $<M_y>$"),
            Line2D([0], [0], marker="none", color="b", label=r"$\mu{}mag$ $<M_z>$"),
        ]
        ax1.legend(handles=legend_elements)
        plt.setp(plt.gca().get_legend().get_texts(), fontsize="14")
        plt.xlabel("Time [s]", fontsize="14")
        plt.ylabel(r"$M_i$" + " [-]", fontsize="14")
        plt.title(f"Standard problem 4, Field {mumag_field}, mesh: {mesh_type}")
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / f"4_field_{mumag_field}_mt_vs_mumag.png")

        if mesh_type == "uniform":
            plot_M_thin_film(M_sq_dym[0], res, "Start_state", figpath=figpath)
            plot_M_thin_film(M_sq_dym[-1], res, "Final_state", figpath=figpath)

    print("int_error: ", int_error)
    print("rel_int_error [%]: ", rel_int_error)
    return int_error, rel_int_error


if __name__ == "__main__":
    int_error, rel_int_error = std_prob_4(
        mumag_field=1,
        mesh_type="uniform",   # 'uniform' or 'unstructuredPrisms'
        cuda=True,
        cvode=False,
        use_avgn=True,
        plotting=True,
        figpath=None,
    )
