from __future__ import annotations

import time
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from magtense.micromag import MicromagProblem
from scipy.interpolate import interp1d

import sys
sys.path.append("../../src/utils")
from utils import setup_grain_problem_from_matfile, load_matlab_struct

def max_energy_product(
    H: np.ndarray,
    M: np.ndarray,
    mu0: float = 4 * np.pi * 1e-7,
) -> tuple[float, float, float, float]:
    """
    Compute maximum energy product (BH)_max from hysteresis data
    using the discrete sample points (no interpolation).

    Returns
    -------
    BH_max : float
        Maximum energy product [J/m^3].
    H_at_max : float
        Field value at (BH)_max [A/m].
    B_at_max : float
        Flux density at (BH)_max [T].
    M_at_max : float
        Magnetization at (BH)_max [A/m].
    """
    H = np.asarray(H, dtype=float).ravel()
    M = np.asarray(M, dtype=float).ravel()

    B = mu0 * (H + M)

    # Second quadrant: H < 0, B > 0
    mask = (H < 0.0) & (B > 0.0)
    if not np.any(mask):
        raise ValueError("No points in the second quadrant (H<0, B>0) found.")

    Hq = H[mask]
    Bq = B[mask]
    Mq = M[mask]

    BH = -Bq * Hq  # J/m^3

    idx = int(np.argmax(BH))

    BH_max = BH[idx]
    H_at_max = Hq[idx]
    B_at_max = Bq[idx]
    M_at_max = Mq[idx]

    return BH_max, H_at_max, B_at_max, M_at_max


def grain_hysteresis(
    use_fmm: bool = True,
    fmm_eps: float = 1e-4,
    fmm_cells_per_node: int = 50,
    cuda: bool = False,
    cvode: bool = False,
    mesh_file: str = "Grid_rasBase_5_nGrains_5_nRef_3_dG_5e-09.mat",
    plotting: bool = True,
    figpath: Path | None = Path(__file__).parent.absolute().joinpath("..", "figs"),
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
) -> None:

    """
    Run a hysteresis simulation for a given grain mesh, compute diagnostics,
    save the full result object `res` plus final scalars, and optionally plot.

    Files written
    -------------
    <mesh_stem>_cuda.npz or <mesh_stem>_fmm.npz   # contains res + final scalars
    <mesh_stem>_cuda.png or <mesh_stem>_fmm.png   # hysteresis figure (if plotting=True)
    """

    print("Loading data from ", mesh_file)

    mu0 = 4 * np.pi * 1e-7
    Hyst_dir = np.array([0.0, 0.0, 1.0])

    Ms = 1.61 / mu0   # [A/m]
    K0 = 4.3e6        # [J/m^3]


    h_ext_base = Hyst_dir / mu0
    #steps = np.arange(1.0, -7.1, -0.1)
    steps = np.arange(1.0, 0.8, -0.1)

    H_ext = np.zeros((len(steps), 4))
    H_ext[:, 0] = steps
    H_ext[:, 1:4] = steps[:, np.newaxis] * h_ext_base





    exma_file = mesh_file[:-4] + "_exMa.npz"
    demag_file = mesh_file[:-4] + "_demag.bin"
    if not Path(exma_file).exists():
        print("File {} does not exist.".format(exma_file))
        print("Generating exchange matrix and saving to file...")
        rng = np.random.default_rng(42)
        if cuda:
            problem_ini = setup_grain_problem_from_matfile(
                fname=mesh_file,
                cuda=cuda,
                cvode=cvode,
                t_end=40e-9,
                nt=2,
                Hyst_dir=Hyst_dir,
                rng=rng,
                input_file_name=demag_file,
            )
        else:
            problem_ini = setup_grain_problem_from_matfile(
                fname=mesh_file,
                cuda=cuda,
                cvode=cvode,
                t_end=40e-9,
                nt=2,
                Hyst_dir=Hyst_dir,
                rng=rng,
            )



        #problem_ini.exch_presize = 10 * problem_ini.nt * len(steps)
        nnbor = 27 # assmue 27 nearest neighbors for exchange
        problem_ini.exch_presize = problem_ini.nt * (nnbor+1)
        problem_ini.dummy_run = 1
        problem_ini.fmm_cells_per_node = fmm_cells_per_node if use_fmm else 0
        problem_ini.fmm_eps = fmm_eps
        problem_ini.ifunif = ifunif
        problem_ini.nlmin = nlmin
        problem_ini.nlmax = nlmax
        problem_ini.allow_fmm_short_circuit = allow_fmm_short_circuit
        problem_ini.fmm_min_n = fmm_min_n

        problem_ini.log_dir = "timer_logs"
        problem_ini.timer_log_file = mesh_file[:-4] + "_timer_ini.log"
        problem_ini.trace_log_file = mesh_file[:-4] + "_trace_ini.log"
        problem_ini.window_enabled = 1
        problem_ini.window_interval = 30.0
        problem_ini.trace_enabled = 0
        problem_ini.flush_each = 1
        problem_ini.trace_verbose = 1


        res_ini = problem_ini.run_hysteresis(H_ext=H_ext)
   
        exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = res_ini[7:]
        print("storing setup...")
        print("exch_nval = ", exch_nval)
        print("exch_nrow = ", exch_nrow)
        print("exch_ncols = ", exch_ncols)
        print("shape ExchMat_r", ExchMat_r.shape)
        print("shape ExchMat_c", ExchMat_c.shape)
        print("shape ExchMat_v", ExchMat_v.shape)
        print("ExchMat_v[:10] = ", ExchMat_v[:10])
        np.savez(
            exma_file,
            exch_nval=exch_nval,
            ExchMat_r=ExchMat_r,
            ExchMat_c=ExchMat_c,
            ExchMat_v=ExchMat_v,
            exch_nrow=exch_nrow,
            exch_ncols=exch_ncols
        )   
        print(f"Saved exchange matrix to {exma_file}")
    else:
        print(f"Loading exchange matrix from {exma_file}...")

    
    exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = np.load(
        exma_file, allow_pickle=True
    ).values()
        


    rng = np.random.default_rng(42)

    if cuda:
        problem = setup_grain_problem_from_matfile(
            fname=mesh_file,
            cuda=cuda,
            cvode=cvode,
            t_end=40e-9,
            nt=2,
            Hyst_dir=Hyst_dir,
            rng=rng,
            exch_val=ExchMat_v,
            exch_rows=ExchMat_r,
            exch_col=ExchMat_c,
            exch_nval=exch_nval,
            exch_nrow=exch_nrow,
            exch_ncols=exch_ncols,
            passexch=1,
            input_file_name=demag_file,
        )
    else:
        problem = setup_grain_problem_from_matfile(
            fname=mesh_file,
            cuda=cuda,
            cvode=cvode,
            t_end=40e-9,
            nt=2,
            Hyst_dir=Hyst_dir,
            rng=rng,
            exch_val=ExchMat_v,
            exch_rows=ExchMat_r,
            exch_col=ExchMat_c,
            exch_nval=exch_nval,
            exch_nrow=exch_nrow,
            exch_ncols=exch_ncols,
            passexch=0,
        )

    #-------- for testing
    #problem.dem_appr = "threshold"
    #problem.thres = 1e-10
    #-------- modify material parameters

    #problem.exch_presize = 10 * problem.nt * len(steps)
    nnbor = 27 # assmue 27 nearest neighbors for exchange
    problem.exch_presize = problem.nt * (nnbor+1)

    problem.dummy_run = 0
    problem.fmm_cells_per_node = fmm_cells_per_node if use_fmm else 0
    problem.fmm_eps = fmm_eps
    problem.ifunif = ifunif
    problem.nlmin = nlmin
    problem.nlmax = nlmax
    problem.allow_fmm_short_circuit = allow_fmm_short_circuit
    problem.fmm_min_n = fmm_min_n

    problem.log_dir = "timer_logs"
    problem.timer_log_file = mesh_file[:-4] + "_timer.log"
    problem.trace_log_file = mesh_file[:-4] + "_trace.log"
    problem.window_enabled = 1
    problem.window_interval = 30.0
    problem.trace_enabled = 0
    problem.flush_each = 1
    problem.trace_verbose = 1


    start_time = time.time()
    res = problem.run_hysteresis(H_ext=H_ext)
    end_time = time.time()
    runtime = end_time - start_time
    print(f"Hysteresis simulation took {runtime:.3f} seconds")

    M_out = res[1][1, :, :, :]  # (N_cell, N_steps, 3)

    mesh_cart, GridInfo, mesh_params, iIn = load_matlab_struct(mesh_file)

    volumes = np.asarray(GridInfo["Volumes"], dtype=float).ravel()
    vol_sum = np.sum(volumes)
    Ms_vol = problem.Ms.ravel() * volumes

    MxMean = np.zeros(len(steps))
    MyMean = np.zeros(len(steps))
    MzMean = np.zeros(len(steps))

    for i in range(len(steps)):
        MxMean[i] = np.sum(Ms_vol * M_out[:, i, 0]) / vol_sum
        MyMean[i] = np.sum(Ms_vol * M_out[:, i, 1]) / vol_sum
        MzMean[i] = np.sum(Ms_vol * M_out[:, i, 2]) / vol_sum

    H_N = 2 * K0 / Ms

    M = Hyst_dir[0] * MxMean + Hyst_dir[1] * MyMean + Hyst_dir[2] * MzMean

    H = (
        np.sign(H_ext[:, 0])
        * np.sqrt(H_ext[:, 1] ** 2 + H_ext[:, 2] ** 2 + H_ext[:, 3] ** 2)
        * mu0
    )

    try:
        f = interp1d(M, H)     # H as function of M
        Hc = float(f(0.0))
    except Exception:
        Hc = np.nan

    # try:
    #     Hc = float(np.interp(0.0, M, H))
    # except Exception:
    #     Hc = np.nan

    try:
        BH_max, H_star, B_star, M_star = max_energy_product(H, M, mu0=mu0)
    except ValueError:
        BH_max = np.nan
        H_star = np.nan
        B_star = np.nan
        M_star = np.nan

    print("BH_max = {:.2e} kJ/m^3".format(BH_max / 1e3))
    print("H*     = {:.3e} A/m".format(H_star))
    print("B*     = {:.3f} T".format(B_star))
    print("M*     = {:.3e} A/m".format(M_star))
    print("Hc     = {:.3e} A/m".format(Hc))
    print("H_N    = {:.3e} A/m".format(H_N))

    mesh_path = Path(mesh_file)
    stem = mesh_path.stem
    #backend_tag = "_cuda" if cuda else "_fmm"


    if use_fmm and allow_fmm_short_circuit == 1 and len(problem.Ms) < fmm_min_n:
        backend_tag = "_cuda_fmmsc"
    elif use_fmm:
        backend_tag = "_fmm"
        backend_tag += f"_eps{fmm_eps:.2e}_cpn{fmm_cells_per_node}"
    else:
        backend_tag = "_cuda"

    #backend_tag = "_fmm" if use_fmm else "_cuda"
    #if use_fmm: 
    #    backend_tag += f"_eps{fmm_eps:.2e}_cpn{fmm_cells_per_node}"

    root_name = stem + backend_tag

    # Use figpath (default ../figs relative to this file)
    if figpath is None:
        out_dir = Path.cwd()
    else:
        out_dir = figpath
        out_dir.mkdir(parents=True, exist_ok=True)

    res_path = out_dir / f"{root_name}.npz"
    fig_path = out_dir / f"{root_name}.png"


    res_obj = np.array(res, dtype=object)  # <-- make `res` a 0-D object array

    np.savez(
        res_path,
        res=res_obj,
        BH_max=BH_max,
        H_star=H_star,
        B_star=B_star,
        M_star=M_star,
        Hc=Hc,
        H_N=H_N,
        runtime=runtime,
        H_array=H,
        M_array=M,
    )
    print(f"Saved result object and scalars to: {res_path}")

    if plotting:
        fig, ax = plt.subplots(figsize=(8, 4))

        ax.plot(H, M / Ms, ".-k", linewidth=1.5, markersize=4, label="M/Ms")
        ax.plot(Hc, 0.0, "r*", markersize=10, label="Hc")
        ax.plot(H_star, M_star / Ms, "b*", markersize=10, label=r"$BH_{\max}$")

        ax.set_xlabel("H [A/m]")
        ax.set_ylabel("M/Ms [-]")
        ax.set_title("Hysteresis Loop")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.legend()

        fig.tight_layout()
        fig.savefig(fig_path, dpi=300)
        print(f"Saved figure to: {fig_path}")

if __name__ == "__main__":
    default_figpath = Path(__file__).parent.absolute().joinpath("..", "figs_perf")

    parser = argparse.ArgumentParser(
        description="Run grain hysteresis experiment"
    )

    parser.add_argument("--use-fmm", action="store_true",
                        help="Enable FMM demagnetisation")

    parser.add_argument("--fmm-eps", type=float, default=1e-4,
                        help="FMM accuracy parameter")

    parser.add_argument("--fmm-cpn", type=int, default=660,
                        help="Cells per FMM tree node")

    parser.add_argument("--cuda", action="store_true",
                        help="Enable CUDA acceleration")

    parser.add_argument("--cvode", action="store_true",
                        help="Enable CVODE time integration")

    parser.add_argument("--mesh-file", type=str, required=True,
                        help="Input mesh .mat file")

    parser.add_argument("--do-plot", action="store_true",
                        help="Enable plotting")

    parser.add_argument("--figpath", type=Path, default=default_figpath,
                        help="Output directory for figures")

    parser.add_argument("--ifunif", type=int, default=1,
                        help="Uniform FMM tree if 1, else adaptive")

    parser.add_argument("--nlmin", type=int, default=1,
                        help="Minimum level for adaptive FMM tree")

    parser.add_argument("--nlmax", type=int, default=5,
                        help="Maximum level for adaptive FMM tree")

    parser.add_argument("--allow-fmm-short-circuit", type=int, default=0,
                        help="Allow FMM short circuit (0/1)")

    parser.add_argument("--fmm-min-n", type=int, default=20000,
                        help="Minimum number of elements for FMM short circuit")


    args = parser.parse_args()
    

    print("Running grain hysteresis experiment with parameters:")
    print(f"  use_fmm = {args.use_fmm}")
    print(f"  fmm_eps = {args.fmm_eps}")
    print(f"  fmm_cpn = {args.fmm_cpn}")
    print(f"  cuda = {args.cuda}")
    print(f"  cvode = {args.cvode}")
    print(f"  mesh_file = {args.mesh_file}")
    print(f"  do_plot = {args.do_plot}")
    print(f"  figpath = {args.figpath}")
    print(f"  ifunif = {args.ifunif}")
    print(f"  nlmin = {args.nlmin}")
    print(f"  nlmax = {args.nlmax}")
    print(f"  allow_fmm_short_circuit = {args.allow_fmm_short_circuit}")
    print(f"  fmm_min_n = {args.fmm_min_n}")

    grain_hysteresis(
        use_fmm=args.use_fmm,
        fmm_eps=args.fmm_eps,
        fmm_cells_per_node=args.fmm_cpn,
        cuda=args.cuda,
        cvode=args.cvode,
        mesh_file=args.mesh_file,
        plotting=args.do_plot,
        figpath=args.figpath,
        ifunif=args.ifunif,
        nlmin=args.nlmin,
        nlmax=args.nlmax,
        allow_fmm_short_circuit=args.allow_fmm_short_circuit,
        fmm_min_n=args.fmm_min_n,
    )


