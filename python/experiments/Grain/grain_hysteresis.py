from __future__ import annotations

import time
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
    cuda: bool = False,
    cvode: bool = False,
    mesh_file: str = "Grid_rasBase_5_nGrains_5_nRef_3_dG_5e-09.mat",
    plotting: bool = True,
    figpath: Path | None = Path(__file__).parent.absolute().joinpath("..", "figs"),
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

    rng = np.random.default_rng(42)
    problem = setup_grain_problem_from_matfile(
        fname=mesh_file,
        cuda=cuda,
        cvode=cvode,
        t_end=40e-9,
        nt=2,
        Hyst_dir=Hyst_dir,
        rng=rng,
    )

    h_ext_base = Hyst_dir / mu0
    steps = np.arange(1.0, -7.1, -0.1)

    H_ext = np.zeros((len(steps), 4))
    H_ext[:, 0] = steps
    H_ext[:, 1:4] = steps[:, np.newaxis] * h_ext_base

    problem.exch_presize = 3 * problem.nt * len(steps)

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

    f = interp1d(M, H)     # H as function of M
    Hc = float(f(0.0))

    # try:
    #     Hc = float(np.interp(0.0, M, H))
    # except Exception:
    #     Hc = np.nan

    BH_max, H_star, B_star, M_star = max_energy_product(H, M, mu0=mu0)

    print("BH_max = {:.2e} kJ/m^3".format(BH_max / 1e3))
    print("H*     = {:.3e} A/m".format(H_star))
    print("B*     = {:.3f} T".format(B_star))
    print("M*     = {:.3e} A/m".format(M_star))
    print("Hc     = {:.3e} A/m".format(Hc))
    print("H_N    = {:.3e} A/m".format(H_N))

    mesh_path = Path(mesh_file)
    stem = mesh_path.stem
    backend_tag = "_cuda" if cuda else "_fmm"
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
    default_figpath = Path(__file__).parent.absolute().joinpath("..", "figs")

    grain_hysteresis(
        cuda=False,
        cvode=False,
        #mesh_file="Grid_rasBase_5_nGrains_5_nRef_3_dG_5e-09.mat",
        mesh_file="Grid_rasBase_5_nGrains_5_nRef_4_dG_3.75e-09.mat",
        plotting=True,
        figpath=default_figpath,
    )
