import os
import shutil
from collections.abc import Callable
from pathlib import Path

import numpy as np

# Windows only
if hasattr(os, "add_dll_directory"):
    mkl_path = Path(__file__).parent / ".." / ".." / ".." / "Library" / "bin"
    if Path.is_dir(mkl_path):
        os.add_dll_directory(mkl_path)

    nvidia_path = Path(__file__).parent / ".." / "nvidia"
    for lib in ["cublas", "cuda_runtime", "cusparse", "nvjitlink"]:
        if Path.is_dir(nvidia_path / lib / "bin"):
            os.add_dll_directory(nvidia_path / lib / "bin")

from magtense.lib import magtensesource


class MicromagProblem:
    """
    Micromagnetic problem using the Fortran implementation of MagTense.

    Args:
        res: Resolution of grid.
        grid_L: Spatial extensions.
        grid_nnod:
        grid_type: Curently only 'uniform' is supported.
        prob_mode:
        solver:
        A0: Anisotropy constant.
        Ms: Saturation magnetization [A/m].
        K0: Exchange constant.
        alpha_mm: Dampening constant.
        gamma: Gyromagnetic factor.
        max_T0:
        nt_conv:
        conv_tol: The convergence tolerence, which is the maximum change in
                  magnetization between two timesteps.
        tol: Relative tolerance for the Fortran ODE solver.
        thres: Fortran ODE solver, when a solution component Y(L) is less in
               magnitude than thres_value its set to zero.
        setTimeDis:
        dem_thres: Initial value of the demag threshold is zero, i.e. it is not used.
        demag_approx:
        exch_nval: Number of values in the exchange operator matrix.
        exch_nrow: Number of rows in the exchange operator matrix.
        filename:
        cuda: Optional GPU support via CUDA.
        cvode:
        precision: Precision for the demag tensor. Only SP is supported.
        n_threads: Number of threads used by OpenMP for building the demag tensor.
        N_ave:
        t_alpha:
        alpha_fct:
    """

    def __init__(
        self,
        res: list[int],
        grid_L: list[float] = (500e-9, 125e-9, 3e-9),
        grid_nnod: int = 0,
        grid_type: str | None = "uniform",
        prob_mode: str | None = "new",
        solver: str | None = "dynamic",
        m0: int | float | list | np.ndarray | None = None,
        A0: int | float | list | np.ndarray | None = None,
        Ms: int | float | list | np.ndarray | None = None,
        K0: int | float | list | np.ndarray | None = None,
        K1: int | float | list | np.ndarray | None = None,
        K2: int | float | list | np.ndarray | None = None,
        K0_arr: np.ndarray | None = None,
        CrysAxis: np.ndarray | None = None,
        alpha: float = 0.02,
        gamma: float = 0.0,
        max_T0: float = 2.0,
        nt_conv: int = 1,
        conv_tol: float = 1e-4,
        tol: float = 1e-4,
        thres: float = 1e-6,
        setTimeDis: int = 10,
        dem_thres: float = 0.0,
        demag_approx: str | None = None,
        cv: float = 0.0,
        grid_pts: list | np.ndarray | None = None,
        grid_abc: list | np.ndarray | None = None,
        exch_val: list | np.ndarray | None = None,
        exch_rows: list | np.ndarray | None = None,
        exch_col: list | np.ndarray | None = None,
        exch_nval: int = 1,
        exch_nrow: int = 1,
        exch_ncols: int = 1,
        exch_intpn: str | None = "extended",
        exch_meth: str | None = "directlaplacianneumann",
        exch_weigh: float = 8,
        exch_presize: int = 12,
        demigstp: int = 0,
        passexch: int = 0,
        usereturnhall: int = 0,
        filename: str = "t",
        cuda: bool = False,
        cvode: bool = False,
        precision: bool = False,
        n_threads: int = 1,
        N_ave: tuple[int] = (1, 1, 1),
        t_alpha: np.ndarray = np.zeros(1),  # noqa: B008
        alpha_fct=lambda t: np.atleast_2d(t).T * 0,
    ) -> None:
        ntot = np.prod(res)
        self.ntot = ntot
        self.grid_nnod = grid_nnod
        self.nt_conv = nt_conv
        self.exch_nval = exch_nval
        self.exch_nrow = exch_nrow
        self.exch_ncols = exch_ncols
        self.exch_intpn = exch_intpn
        self.exch_meth = exch_meth
        self.exch_weigh = exch_weigh
        self.passexch = passexch
        self.demigstp = demigstp
        self.usereturnhall = usereturnhall
        self.exch_presize = exch_presize

        self.grid_n = np.array(res, dtype=np.int32, order="F")
        self.grid_L = np.array(grid_L, dtype=np.float64, order="F")
        self.grid_pts = grid_pts
        self.grid_ele = np.zeros(shape=(4, ntot), dtype=np.float64, order="F")
        self.grid_nod = np.zeros(shape=(grid_nnod, 3), dtype=np.float64, order="F")
        self.grid_abc = grid_abc
        self.u_ea = np.zeros(shape=(ntot, 3), dtype=np.float64, order="F")

        self.grid_type = grid_type
        self.prob_mode = prob_mode
        self.solver = solver

        self.m0 = m0
        self.A0 = A0
        self.Ms = Ms
        self.K0 = K0
        self.K1 = K1
        self.K2 = K2
        self.K0_arr = K0_arr

        # --- Set the local crystal coordinates to the three Cartesian axis
        self.CrysAxis = CrysAxis

        self.alpha_mm = alpha
        self.gamma = gamma
        self.max_T0 = max_T0

        self.t_conv = np.zeros(shape=(nt_conv), dtype=np.float64, order="F")
        self.conv_tol = np.array(
            np.repeat(conv_tol, nt_conv), dtype=np.float64, order="F"
        )
        self.tol = tol
        self.thres = thres
        self.setTimeDis = setTimeDis

        self.dem_thres = dem_thres
        self.dem_appr = demag_approx
        self.cv = cv

        self.nt_alpha = len(t_alpha)
        self.alphat = np.zeros(shape=(self.nt_alpha, 2), dtype=np.float64, order="F")
        self.alphat[:, 0] = t_alpha
        self.alphat[:, 1] = alpha_fct(t_alpha)

        self.exch_val = exch_val
        self.exch_rows = exch_rows
        self.exch_rowe = np.zeros(shape=(self.exch_nrow), dtype=np.int32, order="F")
        self.exch_col = exch_col

        self.N_load = len(filename)
        self.N_file_in = filename
        self.N_ret = len(filename)
        self.N_file_out = filename

        if cuda and not shutil.which("nvidia-smi"):
            print("[WARNING] No GPU available! Falling back to `cuda=false`!")
            self.cuda = 0
        else:
            self.cuda = int(cuda)

        self.cvode = int(cvode)
        self.precision = int(precision)
        self.n_threads = n_threads
        self.N_ave = np.array(N_ave, dtype=np.int32, order="F")

    @property
    def passexch(self) -> int | None:
        return self._passexch

    @passexch.setter
    def passexch(self, val: int) -> None:
        self._passexch = val

    @property
    def exch_nval(self) -> int | None:
        return self._exch_nval

    @exch_nval.setter
    def exch_nval(self, val: int) -> None:
        self._exch_nval = val

    @property
    def exch_ncols(self) -> int | None:
        return self._exch_ncols

    @exch_ncols.setter
    def exch_ncols(self, val: int) -> None:
        self._exch_ncols = val

    @property
    def exch_nrow(self) -> int | None:
        return self._exch_nrow

    @exch_nrow.setter
    def exch_nrow(self, val: int) -> None:
        self._exch_nrow = val

    @property
    def exch_val(self) -> list | np.ndarray | None:
        return self._exch_val

    @exch_val.setter
    def exch_val(self, val: int | None) -> None:
        if val is None:
            self._exch_val = np.zeros(
                shape=(self.exch_nval,), dtype=np.float64, order="F"
            )
        else:
            assert np.asarray(val).shape == (self.exch_nval,)
            self._exch_val = np.asarray(val, dtype=np.float64, order="F")

    @property
    def exch_rows(self) -> list | np.ndarray | None:
        return self._exch_rows

    @exch_rows.setter
    def exch_rows(self, val: list | np.ndarray | None) -> None:
        if val is None:
            self._exch_rows = np.zeros(
                shape=(self.exch_nval,), dtype=np.int32, order="F"
            )
        else:
            assert np.asarray(val).shape == (self.exch_nval,)
            self._exch_rows = np.asarray(val, dtype=np.int32, order="F")

    @property
    def exch_col(self) -> list | np.ndarray | None:
        return self._exch_col

    @exch_col.setter
    def exch_col(self, val: list | np.ndarray | None) -> None:
        if val is None:
            self._exch_col = np.zeros(
                shape=(self.exch_nval,), dtype=np.int32, order="F"
            )
        else:
            assert np.asarray(val).shape == (self.exch_nval,)
            self._exch_col = np.asarray(val, dtype=np.int32, order="F")

    @property
    def grid_pts(self) -> list | np.ndarray | None:
        return self._grid_pts

    @grid_pts.setter
    def grid_pts(self, val: list | np.ndarray | None) -> None:
        if val is None:
            self._grid_pts = np.zeros(shape=(self.ntot, 3), dtype=np.float64, order="F")
        else:
            assert np.asarray(val).shape == (self.ntot, 3)
            self._grid_pts = np.asarray(val, dtype=np.float64, order="F")

    @property
    def grid_abc(self) -> list | np.ndarray | None:
        return self._grid_abc

    @grid_abc.setter
    def grid_abc(self, val: list | np.ndarray | None) -> None:
        if val is None:
            self._grid_abc = np.zeros(shape=(self.ntot, 3), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 3)
            self._grid_abc = np.asarray(val, dtype=np.float64, order="F")

    @property
    def A0(self) -> int | float | list | np.ndarray | None:
        return self._A0

    @A0.setter
    def A0(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._A0 = 1.3e-11 + np.zeros(
                shape=(self.ntot, 1), dtype=np.float64, order="F"
            )

        elif isinstance(val, (int, float)):
            self._A0 = val + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 1)
            self._A0 = np.asarray(val, dtype=np.float64, order="F")

    @property
    def Ms(self) -> int | float | list | np.ndarray | None:
        return self._Ms

    @Ms.setter
    def Ms(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._Ms = 8e5 + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        elif isinstance(val, (int, float)):
            self._Ms = val + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 1)
            self._Ms = np.asarray(val, dtype=np.float64, order="F")

    @property
    def K0(self) -> int | float | list | np.ndarray | None:
        return self._K0

    @K0.setter
    def K0(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._K0 = 0.0 + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        elif isinstance(val, (int, float)):
            self._K0 = val + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 1)
            self._K0 = np.asarray(val, dtype=np.float64, order="F")

    @property
    def K1(self) -> int | float | list | np.ndarray | None:
        return self._K1

    @K1.setter
    def K1(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._K1 = 0.0 + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        elif isinstance(val, (int, float)):
            self._K1 = val + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 1)
            self._K1 = np.asarray(val, dtype=np.float64, order="F")

    @property
    def K2(self) -> int | float | list | np.ndarray | None:
        return self._K2

    @K2.setter
    def K2(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._K2 = 0.0 + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        elif isinstance(val, (int, float)):
            self._K2 = val + np.zeros(shape=(self.ntot, 1), dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == (self.ntot, 1)
            self._K2 = np.asarray(val, dtype=np.float64, order="F")

    @property
    def K0_arr(self) -> int | float | list | np.ndarray | None:
        return self._K0_arr

    @K0_arr.setter
    def K0_arr(self, val: np.ndarray | None) -> None:
        if val is None:
            self._K0_arr = 0.0 + np.zeros(
                shape=(self.ntot, 6, 3), dtype=np.float64, order="F"
            )

        else:
            assert np.asarray(val).shape == (self.ntot, 6, 3)
            self._K0_arr = np.asarray(val, dtype=np.float64, order="F")

    @property
    def CrysAxis(self) -> np.ndarray | None:
        return self._CrysAxis

    @CrysAxis.setter
    def CrysAxis(self, val: np.ndarray | None) -> None:
        if val is None:
            self._CrysAxis = 0.0 + np.zeros(
                shape=(self.ntot, 3, 3), dtype=np.float64, order="F"
            )
            self._CrysAxis[:, 0, 0] = 1
            self._CrysAxis[:, 1, 1] = 1
            self._CrysAxis[:, 2, 2] = 1
        else:
            assert np.asarray(val).shape == (self.ntot, 3, 3)
            self._CrysAxis = np.asarray(val, dtype=np.float64, order="F")

    @property
    def m0(self) -> int | float | list | np.ndarray | None:
        return self._m0

    @m0.setter
    def m0(self, val: int | float | list | np.ndarray | None, seed: int = 0) -> None:
        self._m0 = np.zeros(shape=(self.ntot, 3), dtype=np.float64, order="F")

        if val is None:
            rng = np.random.default_rng(seed)
            theta = np.pi * rng.random(self.ntot)
            phi = 2 * np.pi * rng.random(self.ntot)
            self._m0[:, 0] = np.sin(theta) * np.cos(phi)
            self._m0[:, 1] = np.sin(theta) * np.sin(phi)
            self._m0[:, 2] = np.cos(theta)

        elif isinstance(val, (int, float)):
            self._m0[:] = val

        else:
            assert np.asarray(val).shape == (self.ntot, 3)
            self._m0 = np.asarray(val, dtype=np.float64, order="F")

    @property
    def dem_appr(self) -> int:
        return self._dem_appr

    @dem_appr.setter
    def dem_appr(self, val: str | None = None) -> None:
        self._dem_appr = {
            None: 1,
            "threshold": 2,
            "fft_thres": 3,
            "threshold_fraction": 4,
            "fft_threshold_fraction": 5,
        }[val]

    @property
    def grid_type(self) -> int:
        return self._grid_type

    @grid_type.setter
    def grid_type(self, val: str | None = None) -> None:
        self._grid_type = {
            None: -1,
            "uniform": 1,
            "tetrahedron": 2,
            "unstructuredPrisms": 3,
        }[val]

    @property
    def exch_intpn(self) -> int:
        return self._exch_intpn

    @exch_intpn.setter
    def exch_intpn(self, val: str | None = None) -> None:
        self._exch_intpn = {
            None: -1,
            "extended": 1,
            "compact": 2,
        }[val]

    @property
    def exch_meth(self) -> int:
        return self._exch_meth

    @exch_meth.setter
    def exch_meth(self, val: str | None = None) -> None:
        self._exch_meth = {
            None: -1,
            "directlaplacianneumann": 1,
            "ggneumann": 2,
        }[val]

    @property
    def prob_mode(self) -> int:
        return self._prob_mode

    @prob_mode.setter
    def prob_mode(self, val: str | None = None) -> None:
        self._prob_mode = {None: -1, "new": 1, "old": 2}[val]

    @property
    def solver(self) -> int:
        return self._solver

    @solver.setter
    def solver(self, val: str | None = None) -> None:
        self._solver = {None: -1, "explicit": 1, "dynamic": 2, "implicit": 3}[val]

    def run_simulation(
        self, t_end: float, nt: int, fct_h_ext: Callable, nt_h_ext: int
    ) -> list[np.ndarray | int]:
        """
        Run the micromagnetic simulation.

        Params:
            t_end: End time of the simulation.
            nt: Number of timesteps.
            fct_h_ext: Function to calculate the external field at each timestep.
            nt_h_ext: Number of timesteps for the external field.

        Outputs:
            A list containing the simulation results.
            [0] t_out (np.ndarray): Time output array.
            [1] M_out (np.ndarray): Magnetization output array.
            [2] pts (np.ndarray): Points output array.
            [3] H_exc (np.ndarray): Exchange field output array.
            [4] H_ext (np.ndarray): External field output array.
            [5] H_dem (np.ndarray): Demagnetizing field output array.
            [6] H_ani (np.ndarray): Anisotropy field output array.
            [7] Exch_mat_ntot (int): Total number of exchange matrix elements.
            [8] Exch_mat_r (np.ndarray): Exchange matrix row indices.
            [9] Exch_mat_c (np.ndarray): Exchange matrix column indices.
            [10] Exch_mat_v (np.ndarray): Exchange matrix values.
            [11] Exch_mat_nr (int): Number of rows in the exchange matrix.
            [12] Exch_mat_nc (int): Number of columns in the exchange matrix.

        """
        h_ext = np.zeros(shape=(nt_h_ext, 4), dtype=np.float64, order="F")
        h_ext[:, 0] = np.linspace(0, t_end, nt_h_ext)
        h_ext[:, 1:4] = fct_h_ext(np.linspace(0, t_end, nt_h_ext))

        result = magtensesource.fortrantopythonio.runmicromagsimulation(
            ntot=self.ntot,
            grid_n=self.grid_n,
            grid_l=self.grid_L,
            grid_type=self.grid_type,
            u_ea=self.u_ea,
            problemmode=self.prob_mode,
            solver=self.solver,
            a0=self.A0,
            ms=self.Ms,
            k0=self.K0,
            k1=self.K1,
            k2=self.K2,
            k0_arr=self.K0_arr,
            crysaxis=self.CrysAxis,
            gamma=self.gamma,
            alpha_mm=self.alpha_mm,
            maxt0=self.max_T0,
            nt_hext=nt_h_ext,
            hext=h_ext,
            nt=nt,
            t=np.linspace(0, t_end, nt),
            m0=np.concatenate((self.m0[:, 0], self.m0[:, 1], self.m0[:, 2]), axis=None),
            dem_thres=self.dem_thres,
            usecuda=self.cuda,
            dem_appr=self.dem_appr,
            n_ret=self.N_ret,
            n_file_out=self.N_file_out,
            n_load=self.N_load,
            n_file_in=self.N_file_in,
            settimedis=self.setTimeDis,
            nt_alpha=self.nt_alpha,
            alphat=self.alphat,
            tol=self.tol,
            thres=self.thres,
            usecvode=self.cvode,
            nt_conv=self.nt_conv,
            t_conv=self.t_conv,
            conv_tol=self.conv_tol,
            grid_pts=self.grid_pts,
            grid_ele=self.grid_ele,
            grid_nod=self.grid_nod,
            grid_nnod=self.grid_nnod,
            exch_nval=self.exch_nval,
            exch_nrow=self.exch_nrow,
            exch_val=self.exch_val,
            exch_rows=self.exch_rows,
            exch_rowe=self.exch_rowe,
            exch_col=self.exch_col,
            grid_abc=self.grid_abc,
            useprecision=self.precision,
            nthreadsmatlab=self.n_threads,
            n_ave=self.N_ave,
            cv=self.cv,
            usereturnhall=self.usereturnhall,
            demigstp=self.demigstp,
            exch_weigh=self.exch_weigh,
            exch_meth=self.exch_meth,
            exch_intpn=self.exch_intpn,
            passexch=self.passexch,
            exch_ncols=self.exch_ncols,
            exch_presize=self.exch_presize,
        )

        n_tot_Exch = result[7]
        result = list(result)
        result[8] = result[8][:n_tot_Exch]  # ExchMat_r
        result[9] = result[9][:n_tot_Exch]  # ExchMat_c
        result[10] = result[10][:n_tot_Exch]  # ExchMat_v

        return result
