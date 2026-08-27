import os
import shutil
import sys
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np

# Windows only
if hasattr(os, "add_dll_directory"):
    # First entry is the installed layout
    # (<prefix>/Lib/site-packages/magtense/../../../Library/bin), the second the
    # active environment prefix, which is what applies when running from source.
    dll_paths = [
        Path(__file__).parent / ".." / ".." / ".." / "Library" / "bin",
        Path(sys.prefix) / "Library" / "bin",
    ]

    nvidia_path = Path(__file__).parent / ".." / "nvidia"
    dll_paths += [
        nvidia_path / lib / "bin"
        for lib in ["cublas", "cuda_runtime", "cusparse", "nvjitlink"]
    ]

    for dll_path in dll_paths:
        if Path.is_dir(dll_path):
            os.add_dll_directory(dll_path)
            # libifcoremd.dll resolves its own dependency on libmmd.dll through
            # PATH, which add_dll_directory does not cover.
            os.environ["PATH"] = f"{dll_path.resolve()}{os.pathsep}{os.environ['PATH']}"

from magtense.lib import magtensesource


class MicromagProblem:
    """
    Micromagnetic problem using the Fortran implementation of MagTense.

    Notes:
        All inputs are given in SI units. Magnetisation and field strengths are in A/m rather than T.

    Args:
        grid_type: Currently supports 'uniform', 'tetrahedron' and 'unstructuredPrisms'.
            If 'uniform', grid is inferred from res and grid_L.
            If 'tetrahedron', grid is specified by grid_pts, grid_nnod and grid_ele.
            If 'unstructuredPrisms', grid is specified by grid_pts and grid_abc.
        res: Resolution of grid, i.e. number of micromagnetic tiles along x, y and z.
        grid_L: Spatial extension of simulated domain along x, y and z.
        grid_nnod: Number of nodes in the tetrahedron mesh
        grid_pts: xyz coordinates of micromagnetic tiles
        grid_abc: sidelengths of prism tiles
        prob_mode:
        solver: Options are 'explicit', 'dynamic' and 'implicit'.
            If solver = 'dynamic', a single time-varying magnetic field is constructed
            If solver = 'explicit', the equilibrium configuration is computed at several constant fields
            solver = 'implicit' has not been implemented.
            See documentation under run_simulation for details.
        hysteresis_solver: External-field stepping mode. Options are 'static'
            and 'adaptive'. The default 'static' mode preserves the predefined
            field sequence used by existing simulations.
        A0: Exchange stiffness constant [J/m].
        Ms: Saturation magnetization [A/m].
        K0: Uniaxial anisotropy constant [J/m^3], used together with u_ea.
            Mutually exclusive with K1, K2 and K0_arr.
        K1, K2: Cubic anisotropy constants [J/m^3], used together with CrysAxis.
        K0_arr: General anisotropy expansion [J/m^3], an (ntot, 6, 3) array.
        CrysAxis: Local crystal coordinate system, an (ntot, 3, 3) array.
        alpha: Dampening constant [m/(A*s)]. Product of Gilbert damping and precession parameter.
        T: Temperature [K] ('temp' in fortran part)
        gamma: Gyromagnetic factor [m/(A*s)].
        max_T0:
        nt_conv:
        conv_tol: The convergence tolerance, which is the maximum change in
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
        n_macro: How many copies of the simulated domain to have along x, y and z to represent the macrogeometry
        shiftVec: How far to shift domain copies along x, y and z when constructing the macrogeometry
        macroShape: Sidelengths of a prism representing the shape of the macrogeometry.
        sampleShape: Sidelengths of a prism representing the sample shape.
        use_cdfmm: Use dip-fmm instead of the regular dense demagnetisation
            calculation. Currently supported for uniform cuboid grids only.
            Backend selection follows MagTense: CUDA requires both a
            CUDA-enabled build and ``cuda=True``; otherwise dip-fmm prefers
            oneMKL and falls back to portable CPU execution.
        cdfmm_order: Expansion order used by dip-fmm.
        cdfmm_depth: Fixed octree depth used by dip-fmm.
        cdfmm_basis: Expansion basis, either ``"spherical"`` or ``"cartesian"``.
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
            T: int | float | list | np.ndarray | None = None,
            K0_arr: np.ndarray | None = None,
            CrysAxis: np.ndarray | None = None,
            alpha: float = 4.42e3,
            gamma: float = 2.21e5,
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
            exch_cols: list | np.ndarray | None = None,
            exch_nval: int = 1,
            exch_nrow: int = 1,
            exch_ncols: int = 1,
            exch_intpn: str | None = "extended",
            exch_meth: str | None = "directlaplacianneumann",
            exch_weigh: float = 8,
            exch_presize: int = 12,
            demigstp: int = 0,
            passexch: int = 0,
            filename: str = "t",
            cuda: bool = False,
            cvode: bool = False,
            usedemag: bool = True,
            # True to match the MATLAB DefaultMicroMagProblem, which sets useAvgN = 1
            useavgn: bool = True,
            # False to match the MATLAB DefaultMicroMagProblem, which sets ReturnHall = 0.
            # Set it to True to have run_simulation return H_exc/H_ext/H_dem/H_ani; leaving it
            # off also avoids allocating four more (nt, ntot, nt_h_ext, 3) arrays.
            usereturnhall: bool = False,
            precision: bool = False,
            n_threads: int = 1,
            N_ave: tuple[int] = (1, 1, 1),
            t_alpha: np.ndarray = np.zeros(1),  # noqa: B008
            alpha_fct=lambda t: np.atleast_2d(t).T * 0,
            sampleShape: list | np.ndarray | None = None,
            exchPBC: list | np.ndarray | None = None,
            macroShape: list | np.ndarray | None = None,
            shiftVec: list | np.ndarray | None = np.zeros(3, dtype=np.float64),
            n_macro: list | np.ndarray | None = np.zeros(3),
            hysteresis_solver: str = "static",
            rng_seed: int = 0,
            use_cdfmm: bool = False,
            cdfmm_order: int = 6,
            cdfmm_depth: int = 4,
            cdfmm_basis: str | int = "spherical",
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
        self.useavgn = useavgn
        self.exch_presize = exch_presize

        self.nt_conv = nt_conv

        self.usereturnhall = usereturnhall
        self.demigstp = demigstp
        self.usereturnhall = usereturnhall
        self.exch_presize = exch_presize

        # Set grid, solver and problem type
        self.grid_type = grid_type
        self.prob_mode = prob_mode
        self.solver = solver
        self.hysteresis_solver = hysteresis_solver

        ### Define grid ###
        # Uniform grid
        self.grid_n = np.array(res, dtype=np.int32, order="F")
        self.grid_L = np.array(grid_L, dtype=np.float64, order="F")
        # Tetrahedron and prism grids
        self.grid_pts = grid_pts
        # Prism grid
        self.grid_abc = grid_abc
        # Tetrahedron grid
        self.grid_nnod = grid_nnod
        self.grid_ele = np.zeros(shape=(4, ntot), dtype=np.float64, order="F")
        self.grid_nod = np.zeros(shape=(grid_nnod, 3), dtype=np.float64, order="F")

        # Set macrogeometry
        self.n_macro = n_macro
        self.shiftVec = shiftVec
        self.macroShape = macroShape
        self.sampleShape = sampleShape
        self.exchPBC = exchPBC

        # Set material parameters
        self.A0 = A0
        self.Ms = Ms
        self.K0 = K0
        self.K1 = K1
        self.K2 = K2
        self.K0_arr = K0_arr
        self.T = T

        # Set initial state
        self.m0 = m0

        # --- Set the local crystal coordinates to the three Cartesian axis
        self.CrysAxis = CrysAxis
        self.u_ea = np.zeros(shape=(ntot, 3), dtype=np.float64, order="F")  # Anisotropy axes

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

        # Seed for the stochastic thermal field, drawn inside Fortran. 0 keeps the compiler
        # default sequence, which is identical on every run; a positive value seeds
        # deterministically with that value; a negative value seeds from the clock, which is
        # what independent Monte-Carlo runs need. NumPy's own seed has no effect on it.
        self.rng_seed = int(rng_seed)

        self.nt_alpha = len(t_alpha)
        self.alphat = np.zeros(shape=(self.nt_alpha, 2), dtype=np.float64, order="F")
        self.alphat[:, 0] = t_alpha
        self.alphat[:, 1] = alpha_fct(t_alpha)

        self.exch_nval = exch_nval
        self.exch_nrow = exch_nrow
        self.exch_ncols = exch_ncols
        self.exch_val = exch_val
        self.exch_rows = exch_rows
        self.exch_cols = exch_cols

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
        self.usedemag = int(usedemag)
        self.useavgn = int(useavgn)
        self.usereturnhall = int(usereturnhall)
        self.precision = int(precision)
        self.n_threads = n_threads
        self.N_ave = np.array(N_ave, dtype=np.int32, order="F")


        #--------- default FMM parameters -----------------
        self.dummy_run = 0
        self.fmm_cells_per_node = 100
        self.fmm_eps = 1e-4
        self.ifunif = 1
        self.nlmin = 1
        self.nlmax = 5
        self.allow_fmm_short_circuit = 1
        self.fmm_min_n = 20000
        self.fmm_nterms = -1
        # Off by default: the FMM demag path only pays off for large problems, and it is
        # ignored altogether unless the library is built with USE_FMM3D=1. Opting in
        # explicitly means a problem does not silently change its demag path when the
        # library is rebuilt with FMM enabled. Set `problem.use_fmm = 1` to use it.
        self.use_fmm = 0
        #--------------------------------------------------

        # dip-fmm is independent of the legacy FMM3D backend above. Its current
        # MagTense adapter supports uniform cuboids and FP32 field arrays.
        basis_values = {"spherical": 0, "cartesian": 1}
        if isinstance(cdfmm_basis, str):
            if cdfmm_basis not in basis_values:
                raise ValueError(f"Unknown dip-fmm basis: {cdfmm_basis}")
            cdfmm_basis = basis_values[cdfmm_basis]
        if use_cdfmm and self.grid_type != 1:
            raise ValueError("dip-fmm currently requires grid_type='uniform'")
        if cdfmm_order < 1 or cdfmm_depth < 1:
            raise ValueError("cdfmm_order and cdfmm_depth must be positive")
        if int(cdfmm_basis) not in basis_values.values():
            raise ValueError("cdfmm_basis must be 0 (spherical) or 1 (cartesian)")

        self.use_cdfmm = int(use_cdfmm)
        self.cdfmm_order = int(cdfmm_order)
        self.cdfmm_depth = int(cdfmm_depth)
        self.cdfmm_basis = int(cdfmm_basis)
        self.cdfmm_precision = 0  # MagTense stores the demag field in FP32.

        #---------- timer and trace parameters ----------
        self.log_dir = "logs"
        self.timer_log_file = "timing.log"
        self.trace_log_file = "trace.log"
        self.window_enabled = 1
        self.window_interval = 30.0
        self.trace_enabled = 0
        self.flush_each = 1
        self.trace_verbose = 1
        #-----------------------------------------------

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
    def exch_cols(self) -> list | np.ndarray | None:
        return self._exch_cols

    @exch_cols.setter
    def exch_cols(self, val: list | np.ndarray | None) -> None:
        if val is None:
            self._exch_cols = np.zeros(
                shape=(self.exch_nval,), dtype=np.int32, order="F"
            )
        else:
            assert np.asarray(val).shape == (self.exch_nval,)
            self._exch_cols = np.asarray(val, dtype=np.int32, order="F")

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
    def T(self) -> int | float | list | np.ndarray | None:
        return self._T

    @T.setter
    def T(self, val: int | float | list | np.ndarray | None) -> None:
        if val is None:
            self._T = 0.0 + np.zeros(self.ntot, dtype=np.float64, order="F")

        elif isinstance(val, (int, float)):
            self._T = val + np.zeros(self.ntot, dtype=np.float64, order="F")

        else:
            assert np.asarray(val).shape == self.ntot
            self._T = np.asarray(val, dtype=np.float64, order="F")

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

        # The LLG equation is written for a unit vector, so an un-normalised m0 silently
        # rescales the effective field. MATLAB's DefaultMicroMagProblem.struct() has always
        # checked and corrected this; do the same here rather than passing it through.
        # Rows of zero length are left alone - there is no direction to normalise them to,
        # and dividing them would turn them into NaN.
        norm = np.linalg.norm(self._m0, axis=1)
        nonzero = norm > 0
        if np.any(~nonzero):
            warnings.warn(
                "Zero magnetization in initial array", RuntimeWarning, stacklevel=2
            )
        if np.any(np.abs(norm[nonzero] - 1.0) >= 1e-4):
            warnings.warn(
                "Initial magnetization array not normalized -- normalizing",
                RuntimeWarning,
                stacklevel=2,
            )
        self._m0[nonzero] /= norm[nonzero, np.newaxis]

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

    @property
    def hysteresis_solver(self) -> int:
        return self._hysteresis_solver

    @hysteresis_solver.setter
    def hysteresis_solver(self, val: str = "static") -> None:
        try:
            self._hysteresis_solver = {"static": 1, "adaptive": 2}[val]
        except KeyError as error:
            raise ValueError(
                "hysteresis_solver must be 'static' or 'adaptive'"
            ) from error

    @property
    def n_macro(self) -> list | np.ndarray | None:
        return self._n_macro

    @n_macro.setter
    def n_macro(self, val: list | np.ndarray | None) -> None:
        self._n_macro = np.zeros(shape=3, dtype=np.int32, order="F")

        if val is None:
            pass

        elif isinstance(val, (list, np.ndarray)):
            self._n_macro[:] = val

    @property
    def shiftVec(self) -> list | np.ndarray | None:
        return self._shiftVec

    @shiftVec.setter
    def shiftVec(self, val: list | np.ndarray | None) -> None:
        self._shiftVec = np.zeros(shape=3, dtype=np.float64, order="F")

        if val is None:
            pass

        elif isinstance(val, (list, np.ndarray)):
            self._shiftVec[:] = val

    @property
    def macroShape(self) -> list | np.ndarray | None:
        return self._macroShape

    @macroShape.setter
    def macroShape(self, val: list | np.ndarray | None) -> None:
        self._macroShape = np.ones(shape=3, dtype=np.float64, order="F")

        if val is None:
            pass

        elif isinstance(val, (list, np.ndarray)):
            self._macroShape[:] = val

    @property
    def sampleShape(self) -> list | np.ndarray | None:
        return self._sampleShape

    @sampleShape.setter
    def sampleShape(self, val: list | np.ndarray | None) -> None:
        self._sampleShape = np.ones(shape=3, dtype=np.float64, order="F")

        if val is None:
            pass

        elif isinstance(val, (list, np.ndarray)):
            self._sampleShape[:] = val

    @property
    def exchPBC(self) -> np.ndarray:
        return self._exchPBC

    @exchPBC.setter
    def exchPBC(self, val: list | tuple | np.ndarray | None) -> None:
        self._exchPBC = np.zeros(shape=3, dtype=np.int32)

        if val is None:
            pass

        elif isinstance(val, (list, tuple, np.ndarray)):
            val_arr = np.asarray(val, dtype=np.int32)
            assert val_arr.shape == (3,), "exchPBC must be an array of length 3"
            assert np.all((val_arr == 0) | (val_arr == 1)), "exchPBC elements must be 0 or 1"
            self._exchPBC[:] = val_arr

    @property
    def usereturnhall(self) -> int | None:
        return self._usereturnhall

    @usereturnhall.setter
    def usereturnhall(self, val: bool | None = None) -> None:
        self._usereturnhall = 1 if val else 0

    def run_simulation(
            self, t_end: float, nt: int, fct_h_ext: Callable, nt_h_ext: int
    ) -> list[np.ndarray | int]:
        """
        Run the micromagnetic simulation.

        Params:
            t_end: End time of the simulation.
            nt: Number of timesteps.
            fct_h_ext: Function to calculate the external field at each timestep.
            nt_h_ext: Number of time coordinates at which to evaluate fct_h_ext.
                      Evaluation times are uniformly distributed from t=0 to t=t_end
                      If solver = "dynamic", a single time-varying magnetic field is constructed
                       by linear interpolation of the evaluation points.
                      If solver = "explicit", then each of the nt_h_ext field evaluations are treated
                       as a distinct, constant field and the equilibrium solution is computed for each field

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

        if self.solver == 2: # dynamic solver
            nt_h_ext_out = 1
        else:
            nt_h_ext_out = nt_h_ext 

        if self.solver not in (1, 2):
            print("Only 'dynamic' and 'explicit' solvers are implemented")

        result = magtensesource.fortrantopythonio.runmicromagsimulation(
            ntot=self.ntot,
            grid_type=self.grid_type,
            grid_n=self.grid_n,
            grid_l=self.grid_L,
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
            temperature=self.T,
            maxt0=self.max_T0,
            nt_hext=nt_h_ext,
            nt_hext_out = nt_h_ext_out,
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
            usedemag=self.usedemag,
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
            exch_cols=self.exch_cols,
            grid_abc=self.grid_abc,
            useprecision=self.precision,
            nthreadsmatlab=self.n_threads,
            n_ave=self.N_ave,
            cv=self.cv,
            usereturnhall=self.usereturnhall,
            useavgn=self.useavgn,
            demigstp=self.demigstp,
            exch_weigh=self.exch_weigh,
            exch_meth=self.exch_meth,
            exch_intpn=self.exch_intpn,
            passexch=self.passexch,
            exch_ncols=self.exch_ncols,
            exch_presize=self.exch_presize,
            n_macro=self.n_macro,
            shiftvec=self.shiftVec,
            macroshape=self.macroShape,
            sampleshape=self.sampleShape,
            exchpbc=self.exchPBC,
            hysteresis_solver=1,
            h_start=np.zeros(3, dtype=np.float64),
            h_end=np.zeros(3, dtype=np.float64),
            dh_initial=0.0,
            dh_min=0.0,
            dh_max=0.0,
            maxhextsteps=0,
            dm_min=0.0,
            dm_target=0.0,
            dm_reject=0.0,
            dh_grow=0.0,
            dh_shrink=0.0,
            switch_refine_dh=0.0,
            use_switch_refine=0,
            dummy_run=self.dummy_run,
            fmm_cells_per_node=self.fmm_cells_per_node,
            eps_fmm=self.fmm_eps,
            ifunif=self.ifunif,
            nlmin=self.nlmin,
            nlmax=self.nlmax,
            allow_fmm_short_circuit=self.allow_fmm_short_circuit,
            fmm_min_n=self.fmm_min_n,
            fmm_nterms=self.fmm_nterms,
            usefmm=self.use_fmm,
            usecdfmm=self.use_cdfmm,
            cdfmm_order=self.cdfmm_order,
            cdfmm_depth=self.cdfmm_depth,
            cdfmm_basis=self.cdfmm_basis,
            cdfmm_precision=self.cdfmm_precision,
            log_dir=self.log_dir,
            timer_log_file=self.timer_log_file,
            trace_log_file=self.trace_log_file,
            window_enabled=self.window_enabled,
            window_interval=self.window_interval,
            trace_enabled=self.trace_enabled,
            flush_each=self.flush_each,
            trace_verbose=self.trace_verbose,
            rng_seed=self.rng_seed
        )

        result = list(result)
        result.pop(7)  # n_Hext_accepted is only public for adaptive hysteresis.
        n_tot_Exch = result[7]
        result[8] = result[8][:n_tot_Exch]  # ExchMat_r
        result[9] = result[9][:n_tot_Exch]  # ExchMat_c
        result[10] = result[10][:n_tot_Exch]  # ExchMat_v

        return result

    def run_hysteresis(
            self, H_ext: np.ndarray
    ) -> list[np.ndarray | int]:
        """
        Run the micromagnetic hysteresis simulation.

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

        if self.hysteresis_solver != 1:
            raise ValueError(
                "run_hysteresis requires hysteresis_solver='static'"
            )

        #----------- we always run hysteresis with explicit solver so no need to check here
        nt_h_ext = H_ext.shape[0]
        nt_h_ext_out = nt_h_ext

        if self.solver not in (1, 2):
            print("Only 'dynamic' and 'explicit' solvers are implemented")


        result = magtensesource.fortrantopythonio.runmicromagsimulation(
            ntot=self.ntot,
            grid_type=self.grid_type,
            grid_n=self.grid_n,
            grid_l=self.grid_L,
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
            temperature=self.T,
            maxt0=self.max_T0,
            nt_hext=nt_h_ext,
            nt_hext_out = nt_h_ext_out,
            hext=H_ext,
            nt=self.nt,
            t=self.t,
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
            usedemag=self.usedemag,
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
            exch_cols=self.exch_cols,
            grid_abc=self.grid_abc,
            useprecision=self.precision,
            nthreadsmatlab=self.n_threads,
            n_ave=self.N_ave,
            cv=self.cv,
            usereturnhall=self.usereturnhall,
            useavgn=self.useavgn,
            demigstp=self.demigstp,
            exch_weigh=self.exch_weigh,
            exch_meth=self.exch_meth,
            exch_intpn=self.exch_intpn,
            passexch=self.passexch,
            exch_ncols=self.exch_ncols,
            exch_presize=self.exch_presize,
            n_macro=self.n_macro,
            shiftvec=self.shiftVec,
            macroshape=self.macroShape,
            sampleshape=self.sampleShape,
            exchpbc=self.exchPBC,
            hysteresis_solver=self.hysteresis_solver,
            h_start=np.zeros(3, dtype=np.float64),
            h_end=np.zeros(3, dtype=np.float64),
            dh_initial=0.0,
            dh_min=0.0,
            dh_max=0.0,
            maxhextsteps=0,
            dm_min=0.0,
            dm_target=0.0,
            dm_reject=0.0,
            dh_grow=0.0,
            dh_shrink=0.0,
            switch_refine_dh=0.0,
            use_switch_refine=0,
            dummy_run=self.dummy_run,
            fmm_cells_per_node=self.fmm_cells_per_node,
            eps_fmm=self.fmm_eps,
            ifunif=self.ifunif,
            nlmin=self.nlmin,
            nlmax=self.nlmax,
            allow_fmm_short_circuit=self.allow_fmm_short_circuit,
            fmm_min_n=self.fmm_min_n,
            fmm_nterms=self.fmm_nterms,
            usefmm=self.use_fmm,
            usecdfmm=self.use_cdfmm,
            cdfmm_order=self.cdfmm_order,
            cdfmm_depth=self.cdfmm_depth,
            cdfmm_basis=self.cdfmm_basis,
            cdfmm_precision=self.cdfmm_precision,
            log_dir=self.log_dir,
            timer_log_file=self.timer_log_file,
            trace_log_file=self.trace_log_file,
            window_enabled=self.window_enabled,
            window_interval=self.window_interval,
            trace_enabled=self.trace_enabled,
            flush_each=self.flush_each,
            trace_verbose=self.trace_verbose,
            rng_seed=self.rng_seed
        )

        result = list(result)
        result.pop(7)  # n_Hext_accepted is only public for adaptive hysteresis.
        n_tot_Exch = result[7]
        result[8] = result[8][:n_tot_Exch]  # ExchMat_r
        result[9] = result[9][:n_tot_Exch]  # ExchMat_c
        result[10] = result[10][:n_tot_Exch]  # ExchMat_v

        return result

    def run_hysteresis_adaptive(
            self,
            H_start: np.ndarray,
            H_end: np.ndarray,
            dH_initial: float,
            dH_min: float,
            dH_max: float,
            max_steps: int,
            dM_min: float = 1e-3,
            dM_target: float = 1e-2,
            dM_reject: float = 5e-2,
            dH_grow: float = 1.5,
            dH_shrink: float = 0.75,
            switch_refine_dH: float | None = None,
    ) -> list[np.ndarray | int]:
        """
        Run a micromagnetic hysteresis simulation with adaptive external-field steps.

        The adaptive accept/reject loop is executed by the Fortran backend in a
        single call. Output arrays are preallocated to ``max_steps`` in Fortran
        and sliced here to include only accepted field steps.
        """

        if self.hysteresis_solver != 2:
            raise ValueError(
                "run_hysteresis_adaptive requires hysteresis_solver='adaptive'"
            )
        if self.solver != 1:
            raise ValueError("Adaptive hysteresis requires the explicit solver")

        H_start = np.asarray(H_start, dtype=np.float64)
        H_end = np.asarray(H_end, dtype=np.float64)
        if H_start.shape != (3,) or H_end.shape != (3,):
            raise ValueError("H_start and H_end must be 3-vectors")
        if max_steps <= 0:
            raise ValueError("max_steps must be positive")
        if dH_initial <= 0 or dH_min <= 0 or dH_max <= 0:
            raise ValueError("dH_initial, dH_min and dH_max must be positive")
        if dH_min > dH_max:
            raise ValueError("dH_min must be smaller than or equal to dH_max")
        if np.allclose(H_start, H_end, rtol=0.0, atol=0.0):
            raise ValueError("H_start and H_end must differ")
        if dH_grow <= 1.0:
            raise ValueError("dH_grow must be larger than 1")
        if not (0.0 < dH_shrink < 1.0):
            raise ValueError("dH_shrink must be between 0 and 1")

        # The Fortran loop records the starting field in slot 1 and then one slot per accepted
        # step, so it needs max_steps + 1 slots. Sizing these to max_steps overran problem%Hext
        # by one row on the last step and left the returned M/H arrays one field short of
        # solution%M_out, which made the copy-out a non-conforming array assignment.
        nt_h_ext = int(max_steps) + 1
        nt_h_ext_out = int(max_steps) + 1
        use_switch_refine = switch_refine_dH is not None
        switch_refine_value = float(switch_refine_dH) if use_switch_refine else 0.0
        if use_switch_refine and switch_refine_value <= 0:
            raise ValueError("switch_refine_dH must be positive when provided")

        h_ext = np.zeros(shape=(nt_h_ext, 4), dtype=np.float64, order="F")
        h_ext[:, 1:4] = H_start

        result = magtensesource.fortrantopythonio.runmicromagsimulation(
            ntot=self.ntot,
            grid_type=self.grid_type,
            grid_n=self.grid_n,
            grid_l=self.grid_L,
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
            temperature=self.T,
            maxt0=self.max_T0,
            nt_hext=nt_h_ext,
            nt_hext_out=nt_h_ext_out,
            hext=h_ext,
            nt=self.nt,
            t=self.t,
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
            usedemag=self.usedemag,
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
            exch_cols=self.exch_cols,
            grid_abc=self.grid_abc,
            useprecision=self.precision,
            nthreadsmatlab=self.n_threads,
            n_ave=self.N_ave,
            cv=self.cv,
            usereturnhall=self.usereturnhall,
            useavgn=self.useavgn,
            demigstp=self.demigstp,
            exch_weigh=self.exch_weigh,
            exch_meth=self.exch_meth,
            exch_intpn=self.exch_intpn,
            passexch=self.passexch,
            exch_ncols=self.exch_ncols,
            exch_presize=self.exch_presize,
            n_macro=self.n_macro,
            shiftvec=self.shiftVec,
            macroshape=self.macroShape,
            sampleshape=self.sampleShape,
            exchpbc=self.exchPBC,
            hysteresis_solver=self.hysteresis_solver,
            h_start=H_start,
            h_end=H_end,
            dh_initial=float(dH_initial),
            dh_min=float(dH_min),
            dh_max=float(dH_max),
            maxhextsteps=int(max_steps),
            dm_min=float(dM_min),
            dm_target=float(dM_target),
            dm_reject=float(dM_reject),
            dh_grow=float(dH_grow),
            dh_shrink=float(dH_shrink),
            switch_refine_dh=switch_refine_value,
            use_switch_refine=int(use_switch_refine),
            dummy_run=self.dummy_run,
            fmm_cells_per_node=self.fmm_cells_per_node,
            eps_fmm=self.fmm_eps,
            ifunif=self.ifunif,
            nlmin=self.nlmin,
            nlmax=self.nlmax,
            allow_fmm_short_circuit=self.allow_fmm_short_circuit,
            fmm_min_n=self.fmm_min_n,
            fmm_nterms=self.fmm_nterms,
            usefmm=self.use_fmm,
            usecdfmm=self.use_cdfmm,
            cdfmm_order=self.cdfmm_order,
            cdfmm_depth=self.cdfmm_depth,
            cdfmm_basis=self.cdfmm_basis,
            cdfmm_precision=self.cdfmm_precision,
            log_dir=self.log_dir,
            timer_log_file=self.timer_log_file,
            trace_log_file=self.trace_log_file,
            window_enabled=self.window_enabled,
            window_interval=self.window_interval,
            trace_enabled=self.trace_enabled,
            flush_each=self.flush_each,
            trace_verbose=self.trace_verbose,
            rng_seed=self.rng_seed
        )

        result = list(result)
        n_accepted = int(result[7])
        n_tot_Exch = result[8]

        fixed_order_result = result[0:7] + result[8:14] + [n_accepted]
        for idx in (1, 3, 4, 5, 6):
            fixed_order_result[idx] = fixed_order_result[idx][:, :, :n_accepted, :]
        fixed_order_result[8] = fixed_order_result[8][:n_tot_Exch]  # ExchMat_r
        fixed_order_result[9] = fixed_order_result[9][:n_tot_Exch]  # ExchMat_c
        fixed_order_result[10] = fixed_order_result[10][:n_tot_Exch]  # ExchMat_v

        return fixed_order_result
