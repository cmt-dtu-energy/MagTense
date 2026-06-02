import numpy as np
from magtense.micromag import MicromagProblem
from grain_generator import AdaptiveGrainGenerator
import sys
sys.path.append("../../src")
from utils import *
from typing import Optional

def fmm_test_setup(
                gen: AdaptiveGrainGenerator,
                cuda: bool = True,
                cvode: bool = False,
                cone_angle_degree: float = 20,
                easy_axis_cone_degree: float = 20,
                H_ext_dir: np.ndarray = np.array([0., 0., 1.]),
                passexch: int = 0,
                input_file_name : str = "t",
                timer_log_file: str = "timer_test.log",
                trace_log_file: str = "trace_test.log",
                t_end: float = 40e-9,
                nt: int = 10,
                dummy_run: int = 0,
                nlmax: int = 2,
                allow_fmm_short_circuit: int = 0,
                fmm_min_n: int = 20000,
                fmm_nterms: int = 6, 
                exch_val: list | np.ndarray | None = None,
                exch_rows: list | np.ndarray | None = None,
                exch_col: list | np.ndarray | None = None,
                exch_nval: int = 1,
                exch_nrow: int = 1,
                exch_ncols: int = 1,
                rng: Optional[np.random.Generator] = None,  #
                ) -> MicromagProblem:
        #------------ extract mesh and grid info from generator ------------
        mesh_cart = gen.mesh_cart
        grid_info = gen.grid_info
        n_grains = gen.n_grains
        grid_pts = mesh_cart["pos_out"]
        grid_abc = mesh_cart["dims_out"]
        grain_indices = mesh_cart["GrainIndex"]  # Expected 1D array
        n_cell = grid_pts.shape[0]
        grid_L = gen.grid_l
        res = gen.mesh_params["resolution"]
        #--------------------------------------------------------------
        # --- Material Constants (Based on  standard grain problem) ---
        mu0 = 4 * np.pi * 1e-7
        alpha = 4000.0
        gamma = 0.0
        grid_type = "unstructuredPrisms"
        #----------- grain parameters -----------
        Ms_val = 1.61 / mu0      # 1.61 T
        K0_val = 4.3e6           # J/m3
        A0_val = 7.7e-12         # J/m
        #----------- intergrain parameters -----------
        Ms_ig = 0.5 / mu0        # Intergrain Ms
        A0_ig_factor = 1.0       # Intergrain exchange factor
        #--------------------------------------------------------------

        # --- Generate Easy Axis directions (mag_dir) ---
        # Using 180 degrees for fully random grain orientations on the sphere
        # centered arbitrarily on Z [source: 1]
        grain_mag_dirs = rand_spherical_cap(
                cone_angle_degree=easy_axis_cone_degree, 
                N=n_grains, 
                cone_dir=H_ext_dir, 
                rng=rng
        )

        # --- Initialize MicromagProblem ---
        problem = MicromagProblem(
                res=res,
                grid_L=grid_L,
                alpha=alpha,
                gamma=gamma,
                grid_pts=grid_pts,
                grid_abc=grid_abc,
                grid_type=grid_type,
                cuda=cuda,
                cvode=cvode,
                exch_val=exch_val,
                exch_rows=exch_rows,
                exch_cols=exch_col,
                exch_nval=exch_nval,
                exch_nrow=exch_nrow,
                exch_ncols=exch_ncols,
                passexch=passexch,
                filename=input_file_name,
                useavgn=True,
                usereturnhall=True,
        )

        # --- Assign Material Properties and Easy Axes per Grain ---
        problem.Ms = Ms_val * np.ones((n_cell, 1))
        problem.K0 = K0_val * np.ones((n_cell, 1))
        problem.A0 = A0_val * np.ones((n_cell, 1))

        # Grain 1 to N: Set Easy Axis
        for g_idx in range(1, n_grains + 1):
                cell_mask = (grain_indices == g_idx)
                problem.u_ea[cell_mask, :] = grain_mag_dirs[g_idx - 1]


        # Intergrain (ID = N+1): Set Properties
        ig_mask = (grain_indices == n_grains + 1)
        problem.Ms[ig_mask] = Ms_ig
        problem.A0[ig_mask] = A0_val * A0_ig_factor
        # u_ea for intergrain is usually 0 (already default)


        # --- Set Initial Magnetization (m0) ---------------
        # 1. Grains: Along the Easy Axis
        for g_idx in range(1, n_grains + 1):
                cell_mask = (grain_indices == g_idx)
                problem.m0[cell_mask, :] = grain_mag_dirs[g_idx - 1]
        # 2. Intergrain: Random within a cone [source: 1]
        n_ig = np.sum(ig_mask)
        if n_ig > 0:
                # Standard cone randomization around the Z-axis for the intergrain phase
                ig_m0 = rand_spherical_cap(
                        cone_angle_degree=cone_angle_degree, 
                        N=n_ig, 
                        cone_dir=H_ext_dir, 
                        rng=rng
                )
                problem.m0[ig_mask, :] = ig_m0
        #-------------------------------------------------------------

        problem.solver = "dynamic"

        t = np.linspace(0, t_end, nt)
        problem.t =  t
        problem.nt = nt


        nnbor = 27 # assmue 27 nearest neighbors for exchange
        problem.exch_presize = problem.nt * (nnbor+1)

        problem.dummy_run = dummy_run



        problem.fmm_cells_per_node = 10
        problem.fmm_eps = 1e-4
        problem.ifunif = 1
        problem.nlmin = 1
        problem.nlmax = nlmax
        problem.allow_fmm_short_circuit = allow_fmm_short_circuit
        problem.fmm_min_n = fmm_min_n
        problem.fmm_nterms = fmm_nterms



        problem.log_dir = "timer_logs"
        problem.timer_log_file = timer_log_file
        problem.trace_log_file = trace_log_file
        problem.window_enabled = 1
        problem.window_interval = 5.0
        problem.trace_enabled = 0
        problem.flush_each = 1
        problem.trace_verbose = 2

        
        #problem.usereturnhall = True

        return problem
