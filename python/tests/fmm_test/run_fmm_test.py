
import numpy as np
from grain_generator import AdaptiveGrainGenerator
from magtense.micromag import MicromagProblem
from fmm_test_setup import fmm_test_setup
from experiment_io import save_experiment_data
from auxiliary_functions import calculate_error_metrics, plot_all_comparisons, validate_results
import sys
import argparse
from pathlib import Path
from itertools import product

def fmm_test_main(
                quick_test: bool = False,
                H_ext_dir_type: int = 1,
                seed: int = 100
                ) :
    rng_seed = seed
    rng = np.random.default_rng(rng_seed)
    gen = AdaptiveGrainGenerator(n_grains=25, offset_d=4e-9, rng=rng)
    gen.generate(res_base=6, num_refinements=2)

    cuda = True

    exma_file =  "seed_" + str(rng_seed) + "_exMa.npz"
    demag_file = "seed_" + str(rng_seed) + "_demag.bin"


    if H_ext_dir_type == 1:
        H_ext_dir = np.array([0., 0., 1.])     # For fixed direction along Z
    elif H_ext_dir_type == 2:
        H_ext_dir = rng.random(3)          # For random direction, but "fixed" for reproducibility
    elif H_ext_dir_type == 3:
        H_ext_dir = np.random.random(3)    # For random direction, without reproducibility
    else: 
        raise ValueError("Invalid H_ext_dir_type. Must be 1, 2, or 3.")

    def h_ext_fct(t: np.ndarray) -> np.ndarray:
        return H_ext_dir


    #------------- dummy run to store tensors -------------
    t_end = 40e-9

    if not Path(exma_file).exists():
        print("Running problem to generate exMa...")

        rng = np.random.default_rng(rng_seed)
        problem = fmm_test_setup(
            gen=gen,
            allow_fmm_short_circuit=1,
            timer_log_file="dummy_" + str(rng_seed) + "_timer.log",
            trace_log_file="dummy_" + str(rng_seed) + "_trace.log",
            input_file_name=demag_file,
            rng=np.random.default_rng(rng_seed),
            H_ext_dir=H_ext_dir,
            dummy_run=1,
            cuda = cuda
            )

        res = problem.run_simulation(
            t_end=t_end,
            nt=problem.nt,
            fct_h_ext=h_ext_fct,
            nt_h_ext=problem.nt)
        exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = res[7:]
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
    exch_nval, ExchMat_r, ExchMat_c, ExchMat_v, exch_nrow, exch_ncols = np.load(
        exma_file, allow_pickle=True
    ).values()
    #------------- end of dummy run -------------







    #--------------- run actual tests   -------------
    # Common configuration to avoid repetition
    common_setup = {
        "gen": gen,
        "passexch": 1,
        "input_file_name": demag_file,
        "exch_val": ExchMat_v,
        "exch_rows": ExchMat_r,
        "exch_col": ExchMat_c,
        "exch_nval": exch_nval,
        "exch_nrow": exch_nrow,
        "exch_ncols": exch_ncols,
        "H_ext_dir": H_ext_dir,
    }

    # Define the parameter space
    if quick_test:
        nlmax_list = [2]
    else:
        nlmax_list = [2, 3]
    nterms_list = [6, 10, 15]

    results = {}

    # --- 1. Run CUDA Baseline ---
    print("Running CUDA Baseline...")
    problem_cuda = fmm_test_setup(
        **common_setup,
        allow_fmm_short_circuit=1,
        timer_log_file=f"cuda_{rng_seed}_timer.log",
        trace_log_file=f"cuda_{rng_seed}_trace.log",
        rng = np.random.default_rng(rng_seed),
        cuda = cuda
    )
    res_cuda = problem_cuda.run_simulation(
        t_end=t_end, nt=problem_cuda.nt, 
        fct_h_ext=h_ext_fct, nt_h_ext=problem_cuda.nt
    )
    results["CUDA"] = {"Mout": res_cuda[1], "Hdem": res_cuda[5]}

    # --- 2. Run FMM Iterations ---
    for nl, nt_fmm in product(nlmax_list, nterms_list):
        config_label = f"FMM_L{nl}_N{nt_fmm}"
        print(f"Running {config_label}...")
        
        prob = fmm_test_setup(
            **common_setup,
            allow_fmm_short_circuit=0,
            nlmax=nl,
            fmm_nterms=nt_fmm,
            timer_log_file=f"{config_label}_{rng_seed}_timer.log",
            trace_log_file=f"{config_label}_{rng_seed}_trace.log",
            rng = np.random.default_rng(rng_seed),
            cuda = cuda
        )
        
        res = prob.run_simulation(
            t_end=t_end, nt=prob.nt, 
            fct_h_ext=h_ext_fct, nt_h_ext=prob.nt
        )
        results[config_label] = {"Mout": res[1], "Hdem": res[5]}
    #----------------------------------------------------------------






    #---------------- Analyze and Save Results ----------------
    # --- 1.  Calculate Metrics ---
    # (Assuming your loop to populate 'results' has finished)
    metrics = calculate_error_metrics(results)

    # --- 2. Generate and Capture the Plot ---
    comparison_fig = plot_all_comparisons(results)
    #plt.show() # Display it in your notebook/UI

    # --- 3. Save Everything ---
    save_experiment_data(
        seed=rng_seed, 
        results=results, 
        metrics=metrics, 
        fig=comparison_fig
    )
    #----------------------------------------------------


    #------------------ validate results against thresholds ------------------
    fmm_thresholds = {
        "L2": {
            "Hdem_rel_avg": 1e-3,  
            "Hdem_rel_max": 0.05,
            "Mout_max": 0.03
        },
        "L3": {
            "Hdem_rel_avg": 1e-2,  
            "Hdem_rel_max": 0.4,
            "Mout_max": 0.25
        }
    }

    # Run the adaptive validation
    all_pass = validate_results(metrics, fmm_thresholds)
    #----------------------------------------------------------------


    sys.exit(0 if all_pass else 1) #





if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run grain hysteresis experiment"
    )


    parser.add_argument("--quick-test", action="store_true",
                        help="Only run fast FMM configurations (nlmax=2) for a quick sanity check.")

    parser.add_argument("--Hext-dir-type", type=int, default=2,
                        help="Type of external magnetic field direction (1: fixed along Z, 2: random, 3: random without reproducibility)")
    

    parser.add_argument("--seed", type=int, default=100,
                        help="Random seed for reproducibility (used for grain generation and random H_ext_dir if applicable)")
    
    args = parser.parse_args()

    print("Running test with parameters:")
    print(f"  --quick-test: {args.quick_test}")
    print(f"  --Hext-dir-type: {args.Hext_dir_type}")
    print(f"  --seed: {args.seed}")

    fmm_test_main(
        quick_test=args.quick_test,
        H_ext_dir_type=args.Hext_dir_type,
        seed=args.seed
    )