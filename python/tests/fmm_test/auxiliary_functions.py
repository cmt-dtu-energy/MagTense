import numpy as np
import matplotlib.pyplot as plt

def plot_all_comparisons(results, target_cell_idx=9):
    components = ['x', 'y', 'z']
    fields = ['Hdem', 'Mout']
    cuda = results["CUDA"]
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    fig.suptitle(f"FMM vs CUDA Error Analysis (Cell {target_cell_idx})", fontsize=16, fontweight='bold')

    for col, field in enumerate(fields):
        for row, comp in enumerate(components):
            ax = axes[row, col]
            cuda_data = cuda[field][target_cell_idx, :, 0, row]
            norm = np.mean(np.abs(cuda_data)) if field == "Hdem" else 1.0
            
            for label, data in results.items():
                if label == "CUDA": continue
                fmm_data = data[field][target_cell_idx, :, 0, row]
                error = (cuda_data - fmm_data) / norm
                ax.plot(error, label=label)

            ax.set_title(f"{field} - {comp.upper()}")
            ax.grid(True, linestyle='--', alpha=0.6)
            if row == 0 and col == 1:
                ax.legend(loc='upper right', fontsize='8')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # Return the figure object instead of just calling plt.show()
    return fig

def calculate_error_metrics(results, target_cell_idx=9):
    metrics_list = []
    cuda = results["CUDA"]
    
    for label, data in results.items():
        if label == "CUDA": continue
        
        entry = {"Config": label}
        
        for field in ["Hdem", "Mout"]:
            # Combine x, y, z for a total field error [time, 3]
            c_data = cuda[field][target_cell_idx, :, 0, :]
            f_data = data[field][target_cell_idx, :, 0, :]
            
            diff = np.abs(c_data - f_data)
            
            # Standard Errors
            entry[f"{field}_avg"] = np.mean(diff)
            entry[f"{field}_max"] = np.max(diff)
            
            # Relative Errors (Hdem only)
            if field == "Hdem":
                norm = np.mean(np.abs(c_data))
                entry["Hdem_rel_avg"] = entry["Hdem_avg"] / norm
                entry["Hdem_rel_max"] = entry["Hdem_max"] / norm
                
        metrics_list.append(entry)
        
    return metrics_list

def validate_results(metrics, threshold_map, default_key="L2"):
    """
    Validates results based on config-specific thresholds.
    threshold_map: {
        "L2": {"Hdem_rel_avg": 0.01, ...},
        "L3": {"Hdem_rel_avg": 0.02, ...}
    }
    """
    header = f"{'Configuration':<15} | {'Status':<8} | {'Limit Used':<10} | Errors Found"
    print(f"\n{header}")
    print("-" * 85)
    
    overall_pass = True
    
    for m in metrics:
        # Determine which threshold set to use based on the label
        # Checks if "L2", "L3", etc., is in the config string
        active_thresholds = None
        limit_label = "Default"
        
        for key in threshold_map:
            if key in m['Config']:
                active_thresholds = threshold_map[key]
                limit_label = key
                break
        
        # Fallback to default if no match found
        if active_thresholds is None:
            active_thresholds = threshold_map.get(default_key, {})
            limit_label = default_key

        failures = []
        for metric_name, limit in active_thresholds.items():
            val = m.get(metric_name, 0)
            if val > limit:
                failures.append(f"{metric_name}({val:.2e} > {limit:.2e})")
        
        status = "✅ PASS" if not failures else "❌ FAIL"

        overall_pass = overall_pass and (not failures)


        error_str = ", ".join(failures) if failures else "Within tolerance"
        print(f"{m['Config']:<15} | {status:<8} | {limit_label:<10} | {error_str}")
        
    print("-" * 85)
    if overall_pass:
        print("OVERALL TEST SUITE: SUCCESS 🚀")
    else:
        print("OVERALL TEST SUITE: FAILURE 💀")
    
    return overall_pass