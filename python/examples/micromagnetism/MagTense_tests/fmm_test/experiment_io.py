import os
import pickle
import matplotlib.pyplot as plt

def save_experiment_data(seed, results, metrics, fig=None):
    """
    Saves simulation results, error metrics, and the comparison plot
    into a seed-specific subfolder.
    """
    # Create the base results directory
    base_dir = "results"
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
        
    # Create the seed-specific subfolder
    subfolder = os.path.join(base_dir, "seed_" + str(seed))
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    
    # Save the results variable (binary)
    results_path = os.path.join(subfolder, "results.pkl")
    with open(results_path, 'wb') as f:
        pickle.dump(results, f)
        
    # Save the metrics variable (binary)
    metrics_path = os.path.join(subfolder, "metrics.pkl")
    with open(metrics_path, 'wb') as f:
        pickle.dump(metrics, f)
        
    # Save the figure if provided
    if fig is not None:
        fig_path = os.path.join(subfolder, "comparison_plot.png")
        fig.savefig(fig_path)
        print(f"Figure saved to: {fig_path}")

    print(f"Experiment data for seed {seed} saved successfully in '{subfolder}'.")

def load_experiment_data(seed):
    """
    Loads simulation results and metrics from a seed-specific subfolder.
    """
    subfolder = os.path.join("results", "seed_" + str(seed))
    
    results_path = os.path.join(subfolder, "results.pkl")
    metrics_path = os.path.join(subfolder, "metrics.pkl")
    
    if not os.path.exists(results_path) or not os.path.exists(metrics_path):
        raise FileNotFoundError(f"No saved data found for seed {seed} in {subfolder}")

    with open(results_path, 'rb') as f:
        results = pickle.load(f)
        
    with open(metrics_path, 'rb') as f:
        metrics = pickle.load(f)
        
    print(f"Data for seed {seed} loaded successfully.")
    return results, metrics