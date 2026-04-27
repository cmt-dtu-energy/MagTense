#!/bin/bash

# Function to run script and copy results
run_and_copy() {
    local script_name=$1
    local result_folder=$2
    
    echo "Running $script_name..."
    bash "$script_name"
    
    echo "Creating result folder: $result_folder"
    mkdir -p "$result_folder"
    
    # Copy logs_perf if it exists
    if [ -d "logs_perf" ]; then
        mv logs_perf "$result_folder/"
        echo "Copied logs_perf to $result_folder/"
    fi
    
    # Move timer_logs if it exists
    if [ -d "timer_logs" ]; then
        mv timer_logs "$result_folder/"
        echo "Moved timer_logs to $result_folder/"
    fi
    
    # Move results if it exists
    if [ -d "results" ]; then
        mv results "$result_folder/"
        echo "Moved results to $result_folder/"
    fi
    
    if [ $? -eq 0 ]; then
        echo "$script_name completed successfully!"
    else
        echo "Error: $script_name failed to complete!"
        return 1
    fi
    
    echo ""
}

# Run all scripts
#run_and_copy "run_script_cuda.sh" "cuda_res"
#run_and_copy "run_script_n4_l2.sh" "n4_l2_res"
#run_and_copy "run_script_n4_l3.sh" "n4_l3_res"
#run_and_copy "run_script_n6_l2.sh" "n6_l2_res"
#run_and_copy "run_script_n6_l3.sh" "n6_l3_res"
run_and_copy "run_script_n8_l2.sh" "n8_l2_res"
run_and_copy "run_script_n8_l3.sh" "n8_l3_res"
run_and_copy "run_script_n10_l2.sh" "n10_l2_res"
run_and_copy "run_script_n10_l3.sh" "n10_l3_res"

echo "All scripts completed!"
