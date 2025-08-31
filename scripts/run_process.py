#!/usr/bin/env python3
import yaml
import argparse
import os
import subprocess
import math
import json
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Run the AMS Isotope Framework IsoProcess stage.")
    parser.add_argument("config", help="Path to the master YAML configuration file for physics settings.")
    parser.add_argument("--mode", choices=['local', 'condor'], default='local', help="Execution mode.")
    parser.add_argument("--samples", default="all", help="Comma-separated list of samples to process (e.g., ISS_Beryllium,MC_Be10).")
    parser.add_argument("--files-per-job", type=int, default=20, help="Number of input files per job.")
    args = parser.parse_args()

    # --- 1. Load Physics Configuration ---
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # --- 2. Setup Paths ---
    framework_base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    executable_path = os.path.join(framework_base, "build/apps/process/process")
    
    # Create directories for tasks and results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    task_dir = os.path.join(framework_base, f"tasks_{timestamp}")
    results_dir = os.path.join(framework_base, f"results_{timestamp}")
    os.makedirs(task_dir, exist_ok=True)
    
    if not os.path.exists(executable_path):
        print(f"Error: Executable not found at {executable_path}. Please build the project.")
        return

    # --- 3. Generate and Execute Tasks ---
    all_samples = {s['name']: s for s in config['samples_to_process']}
    samples_to_run = all_samples.keys() if args.samples == "all" else args.samples.split(',')

    for sample_name in samples_to_run:
        if sample_name not in all_samples:
            print(f"Warning: Sample '{sample_name}' not found in config. Skipping.")
            continue
        
        sample = all_samples[sample_name]
        print(f"--- Preparing tasks for sample: {sample['name']} ---")

        sample_results_dir = os.path.join(results_dir, sample_name)
        os.makedirs(sample_results_dir, exist_ok=True)
        
        input_files = sample['files']
        num_jobs = math.ceil(len(input_files) / args.files_per_job)

        for i in range(num_jobs):
            start_idx = i * args.files_per_job
            end_idx = start_idx + args.files_per_job
            file_chunk = input_files[start_idx:end_idx]
            
            task_data = {
                "sample_name": sample['name'],
                "input_files": file_chunk,
                "output_file": os.path.join(sample_results_dir, f"chunk_{i:04d}.root")
            }
            
            task_file_path = os.path.join(task_dir, f"task_{sample_name}_{i:04d}.json")
            with open(task_file_path, 'w') as f:
                json.dump(task_data, f, indent=4)

            # --- 4. Execute Task based on mode ---
            cmd = [executable_path, args.config, "--task", task_file_path]
            
            if args.mode == 'local':
                print(f"Running task locally: {task_file_path}")
                subprocess.run(cmd, check=True)
            
            elif args.mode == 'condor':
                # This part would generate a .sub file and submit it.
                # For now, it just prints the command that would go into the condor script.
                print(f"Generated Condor task: {task_file_path}")
                print(f"  Command to run: {' '.join(cmd)}")

if __name__ == "__main__":
    main()