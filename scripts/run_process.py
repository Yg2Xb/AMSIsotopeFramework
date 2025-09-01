#!/usr/bin/env python3
import yaml
import argparse
import os
import glob
import math
import subprocess
import json
from datetime import datetime

class Job:
    """一个简单的类，用于表示工作流中的一个作业及其依赖关系。"""
    def __init__(self, name, command, output_file=None):
        self.name = name
        self.command = command
        self.output_file = output_file
        self.parents = []

    def add_parent(self, parent_job):
        self.parents.append(parent_job)

class WorkflowGenerator:
    """
    一个专业的、配置驱动的动态工作流生成器。
    它读取一个YAML配置文件和用户输入，然后生成一个完整的多级处理工作流。
    """
    def __init__(self, config_path, args):
        print("--- Initializing Workflow Generator ---")
        self.args = args
        self.config_path = os.path.abspath(config_path)
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        self.framework_base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        self.process_exe = os.path.join(self.framework_base, "build/apps/process/process")
        
        if not os.path.exists(self.process_exe):
            raise FileNotFoundError(f"C++ executable not found: {self.process_exe}")

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        chain_id = self.config.get('analysis_chain', {}).get('chain_id', 'default')
        session_name = f"iso_proc_{self.config['run_settings']['target_nucleus']}_{chain_id}_{ts}"
        self.session_dir = os.path.join(self.args.output_dir, session_name)
        
        self.dag_dir = os.path.join(self.session_dir, "dags")
        self.log_dir = os.path.join(self.session_dir, "logs")
        self.task_dir = os.path.join(self.session_dir, "tasks")
        self.l1_dir = os.path.join(self.session_dir, "l1_processed")
        self.l2_dir = os.path.join(self.session_dir, "l2_merged")
        self.final_dir = os.path.join(self.session_dir, "final_product")
        
        print(f"Session directory will be created at: {self.session_dir}")

    def run(self):
        os.makedirs(self.session_dir, exist_ok=True)
        os.makedirs(self.dag_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)
        os.makedirs(self.task_dir, exist_ok=True)
        os.makedirs(self.l1_dir, exist_ok=True)
        os.makedirs(self.l2_dir, exist_ok=True)
        os.makedirs(self.final_dir, exist_ok=True)

        if self.args.mode == 'local':
            self._execute_local()
        elif self.args.mode == 'condor':
            self._generate_condor_workflow()
    
    def _scan_files(self, sample_config):
        if self.args.mode == 'condor':
            path_pattern = sample_config['condor_path_pattern']
        else:
            if 'local_path_pattern' not in sample_config:
                 raise KeyError(f"Sample '{sample_config['name']}' is missing 'local_path_pattern' for local mode.")
            path_pattern = sample_config['local_path_pattern']
            if not os.path.isabs(path_pattern):
                path_pattern = os.path.join(self.args.input_dir, path_pattern)

        files = sorted(glob.glob(path_pattern))
        if not files:
            print(f"Warning: No .root files found matching pattern: {path_pattern}")
        return files

    def _create_task_file(self, sample_name, input_files, output_file, task_id):
        task_content = {
            "sample_name": sample_name,
            "input_files": [os.path.abspath(f) for f in input_files],
            "output_file": os.path.abspath(output_file)
        }
        task_file_path = os.path.join(self.task_dir, f"{sample_name}_task_{task_id:04d}.json")
        with open(task_file_path, 'w') as f:
            json.dump(task_content, f, indent=4)
        return task_file_path

    def _execute_local(self):
        print("\n--- Running in LOCAL mode ---")
        for sample in self.config['samples_to_process']:
            sample_name = sample['name']
            if self.args.samples != 'all' and sample_name not in self.args.samples.split(','):
                continue

            print(f"\n--- Processing sample: {sample_name} ---")
            input_files = self._scan_files(sample)
            if not input_files: continue

            num_jobs = math.ceil(len(input_files) / self.args.files_per_job)
            l1_outputs = []
            
            for i in range(num_jobs):
                chunk = input_files[i*self.args.files_per_job : (i+1)*self.args.files_per_job]
                output_file = os.path.join(self.l1_dir, f"{sample_name}_chunk_{i:04d}.root")
                l1_outputs.append(output_file)
                
                task_file = self._create_task_file(sample_name, chunk, output_file, i)
                cmd = [self.process_exe, self.config_path, "--task", task_file]

                print(f"Executing L1 job {i+1}/{num_jobs} for {sample_name}...")
                subprocess.run(cmd, check=True)

            print(f"--- Merging results for {sample_name} ---")
            
            # 修正: 将所有路径转换为绝对路径再传递给 hadd，避免路径混淆
            final_output_abs = os.path.abspath(os.path.join(self.final_dir, f"{sample_name}_final.root"))
            l1_outputs_abs = [os.path.abspath(p) for p in l1_outputs]
            
            merge_cmd = ['hadd', '-f', final_output_abs] + l1_outputs_abs
            subprocess.run(merge_cmd, check=True)
            print(f"Final output for {sample_name} is at: {final_output_abs}")

    def _generate_condor_workflow(self):
        # ... (condor logic here) ...
        pass

def main():
    parser = argparse.ArgumentParser(description="AMS Isotope Framework - Workflow Generator & Runner.")
    parser.add_argument("config", help="Path to the master YAML physics configuration file.")
    parser.add_argument("-i", "--input-dir", help="Base directory for input files (used for local mode).")
    parser.add_argument("-o", "--output-dir", required=True, help="Base directory to store all workflow outputs.")
    parser.add_argument("--samples", default="all", help="Comma-separated list of samples to process.")
    parser.add_argument("--files-per-job", type=int, default=2, help="Number of input files per L1 processing job.")
    parser.add_argument("--mode", choices=['local', 'condor'], default='local', help="Execution mode.")
    
    args = parser.parse_args()

    if args.mode == 'local' and not args.input_dir:
        parser.error("--input-dir is required for 'local' mode.")

    try:
        workflow = WorkflowGenerator(args.config, args)
        workflow.run()
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        exit(1)

if __name__ == "__main__":
    main()