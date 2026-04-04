import os
import sys
import subprocess

# Get the absolute path of the project root directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Define standard directory paths
CONFIG = {
    "INPUT_DIR": os.path.join(BASE_DIR, "data", "input"),
    "REF_DIR": os.path.join(BASE_DIR, "data", "reference"),
    "OUTPUT_DIR": os.path.join(BASE_DIR, "output")
}

def check_env():
    """Check and create necessary directories if they don't exist."""
    for path in CONFIG.values():
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Created directory: {path}")

def run_stage(script_name):
    """Run a specific stage script located in the src folder."""
    script_path = os.path.join(BASE_DIR, "src", script_name)
    print(f"\n>>> Executing: {script_name}...")
    
    # Pass the BASE_DIR as an environment variable to the subprocess
    env = os.environ.copy()
    env["PROJ_BASE"] = BASE_DIR
    
    result = subprocess.run([sys.executable, script_path], env=env)
    
    if result.returncode != 0:
        print(f"Error: {script_name} failed to execute.")
        sys.exit(1)

if __name__ == "__main__":
    print("=== PanOperon-Tool Automated Pipeline ===")
    check_env()
    
    # The 5-stage pipeline
    stages = [
        "stage_0_input.py",
        "stage_1_preprocess.py",
        "stage_2_popgid.py",
        "stage_3_network.py",
        # "stage_4_finalize.py"
    ]
    
    for stage in stages:
        run_stage(stage)
        
    print("\n[Success] Pipeline finished. Intermediate files cleaned.")
    print(f"Final result: {os.path.join(BASE_DIR, 'output/network_results/network_final.csv')}")
