import sys
from pathlib import Path
import json
import subprocess

# STEP 1: get path for current user input json file------------------

if len(sys.argv) < 2:
    raise ValueError("Usage: python run_simulation_with_extraction.py </scenarios_save/config.json>")

json_file_path = Path(sys.argv[1]).resolve()

# STEP 2: read model path from the loaded json file ---------------------------

with open(f"{json_file_path}", 'r') as f:
    config = json.load(f)

print(f"<<<< Loaded user inputs from: {json_file_path} >>>>")

model_path = Path(config["SIMULATION_MODEL"]["kivu_simstrat_path"]).resolve()

# --- validate the path --------------
if not model_path.exists():
    raise FileNotFoundError(f"Model path does not exist: {model_path}")

# STEP 3: prepare to run extraction simulation launcher file -------------
launcher_file_path = model_path / "simulation_with_extraction_launcher.py" 

if not launcher_file_path.exists():
    raise FileNotFoundError(f"Target script not found: {launcher_file_path}")

# STEP 4: Run launcher file inside the model directory
print(f"<<<Running secondary script: {launcher_file_path}>>>>")
subprocess.run(
    [sys.executable, str(launcher_file_path), str(json_file_path)],
    cwd=model_path,
    check=True
)



