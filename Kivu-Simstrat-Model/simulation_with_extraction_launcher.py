import sys
from pathlib import Path
import subprocess

from src_extraction.kivu_simstrat_processor import create_scenarios_extraction, delete_existing_extraction_parfile, copy_original_parfile, update_simstart_file_par
from src_extraction.write_inflows import write_simstrat_inflows_file, write_aed_inflows_file


#STEP 1: Get Json file path -----------------------------------------

# Ensure a JSON file path was passed
#if len(sys.argv) < 2:
    #raise ValueError("Usage: python simulation_with_extraction_launcher.py <config.json>")

# Get the JSON file path passed from run_file.py
#json_file_path = Path(sys.argv[1]).resolve()


# STEP 2: load the current location / model path
model_path = Path(__file__).resolve().parent

# Get the JSON file is one level out, inside scenarios_save
scenario_folder = model_path.parent / "scenarios_save"
if not scenario_folder.exists():
    raise FileNotFoundError(f"Scenario folder not found: {scenario_folder}")
json_files = list(scenario_folder.glob("*.json"))
if not json_files:
    raise FileNotFoundError(f"No JSON scenario file found in: {scenario_folder}")

# Use the first (and only) JSON file
json_file_path = json_files[0]

if not json_file_path.exists():
    raise FileNotFoundError(f"Could not find scenario file at: {json_file_path}")

#---- CREATE SCENARIOS_EXTRACTION OR DESTROY IN CASE ----------
scenarios_extraction_path = create_scenarios_extraction(model_path) #------check

# STEP 3: Load intrinsic paths
def ensure_single_file(path_dir, list_file_paths):
    if len(list_file_paths) == 0:
        raise FileNotFoundError(f"No .dat file found in {path_dir}")
    elif len(list_file_paths) > 1:
        raise RuntimeError(f"Multiple .dat files found in {path_dir}: {[f.name for f in list_file_paths]}")
    # If exactly one, use it
    return list_file_paths[0].resolve()
    

# ------3.1: path for simstrat initial conditions-----------------
# Define the directory containing the expected .dat file
sim_initcond_dir = Path(model_path) / "scenarios" / "Initcond"
# Search for all .dat files
sim_dat_files = list(sim_initcond_dir.glob("*.dat"))
simstrat_initcond_file = ensure_single_file(sim_initcond_dir, sim_dat_files)

# ------3.2: path for simstrat inflows data-----------------
simstrat_inflows_dir = Path(model_path) / "scenarios" / "Inflow"

# ------3.3: path for simstrat config file -----------------
sim_config_dir = Path(model_path) / "config_files" 

# Delete any existing extraction config file first
new_config_name = "simstrat_config_steady_ch4inflow_with_extraction.par"
delete_existing_extraction_parfile(sim_config_dir, new_config_name) 

# Now Search for all .dat files to ensure only one exist
sim_config_files = list(sim_config_dir.glob("*.par"))
simstrat_config_file = ensure_single_file(sim_config_dir, sim_config_files) 


new_simstrat_config_file = copy_original_parfile(simstrat_config_file, new_config_name) 

# ----- WRITE SIMSTRAT INFLOWS -----------
write_simstrat_inflows_file(simstrat_inflows_dir, json_file_path, simstrat_initcond_file, new_simstrat_config_file, scenarios_extraction_path)


# ------3.4: path for aed2 inflows and IC data-----------------
aed2_inflows_dir = Path(model_path) / "scenarios" / "AED2_inflow_ch4inflow"
aed2_initcond_dir = Path(model_path) / "scenarios" / "AED2_initcond"

# ----- WRITE AED2 INFLOWS -----------
write_aed_inflows_file(aed2_inflows_dir, json_file_path, aed2_initcond_dir, new_simstrat_config_file, scenarios_extraction_path)


# ----- UPDATE CONFIG FILE -----------
update_simstart_file_par(new_simstrat_config_file, json_file_path, simstrat_inflows_dir)


# ----- RUN THE MODEL/ SIMULATION --------------
script_path = model_path / "run_simstrat_aed2_with_extraction.sh"

# This works on Linux and WSL (Windows Subsystem for Linux)
subprocess.run(["bash", str(script_path)], cwd=model_path, check=True)