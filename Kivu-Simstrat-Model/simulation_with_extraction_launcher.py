import kivu_simstrat_processor
import write_inflows

import sys
from pathlib import Path
import subprocess


#STEP 1: Get Json file path -----------------------------------------

# Ensure a JSON file path was passed
if len(sys.argv) < 2:
    raise ValueError("Usage: python simulation_with_extraction_launcher.py <config.json>")

# Get the JSON file path passed from run_file.py
json_file_path = Path(sys.argv[1]).resolve()


# STEP 2: load the current location / model path
model_path = current_loc = Path(__file__).resolve().parent

#---- CREATE SCENARIOS_EXTRACTION OR DESTROY IN CASE ----------
scenarios_extraction_path = kivu_simstrat_processor.create_scenarios_extraction(model_path) #------check

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
kivu_simstrat_processor.delete_existing_extraction_parfile(sim_config_dir, new_config_name) 

# Now Search for all .dat files to ensure only one exist
sim_config_files = list(sim_config_dir.glob("*.par"))
simstrat_config_file = ensure_single_file(sim_config_dir, sim_config_files) 


new_simstrat_config_file = kivu_simstrat_processor.copy_original_parfile(simstrat_config_file, new_config_name) 

# ----- WRITE SIMSTRAT INFLOWS -----------
write_inflows.write_simstrat_inflows_file(simstrat_inflows_dir, json_file_path, simstrat_initcond_file, new_simstrat_config_file, scenarios_extraction_path)


# ------3.4: path for aed2 inflows and IC data-----------------
aed2_inflows_dir = Path(model_path) / "scenarios" / "AED2_inflows_ch4inflows"
aed2_initcond_dir = Path(model_path) / "scenarios" / "AED2_initcond"

# ----- WRITE AED2 INFLOWS -----------
write_inflows.write_aed_inflows_file(aed2_inflows_dir, json_file_path, aed2_initcond_dir, new_simstrat_config_file, scenarios_extraction_path)


# ----- UPDATE CONFIG FILE -----------
kivu_simstrat_processor.update_simstart_file_par(new_simstrat_config_file, json_file_path, simstrat_inflows_dir)


# ----- RUN THE MODEL/ SIMULATION --------------
script_path = model_path / "run_simstrat_aed2_with_extraction.sh"

# This works on Linux and WSL (Windows Subsystem for Linux)
subprocess.run(["bash", str(script_path)], cwd=model_path, check=True)