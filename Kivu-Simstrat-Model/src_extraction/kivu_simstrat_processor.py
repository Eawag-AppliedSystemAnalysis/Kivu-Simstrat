#--------This script contains all necessary functions to process inflows for Kivu-Simstrat model V.1.1 to simulate improved methane extraction operations---------
#=======<< General packages >>======
from datetime import datetime, date
from pathlib import Path
import os
import shutil
import re

#=======<< Import modules >>======
import load_user_inputs
import inflow_processor

#------ Function: main Function to create new scenarios with extraction directory -----------
def create_scenarios_extraction(model_dir):
    """
    Create a lightweight 'scenarios_extraction' directory structure inside the model directory.
    - No large files are copied (only folders and small configuration files if needed).
    - 'Inflow' and 'Aed_Inflow' subdirectories are created empty.
    - Any 'Results_*' folders in the source 'scenarios' are deleted before processing.
    - If 'scenarios_extraction' exists, it is deleted before recreating.
    
    Parameters
    ----------
    model_dir : str or Path
        Path to the model directory containing 'scenarios'.
    """

    model_path = Path(model_dir).resolve()
    scenarios_path = model_path / "scenarios"
    new_scenarios_path = model_path / "scenarios_extraction"

    # --- Safety check ---
    if not scenarios_path.exists() or not scenarios_path.is_dir():
        raise FileNotFoundError(f"Scenarios directory not found: {scenarios_path}")

    # --- Remove the old extraction folder if it exists ---
    if new_scenarios_path.exists():
        shutil.rmtree(new_scenarios_path)
        print(f"Removed existing directory: {new_scenarios_path}")

    # --- Recreate a lightweight structure ---
    new_scenarios_path.mkdir()

    for item in scenarios_path.iterdir():

        if item.is_dir():
            if item.name in ["Inflow", "AED2_inflow_ch4inflow"]: # the choice of inflows for aed2 can be used to allow aed_inflow
                # Create empty inflow directories
                (new_scenarios_path / item.name).mkdir()
            elif not item.name.startswith("Results_"):
                # Create symbolic link (not a full copy)
                shutil.copytree(item, new_scenarios_path / item.name, dirs_exist_ok=True)
    
        elif item.is_file():
            # Copy files
            shutil.copy2(item, new_scenarios_path / item.name)

    print("-----'scenarios_extraction' structure prepared successfully!------")
    return new_scenarios_path

#------ Function: function to delete any existing config file with extraction -----------
def delete_existing_extraction_parfile(simstrat_config_path, new_filename):
    """
    If the new file already exists, it is deleted first.
    """

    # Convert to Path object for easy manipulation
    config_path = Path(simstrat_config_path)
    new_path = config_path / new_filename

    # Check existence of source
    if not config_path.exists():
        raise FileNotFoundError(f"The source file does not exist: {config_path}!!")

    # Remove destination if it already exists
    if new_path.exists():
        os.remove(new_path)
        print(f"Existing file '{new_path.name}' deleted.")

#------ Function: function to copy an original config fileto new one with extraction -----------
def copy_original_parfile(simstrat_config_path, new_filename):
    """
    Copies an existing .par file to the same directory with a new name.
    If the new file already exists, it is deleted first.
    """

    # Convert to Path object for easy manipulation
    config_path = Path(simstrat_config_path)
    new_path = config_path.with_name(new_filename)

    # Check existence of source
    if not config_path.exists():
        raise FileNotFoundError(f"The source file does not exist: {config_path}!!")

    # Copy file
    shutil.copy2(config_path, new_path)
    print(f"File copied: '{config_path.name}' to '{new_path.name}'!!")

    return str(new_path)


#------ 24. main Function to run and update simstrat configuration file -----------
def update_simstart_file_par(new_simstrat_config_path, json_file_path, inflow_file_path):
    update_scenario_paths(new_simstrat_config_path)
    inflow_file_path = sorted(Path(inflow_file_path).glob("*.dat"))
    update_simulation_days(json_file_path, new_simstrat_config_path)
    add_extraction_inputs_config(new_simstrat_config_path, json_file_path, inflow_file_path[0]) 

#------ 24. main Function to update scenarios to scenarios_extraction in path provided in config file -----------
def update_scenario_paths(simstrat_config_path): #==== Should be adapted for windows config as well ===========????????????????
    """
    Updates all file paths containing 'scenarios//' to 'scenarios_extraction//'
    inside the Simstrat configuration file.
    """

    # Load the file as text because the file might include comments (invalid JSON)
    with open(simstrat_config_path, 'r') as f:
        content = f.read()

    # Replace all occurrences robustly, preserving path structure
    updated_content = re.sub(
        r'(?<!\w)scenarios//',  # ensures we match exact 'scenarios//' not e.g. 'subscenarios'
        'scenarios_extraction//',
        content
    )

    # Write the modified file back
    with open(simstrat_config_path, 'w') as f:
        f.write(updated_content)

#=========== Fuction 18: convert datesto referenced year days ================================
def convert_date_to_days(json_file_path, simstrat_config_path):
    referenced_start_dates = list()
    referenced_end_dates = list()

    start_ext_periods, end_ext_periods = load_user_inputs.process_extraction_periods(json_file_path)

    # Extract reference year from par file
    reference_year = None
    with open(simstrat_config_path, 'r') as file:
        for line in file:
            if "Reference year" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    value = parts[1].strip().rstrip(',').strip()
                    reference_year = int(value)
                    break  # Stop once we find it

    if reference_year is None:
        raise ValueError("Reference year not found in the provided .par file.")
    
    def reference_date(reference_year, date_j):
        # convert it to date
        date_obj = datetime.strptime(date_j, "%d-%m-%Y")
        # Define the start and end dates (year, month, day)
        start_date = date(reference_year, 1, 1)
        end_date = date(date_obj.year, date_obj.month, date_obj.day)
        # Calculate the difference in days between the two dates
        delta_date = end_date - start_date
        return delta_date.days
    
    if reference_year is not None:
        for i in range(len(start_ext_periods)):
            #-- for start extraction dates
            referenced_start_dates.append(reference_date(reference_year, start_ext_periods[i]))
            #-- for end extraction dates
            referenced_end_dates.append(reference_date(reference_year, end_ext_periods[i]))

    return referenced_start_dates, referenced_end_dates

#=========== Fuction 22: update the starting and ending days of simulations ================================
def update_simulation_days(json_file_path, simstrat_config_path):
    start_day = convert_date_to_days(json_file_path, simstrat_config_path)[0][0] # minimum date is the first elent in the list
    end_day = max(convert_date_to_days(json_file_path, simstrat_config_path)[1])
    updated_lines = []
    found = [False, False]  # Flags for "Start d" and "End d"

    with open(simstrat_config_path, 'r') as file:
        for line in file:
            if "Start d" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    #indentation = line[:line.index(parts[0])]  # preserve leading spaces
                    indentation = updated_lines[-1][:len(updated_lines[-1]) - len(updated_lines[-1].lstrip())]
                    key = parts[0].strip()
                    trailing_comma = ',' if line.strip().endswith(',') else ''
                    updated_line = f"{indentation}{key}                :{start_day}{trailing_comma}\n"
                    updated_lines.append(updated_line)
                    found[0] = True
                else:
                    updated_lines.append(line)

            elif "End d" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    #indentation = line[:line.index(parts[0])]  # preserve leading spaces
                    indentation = updated_lines[-1][:len(updated_lines[-1]) - len(updated_lines[-1].lstrip())]
                    key = parts[0].strip()
                    trailing_comma = ',' if line.strip().endswith(',') else ''
                    updated_line = f"{indentation}{key}                :{end_day}{trailing_comma}\n"
                    updated_lines.append(updated_line)
                    found[1] = True
                else:
                    updated_lines.append(line)

            else:
                updated_lines.append(line)

    if not all(found):
        missing = []
        if not found[0]:
            missing.append('"Start d"')
        if not found[1]:
            missing.append('"End d"')
        raise ValueError(f"Missing key(s) in configuration file: {', '.join(missing)}")

    with open(simstrat_config_path, 'w') as file:
        file.writelines(updated_lines)

#=========== Fuction 23: add n_extractions, extraction_depths and dates ================================
def add_extraction_inputs_config(simstrat_config_path, json_file_path, inflow_file_path):
    
    ext_and_wash_z = load_user_inputs.process_extraction_depths(json_file_path)
    rei_and_wash_z = load_user_inputs.process_reinjection_depths(json_file_path)
    _, n_deep_z, _, depths, _ = inflow_processor.read_inflow_file(inflow_file_path)
    _, ext_new_indices, wash_ext_new_indices, rei_new_indices, wash_rei_new_indices, _ = inflow_processor.combine_depths_and_track_indices(ext_and_wash_z, rei_and_wash_z, depths[0:n_deep_z])
    extraction_start_dates, extraction_end_dates = convert_date_to_days(json_file_path, simstrat_config_path)
    ch4_rei_percents, co2_rei_percents = load_user_inputs.process_reinjection_efficieny(json_file_path)
    
    ext_new_indices = [ext_index + 1 for ext_index in ext_new_indices]
    wash_ext_new_indices = [wash_ext_index + 1 for wash_ext_index in wash_ext_new_indices]
    rei_new_indices = [rei_index + 1 for rei_index in rei_new_indices] # convert from python to fortran indices format
    wash_rei_new_indices = [wash_rei_index + 1 for wash_rei_index in wash_rei_new_indices] # convert from python to fortran indices format
    
    updated_lines = []
    in_simulation_block = False
    last_key_index = None

    extraction_number_str = f'"Extraction number"                : {int(len(ext_new_indices)/4)}'
    extraction_depths_str = f'"Extraction depths"                : {ext_new_indices}'
    wash_extraction_depths_str = f'"Wash extraction depths"      : {wash_ext_new_indices}'
    reinjection_depths_str = f'"Reinjection depths"              : {rei_new_indices}'
    wash_reinjection_depths_str = f'"Wash reinjection depths"    : {wash_rei_new_indices}'
    co2_reinjection_percents_str = f'"DIC reinjection percents"    : {list(co2_rei_percents[0,:])}'
    co2_wash_reinjection_percents_str = f'"DIC wash reinjection percents"    : {list(co2_rei_percents[1,:])}'
    ch4_reinjection_percents_str = f'"CH4 reinjection percents"    : {list(ch4_rei_percents[0,:])}'
    ch4_wash_reinjection_percents_str = f'"CH4 wash reinjection percents"    : {list(ch4_rei_percents[1,:])}'
    extraction_start_dates_str = f'"Start extraction dates"                  : {extraction_start_dates}'
    extraction_end_dates_str = f'"End extraction dates"                  : {extraction_end_dates}'

    with open(simstrat_config_path, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        # Detect start of Simulation block
        if '"Simulation"' in line and '{' in line:
            in_simulation_block = True

        if in_simulation_block and '}' in line:
            #if '}' in line :  # End of block found
            # Ensure previous line ends with a comma
            if not lines[i - 1].rstrip().endswith(','): # THIS COMMA CANNOT BE ADDED---------------
                updated_lines[i - 1] = updated_lines[i - 1].rstrip() + ',\n'

            # Detect indentation from previous key line
            prev_indent = lines[i - 1][:len(lines[i - 1]) - len(lines[i - 1].lstrip())]

            # for extraction number
            updated_lines.append(f"{prev_indent}{extraction_number_str},\n")
            # for extraction depths
            updated_lines.append(f"{prev_indent}{extraction_depths_str},\n")
            updated_lines.append(f"{prev_indent}{wash_extraction_depths_str},\n")
            # for reinjection depths
            updated_lines.append(f"{prev_indent}{reinjection_depths_str},\n")
            updated_lines.append(f"{prev_indent}{wash_reinjection_depths_str},\n")
            # for co2 reinjection percents
            updated_lines.append(f"{prev_indent}{co2_reinjection_percents_str},\n")
            updated_lines.append(f"{prev_indent}{co2_wash_reinjection_percents_str},\n")     
            # for ch4 reinjection percents
            updated_lines.append(f"{prev_indent}{ch4_reinjection_percents_str},\n")
            updated_lines.append(f"{prev_indent}{ch4_wash_reinjection_percents_str},\n")        
            # for extraction start dates
            updated_lines.append(f"{prev_indent}{extraction_start_dates_str},\n")
            # for extraction end dates
            updated_lines.append(f"{prev_indent}{extraction_end_dates_str}\n")

            in_simulation_block = False  # Exit Simulation block

        updated_lines.append(line)

    with open(simstrat_config_path, 'w') as file:
        file.writelines(updated_lines)

