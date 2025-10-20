#--------This script contains all necessary functions to process inflows for Kivu-Simstrat model V.1.1 to simulate improved methane extraction operations---------
import json
import numpy as np
from datetime import datetime, date
import pandas as pd
from pathlib import Path

#=========== Fuction 1: validate first hand data consistence ====================================
def check_inconsistence_data_length(config, section, required_keys=None):
    """
    Validate that the required keys exist in the config[section] and have the same length.
    """
    if section == "EXTRACTION":
        if required_keys is None:
            required_keys = [
                "power_production_MW",
                "extraction_depth_m",
                "extraction_range_m",
                "extraction_water_m3s",
                "ch4_extraction_efficiency_percent"
            ]
    
    if section == "REINJECTION":
        if required_keys is None:
            required_keys = [
                "reinjection_depth_m",
                "reinjection_water_m3s",
                "ch4_reinjection_percent",
                "co2_reinjection_percent"
            ]

    if section == "WASHING_EXTRACTION":
        if required_keys is None:
            required_keys = [
                "washing_extraction_depth_m",
                "washing_extraction_range_m",
                "washing_extraction_water_m3s"
            ]

    if section == "WASHING_REINJECTION":
        if required_keys is None:
            required_keys = [
                "washing_reinjection_depth_m",
                "washing_reinjection_water_m3s",
                "ch4_washing_reinjection_percent",
                "co2_washing_reinjection_percent"
            ]

    if section == "EXTRACTION_PERIOD":
        if required_keys is None:
            required_keys = [
                "extraction_start_date",
                "extraction_end_date"
            ]
    
    if section not in config:
        raise KeyError(f"Section '{section}' not found in configuration.")

    data = config[section]
    
    # Check all keys
    for key in required_keys:
        # do they all exist
        if key not in data:
            raise KeyError(f"Missing required key: '{key}' in section '{section}'")
        
        # Ensure values for range keys are strictly > 0
        if key in ("extraction_range_m", "washing_extraction_range_m"):
            if not all(val > 0 for val in data[key]):
                raise ValueError(f"All values for '{key}' in section '{section}' must be greater than 0.")
        
        # Ensure values for range keys are strictly aligning with max-depth of the lake
        if key in ("extraction_depth_m", "reinjection_depth_m", "washing_extraction_depth_m", "washing_reinjection_depth_m"):
            if not all(val > 0 and val <485 for val in data[key]):
                raise ValueError(f"All values for '{key}' in section '{section}' must be in the range of lake depths.")
        
        # Ensure user date periods are validated
        if section == "EXTRACTION_PERIOD" and key in ("extraction_start_date", "extraction_end_date"):
            start_dates = data.get("extraction_start_date", [])
            end_dates = data.get("extraction_end_date", [])

            # Parse all date strings first
            parsed_starts = []
            parsed_ends = []

            for i, (start, end) in enumerate(zip(start_dates, end_dates)):
                # Conditioning for the required date format, otherwise raise an error
                try:
                    start_dt = datetime.strptime(start, "%d-%m-%Y")
                    end_dt = datetime.strptime(end, "%d-%m-%Y")
                except ValueError:
                    raise ValueError(
                        f"Invalid date format at index {i+1} in section '{section}'. "
                        f"Expected 'dd-mm-yyyy'. Got start='{start}', end='{end}'"
                    )
                # Always start date should be less than end date
                if end_dt <= start_dt:
                    raise ValueError(
                        f"End date must be greater than start date at index {i+1} in section '{section}'. "
                        f"Got start='{start}', end='{end}'"
                    )
                
                # append all validated dates in the lists for subsquent validation
                parsed_starts.append(start_dt)
                parsed_ends.append(end_dt) 

            # Now do continuity check (your original logic but with datetime objects)
            min_date = parsed_starts[0]
            max_date = parsed_ends[0]

            for i in range(len(parsed_starts) - 1):
                if not (min_date <= parsed_starts[i+1] <= max_date):
                    raise ValueError(
                        f"Continuity in the given dates is required, there should be no gap between {max_date} and {start_dates[i+1]}."
                    )
                if parsed_ends[i+1] > max_date:
                    max_date = parsed_ends[i+1]  

    # Check all keys have same length
    lengths = [len(data[key]) for key in required_keys]
    if len(set(lengths)) != 1:
        raise ValueError(f"Inconsistent lengths among keys in '{section}': {dict(zip(required_keys, lengths))}")
    section_length = lengths[0]
    return section_length

#=========== Fuction 2: validate extraction user input data =====================================
def validate_data_inputs(config):
    # Required data sections
    sections = ["EXTRACTION", "REINJECTION", "WASHING_EXTRACTION", "WASHING_REINJECTION", "EXTRACTION_PERIOD"]
    #----data length validation across sections
    data_lengths = [check_inconsistence_data_length(config, section=section) for section in sections]
    
    #---- revalidate the validated length data
    if len(set(data_lengths)) != 1:
        raise ValueError("Inconsistent lengths among the data sections")
    else:
        valid_data = [
            config[section] for section in sections
            ]
        n_extractions = data_lengths[0]
        #print("Data Inputs Validation Done!")
        
    return valid_data, n_extractions

#=========== Fuction 3: load the validated json file data ========================================
def load_inputs_data(json_file_path):

    with open(f"{json_file_path}", "r") as f:
        config = json.load(f)
    loaded_data, n_extractions = validate_data_inputs(config)
    
    extraction, reinjection, washing_ext, washing_rei, ext_period = loaded_data

    return extraction, reinjection, washing_ext, washing_rei, ext_period, n_extractions

#=========== Fuction 4: to process the extraction depths ====================================
def process_extraction_depths(json_file_path): 
    
    extraction = load_inputs_data(json_file_path)[0]
    washing_ext = load_inputs_data(json_file_path)[2]
    n_extractions = load_inputs_data(json_file_path)[5]

    # lists to save all extracted data
    ext_depths = np.zeros((2,4*n_extractions))

    # Extract individual values (first item of each list)
    for i in range(n_extractions):

        #------ Extraction depths -----------------------
        ext_depth = extraction.get("extraction_depth_m", [None])[i]
        ext_range = extraction.get("extraction_range_m", [None])[i]
        d_ext_range = ext_range / 2
        ext_depth_lst = [-ext_depth-d_ext_range, -ext_depth-d_ext_range, -ext_depth+d_ext_range, -ext_depth+d_ext_range]

        #----- Washing extraction depths -----------------
        wash_ext_depth = washing_ext.get("washing_extraction_depth_m", [None])[i]
        wash_ext_range = washing_ext.get("washing_extraction_range_m", [None])[i]
        d_wash_range = wash_ext_range / 2
        wash_ext_depth_lst = [-wash_ext_depth-d_wash_range, -wash_ext_depth-d_wash_range, -wash_ext_depth+d_wash_range, -wash_ext_depth+d_wash_range]
        #####--loop out the extraction depth format(eg, at 450: -451, -451, -449, -449)
        for j in range(4):
            idx = 4 * i + j
            ext_depths[0][idx] = ext_depth_lst[j]
            ext_depths[1][idx] = wash_ext_depth_lst[j]

    return ext_depths

#=========== Fuction 5: to process the extraction discharges ====================================
def process_extraction_discharges(json_file_path): 
    
    extraction = load_inputs_data(json_file_path)[0]
    washing_ext = load_inputs_data(json_file_path)[2]
    n_extractions = load_inputs_data(json_file_path)[5]

    # lists to save all extracted data
    ext_discharges = np.zeros((2,4*n_extractions))

    # Extract individual values (first item of each list)
    for i in range(n_extractions):

        #------ Extraction depths -----------------------
        ext_discharge = extraction.get("extraction_water_m3s", [None])[i]
        ext_discharge_lst = [0, -ext_discharge/2, -ext_discharge/2, 0]

        #----- Washing extraction depths -----------------
        wash_ext_discharge = washing_ext.get("washing_extraction_water_m3s", [None])[i]
        wash_ext_discharge_lst = [0, -wash_ext_discharge/2, -wash_ext_discharge/2, 0]
        #####--loop out the extraction depth format(eg, at 450: -451, -451, -449, -449)
        for j in range(4):
            idx = 4 * i + j
            ext_discharges[0][idx] = ext_discharge_lst[j]
            ext_discharges[1][idx] = wash_ext_discharge_lst[j]

    return ext_discharges

#=========== Fuction 6: to process the reinjection depths ====================================
def process_reinjection_depths(json_file_path):
    
    reinjection = load_inputs_data(json_file_path)[1]
    washing_rei = load_inputs_data(json_file_path)[3]
    n_extractions = load_inputs_data(json_file_path)[5]

    # declare matrix to save data
    rei_depths = np.zeros((2,n_extractions))

    # Extract individual values
    for i in range(n_extractions):

        #----- Reinjection depths -----------------------
        rei_depth = reinjection.get("reinjection_depth_m", [None])[i]
        rei_depths[0][i] = -rei_depth

        #----- Washing reinjection depths -----------------------
        wash_rei_depth = washing_rei.get("washing_reinjection_depth_m", [None])[i]
        rei_depths[1][i] = -wash_rei_depth

    return rei_depths

#=========== Fuction 7: to process the reinjection discharges ====================================
def process_reinjection_discharges(json_file_path):
    
    reinjection = load_inputs_data(json_file_path)[1]
    washing_rei = load_inputs_data(json_file_path)[3]
    n_extractions = load_inputs_data(json_file_path)[5]

    # declare matrix to save data
    rei_discharges = np.zeros((2,n_extractions))

    # Extract individual values
    for i in range(n_extractions):

        #----- Reinjection discharges -----------------------
        rei_discharges[0][i] = reinjection.get("reinjection_water_m3s", [None])[i]

        #----- Washing reinjection discharges -----------------------
        rei_discharges[1][i] = washing_rei.get("washing_reinjection_water_m3s", [None])[i]

    return rei_discharges

#=========== Fuction 8: to process the extraction period dates ================================
def process_extraction_periods(json_file_path):
    
    ext_periods = load_inputs_data(json_file_path)[4]
    n_extractions = load_inputs_data(json_file_path)[5]

    # declare list to save data
    start_periods = list()
    end_periods = list()

    # Extract individual values
    for i in range(n_extractions):

        #----- Extraction periods -----------------------
        ext_period_start = ext_periods.get("extraction_start_date", [None])[i]
        ext_period_end = ext_periods.get("extraction_end_date", [None])[i]
        start_periods.append(ext_period_start)
        end_periods.append(ext_period_end)

    return start_periods, end_periods

#=========== Fuction 9: to process the reinjection efficiency fraction ====================================
def process_reinjection_efficieny(json_file_path):
    
    reinjection = load_inputs_data(json_file_path)[1]
    washing_rei = load_inputs_data(json_file_path)[3]
    n_extractions = load_inputs_data(json_file_path)[5]

    # declare matrix to save data
    ch4_rei_efficiencies = np.zeros((2,n_extractions))
    co2_rei_efficiencies = np.zeros((2,n_extractions))

    # Extract individual values
    for i in range(n_extractions):

        #----- Reinjection eff -----------------------
        ch4_rei_eff = reinjection.get("ch4_reinjection_percent", [None])[i]
        ch4_rei_efficiencies[0][i] = ch4_rei_eff / 100 
        co2_rei_eff = reinjection.get("co2_reinjection_percent", [None])[i]
        co2_rei_efficiencies[0][i] = co2_rei_eff / 100

        #----- Washing reinjection eff -----------------------
        ch4_wash_rei_eff = washing_rei.get("ch4_washing_reinjection_percent", [None])[i]
        ch4_rei_efficiencies[1][i] = ch4_wash_rei_eff / 100
        co2_wash_rei_eff = washing_rei.get("co2_washing_reinjection_percent", [None])[i]
        co2_rei_efficiencies[1][i] = co2_wash_rei_eff / 100

    return ch4_rei_efficiencies, co2_rei_efficiencies

#====== Function 9: load the model path ======================
def load_model_path(json_file_path):

    with open(f"{json_file_path}", 'r') as f:
        config = json.load(f)
    print(f"<<<< Loaded user inputs from: {json_file_path} >>>>")

    model_path = Path(config["SIMULATION_MODEL"]["kivu_simstrat_path"])
    #--- validate the path--------------
    if not model_path.exists():
        raise FileNotFoundError(f"Model path does not exist: {model_path}")
    return model_path