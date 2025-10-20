#--------This script contains all necessary functions to process inflows for Kivu-Simstrat model V.1.1 to simulate improved methane extraction operations---------
import json
import numpy as np
from datetime import datetime, date
import pandas as pd
from pathlib import Path
import load_user_inputs
import inflow_processor
import kivu_simstrat_processor

#=========== Fuction 10: to interpolate simstrat inflows values ================================
def interpolate_simstrat_initial_condition(simstrat_ic_file_path, intrp_z): # intrp_z has to be nedative
    # Load Simstrat initial conditions data
    simstrat_ic_df = pd.read_csv(simstrat_ic_file_path, sep='\t')
    #depths = simstrat_ic_df["Depth [m]"].values
    depths = simstrat_ic_df.iloc[:, 0].values
    #temperature = simstrat_ic_df["T [°C]"].values
    temperature = simstrat_ic_df.iloc[:, 3].values
    #salinity = simstrat_ic_df["S [permil]"].values
    salinity = simstrat_ic_df.iloc[:, 4].values

    # Ensure depths are sorted increasing (required for np.interp)
    depths_sorted = depths[::-1]
    temp_sorted = temperature[::-1]
    sal_sorted = salinity[::-1]

    # Interpolate (or extrapolate) using np.interp default behavior
    ic_intrp_temp_value = np.interp(intrp_z, depths_sorted, temp_sorted)
    ic_intrp_sal_value = np.interp(intrp_z, depths_sorted, sal_sorted)

    return ic_intrp_temp_value, ic_intrp_sal_value

#=========== Fuction 11: to interpolate aed inflows values ================================
def interpolate_aed_initial_condition(aed_ic_file_path, intrp_z): # intrp_z has to be nedative

    aed_ic = pd.read_csv(aed_ic_file_path, sep='\t', header=None, skiprows=1) # skip header otherwise themay not tab separated
    if pd.isna(aed_ic.iloc[0, 1]): # some times first depth value (zero) is way and considered as NAN
        aed_ic.iloc[0, 1] = 0.0

    depths = aed_ic.iloc[:, 0].values
    par_ic = aed_ic.iloc[:, 1].values


    # Ensure depths are sorted increasing (required for np.interp)
    depths_sorted = depths[::-1]
    par_sorted = par_ic[::-1]

    # Interpolate (or extrapolate) using np.interp default behavior
    ic_intrp_par_value = np.interp(intrp_z, depths_sorted, par_sorted)

    return ic_intrp_par_value


#=========== Fuction 13: formulate simstrat inflows: Tin and Sin values ================================
def formulate_given_simstrat_inflows(simstrat_ic_file_path, ext_and_wash_z, rei_and_wash_z):

    #zeros matrix for interpolated values: temperature
    ext_and_wash_temp = np.zeros_like(ext_and_wash_z)
    rei_and_wash_temp = np.zeros_like(rei_and_wash_z)

    #zeros matrix for interpolated values: salinity
    ext_and_wash_sal = np.zeros_like(ext_and_wash_z)
    rei_and_wash_sal = np.zeros_like(rei_and_wash_z)

    #corresponding initial conditions
    #reinjection depths
    #for j in range(len(rei_and_wash_z)):
    #for k in range (len(rei_and_wash_z[0])):
        #reinjection values
        #rei_and_wash_temp[0][k], rei_and_wash_sal[0][k] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, rei_and_wash_z[0][k])
        # washing reinjection values
        #rei_and_wash_temp[1][k], rei_and_wash_sal[1][k] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, rei_and_wash_z[1][k])

    #extraction depths
    #for n in range(len(ext_and_wash_z)):
    for m in range(0, len(ext_and_wash_z[0])-3, 4):
        #extraction values: only fill in the two in the middle out of evry 4 chunk
        ext_and_wash_temp[0][m+1], ext_and_wash_sal[0][m+1] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, ext_and_wash_z[0][m+1])
        ext_and_wash_temp[0][m+2], ext_and_wash_sal[0][m+2] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, ext_and_wash_z[0][m+2])
        # washing extraction values:  only fill in the two in the middle out of evry 4 chunk
        ext_and_wash_temp[1][m+1], ext_and_wash_sal[1][m+1] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, ext_and_wash_z[1][m+1])
        ext_and_wash_temp[1][m+2], ext_and_wash_sal[1][m+2] = interpolate_simstrat_initial_condition(simstrat_ic_file_path, ext_and_wash_z[1][m+2])
        
        #--for simstrat parameters reinjections remain the same as averaged extraction values
        k = m // 4 # interger division works perfect from m

        #reinjection averaged values
        rei_and_wash_temp[0][k] = (ext_and_wash_temp[0][m+1] + ext_and_wash_temp[0][m+2]) / 2 
        rei_and_wash_sal[0][k] = (ext_and_wash_sal[0][m+1] + ext_and_wash_sal[0][m+2]) / 2 

        #washing reinjection averaged values
        rei_and_wash_temp[1][k] = (ext_and_wash_temp[1][m+1] + ext_and_wash_temp[1][m+2]) / 2 
        rei_and_wash_sal[1][k] = (ext_and_wash_sal[1][m+1] + ext_and_wash_sal[1][m+2]) / 2 

    return ext_and_wash_temp, rei_and_wash_temp, ext_and_wash_sal, rei_and_wash_sal

#=========== Fuction 14: formulate aed inflows: each by each ================================
def formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z):

    #zeros matrix for interpolated values: temperature
    ext_and_wash_aed_par = np.zeros_like(ext_and_wash_z)
    rei_and_wash_aed_par = np.zeros_like(rei_and_wash_z)

    #corresponding initial conditions
    #reinjection depths
    #for j in range(len(rei_and_wash_z)):
        #for k in range (len(rei_and_wash_z[0])):
            #reinjection values
            #rei_and_wash_aed_par[0][k] = interpolate_aed_initial_condition(aed_ic_file_path, rei_and_wash_z[0][k])
            # washing reinjection values
            #rei_and_wash_aed_par[1][k] = interpolate_aed_initial_condition(aed_ic_file_path, rei_and_wash_z[1][k])

    #extraction depths
    #for n in range(len(ext_and_wash_z)):
    for m in range(0, len(ext_and_wash_z[0])-3, 4):
        #extraction values: only fill in the two in the middle out of evry 4 chunk
        ext_and_wash_aed_par[0][m+1] = interpolate_aed_initial_condition(aed_ic_file_path, ext_and_wash_z[0][m+1])
        ext_and_wash_aed_par[0][m+2] = interpolate_aed_initial_condition(aed_ic_file_path, ext_and_wash_z[0][m+2])
        # washing extraction values:  only fill in the two in the middle out of evry 4 chunk
        ext_and_wash_aed_par[1][m+1] = interpolate_aed_initial_condition(aed_ic_file_path, ext_and_wash_z[1][m+1])
        ext_and_wash_aed_par[1][m+2] = interpolate_aed_initial_condition(aed_ic_file_path, ext_and_wash_z[1][m+2])

        #--for simstrat parameters reinjections remain the same as averaged extraction values
        k = m // 4 # interger division works perfect from m

        #reinjection averaged values
        rei_and_wash_aed_par[0][k] = (ext_and_wash_aed_par[0][m+1] + ext_and_wash_aed_par[0][m+2]) / 2
        # washing reinjection averaged values
        rei_and_wash_aed_par[1][k] = (ext_and_wash_aed_par[1][m+1] + ext_and_wash_aed_par[1][m+2]) / 2

    return ext_and_wash_aed_par, rei_and_wash_aed_par

#=========== Fuction 20: write simstrat state and constant inflow files and save them ================================
def write_simstrat_inflows_file(data_dir, json_file_path, simstrat_ic_file_path, simstrat_config_path, scenarios_extraction_path):
    # Get all .dat files in the directory
    dat_files = sorted(Path(data_dir).glob("*.dat"))
    ext_and_wash_z = load_user_inputs.process_extraction_depths(json_file_path)
    rei_and_wash_z = load_user_inputs.process_reinjection_depths(json_file_path)
    ext_and_wash_temp, rei_and_wash_temp, ext_and_wash_sal, rei_and_wash_sal = formulate_given_simstrat_inflows(simstrat_ic_file_path, ext_and_wash_z, rei_and_wash_z) # This will be called accordingly for Simstrat and Aed2
    #print(ext_and_wash_temp)
    
    # Iterate over each file
    for inflow_file_path in dat_files:

        # For temperature inflows
        if inflow_file_path.name == "Tin.dat":
            # original inflows data
            state_inflow = True # identify wehter the inflow is state or constant
            inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
            Tin_depths, Tin_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_temp, rei_and_wash_temp, n_deep_z, n_surface_z, depths, inflows_matrix)
            #new_file_path = inflow_file_path.with_name("Tin_new_1.dat")
            new_file_path = Path(scenarios_extraction_path) / "Inflow" / inflow_file_path.name
            extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
            inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(Tin_depths)-n_surface_z, n_surface_z, Tin_depths, Tin_values, extraction_start_dates, state_inflow)
        
        # For salinity inflows
        if inflow_file_path.name == "Sin.dat":
            # original inflows data
            state_inflow = True # identify wehter the inflow is state or constant
            inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
            Sin_depths, Sin_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_sal, rei_and_wash_sal, n_deep_z, n_surface_z, depths, inflows_matrix)
            new_file_path = Path(scenarios_extraction_path) / "Inflow" / inflow_file_path.name
            extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
            inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(Sin_depths)-n_surface_z, n_surface_z, Sin_depths, Sin_values, extraction_start_dates, state_inflow)

        # For salinity inflows
        if inflow_file_path.name == "Qin.dat":
            # original inflows data
            state_inflow = False # identify wehter the inflow is state or constant
            inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
            ext_and_wash_qins = np.zeros_like(load_user_inputs.process_extraction_discharges(json_file_path)) # zeros because its q_out
            rei_and_wash_qins = load_user_inputs.process_reinjection_discharges(json_file_path)
            extraction_start_dates, extraction_end_dates = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for constant inflow get all the extraction starting dates
            Qin_depths, Qin_values, iterated_days = inflow_processor.get_constant_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_qins, rei_and_wash_qins, n_deep_z, n_surface_z, depths, inflows_matrix, extraction_start_dates, extraction_end_dates)
            new_file_path = Path(scenarios_extraction_path) / "Inflow" / inflow_file_path.name
            inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(Qin_depths)-n_surface_z, n_surface_z, Qin_depths, Qin_values, iterated_days, state_inflow)

        # For salinity inflows
        if inflow_file_path.name == "Qout.dat":
            # original inflows data
            state_inflow = False # identify wehter the inflow is state or constant
            inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
            ext_and_wash_qouts = load_user_inputs.process_extraction_discharges(json_file_path)
            rei_and_wash_qouts = np.zeros_like(load_user_inputs.process_reinjection_discharges(json_file_path)) # zeros because its q_in
            extraction_start_dates, extraction_end_dates = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for constant inflow get all the extraction starting dates
            Qout_depths, Qout_values, iterated_days = inflow_processor.get_constant_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_qouts, rei_and_wash_qouts, n_deep_z, n_surface_z, depths, inflows_matrix, extraction_start_dates, extraction_end_dates)
            new_file_path = Path(scenarios_extraction_path) / "Inflow" / inflow_file_path.name
            inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(Qout_depths)-n_surface_z, n_surface_z, Qout_depths, Qout_values, iterated_days, state_inflow)




#=========== Fuction 21: write aed state and constant inflow files and save them ================================
def write_aed_inflows_file(aed_inflow_dir, json_file_path, aed_ic_files_dir, simstrat_config_path, scenarios_extraction_path):
    # Get all .dat files in the directory
    inflow_dat_files = sorted(Path(aed_inflow_dir).glob("*.dat"))
    ic_dat_files = sorted(Path(aed_ic_files_dir).glob("*.dat"))
    ext_and_wash_z = load_user_inputs.process_extraction_depths(json_file_path)
    rei_and_wash_z = load_user_inputs.process_reinjection_depths(json_file_path)
    
    # Iterate over each inflow file
    for inflow_file_path in inflow_dat_files:

        # Iterate over each aed initial conditions files
        for aed_ic_file_path in ic_dat_files:

            # For CAR_ch4_bub inflows
            if (inflow_file_path.name == "CAR_ch4_bub_inflow.dat") and (aed_ic_file_path.name == "CAR_ch4_bub_ini.dat"):
                ext_and_wash_ch4_bub, rei_and_wash_ch4_bub = formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z) # for only Aed2
                # original inflows data
                state_inflow = True # identify wehter the inflow is state or constant
                inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
                ch4_bub_in_depths, ch4_bub_in_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z, ext_and_wash_ch4_bub, rei_and_wash_ch4_bub, n_deep_z, n_surface_z, depths, inflows_matrix)
                #new_file_path = inflow_file_path.with_name("CAR_ch4_bub_inflow_new_1.dat")
                new_file_path = Path(scenarios_extraction_path) / "AED2_inflow_ch4inflow" / inflow_file_path.name
                extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
                inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(ch4_bub_in_depths)-n_surface_z, n_surface_z, ch4_bub_in_depths, ch4_bub_in_values, extraction_start_dates, state_inflow)

            # For CAR_ch4 inflows
            if (inflow_file_path.name == "CAR_ch4_inflow.dat") and (aed_ic_file_path.name == "CAR_ch4_ini.dat"):
                ext_and_wash_ch4, rei_and_wash_ch4 = formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z) # for only Aed2
                ch4_rei_eff, _ = load_user_inputs.process_reinjection_efficieny(json_file_path)
                rei_and_wash_ch4 = ch4_rei_eff * rei_and_wash_ch4 #--element-wise multiplication
                # original inflows data
                state_inflow = True # identify wehter the inflow is state or constant
                inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
                ch4_in_depths, ch4_in_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z, ext_and_wash_ch4, rei_and_wash_ch4, n_deep_z, n_surface_z, depths, inflows_matrix)
                new_file_path = Path(scenarios_extraction_path) / "AED2_inflow_ch4inflow" / inflow_file_path.name
                extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
                inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(ch4_in_depths)-n_surface_z, n_surface_z, ch4_in_depths, ch4_in_values, extraction_start_dates, state_inflow)

            # For CAR_dic inflows
            if (inflow_file_path.name == "CAR_dic_inflow.dat") and (aed_ic_file_path.name == "CAR_dic_ini.dat"):
                ext_and_wash_dic, rei_and_wash_dic = formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z) # for only Aed2
                dic_rei_eff, _ = load_user_inputs.process_reinjection_efficieny(json_file_path)
                rei_and_wash_dic = dic_rei_eff * rei_and_wash_dic #--element-wise multiplication
                # original inflows data
                state_inflow = True # identify wehter the inflow is state or constant
                inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
                dic_in_depths, dic_in_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z, ext_and_wash_dic, rei_and_wash_dic, n_deep_z, n_surface_z, depths, inflows_matrix)
                new_file_path = Path(scenarios_extraction_path) / "AED2_inflow_ch4inflow" / inflow_file_path.name
                extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
                inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(dic_in_depths)-n_surface_z, n_surface_z, dic_in_depths, dic_in_values, extraction_start_dates, state_inflow)

            # For CAR_ph inflows
            if (inflow_file_path.name == "CAR_pH_inflow.dat") and (aed_ic_file_path.name == "CAR_pH_ini.dat"): #----------naming concern --------------
                ext_and_wash_ph, rei_and_wash_ph = formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z) # for only Aed2
                # original inflows data
                state_inflow = True # identify wehter the inflow is state or constant
                inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
                ph_in_depths, ph_in_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z, ext_and_wash_ph, rei_and_wash_ph, n_deep_z, n_surface_z, depths, inflows_matrix)
                new_file_path = Path(scenarios_extraction_path) / "AED2_inflow_ch4inflow" / inflow_file_path.name
                extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
                inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(ph_in_depths)-n_surface_z, n_surface_z, ph_in_depths, ph_in_values, extraction_start_dates, state_inflow,)

            # For OXY_oxy inflows
            if (inflow_file_path.name == "OXY_oxy_inflow.dat") and (aed_ic_file_path.name == "OXY_oxy_ini.dat"):
                ext_and_wash_oxy, rei_and_wash_oxy = formulate_given_aed_inflows(aed_ic_file_path, ext_and_wash_z, rei_and_wash_z) # for only Aed2
                # original inflows data
                state_inflow = True # identify wehter the inflow is state or constant
                inflow_header, n_deep_z, n_surface_z, depths, inflows_matrix = inflow_processor.read_inflow_file(inflow_file_path)
                oxy_in_depths, oxy_in_values = inflow_processor.get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z, ext_and_wash_oxy, rei_and_wash_oxy, n_deep_z, n_surface_z, depths, inflows_matrix)
                new_file_path = Path(scenarios_extraction_path) / "AED2_inflow_ch4inflow" / inflow_file_path.name
                extraction_start_dates, _ = kivu_simstrat_processor.convert_date_to_days(json_file_path, simstrat_config_path) # for state inflow get only the first extraction starting date
                inflow_processor.write_and_save_inflow_file(new_file_path, inflow_header, len(oxy_in_depths)-n_surface_z, n_surface_z, oxy_in_depths, oxy_in_values, extraction_start_dates, state_inflow)