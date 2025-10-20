#--------This script contains all necessary functions to process inflows for Kivu-Simstrat model V.1.1 to simulate improved methane extraction operations---------
import json
import numpy as np
from datetime import datetime, date
import pandas as pd
from pathlib import Path

#=========== Fuction 12: read inflow file ================================
def read_inflow_file(inflow_file_path):
    with open(inflow_file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Extract inflow counts
    inflow_header = lines[0].strip()
    inflow_counts = list(map(int, lines[1].split()))
    n_deep_z, n_surface_z = inflow_counts[0], inflow_counts[1]
    n_total = n_deep_z + n_surface_z

    # Extract depth values
    depth_line = list(map(float, lines[2].split()))
    depths = depth_line[1:]  # actual depth entries

    assert len(depths) == n_total, "Mismatch in number of inflow depth entries."

    # Extract inflow values: each line after line 2 has 1 time value + n_total inflow values
    inflow_matrix = []
    for line in lines[3:]:
        row = list(map(float, line.split()))
        time = row[0]
        values = row[1:]
        assert len(values) == n_total, f"Inflow data does not match inflow depth count at time {time}."
        inflow_matrix.append(values) #without time data points

    return inflow_header, n_deep_z, n_surface_z, depths, inflow_matrix

#=========== Fuction 15: formulate simstrat inflows: Tin and Sin values ================================
def combine_depths_and_track_indices(ext_and_wash_z, rei_and_wash_z, original_deep_depths): #Confirmed

    # depths' lists with source tags
    ext_z_lst = [(val, 'ext_z', i) for i, val in enumerate(ext_and_wash_z[0])]
    wash_ext_z_lst = [(val, 'wash_ext_z', i) for i, val in enumerate(ext_and_wash_z[1])]
    rei_z_lst = [(val, 'rei_z', i) for i, val in enumerate(rei_and_wash_z[0])]
    wash_rei_z_lst = [(val, 'wash_rei_z', i) for i, val in enumerate(rei_and_wash_z[1])]
    original_z_lst = [(val, 'org_z', i) for i, val in enumerate(original_deep_depths)]

    # Step 1: Combine depth lists 
    combined_depths = ext_z_lst + wash_ext_z_lst + rei_z_lst + wash_rei_z_lst + original_z_lst

    # Step 2: Sort the combined list based on the values
    combined_depths_sorted = sorted(combined_depths, key=lambda x: x[0])

    # Step 3: Track new indices
    ext_new_indices = [None] * len(ext_z_lst)
    wash_ext_new_indices = [None] * len(wash_ext_z_lst)
    rei_new_indices = [None] * len(rei_z_lst)
    wash_rei_new_indices = [None] * len(wash_rei_z_lst)
    org_new_indices = [None] * len(original_z_lst)

    for new_index, (val, origin, original_index) in enumerate(combined_depths_sorted):
        if origin == 'ext_z':
            ext_new_indices[original_index] = new_index
        elif origin == 'wash_ext_z':
            wash_ext_new_indices[original_index] = new_index
        elif origin == 'rei_z':
            rei_new_indices[original_index] = new_index
        elif origin == 'wash_rei_z':
            wash_rei_new_indices[original_index] = new_index
        else:
            org_new_indices[original_index] = new_index

    # Output
    merged_sorted_depths = [item[0] for item in combined_depths_sorted]

    return merged_sorted_depths, ext_new_indices, wash_ext_new_indices, rei_new_indices, wash_rei_new_indices, org_new_indices

#=========== Fuction 16: formulate simstrat inflows: Tin and Sin values ================================
def get_state_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_par, rei_and_wash_par, old_n_deep_z, old_n_surface_z, old_depths, old_inflows_matrix):

    merged_sorted_depths, ext_new_indices, wash_ext_new_indices, rei_new_indices, wash_rei_new_indices, org_new_indices = combine_depths_and_track_indices(ext_and_wash_z, rei_and_wash_z, old_depths[0:old_n_deep_z])
    #combine inflow depths into one list and declare zores for their corresponding values
    new_deep_and_surface_depths = merged_sorted_depths + old_depths[old_n_deep_z:]

    #if state_inflow: # wherever +1 indicates the reserved position for date (which is set to zero here to be jumped later in write_inflow function)
    par_inflow_list = [0] * (len(new_deep_and_surface_depths)) # accomodate space for date +1

    for i in range(len(ext_and_wash_z[0])): # for extraction
        par_inflow_list[ext_new_indices[i]] = ext_and_wash_par[0][i] # extr
        par_inflow_list[wash_ext_new_indices[i]] = ext_and_wash_par[1][i] # wash extr

    for j in range(len(rei_and_wash_z[0])): # for reinjection
        par_inflow_list[rei_new_indices[j]] = rei_and_wash_par[0][j] # rei
        par_inflow_list[wash_rei_new_indices[j]] = rei_and_wash_par[1][j] # wash rei
    
    for k in range(old_n_deep_z):# for old deep par values
        par_inflow_list[org_new_indices[k]] = old_inflows_matrix[0][k] # old inflow values from the first row of the file data without (+1) date

    for l in range(old_n_surface_z): # for old surface par values
        par_inflow_list[len(merged_sorted_depths)+l] = old_inflows_matrix[0][old_n_deep_z+l] # old inflow values from the first row of the file data
    
    return new_deep_and_surface_depths, par_inflow_list

#=========== Fuction 17: formulate simstrat inflows: Tin and Sin values ================================
def get_constant_inflow_depths_and_values(ext_and_wash_z, rei_and_wash_z,ext_and_wash_par, rei_and_wash_par, old_n_deep_z, old_n_surface_z, old_depths, old_inflows_matrix, start_extraction_dates, end_extraction_dates):

    merged_sorted_depths, ext_new_indices, wash_ext_new_indices, rei_new_indices, wash_rei_new_indices, org_new_indices = combine_depths_and_track_indices(ext_and_wash_z, rei_and_wash_z, old_depths[0:old_n_deep_z])
    #combine inflow depths into one list and declare zores for their corresponding values
    new_deep_and_surface_depths = merged_sorted_depths + old_depths[old_n_deep_z:]

    #if state_inflow: # wherever +1 indicates the reserved position for date (which is set to zero here to be jumped later in write_inflow function)
    num_extractions = int(len(ext_and_wash_z[0])/4) # another way to know number of extractions
    min_date = min(start_extraction_dates)
    max_date = max(end_extraction_dates)
    if min_date == 0:
        par_inflow_array = np.zeros((max_date-min_date, len(new_deep_and_surface_depths)))
    else:
        par_inflow_array = np.zeros((max_date-min_date+1, len(new_deep_and_surface_depths)))
    
    iterated_days_lst = list()
    for date_iter in range(min_date, max_date+1): # at some cases this can be a time consuming loop, luckily the computations below are very simple then python can still handle the iterations
        for n_ext in range(num_extractions):
            if (date_iter >= start_extraction_dates[n_ext]) and (date_iter <= end_extraction_dates[n_ext]):
                for i in range(n_ext*4, 4*(n_ext+1)): # for extraction
                    #for i in range(len(ext_and_wash_z[0])): # for extraction
                    par_inflow_array[date_iter-min_date, ext_new_indices[i]] = ext_and_wash_par[0][i] # extr
                    par_inflow_array[date_iter-min_date, wash_ext_new_indices[i]] = ext_and_wash_par[1][i] # wash extr

                #for j in range(len(rei_and_wash_z[0])): # for reinjection --- No need for loop because n_extractions is equal to that
                par_inflow_array[date_iter-min_date, rei_new_indices[n_ext]] = rei_and_wash_par[0][n_ext] # rei
                par_inflow_array[date_iter-min_date, wash_rei_new_indices[n_ext]] = rei_and_wash_par[1][n_ext] # wash rei

                for k in range(old_n_deep_z):# for old deep par values
                    #--here we can add date condition from old inflow dates (if date_iter>=old_date_i and date_iter<= next_old_date_i use the given values in the list, else use only the line below to repeat them)
                    par_inflow_array[date_iter-min_date, org_new_indices[k]] = old_inflows_matrix[0][k] # old inflow values from the first row of the file data without (+1) date

                for m in range(old_n_surface_z): # for old surface par values
                    #-- the same here as well like deep old inflows
                    par_inflow_array[date_iter-min_date, len(merged_sorted_depths)+m] = old_inflows_matrix[0][old_n_deep_z+m] # old inflow values from the first row of the file data
        iterated_days_lst.append(date_iter)
 
    return new_deep_and_surface_depths, par_inflow_array, iterated_days_lst

#=========== Fuction 19: write inflow file and save it ================================
def write_and_save_inflow_file(file_path_dat, header, new_n_deep_z, new_n_surface_z, new_depths, new_inflow_data, extraction_dates, state_inflow=False):
    #print(new_inflow_data)
    depths = [-1] + new_depths
    with open(file_path_dat, "w") as f:
        # Write header
        f.write(header + "\n")
        
        # Write n_deep_z and n_surface_z
        f.write(f"{new_n_deep_z} {new_n_surface_z}\n")
        
        # Write depth line
        depth_line = " ".join(str(d) for d in depths)
        f.write(depth_line + "\n")
        
        # Write inflow data rows
        #for row in new_inflow_data:
            #row_line = " ".join(str(val) for val in row)
        if state_inflow:
            # for state inflow get only the first extraction starting date
            row_line = f"{int(extraction_dates[0])}" + " " + " ".join(str(val) for val in new_inflow_data)
            f.write(row_line + "\n")
        else: # constant inflow (FOR ONLY ONE EXTRATION CASE)
            for i in range(len(extraction_dates)):
                row_line = f"{int(extraction_dates[i])}" + " " + " ".join(str(val) for val in new_inflow_data[i])
                f.write(row_line + "\n")