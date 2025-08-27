import mesa_reader as mr
import numpy as np
import scipy as sp
import scipy.constants as spc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.optimize as optimize
import glob
import os
import csv
import configparser
from pylab import rcParams

import math
import pandas as pd
import importlib

from astropy import constants as const

print_details=0

#==================================================================
# MESA SIMULATION MATCHING FUNCTION
#==================================================================
def _find_nearest_indices(sorted_grid, values):
    """A helper function to find the indices of the nearest values in a sorted grid."""
    # Find insertion points for all values at once using binary search
    indices = np.searchsorted(sorted_grid, values, side='left')

    # Handle edge cases
    indices = np.maximum(indices, 1)
    indices = np.minimum(indices, len(sorted_grid) - 1)

    # Find the two nearest grid points for each value
    left_neighbors = sorted_grid[indices - 1]
    right_neighbors = sorted_grid[indices]

    # Find which of the two neighbors is closer for each value
    is_left_closer = np.abs(values - left_neighbors) <= np.abs(values - right_neighbors)

    # Return the indices of the closer neighbors
    return np.where(is_left_closer, indices - 1, indices)

def match_to_mesa_sims_fast(quantized_df, base_path, mesa_cols_to_extract, WR_cond="any"):
    """
    Matches quantized stars to MESA simulations using a fast, vectorized approach.
    It classifies stars into multiple detailed types (WR, RSG, OB, A, Dwarfs, etc.)
    into a single 'stellar_type' column based on a physically motivated hierarchy.

    Args:
        quantized_df (pd.DataFrame): DataFrame of stars with quantized mass and FeH.
        base_path (str): The base directory path for the MESA simulations.
        mesa_cols_to_extract (list): A list of column names to extract from history.data.
        WR_cond (str): The condition to classify a star as Wolf-Rayet ('WR').
                       Options: "Xsurf", "eta", "gamma", "gamma_eta", "any".

    Returns:
        pd.DataFrame: A new DataFrame containing the matched stars, their
                      extracted MESA parameters, and a new 'stellar_type' column.
    """
    if not mr:
        print("Error: mesa_reader is not installed. Cannot perform MESA matching.")
        return None

    print(f"\nMatching population to MESA simulations (vectorized) with WR_cond='{WR_cond}'...")
    matched_groups = []

    # Define all columns needed for ANY classification
    wr_classification_cols = {'surface_h1', 'eta', 'eta_trans', 'gamma_edd'}
    general_classification_cols = {'log_Teff', 'log_L', 'center_h1', 'center_he4'}
    all_classification_cols = wr_classification_cols.union(general_classification_cols)
    all_cols_to_extract = set(mesa_cols_to_extract).union(all_classification_cols)

    for (mass, feh), group in quantized_df.groupby(['Mzams', 'Fe/H']):
        z_val = 10**(feh*0.977)
        history_path = os.path.join(base_path, f"{z_val:g}Zsun", f"{mass:g}", "LOGS", "history.data")

        if not os.path.exists(history_path):
            continue

        m = mr.MesaData(history_path)
        mesa_ages_yrs = m.star_age
        max_mesa_age = mesa_ages_yrs[-1]

        star_ages_yrs = group['Age [Gyr]'].values * 1e9
        valid_age_mask = star_ages_yrs <= max_mesa_age
        valid_group = group[valid_age_mask].copy()
        if valid_group.empty:
            continue

        valid_star_ages = star_ages_yrs[valid_age_mask]
        closest_mesa_indices = _find_nearest_indices(mesa_ages_yrs, valid_star_ages)

        for col in all_cols_to_extract:
            if col == "log_Rzams": valid_group[col] = m.log_R[0]
            elif col == "log_Lzams": valid_group[col] = m.log_L[0]
            elif col == "log_Teffzams": valid_group[col] = m.log_Teff[0]
            elif col == "max_age": valid_group[col] = max_mesa_age
            elif col in ["log_RMAX", "RMAX_logL", "RMAX_logTeff"]:
                rmax_idx = np.argmax(m.log_R)
                if col == "log_RMAX": valid_group[col] = m.log_R[rmax_idx]
                elif col == "RMAX_logL": valid_group[col] = m.log_L[rmax_idx]
                elif col == "RMAX_logTeff": valid_group[col] = m.log_Teff[rmax_idx]
            else:
                if col in m.bulk_names:
                    valid_group[col] = m.data(col)[closest_mesa_indices]
                else:
                    valid_group[col] = np.nan

        # --- NEW: UNIFIED HIERARCHICAL STELLAR TYPE CLASSIFICATION ---

        # Get all necessary data columns safely
        s_h1 = valid_group.get('surface_h1', pd.Series(np.nan, index=valid_group.index))
        eta = valid_group.get('eta', pd.Series(np.nan, index=valid_group.index))
        eta_trans = valid_group.get('eta_trans', pd.Series(np.nan, index=valid_group.index))
        gamma_edd = valid_group.get('gamma_edd', pd.Series(np.nan, index=valid_group.index))
        log_Teff = valid_group.get('log_Teff', pd.Series(np.nan, index=valid_group.index))
        log_L = valid_group.get('log_L', pd.Series(np.nan, index=valid_group.index))
        center_h1 = valid_group.get('center_h1', pd.Series(np.nan, index=valid_group.index))
        center_he4 = valid_group.get('center_he4', pd.Series(np.nan, index=valid_group.index))
        Teff = 10**log_Teff

        # --- Define boolean masks for each classification condition ---

        # Priority 1: Wolf-Rayet and naked Helium conditions
        wr_mask = pd.Series(False, index=valid_group.index)
        if WR_cond == "Xsurf": wr_mask = (s_h1 < 0.4)
        elif WR_cond == "eta": wr_mask = (eta > eta_trans)
        elif WR_cond == "gamma": wr_mask = (gamma_edd > 0.5) | (s_h1 < 1e-7)
        elif WR_cond == "gamma_eta": wr_mask = (gamma_edd > 0.5) | (s_h1 < 1e-7) | (eta > eta_trans)
        elif WR_cond == "any": wr_mask = (gamma_edd > 0.5) | (s_h1 < 0.4) | (eta > eta_trans)

        naked_he_mask = pd.Series(False, index=valid_group.index)
        if WR_cond not in ["any", "Xsurf"]:
            naked_he_mask = (~wr_mask) & (s_h1 < 0.4)

        # Priority 2: Supergiants (must be evolved AND highly luminous)
        is_evolved = (center_h1 < 1e-3)
        is_luminous = (log_L > 4.0)
        rsg_mask = (Teff < 4000) & is_luminous & is_evolved
        ysg_mask = (Teff >= 4000) & (Teff < 7500) & is_luminous & is_evolved

        # Priority 3: Hertzsprung Gap (post-MS, pre-CHeB, not a supergiant)
        hg_mask = is_evolved & (center_he4 > 0.95) & (~is_luminous)

        # Priority 4: Red Giants (evolved, cool, but not supergiants)
        rg_mask = (Teff < 5200) & (~is_luminous) & is_evolved

        # Priority 5: Main Sequence Stars ("Dwarfs", Luminosity Class V)
        is_main_sequence = (center_h1 >= 1e-3)
        ob_mask = (Teff >= 10000) & is_main_sequence
        a_mask = (Teff >= 7500) & (Teff < 10000) & is_main_sequence
        f_mask = (Teff >= 6000) & (Teff < 7500) & is_main_sequence
        g_mask = (Teff >= 5200) & (Teff < 6000) & is_main_sequence
        k_mask = (Teff >= 3700) & (Teff < 5200) & is_main_sequence
        m_mask = (Teff < 3700) & is_main_sequence

        # Use np.select for efficient, hierarchical assignment.
        conditions = [
            wr_mask, naked_he_mask, rsg_mask, ysg_mask, hg_mask, rg_mask,
            ob_mask, a_mask, f_mask, g_mask, k_mask, m_mask
        ]
        choices = [
            "WR", "naked_He", "RSG", "YSG", "HG_star", "Red_Giant",
            "OB-type (MS)", "A-type (MS)", "F-type (MS)", "G-type (Yellow Dwarf)",
            "K-type (Orange Dwarf)", "M-type (Red Dwarf)"
        ]

        valid_group['stellar_type'] = np.select(conditions, choices, default="Other_Evolved")

        matched_groups.append(valid_group)

    if not matched_groups:
        print("Warning: No matching MESA files found or no stars fell within model age ranges.")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} stars.")
    return final_df

def match_to_mesa_sims(quantized_df, base_path, mesa_cols_to_extract, Zsun = 0.0142):
    """
    Matches quantized stars to MESA simulation outputs and extracts physical parameters.

    Args:
        quantized_df (pd.DataFrame): DataFrame of stars with quantized mass and FeH.
        base_path (str): The base directory path for the MESA simulations.
        mesa_cols_to_extract (list): A list of column names to extract from history.data.

    Returns:
        pd.DataFrame: A new DataFrame containing the matched stars and their
                      extracted MESA parameters. Returns None if mesa_reader is not available.
    """
    if not mr:
        print("Error: mesa_reader is not installed. Cannot perform MESA matching.")
        return None

    print("\nMatching population to MESA simulations...")
    results = []

    # Group by the quantized parameters to process one MESA file at a time
    for (mass, feh), group in quantized_df.groupby(['Mzams', 'Fe/H']):

        # Construct the path to the MESA history.data file
        # Use '%g' for smart formatting that removes trailing zeros.
        z_val = 10**(feh*0.977)
        z_folder = f"{z_val:g}Zsun"      # e.g., 0.1, 0.2, 0.45, 1
        m_folder = f"{mass:g}"          # e.g., 1.25, 2, 10

        print(f"Z = {z_val:g} Zsun; Mzams = {mass:g}")

        history_path = os.path.join(base_path, z_folder, m_folder, "LOGS", "history.data")

        if not os.path.exists(history_path):
            # print(f"Warning: MESA file not found, skipping: {history_path}")
            continue

        # Load MESA data
        m = mr.MesaData(history_path)
        mesa_ages_yrs = m.star_age # MESA ages are in years

        # For each star in the group that matches this MESA file
        for index, star in group.iterrows():
            star_age_yrs = star['Age [Gyr]'] * 1e9 # Convert Gyr to years

            # Check if the star's age is within the MESA model's lifetime
            if star_age_yrs <= mesa_ages_yrs[-1]:
                # Find the index of the closest age in the MESA track
                closest_mesa_idx = np.abs(mesa_ages_yrs - star_age_yrs).argmin()

                # --- DATA ASSEMBLY ---
                # 1. Start with a copy of all the original data for this star.
                #    This dictionary already contains Mass, Age, FeH, Radius, and Z.
                final_star_data = star.to_dict()

                # 2. Extract the special and regular MESA parameters and add them
                #    to the dictionary for our final output.
                for col in mesa_cols_to_extract:
                    if col == "log_Rzams":
                        final_star_data[f"{col}"] = m.log_R[0]
                    elif col == "log_Lzams":
                        final_star_data[f"{col}"] = m.log_L[0]
                    elif col == "log_Teffzams":
                        final_star_data[f"{col}"] = m.log_Teff[0]
                    elif col == "log_RMAX":
                        final_star_data[f"{col}"] = np.max(m.log_R)
                    elif col == "RMAX_logL":
                        final_star_data[f"{col}"] = m.log_L[np.argmax(m.log_R)]
                    elif col == "RMAX_logTeff":
                        final_star_data[f"{col}"] = m.log_Teff[np.argmax(m.log_R)]
                    elif col == "max_age":
                        final_star_data[f"{col}"] = m.star_age[-1]
                    else:
                        # Default case for any other column
                        final_star_data[f"{col}"] = m.data(col)[closest_mesa_idx]

                # 3. Add the completed dictionary for this star to our results list.
                results.append(final_star_data)

    if not results:
        print("Warning: No matching MESA files found or no stars fell within model age ranges.")
        return None

    print(f"Successfully matched {len(results)} stars to MESA tracks.")
    return pd.DataFrame(results)



def extract_stars(stellar_masses,masses_list, directory, Z_age):

    masses_list=np.sort(masses_list)

    mesa_data_dict, evo_stages_dict,Mzams_counter=process_stars(stellar_masses,masses_list,directory)

    star_para = {'Mzams': [],'Age': [],'Mass': [],'Luminosity': [],'Teff': [],'Radius': [], 'Evo_stage': [],
    'tend_MS': [],'tend_giant': [],'Rmax_MS': [],'Rmax_giant': [],'L_Rmax_MS': [],'L_Rmax_giant': [],'L_WD': []}

    Mass_count=0
    for key in mesa_data_dict:

        # Access the MesaData object for the current key
        evo_stages = evo_stages_dict[key]
        mesa_data = mesa_data_dict[key]

        HG_CHeB_i=max([evo_stages["HG"],evo_stages["CHeB"]])
        late_stages_i=max([evo_stages['HeHG'],evo_stages['CCB']])
        Mzams=round(mesa_data["mass"][0],2)

        for MMM in range(Mzams_counter[key]):

            # Generate a random time between 0 and t
            random_time = np.random.uniform(mesa_data["time"][max([0,evo_stages['ZAMS']])], Z_age)

            if random_time <= mesa_data["time"].iloc[-1]:
                # Find the index in star_age that is the closest to random_time
                closest_index = np.argmin(np.abs(mesa_data["time"] - random_time))

                if Mzams<0.8:                                                   # 1/2/10: MS/giant/WD
                    if closest_index<=evo_stages["MS"]:
                        star_para['Evo_stage'].append(1)
                    else:
                        star_para['Evo_stage'].append(10)
                else:
                    if closest_index<=evo_stages["MS"]:
                        star_para['Evo_stage'].append(1)
                    elif closest_index<=max([HG_CHeB_i,late_stages_i]):
                        star_para['Evo_stage'].append(2)
                    else:
                        star_para['Evo_stage'].append(10)


                star_para['Mzams'].append(Mzams)
                star_para['Age'].append(mesa_data["time"][closest_index])
                star_para['Mass'].append(mesa_data["mass"][closest_index])
                star_para['Luminosity'].append(10**mesa_data["log_L"][closest_index])
                star_para['Teff'].append(10**mesa_data["log_Teff"][closest_index])
                star_para['Radius'].append(10**(mesa_data["log_R"][closest_index]))
                star_para['tend_MS'].append(mesa_data["time"][evo_stages['MS']])
                try:
                    star_para['tend_giant'].append((mesa_data["time"][HG_CHeB_i]))
                except:
                    star_para['tend_giant'].append(-1)
                try:
                    star_para['Rmax_MS'].append(10**max(mesa_data["log_R"][max([0,evo_stages["ZAMS"]]):evo_stages["MS"]]))
                    star_para['L_Rmax_MS'].append(10**(mesa_data["log_L"][np.argmax(mesa_data["log_R"][max([0,evo_stages["ZAMS"]]):evo_stages["MS"]])]))

                except:
                    star_para['Rmax_MS'].append(-1)
                    star_para['L_Rmax_MS'].append(-1)



                try:
                    star_para['Rmax_giant'].append(10**max(mesa_data["log_R"][max([0,evo_stages["MS"]]):max([HG_CHeB_i,late_stages_i])]))
                    star_para['L_Rmax_giant'].append(10**(mesa_data["log_L"][np.argmax(mesa_data["log_R"][max([0,evo_stages["MS"]]):max([HG_CHeB_i,late_stages_i])])]))
                except:
                    star_para['Rmax_giant'].append(-1)
                    star_para['L_Rmax_giant'].append(-1)

                # if HG_CHeB_i==-1 or evo_stages["MS"]==len(mesa_data["time"]):
                    # star_para['L_WD'].append(-1)
                    # star_para['rot_WD'].append(-1)
                # else:
                    # star_para['L_WD'].append(10**mesa_data.mesa_data["log_L"][-1])
                    # star_para['rot_WD'].append(mesa_data.surf_avg_v_rot[-1])

        print("All the stars for {} M\u2609 done".format(Mzams))


        Mass_count+=1


    return star_para

def process_stars(stellar_masses,masses_list,directory):

    mesa_data_dict = {}  # Initialize mesa_data_dict here
    evo_stages_dict = {}  # Initialize evo_stages_dict here

    files_list = []
    data_files = []

    # Mzams_counter= []
    unique_elements, counts = np.unique(np.array(masses_list), return_counts=True)
    Mzams_counter = dict(zip(unique_elements, counts))


    evo_stages = {'ZAMS' : -1, 'MS' : -1, 'HG' : -1, 'CHeB' : -1, 'HeHG' : -1, 'CCB' : -1, 'COB': -1}

    for i in range(len(stellar_masses)):

        data_files.append("{0}/{1}.csv".format(directory,stellar_masses[i]))

    for k in range(len(stellar_masses)):

        df = pd.read_csv(data_files[k])

        key_m = str(stellar_masses[k])
        key_evo_stages = "evo_stages_m{}".format(k)

        mesa_data_dict[key_m] = df

        # Iterate through your sorted evo list to find the last index of each string
        last_indices = {}
        current_string = None
        for index, string in enumerate(df["evo_phase"]):
            if string != current_string:
                current_string = string
            last_indices[string] = index

        # Assign the last indices to your strings in the dictionary
        for string, last_index in last_indices.items():
            evo_stages[string] = last_index

        evo_stages_dict[str(stellar_masses[k])] = evo_stages
        # print(stellar_masses[k],evo_stages)

        evo_stages = {'ZAMS' : -1, 'MS' : -1, 'HG' : -1, 'CHeB' : -1, 'HeHG' : -1, 'CCB' : -1, 'COB' : -1}

    print("Stars processed \n")
    return mesa_data_dict,evo_stages_dict,Mzams_counter

def process_stars_MESA(selected_masses, directory):

    mesa_data_dict = {}  # Initialize mesa_data_dict here
    evo_stages_dict = {}  # Initialize evo_stages_dict here

    files_list = []
    data_files = []
    file_str = 'history'

    Mzams_already_done = []                                                           # I am not processing a file twice
    Mzams_counter= []
    evo_stages = {'ZAMS' : -1, 'MS' : -1, 'HG' : -1, 'CHeB' : -1, 'HeHG' : -1, 'CCB' : -1, 'COB': -1}

    for i in range(len(selected_masses)):

        if selected_masses[i] not in Mzams_already_done:
            Mzams_already_done.append(selected_masses[i])
            Mzams_counter.append(1)
        else:
            Mzams_counter[-1]+=1
            continue


        files_list.append(glob.glob("{0}/{1}/LOGS/*.data".format(directory,selected_masses[i])))
        data_files.append("{0}/{1}/LOGS/{2}.data".format(directory,selected_masses[i],file_str))


    for i in range(len(data_files)):
        exec("m{0} = mr.MesaData(data_files[i])".format(i))
        # Mixing regions list
        exec("m{}.read_data()".format(i))

    for k in range(len(Mzams_already_done)):

        key_m = str(Mzams_already_done[k])
        key_evo_stages = "evo_stages_m{}".format(k)

        mesa_data_dict[key_m] = mr.MesaData(data_files[k])


    frac_exp=-6

    for ind in range(len(Mzams_already_done)):

        m=mesa_data_dict[str(Mzams_already_done[ind])]

        in_H1=0
        out_preZAMS=0

        for i in range(len(m.he_core_mass)):
            #         if 10**m.pp[i]<10**(-3):
            if out_preZAMS==0 and (10**(m.log_Lnuc[i]) < 0.999*10**(m.log_L[i]) or 10**(m.log_Lnuc[i]) > 1.001*10**(m.log_L[i])):
                # When L_nucl != L surface
                evo_stages['ZAMS'] = i
            else:
                out_preZAMS=1
                evo_stages['MS'] = i
                if in_H1==0:
                    in_H1=1
                    if print_details:
                        print('MS started!!', i)


                if m.center_h1[i] < 10**(frac_exp):
                        #             if m.he_core_mass[i]>0.0:
                    evo_stages['MS'] = i
                    real_end_MS=1

                    if print_details:
                        print('We reached the end of the MS!!\n')

                    break
                else:
                    evo_stages['MS'] = i

        in_He4=0
        in_HG=0
        out_HG=0
        out_CHeB=0

        for j in np.arange(i,len(m.he_core_mass)):
            if j >= evo_stages['MS'] and "real_end_MS" in locals():
                if out_HG==0 and (math.isclose(10**(m.log_Lnuc[j]),10**(m.log_L[j]),abs_tol=0.001)) and math.isclose(m.center_he4[j],  m.center_he4[i], abs_tol = 0.01):
                    if in_HG==0:
                        in_HG=1
                        if print_details:
                            print('HG started!!', j)

                        evo_stages['HG'] = j

                    else:
                        if m.center_he4[j] >= 10**(frac_exp):
                            out_HG=1
                            evo_stages['CHeB'] = j
                            if in_He4==0:
                                in_He4=1
                                if print_details:
                                    print()
                                    print('CHeB started!!', j)

                            elif m.center_he4[j] < 10**(frac_exp):
                                evo_stages['CHeB'] = j

                                CO_ratio=m.center_c12[j]/m.center_o16[j]
                                out_CHeB=1

                                if print_details:
                                    print('We reached the end of CHeB!!\n')

                                break

        in_C12=0
        in_HeHG=0
        out_HeHG=0
        out_CCB=0

        #     for k in np.arange(j,len(m.o_core_mass)):
        for k in np.arange(j,len(m.co_core_mass)):

            if k >= evo_stages['CHeB'] and out_CHeB==1:

                try:
                    if m.center_c12[k]/m.center_o16[k]==CO_ratio and out_HeHG==0:
                        evo_stages['HeHG'] = k
                        if in_HeHG==0 and print_details:
                            print('HeHG started!!', k)

                        else:
                            if m.center_c12[k] > 10**(frac_exp):
                                evo_stages['CCB'] = k
                                if in_C12==0:
                                    out_HeHG=1
                                    in_C12=1

                                    if print_details:
                                        print()
                                        print('C burning started!!')


                                elif m.center_c12[k] < 10**(frac_exp):
                                    evo_stages['CCB'] = k
                                    out_CCB = 1

                                    if print_details:
                                        print('We reached the end of the C burning!!')

                                    break
                except:
                    break

        in_O16=0

#     for l in np.arange(k,len(m.o_core_mass)):
        for l in np.arange(k,len(m.co_core_mass)):
            if l >= evo_stages['CCB'] and out_CCB==1:
                try:

                    if m.center_o16[l] > 10**(frac_exp):
                        evo_stages['COB'] = l
                        if in_O16==0:
                            in_O16=1
                            if print_details:
                                print()
                                print('O burning started!!')


                    elif m.center_o16[l] <= 10**(frac_exp):
                        evo_stages['COB'] = l
                        outCCB = 1

                        if print_details:
                            print('We reached the end of the O burning!!')

                        break
                except:
                    pass

        exec('evo_stages_m{} = evo_stages'.format(ind))
        evo_stages_dict[str(Mzams_already_done[ind])] = evo_stages

        evo_stages = {'ZAMS' : -1, 'MS' : -1, 'HG' : -1, 'CHeB' : -1, 'HeHG' : -1, 'CCB' : -1, 'COB' : -1}

    print("Stars processed \n")
    return mesa_data_dict, evo_stages_dict,Mzams_counter



def lumHR(R,logT):

    Rsun =6.9599e+10                                 # cm
    Lsun =3.826e+33                                  # [egr s-1]
    sigma=5.6724e-05*Rsun*Rsun/Lsun                   # Stefan constant [L_sun R_sun^{-2} K^{-4}]
    L = np.log10(4*np.pi*R**2*spc.sigma*(10**logT)**4)

    return L

def lum_Lsun(R,logT):

    Rsun =6.9599e+10                                 # cm
    Lsun =3.826e+33                                  # [egr s-1]
    sigma=5.6724e-05*Rsun*Rsun/Lsun                   # Stefan constant [L_sun R_sun^{-2} K^{-4}]
    L = np.log10(R**2*((10**logT)/5778)**4)

    return L

def Teff(R,logL):

    Rsun =6.9599e+10                                 # cm
    Lsun =3.826e+33                                  # [egr s-1]
    sigma=5.6724e-05*Rsun*Rsun/Lsun                   # Stefan constant [L_sun R_sun^{-2} K^{-4}]
    T = np.log10(((10**logL/R**2)**(1/4))*5778)

    return T


def HRD(star_para):

    # Extract data from star_para
    Mzams_values = star_para['Mzams']
    Luminosity_values = star_para['Luminosity']
    Teff_values = star_para['Teff']

    # Create a colormap for different ZAMS masses
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(min(Mzams_values), max(Mzams_values))

    # Create a scatter plot
    plt.scatter(np.log10(Teff_values), np.log10(Luminosity_values), c=Mzams_values, cmap=cmap, norm=norm, marker='o', alpha=0.85,s=25)

    # Add labels and a colorbar
    plt.xlabel('log($T_{eff}$ [K])',fontsize=25)
    plt.ylabel('log($L$ [L$_{\odot}$])',fontsize=25)
    #plt.title('HR Diagram',labelpad=20)
    plt.tick_params(labelsize=25,axis='both', which='both',top=True,right=True,length=10)

    cbar = plt.colorbar(label='ZAMS Mass')
    cbar.set_label(r'$M_\mathrm{ZAMS}}$ [M$_{\odot}$]', rotation=270, labelpad=25,fontsize=25)

    # invert the x-axis
    plt.gca().invert_xaxis()

    # Show the plot
    plt.show()
