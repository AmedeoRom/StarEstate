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
import dask.dataframe as dd

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

        # --- UNIFIED HIERARCHICAL STELLAR TYPE CLASSIFICATION ---

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

        # Priority 3: AGB Stars
        # Criteria: Evolved, Core Helium Exhausted (center_he4 < 0.01), and NOT classified as RSG/YSG above.
        # Note: Some very luminous AGBs might get caught in RSG above.
        agb_mask = is_evolved & (center_he4 < 0.01) & (~is_luminous)

        # Priority 4: Hertzsprung Gap (post-MS, pre-CHeB, not a supergiant)
        hg_mask = is_evolved & (center_he4 > 0.95) & (~is_luminous)

        # Priority 5: Red Giants (evolved, cool, but not supergiants)
        rg_mask = (Teff < 5200) & (~is_luminous) & is_evolved

        # Priority 6: Main Sequence Stars
        is_main_sequence = (center_h1 >= 1e-3)

        o_mask = (Teff >= 30000) & is_main_sequence
        b_mask = (Teff >= 10000) & (Teff < 30000) & is_main_sequence
        a_mask = (Teff >= 7500) & (Teff < 10000) & is_main_sequence
        f_mask = (Teff >= 6000) & (Teff < 7500) & is_main_sequence
        g_mask = (Teff >= 5200) & (Teff < 6000) & is_main_sequence
        k_mask = (Teff >= 3700) & (Teff < 5200) & is_main_sequence
        m_mask = (Teff < 3700) & is_main_sequence

        conditions = [
            wr_mask, naked_he_mask, rsg_mask, ysg_mask, agb_mask, hg_mask, rg_mask,
            o_mask, b_mask, a_mask, f_mask, g_mask, k_mask, m_mask
        ]
        choices = [
            "WR", "naked_He", "RSG", "YSG", "AGB", "HG_star", "Red_Giant",
            "O-type", "B-type", "A-type", "F-type", "G-type",
            "K-type", "M-type"
        ]

        valid_group['stellar_type'] = np.select(conditions, choices, default="Other_Evolved")

        matched_groups.append(valid_group)

    if not matched_groups:
        print("Warning: No matching MESA files found or no stars fell within model age ranges.")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} stars.")
    return final_df

def match_to_legacy_sims_fast(quantized_df, base_path, WR_cond="any", Zsun=0.0142):
    """
    Matches quantized stars to Legacy (Hurley/BSE/SSE style) simulations.
    Parses 'evolution.dat' files and maps 'Ka' (stellar type) integers to
    the string-based classification system used in MESA matching.

    Implements specific WR/Naked Helium conditions based on Merritt et al. (2025).

    Args:
        quantized_df (pd.DataFrame): DataFrame of stars with quantized mass and FeH.
        base_path (str): The base directory path. Assumes structure: base_path/Z/M/evolution.dat
        WR_cond (str): "Xsurf", "eta", "gamma", "gamma_eta", or "any".
                       Controls which conditions distinguish WR from naked_He.
        Zsun (float): Solar metallicity fraction (default 0.0142).

    Returns:
        pd.DataFrame: DataFrame with legacy parameters and 'stellar_type' column.
    """
    print(f"\nMatching population to LEGACY simulations with WR_cond='{WR_cond}'...")
    matched_groups = []

    for (mass, feh), group in quantized_df.groupby(['Mzams', 'Fe/H']):
        # 1. Determine Metallicity Z for Merritt+ (2025) equations
        z_multiplier = 10**(feh*0.977)
        Z_abs = Zsun * z_multiplier

        # 2. Locate evolution.dat
        # Assuming folder structure: base_path / <z_multiplier>Zsun / <mass> / evolution.dat
        # e.g. /data/0.02Zsun/15/evolution.dat
        evo_path = os.path.join(base_path, f"{z_multiplier:g}Zsun", f"{mass:g}", "evolution.dat")

        if not os.path.exists(evo_path):
            continue

        # 3. Read Legacy Data
        try:
            # Reads: t Ka M R logL logTeff
            legacy_data = pd.read_csv(evo_path, delim_whitespace=True)

            if legacy_data.empty:
                continue

            # Filter out rows where logL is -99 (invalid/terminated steps)
            legacy_data = legacy_data[legacy_data['logL'] != -99]

            # Check if data is empty after filtering
            if legacy_data.empty:
                continue
        except Exception as e:
            print(f"Error reading {evo_path}: {e}")
            continue

        # 4. Match Ages
        model_ages_yrs = legacy_data['t'].values * 1e6 # Legacy 't' is usually in Myr
        max_model_age = model_ages_yrs[-1]

        star_ages_yrs = group['Age [Gyr]'].values * 1e9
        valid_age_mask = star_ages_yrs <= max_model_age
        valid_group = group[valid_age_mask].copy()

        if valid_group.empty:
            continue

        valid_star_ages = star_ages_yrs[valid_age_mask]
        # Find nearest indices
        closest_indices = _find_nearest_indices(model_ages_yrs, valid_star_ages)

        # 5. Extract Legacy Parameters
        # Map columns to standard names for consistency
        valid_group['log_L'] = legacy_data['logL'].values[closest_indices]
        valid_group['log_Teff'] = legacy_data['logTeff'].values[closest_indices]
        valid_group['Mass'] = legacy_data['M'].values[closest_indices]
        valid_group['Radius'] = legacy_data['R'].values[closest_indices] # Solar radii
        valid_group['Ka'] = legacy_data['Ka'].values[closest_indices].astype(int)

        # 6. Unified Classification Logic

        # Prepare vectors
        Ka = valid_group['Ka']
        log_L = valid_group['log_L']
        log_Teff = valid_group['log_Teff']
        mass_current = valid_group['Mass']
        Teff = 10**log_Teff
        L = 10**log_L

        # --- WR / Naked Helium Logic (Merritt 2025) ---
        # Equation 1: Gamma_e = 2.49e-5 * L * M^-1 (Solar units)
        gamma_e = 2.49e-5 * L * (mass_current**-1)

        # Equation 5 (Switch Luminosity for Eta condition): L_switch = 10^2.36 * Z^-1.91
        # Note: This uses absolute Z (mass fraction).
        l_switch = (10**2.36) * (Z_abs**-1.91)

        # Evaluation Booleans
        # "For the Xsurf transition... assume that any K value for naked He is a WR"
        # This implies X_surf condition is effectively True for all Ka in [7,8,9]
        cond_xsurf_met = pd.Series(True, index=valid_group.index)
        cond_gamma_met = gamma_e > 0.5
        cond_eta_met = L > l_switch

        if WR_cond == "any":
            # If any condition is met. Since Xsurf is assumed True for Ka 7-9,
            # this makes ALL Ka = 7-9 WRs, as well as some MS VMSs.
            wr_satisfied = cond_xsurf_met | cond_gamma_met | cond_eta_met
        elif WR_cond == "Xsurf":
            wr_satisfied = cond_xsurf_met
        elif WR_cond == "gamma":
            wr_satisfied = cond_gamma_met
        elif WR_cond == "eta":
            wr_satisfied = cond_eta_met
        elif WR_cond == "gamma_eta":
            wr_satisfied = cond_gamma_met | cond_eta_met
        else:
            wr_satisfied = cond_xsurf_met

        # Define Masks

        # 1. Naked Helium / WR (Ka 7, 8, 9)
        is_he_star = Ka.isin([7, 8, 9])
        wr_mask = is_he_star & wr_satisfied
        naked_he_mask = is_he_star & (~wr_satisfied)

        # 2. Supergiants (Luminous & Evolved)
        # Ka 2-6 (HG through AGB) can be supergiants if luminous
        is_post_ms = Ka.isin([2, 3, 4, 5, 6])
        is_luminous = log_L > 4.0
        rsg_mask = is_post_ms & is_luminous & (Teff < 4000)
        ysg_mask = is_post_ms & is_luminous & (Teff >= 4000) & (Teff < 7500)

        # 3. AGB (Ka 5, 6) - Explicit in Legacy
        # Note: We check RSG/YSG first. If not super luminous/cool, it's generic AGB.
        agb_mask = Ka.isin([5, 6]) & (~rsg_mask) & (~ysg_mask)

        # 4. Hertzsprung Gap (Ka 2)
        hg_mask = (Ka == 2) & (~rsg_mask) & (~ysg_mask)

        # 5. Red Giants (Ka 3 = GB, Ka 4 = CHeB)
        # CHeB (Ka 4) can be blue loops. If hot, we shouldn't call it "Red Giant".
        # If Ka=4 and hot (and not YSG), it falls to default "Other_Evolved".
        rg_mask = Ka.isin([3, 4]) & (Teff < 5200) & (~rsg_mask) & (~ysg_mask)

        # 6. Main Sequence (Ka 0, 1)
        is_main_sequence = Ka.isin([0, 1])
        o_mask = (Teff >= 30000) & is_main_sequence
        b_mask = (Teff >= 10000) & (Teff < 30000) & is_main_sequence
        a_mask = (Teff >= 7500) & (Teff < 10000) & is_main_sequence
        f_mask = (Teff >= 6000) & (Teff < 7500) & is_main_sequence
        g_mask = (Teff >= 5200) & (Teff < 6000) & is_main_sequence
        k_mask = (Teff >= 3700) & (Teff < 5200) & is_main_sequence
        m_mask = (Teff < 3700) & is_main_sequence

        # 7. Remnants (Ka 10, 11, 12, 13, 14, 15)
        wd_mask = Ka.isin([10, 11, 12])
        ns_mask = (Ka == 13)
        bh_mask = (Ka == 14)

        conditions = [
            wr_mask, naked_he_mask, rsg_mask, ysg_mask, agb_mask, hg_mask, rg_mask,
            o_mask, b_mask, a_mask, f_mask, g_mask, k_mask, m_mask,
            wd_mask, ns_mask, bh_mask
        ]
        choices = [
            "WR", "naked_He", "RSG", "YSG", "AGB", "HG_star", "Red_Giant",
            "O-type (MS)", "B-type (MS)", "A-type (MS)", "F-type (MS)", "G-type (Yellow Dwarf)",
            "K-type (Orange Dwarf)", "M-type (Red Dwarf)",
            "White_Dwarf", "Neutron_Star", "Black_Hole"
        ]

        valid_group['stellar_type'] = np.select(conditions, choices, default="Other_Evolved")

        matched_groups.append(valid_group)

    if not matched_groups:
        print("Warning: No matching LEGACY files found.")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} stars (Legacy).")
    return final_df

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

def process_matched_population(
    q_df=None,
    mode="MESA",
    what_to_do="generate",
    mesa_dir="./MESA Simulations",
    sse_dir="../SSE/Simulations",
    output_dir_mesa="./mesa_matched_chunks",
    output_dir_sse="./sse_matched_chunks",
    chunk_size=10_000_000,
    wr_cond="any",
    mesa_cols=None
):
    """
    Wrapper function to handle chunking, generating, or extracting MESA/SSE populations.

    Args:
        q_df (pd.DataFrame, optional): Input dataframe of stars (quantized). Defaults to None.
        mode (str): "MESA" or "SSE".
        what_to_do (str): "generate" or "just_extract".
        mesa_dir (str): Path to MESA simulations.
        sse_dir (str): Path to SSE simulations.
        output_dir_mesa (str): Directory to save/load MESA parquets.
        output_dir_sse (str): Directory to save/load SSE parquets.
        chunk_size (int): Number of rows per chunk for large generation.
        wr_cond (str): Wolf-Rayet condition ('any', 'Xsurf', etc.).
        mesa_cols (list): Columns to extract for MESA. Defaults to standard set if None.

    Returns:
        dask.dataframe: The resulting dataframe (lazy loaded).
    """

    # 1. Select the output directory based on mode
    if mode == "MESA":
        output_dir = output_dir_mesa
    elif mode == "SSE":
        output_dir = output_dir_sse
    else:
        print(f"Error: Unknown mode '{mode}'. Choose 'MESA' or 'SSE'.")
        return None

    # Default MESA columns if not provided
    if mesa_cols is None:
        mesa_cols = [
            "star_age", "star_mass", "log_R", "log_Teff", "log_L",
            "log_Rzams", "log_Lzams", "log_Teffzams", "log_RMAX",
            "RMAX_logL", "RMAX_logTeff", "max_age"
        ]

    # Handle numeric limit logic
    # If q_df is provided, we use its length. If not (extract only), we assume large or handle differently.
    num_stars = len(q_df) if q_df is not None else 0

    # Check dependencies
    if what_to_do == 'generate' and q_df is None:
        print("Error: 'generate' mode requires an input DataFrame (q_df).")
        return None

    df_stars = None

    # --- MESA LOGIC ---
    if mode == "MESA":

        # 1. Just Extract
        if what_to_do == 'just_extract':
            if num_stars > 0 and num_stars < 1e7:
                 # Try single file first
                 single_path = os.path.join(output_dir, 'matched_chunk.parquet')
                 if os.path.exists(single_path):
                     df_stars = dd.read_parquet(single_path)
                 else:
                     # Fallback to wildcard
                     df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            else:
                 # Default to wildcard for extract mode if N is unknown or large
                 df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))

            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

        # 2. Generate (Small N)
        elif num_stars < 1e7:
            os.makedirs(output_dir, exist_ok=True)
            df_mesa = match_to_mesa_sims_fast(q_df, mesa_dir, mesa_cols, WR_cond=wr_cond)

            if df_mesa is not None and not df_mesa.empty:
                output_path = os.path.join(output_dir, "matched_chunk.parquet")
                df_mesa.to_parquet(output_path)
                print(f"Chunk successfully processed and saved to {output_path}")
                df_stars = dd.read_parquet(output_path)

        # 3. Generate (Large N - Chunking)
        else:
            os.makedirs(output_dir, exist_ok=True)
            num_chunks = int(np.ceil(num_stars / chunk_size))
            print(f"Processing {num_stars} stars in {num_chunks} chunks of size {chunk_size}...")

            for i in range(num_chunks):
                start_index = i * chunk_size
                end_index = (i + 1) * chunk_size
                df_chunk = q_df.iloc[start_index:end_index]

                print(f"\n--- Processing Chunk {i+1}/{num_chunks} ---")
                matched_chunk_df = match_to_mesa_sims_fast(df_chunk, mesa_dir, mesa_cols, WR_cond=wr_cond)

                if matched_chunk_df is not None and not matched_chunk_df.empty:
                    output_path = os.path.join(output_dir, f"matched_chunk_{i+1}.parquet")
                    matched_chunk_df.to_parquet(output_path)
                    print(f"Chunk {i+1} successfully processed and saved to {output_path}")

            print("\nAll chunks processed!")
            df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

    # --- SSE LOGIC ---
    elif mode == "SSE":

        # 1. Just Extract
        if what_to_do == 'just_extract':
            if num_stars > 0 and num_stars < 1e7:
                 single_path = os.path.join(output_dir, 'matched_chunk.parquet')
                 if os.path.exists(single_path):
                     df_stars = dd.read_parquet(single_path)
                 else:
                     df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            else:
                 df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))

            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

        # 2. Generate (Small N)
        elif num_stars < 1e7:
            os.makedirs(output_dir, exist_ok=True)
            # Pass WR_cond to SSE as well since the function supports it
            df_sse = match_to_legacy_sims_fast(q_df, sse_dir, WR_cond=wr_cond)

            if df_sse is not None and not df_sse.empty:
                output_path = os.path.join(output_dir, "matched_chunk.parquet")
                df_sse.to_parquet(output_path)
                print(f"Processed and saved to {output_path}")
                df_stars = dd.read_parquet(output_path)

        # 3. Generate (Large N - Chunking)
        else:
            os.makedirs(output_dir, exist_ok=True)
            num_chunks = int(np.ceil(num_stars / chunk_size))
            print(f"Processing {num_stars} stars in {num_chunks} chunks of size {chunk_size} using SSE tracks...")

            for i in range(num_chunks):
                start_index = i * chunk_size
                end_index = (i + 1) * chunk_size
                df_chunk = q_df.iloc[start_index:end_index]

                print(f"\n--- Processing Chunk {i+1}/{num_chunks} ---")
                matched_chunk_df = match_to_legacy_sims_fast(df_chunk, sse_dir, WR_cond=wr_cond)

                if matched_chunk_df is not None and not matched_chunk_df.empty:
                    output_path = os.path.join(output_dir, f"matched_chunk_{i+1}.parquet")
                    matched_chunk_df.to_parquet(output_path)
                    print(f"Chunk {i+1} successfully processed and saved to {output_path}")

            print("\nAll chunks processed!")
            df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

    return df_stars
