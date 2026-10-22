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
# HELPER FUNCTIONS
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

def _generate_binary_params(df, q_grid, p_grid, q_dist_type="uniform", p_dist_type="power_law_Sana2025", min_mass=0.11):
    """
    Generates 'q' and 'Period' columns for the dataframe based on specified distributions
    and snaps the values to the nearest available grid points.

    Args:
        df (pd.DataFrame): Input dataframe.
        q_grid (np.array): Array of available q values.
        p_grid (np.array): Array of available logP values.
        q_dist_type (str): 'uniform', 'other'.
        p_dist_type (str): 'uniform', 'power_law','power_law_Sana2025'.

    Returns:
        pd.DataFrame: Dataframe with added 'q' and 'Period' columns.
    """
    n = len(df)

    # --- Generate q (Mass Ratio) ---
    # 1. Calculate the minimum valid q for each star: q_min = min_mass / Mzams
    # 2. Clip it to ensure it stays within the grid bounds.
    #    CRITICAL: If Mzams is very small (e.g. < min_mass), the required q would be > 1.0.
    #    We must clip this to the maximum available q in the grid (e.g. 1.0 or 0.9)
    #    to avoid generating invalid mass ratios.
    min_q_per_star = np.clip(min_mass / df['Mzams'].values, q_grid.min(), q_grid.max())

    # --- Generate q (Mass Ratio) ---
    if q_dist_type == "uniform":
        # Continuous uniform [min, max], then snap to grid
        # We generate a random number 0-1 and scale it to the valid range for each star
        u = np.random.random(n)
        q_max_grid = q_grid.max()
        range_width = q_max_grid - min_q_per_star
        raw_vals = min_q_per_star + (u * range_width)

        # Snap to grid
        idx = _find_nearest_indices(q_grid, raw_vals)
        df['q'] = q_grid[idx]

        # --- SAFETY CHECK ---
        # Snapping to 'nearest' might round down (e.g. 0.16 -> 0.1).
        # If that pushes M2 below min_mass, we must force it up to the next grid point.
        # We use a small tolerance (0.001) for float comparison.
        actual_m2 = df['q'] * df['Mzams']
        invalid_mask = actual_m2 < (min_mass - 0.001)

        if invalid_mask.any():
            # For these rows, we bump the index up by 1 to select the next higher q.
            # np.minimum ensures we don't go out of bounds if we are already at max q.
            # If we are at max q and M2 is still < min_mass (because M1 is tiny),
            # it stays at max q (best effort).
            new_idx = np.minimum(idx[invalid_mask] + 1, len(q_grid) - 1)
            df.loc[invalid_mask, 'q'] = q_grid[new_idx]

    elif q_dist_type == "other":
        # Placeholder for other q distributions (e.g. f(q) ~ q^0)
        # For now, default to uniform logic
        print("Warning: 'other' q_dist_type not implemented. Using uniform selection.")
        raw_vals = np.random.uniform(q_grid.min(), q_grid.max(), n)
        idx = _find_nearest_indices(q_grid, raw_vals)
        df['q'] = q_grid[idx]

    else:
        # Default fallback
        print(f"Warning: Unknown q_dist_type '{q_dist_type}'. Using uniform.")
        raw_vals = np.random.uniform(q_grid.min(), q_grid.max(), n)
        idx = _find_nearest_indices(q_grid, raw_vals)
        df['q'] = q_grid[idx]


    # --- Generate Period (logP) ---
    # Note: The directory structure uses 'per_{val}', which corresponds to logP grid values.

    if p_dist_type == "uniform":
        # Continuous uniform on logP [min, max], then snap
        raw_vals = np.random.uniform(p_grid.min(), p_grid.max(), n)
        idx = _find_nearest_indices(p_grid, raw_vals)
        df['Period'] = p_grid[idx]

    elif "power_law"  in p_dist_type:
        # Sana et al. style: f(logP) ~ (logP)^alpha
        # log(P) in [0.5, 5.5] with exponent alpha = -0.55 (Sana et al. 2012; de Mink & Belczynski 2015)
        if p_dist_type == "power_law":
            alpha = -0.55
        if p_dist_type == "power_law_Sana2025":
            alpha = -0.1
        # Use grid bounds or hardcoded standard? usually better to use grid bounds for consistency
        vmin, vmax = np.min(p_grid), np.max(p_grid)

        # Inverse transform sampling
        # CDF(x) ~ x^(alpha+1)
        u = np.random.random(n)
        norm_factor = (vmax**(alpha+1) - vmin**(alpha+1))
        raw_vals = (u * norm_factor + vmin**(alpha+1))**(1/(alpha+1))

        # Snap to grid
        idx = _find_nearest_indices(p_grid, raw_vals)
        df['Period'] = p_grid[idx]

    else:
        # Default fallback
        print(f"Warning: Unknown p_dist_type '{p_dist_type}'. Using uniform.")
        raw_vals = np.random.uniform(p_grid.min(), p_grid.max(), n)
        idx = _find_nearest_indices(p_grid, raw_vals)
        df['Period'] = p_grid[idx]

    return df

def _classify_legacy_rows(df, WR_cond, Z_abs, use_orig_thick_winds=False):
    """
    Helper function to apply unified stellar type classification to a dataframe
    containing Legacy (BSE/SSE) parameters (Ka, log_L, log_Teff, Mass).

    Args:
        df (pd.DataFrame): Dataframe with 'Ka', 'log_L', 'log_Teff', 'Mass'.
        WR_cond (str): Wolf-Rayet condition string.
        Z_abs (float): Absolute metallicity for the group.
        use_orig_thick_winds (bool): If True, dynamically determines WR condition from 'thick_winds_trans_condition' column.

    Returns:
        pd.DataFrame: The input dataframe with a new 'stellar_type' column.
    """
    # Prepare vectors
    Ka = df['Ka']
    log_L = df['log_L']
    log_Teff = df['log_Teff']
    mass_current = df['Mass']
    Teff = 10**log_Teff
    L = 10**log_L

    # --- WR / Naked Helium Logic ---
    # Check if detailed transition values from BSE/SSE winds are present
    has_detailed_vars = all(col in df.columns for col in ['Gamma_e', 'eta', 'xsurf', 'Gamma_e_trans', 'eta_trans'])

    if has_detailed_vars:
        gamma_e = df['Gamma_e']
        gamma_trans = df['Gamma_e_trans'] if 'Gamma_e_trans' in df.columns else 0.5
        eta = df['eta']
        eta_trans = df['eta_trans']
        xsurf = df['xsurf']

        cond_xsurf_met = xsurf < 0.4
        cond_gamma_met = gamma_e > gamma_trans
        cond_eta_met = eta > eta_trans
    else:
        # Fallback approximation (Merritt 2025)
        # Equation 1: Gamma_e = 2.49e-5 * L * M^-1 (Solar units)
        gamma_e = 2.49e-5 * L * (mass_current**-1)
        cond_gamma_met = gamma_e > 0.5

        # Equation 5 (Switch Luminosity for Eta condition): L_switch = 10^2.36 * Z^-1.91
        l_switch = (10**2.36) * (Z_abs**-1.91)
        cond_eta_met = L > l_switch

        # Without detailed tracking, we cannot strictly evaluate xsurf < 0.4, defaulting to False
        cond_xsurf_met = pd.Series(False, index=df.index)

    # Legacy Option (SSE standard logic for WR / Naked Helium)
    cond_sse_met = (Ka.isin([7, 8, 9])) | (log_Teff > 4.4)

    # Evaluation Booleans
    if use_orig_thick_winds:
        if 'thick_winds_trans_condition' in df.columns:
            trans_cond = df['thick_winds_trans_condition'].astype(str)
            wr_satisfied = (
                ((trans_cond == "gamma") & cond_gamma_met) |
                ((trans_cond == "eta") & cond_eta_met) |
                ((trans_cond.isin(["Xsurf", "xsurf"])) & cond_xsurf_met)
            )
        else:
            wr_satisfied = cond_sse_met
    else:
        if WR_cond == "any":
            wr_satisfied = cond_xsurf_met | cond_gamma_met | cond_eta_met | cond_sse_met
        elif WR_cond == "Xsurf":
            wr_satisfied = cond_xsurf_met
        elif WR_cond == "gamma":
            wr_satisfied = cond_gamma_met
        elif WR_cond == "eta":
            wr_satisfied = cond_eta_met
        elif WR_cond == "gamma_eta":
            wr_satisfied = cond_gamma_met | cond_eta_met
        elif WR_cond == "SSE":
            wr_satisfied = cond_sse_met
        else:
            wr_satisfied = cond_sse_met

    # Define Masks

    # 1. WR Stars: Any star satisfying the physical WR conditions
    # This captures classic WRs (Ka 7-9 + hot) AND potential pre-WR stripped stars
    wr_mask = wr_satisfied

    # 2. Naked Helium: Stars that are explicitly naked helium (Ka 7-9)
    # but did NOT satisfy the WR condition (e.g. too cool/dim)
    naked_he_mask = Ka.isin([7, 8, 9]) & (~wr_mask)

    # 3. Supergiants (Luminous & Evolved & Still H-rich)
    is_post_ms = Ka.isin([2, 3, 4, 5, 6])
    is_luminous = log_L > 4.0
    rsg_mask = is_post_ms & is_luminous & (Teff < 4000)
    ysg_mask = is_post_ms & is_luminous & (Teff >= 4000) & (Teff < 7500)

    # 4. AGB (Ka 5, 6)
    agb_mask = Ka.isin([5, 6]) & (~rsg_mask) & (~ysg_mask)

    # 5. Hertzsprung Gap (Ka 2)
    hg_mask = (Ka == 2) & (~rsg_mask) & (~ysg_mask)

    # 6. Red Giants (Ka 3 = GB, Ka 4 = CHeB)
    rg_mask = Ka.isin([3, 4]) & (Teff < 5200) & (~rsg_mask) & (~ysg_mask)

    # 7. Main Sequence (Ka 0, 1)
    is_main_sequence = Ka.isin([0, 1])
    o_mask = (Teff >= 30000) & is_main_sequence
    b_mask = (Teff >= 10000) & (Teff < 30000) & is_main_sequence
    a_mask = (Teff >= 7500) & (Teff < 10000) & is_main_sequence
    f_mask = (Teff >= 6000) & (Teff < 7500) & is_main_sequence
    g_mask = (Teff >= 5200) & (Teff < 6000) & is_main_sequence
    k_mask = (Teff >= 3700) & (Teff < 5200) & is_main_sequence
    m_mask = (Teff < 3700) & is_main_sequence

    # 8. Remnants
    wd_mask = Ka.isin([10, 11, 12])
    ns_mask = (Ka == 13)
    bh_mask = (Ka == 14)

    conditions = [
        bh_mask, ns_mask, wd_mask,  # <--- Evaluated First
        wr_mask, naked_he_mask, rsg_mask, ysg_mask, agb_mask, hg_mask, rg_mask,
        o_mask, b_mask, a_mask, f_mask, g_mask, k_mask, m_mask
    ]
    choices = [
        "Black_Hole", "Neutron_Star", "White_Dwarf", # BH is index 0, NS is 1, WD is 2
        "WR", "naked_He", "RSG", "YSG", "AGB", "HG_star", "Red_Giant",
        "O-type (MS)", "B-type (MS)", "A-type (MS)", "F-type (MS)", "G-type (Yellow Dwarf)",
        "K-type (Orange Dwarf)", "M-type (Red Dwarf)"
    ]

    df['stellar_type'] = np.select(conditions, choices, default="Other_Evolved")
    return df

#==================================================================
# MESA SIMULATION MATCHING FUNCTION
#==================================================================

def match_to_mesa_sims_fast(quantized_df, base_path, mesa_cols_to_extract, WR_cond="Xsurf", no_WR=False):
    """
    Matches quantized stars to MESA simulations.

    Args include:
        no_WR (bool): If True, filters out stars classified as Wolf-Rayet.
    """
    if not mr:
        print("Error: mesa_reader is not installed. Cannot perform MESA matching.")
        return None

    print(f"\nMatching population to MESA simulations (vectorized) with WR_cond='{WR_cond}'...")
    matched_groups = []

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

        # --- UNIFIED HIERARCHICAL STELLAR TYPE CLASSIFICATION (MESA SPECIFIC) ---
        s_h1 = valid_group.get('surface_h1', pd.Series(np.nan, index=valid_group.index))
        eta = valid_group.get('eta', pd.Series(np.nan, index=valid_group.index))
        eta_trans = valid_group.get('eta_trans', pd.Series(np.nan, index=valid_group.index))
        gamma_edd = valid_group.get('gamma_edd', pd.Series(np.nan, index=valid_group.index))
        log_Teff = valid_group.get('log_Teff', pd.Series(np.nan, index=valid_group.index))
        log_L = valid_group.get('log_L', pd.Series(np.nan, index=valid_group.index))
        center_h1 = valid_group.get('center_h1', pd.Series(np.nan, index=valid_group.index))
        center_he4 = valid_group.get('center_he4', pd.Series(np.nan, index=valid_group.index))
        Teff = 10**log_Teff

        wr_mask = pd.Series(False, index=valid_group.index)
        if WR_cond == "Xsurf": wr_mask = (s_h1 < 0.4)
        elif WR_cond == "eta": wr_mask = (eta > eta_trans)
        elif WR_cond == "gamma": wr_mask = (gamma_edd > 0.5) | (s_h1 < 1e-7)
        elif WR_cond == "gamma_eta": wr_mask = (gamma_edd > 0.5) | (s_h1 < 1e-7) | (eta > eta_trans)
        elif WR_cond == "any": wr_mask = (gamma_edd > 0.5) | (s_h1 < 0.4) | (eta > eta_trans)

        naked_he_mask = pd.Series(False, index=valid_group.index)
        if WR_cond not in ["any", "Xsurf"]:
            naked_he_mask = (~wr_mask) & (s_h1 < 0.4)

        is_evolved = (center_h1 < 1e-3)
        is_luminous = (log_L > 4.0)
        rsg_mask = (Teff < 4000) & is_luminous & is_evolved
        ysg_mask = (Teff >= 4000) & (Teff < 7500) & is_luminous & is_evolved
        agb_mask = is_evolved & (center_he4 < 0.01) & (~is_luminous)
        hg_mask = is_evolved & (center_he4 > 0.95) & (~is_luminous)
        rg_mask = (Teff < 5200) & (~is_luminous) & is_evolved
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

        # --- APPLY FILTERS ---
        if no_WR:
            valid_group = valid_group[valid_group['stellar_type'] != 'WR']

        matched_groups.append(valid_group)

    if not matched_groups:
        print("Warning: No matching MESA files found or all filtered out.")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} stars.")
    return final_df

def match_to_legacy_sims_fast(quantized_df, base_path, WR_cond="any", Zsun=0.0142, no_WR=False, no_compactOBJ=False, use_orig_thick_winds=False):
    """
    Matches quantized stars to Legacy (Hurley/SSE style) simulations.

    Args:
        no_WR (bool): If True, filters out stars classified as Wolf-Rayet.
        no_compactOBJ (bool): If True, filters out stars with K > 9 (WDs, NS, BH).
                              Always filters out K > 14 (massless/SN remnants).
    """
    print(f"\nMatching population to LEGACY simulations with WR_cond='{WR_cond}'...")
    matched_groups = []

    for (mass, feh), group in quantized_df.groupby(['Mzams', 'Fe/H']):
        z_multiplier = 10**(feh*0.977)
        Z_abs = Zsun * z_multiplier
        evo_path = os.path.join(base_path, f"{z_multiplier:g}Zsun", f"{mass:g}", "evolution.dat")

        if not os.path.exists(evo_path):
            continue

        try:
            legacy_data = pd.read_csv(evo_path, sep=r'\s+')
            if legacy_data.empty: continue
            legacy_data = legacy_data[legacy_data['logL'] != -99]
            if legacy_data.empty: continue
        except Exception as e:
            print(f"Error reading {evo_path}: {e}")
            continue

        model_ages_yrs = legacy_data['t'].values * 1e6
        max_model_age = model_ages_yrs[-1]

        star_ages_yrs = group['Age [Gyr]'].values * 1e9
        valid_age_mask = star_ages_yrs <= max_model_age
        valid_group = group[valid_age_mask].copy()

        if valid_group.empty:
            continue

        valid_star_ages = star_ages_yrs[valid_age_mask]
        closest_indices = _find_nearest_indices(model_ages_yrs, valid_star_ages)

        # Extract Parameters
        valid_group['log_L'] = legacy_data['logL'].values[closest_indices]
        valid_group['log_Teff'] = legacy_data['logTeff'].values[closest_indices]
        valid_group['Mass'] = legacy_data['M'].values[closest_indices]
        valid_group['Radius'] = legacy_data['R'].values[closest_indices]

        ka_values = legacy_data['Ka'].values[closest_indices].astype(int)
        valid_group['Ka'] = ka_values

        # --- FILTER 1: Compact Objects ---
        # "Always make it skip with K > 14" -> Filter K <= 14
        # "no_compactOBJ" -> Filter K <= 9

        max_k = 9 if no_compactOBJ else 14

        # Apply filter
        k_mask = (valid_group['Ka'] <= max_k)
        if not k_mask.any():
            continue

        valid_group = valid_group[k_mask]

        # Classify
        valid_group = _classify_legacy_rows(valid_group, WR_cond, Z_abs, use_orig_thick_winds)

        # --- FILTER 2: WR Stars ---
        if no_WR:
            valid_group = valid_group[valid_group['stellar_type'] != 'WR']

        if not valid_group.empty:
            matched_groups.append(valid_group)

    if not matched_groups:
        print("Warning: No matching LEGACY files found or all filtered out.")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} stars (Legacy).")
    return final_df

def match_to_legacy_sims_BSE_fast(
    quantized_df, base_path, WR_cond="any", Zsun=0.0142,
    q_dist_type="uniform", p_dist_type="uniform", q_grid=None, p_grid=None,
    no_WR=False, no_compactOBJ=False, use_orig_thick_winds=False
):
    """
    Matches quantized binary stars to Legacy (Hurley/BSE style) simulations.

    Filters:
    - Always skips systems where Ka > 14 or Kb > 14.
    - If no_compactOBJ=True, skips systems where Ka > 9 or Kb > 9.
    - If no_WR=True, skips systems where EITHER star is classified as 'WR'.
    """
    print(f"\nMatching BINARY population to BSE simulations with WR_cond='{WR_cond}'...")

    if q_grid is None: q_grid = DEFAULT_BSE_GRID_Q
    if p_grid is None: p_grid = DEFAULT_BSE_GRID_LOGP

    # 0. Generate Params if needed
    if 'q' not in quantized_df.columns or 'Period' not in quantized_df.columns:
        print(f"Generating binary parameters (q: {q_dist_type}, Period: {p_dist_type})...")
        quantized_df = _generate_binary_params(quantized_df, q_grid, p_grid, q_dist_type, p_dist_type)

    required_cols = ['Mzams', 'Fe/H', 'q', 'Period']
    for col in required_cols:
        if col not in quantized_df.columns:
            print(f"Error: Missing column '{col}' in quantized_df for BSE matching.")
            return None

    matched_groups = []
    group_cols = ['Mzams', 'Fe/H', 'q', 'Period']

    if quantized_df.empty:
        print("Warning: Input dataframe is empty.")
        return None

    for (mass, feh, q, per), group in quantized_df.groupby(group_cols):

        z_multiplier = 10**(feh*0.977)
        Z_abs = Zsun * z_multiplier

        path_parts = [
            base_path,
            f"{z_multiplier:g}Zsun",
            f"{mass:g}",
            f"q_{q:g}",
            f"per_{per:g}",
            "evolution.dat"
        ]
        evo_path = os.path.join(*path_parts)

        if not os.path.exists(evo_path):
            continue

        try:
            legacy_data = pd.read_csv(evo_path, sep=r'\s+')
            if legacy_data.empty: continue
            is_remnant = legacy_data['Ka'] >= 10
            legacy_data = legacy_data[(legacy_data['logLa'] != -99) | is_remnant]
            if legacy_data.empty: continue
        except Exception as e:
            print(f"Error reading {evo_path}: {e}")
            continue

        model_ages_yrs = legacy_data['t'].values * 1e6
        max_model_age = model_ages_yrs[-1]

        star_ages_yrs = group['Age [Gyr]'].values * 1e9
        valid_age_mask = star_ages_yrs <= max_model_age
        valid_group = group[valid_age_mask].copy()

        if valid_group.empty:
            continue

        valid_star_ages = star_ages_yrs[valid_age_mask]
        closest_indices = _find_nearest_indices(model_ages_yrs, valid_star_ages)

        # --- EXTRACT PARAMETERS ---
        kas_at_age = legacy_data['Ka'].values[closest_indices].astype(int)
        kbs_at_age = legacy_data['Kb'].values[closest_indices].astype(int)

        # --- FILTER 1: COMPACT OBJECTS ---
        max_k = 9 if no_compactOBJ else 14
        k_mask = (kas_at_age <= max_k) & (kbs_at_age <= max_k)

        if not k_mask.any():
            continue

        final_group = valid_group[k_mask].copy()
        final_indices = closest_indices[k_mask]

        final_group['binary_id'] = final_group.index

        # --- CONSTRUCT DATAFRAMES (Before WR Filter) ---

        # Primary (A)
        df_a = final_group.copy()
        df_a['component'] = 'A'
        df_a['Ka'] = legacy_data['Ka'].values[final_indices].astype(int)
        df_a['Mass'] = legacy_data['Ma'].values[final_indices]
        df_a['Radius'] = legacy_data['Ra'].values[final_indices]
        df_a['log_L'] = legacy_data['logLa'].values[final_indices]
        df_a['log_Teff'] = legacy_data['logTeffa'].values[final_indices]
        df_a['separation'] = legacy_data['a'].values[final_indices]

        # Secondary (B)
        df_b = final_group.copy()
        df_b['component'] = 'B'
        df_b['Ka'] = legacy_data['Kb'].values[final_indices].astype(int)
        df_b['Mass'] = legacy_data['Mb'].values[final_indices]
        df_b['Radius'] = legacy_data['Rb'].values[final_indices]
        df_b['log_L'] = legacy_data['logLb'].values[final_indices]
        df_b['log_Teff'] = legacy_data['logTeffb'].values[final_indices]
        df_b['separation'] = legacy_data['a'].values[final_indices]

        if "thick_winds_trans_condition" in legacy_data.columns:
            df_a['xsurf'] = legacy_data["xsurfa"].values[final_indices]
            df_a['eta'] = legacy_data["etaa"].values[final_indices]
            df_a['eta_trans'] = legacy_data["eta_transa"].values[final_indices]
            df_a['Gamma_e'] = legacy_data["Gamma_ea"].values[final_indices]
            df_a['Gamma_e_trans'] = legacy_data["Gamma_e_transa"].values[final_indices]
            df_a['thick_winds_trans_condition'] = legacy_data["thick_winds_trans_condition"].values[final_indices]

            df_b['xsurf'] = legacy_data["xsurfb"].values[final_indices]
            df_b['eta'] = legacy_data["etab"].values[final_indices]
            df_b['eta_trans'] = legacy_data["eta_transb"].values[final_indices]
            df_b['Gamma_e'] = legacy_data["Gamma_eb"].values[final_indices]
            df_b['Gamma_e_trans'] = legacy_data["Gamma_e_transb"].values[final_indices]
            df_b['thick_winds_trans_condition'] = legacy_data["thick_winds_trans_condition"].values[final_indices]

        # Classify
        df_a = _classify_legacy_rows(df_a, WR_cond, Z_abs, use_orig_thick_winds)
        df_b = _classify_legacy_rows(df_b, WR_cond, Z_abs, use_orig_thick_winds)

        # --- FILTER 2: WR STARS ---
        if no_WR:
            # Use .values to avoid index alignment issues with duplicate indices
            a_is_wr = (df_a['stellar_type'].values == 'WR')
            b_is_wr = (df_b['stellar_type'].values == 'WR')

            # We want to KEEP systems where NEITHER is WR
            # i.e. NOT (a_is_wr OR b_is_wr)
            keep_mask = ~ (a_is_wr | b_is_wr)

            if not keep_mask.any():
                continue

            df_a = df_a.iloc[keep_mask]
            df_b = df_b.iloc[keep_mask]

        matched_groups.append(df_a)
        matched_groups.append(df_b)

    if not matched_groups:
        print("Warning: No matching BSE files found (or all filtered out).")
        return None

    final_df = pd.concat(matched_groups, ignore_index=True)
    print(f"Successfully matched and classified {len(final_df)} binary components.")
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
    plt.xlabel(r'log($T_{eff}$ [K])',fontsize=25)
    plt.ylabel(r'log($L$ [L$_{\odot}$])',fontsize=25)
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
    bse_dir="../BSE/Simulations",
    output_dir_mesa="./mesa_matched_chunks",
    output_dir_sse="./sse_matched_chunks",
    output_dir_bse="./bse_matched_chunks",
    chunk_size=10_000_000,
    wr_cond="any",
    mesa_cols=None,
    q_dist_type="uniform",
    p_dist_type="uniform",
    bse_grid_q=None,
    bse_grid_logp=None,
    no_WR=False,
    no_compactOBJ=True
):
    """
    Wrapper function to handle chunking, generating, or extracting MESA/SSE/BSE populations.

    Args:
        no_WR (bool): If True, filters out systems containing WR stars.
        no_compactOBJ (bool): If True, filters out systems with Compact Objects (K>9) for SSE/BSE.
                              Always filters K>14.
    """

    if bse_grid_q is None: bse_grid_q = DEFAULT_BSE_GRID_Q
    if bse_grid_logp is None: bse_grid_logp = DEFAULT_BSE_GRID_LOGP

    if mode == "MESA":
        output_dir = output_dir_mesa
    elif mode == "SSE":
        output_dir = output_dir_sse
    elif mode == "BSE":
        output_dir = output_dir_bse
    else:
        print(f"Error: Unknown mode '{mode}'. Choose 'MESA', 'SSE', or 'BSE'.")
        return None

    if mesa_cols is None:
        mesa_cols = [
            "star_age", "star_mass", "log_R", "log_Teff", "log_L",
            "log_Rzams", "log_Lzams", "log_Teffzams", "log_RMAX",
            "RMAX_logL", "RMAX_logTeff", "max_age"
        ]

    num_stars = len(q_df) if q_df is not None else 0

    if what_to_do == 'generate' and q_df is None:
        print("Error: 'generate' mode requires an input DataFrame (q_df).")
        return None

    df_stars = None

    # --- MESA LOGIC ---
    if mode == "MESA":
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

        elif num_stars < 1e7:
            os.makedirs(output_dir, exist_ok=True)
            df_mesa = match_to_mesa_sims_fast(q_df, mesa_dir, mesa_cols, WR_cond=wr_cond, no_WR=no_WR)

            if df_mesa is not None and not df_mesa.empty:
                output_path = os.path.join(output_dir, "matched_chunk.parquet")
                df_mesa.to_parquet(output_path)
                print(f"Chunk successfully processed and saved to {output_path}")
                df_stars = dd.read_parquet(output_path)

        else:
            os.makedirs(output_dir, exist_ok=True)
            num_chunks = int(np.ceil(num_stars / chunk_size))
            print(f"Processing {num_stars} stars in {num_chunks} chunks of size {chunk_size}...")

            for i in range(num_chunks):
                start_index = i * chunk_size
                end_index = (i + 1) * chunk_size
                df_chunk = q_df.iloc[start_index:end_index]

                print(f"\n--- Processing Chunk {i+1}/{num_chunks} ---")
                matched_chunk_df = match_to_mesa_sims_fast(df_chunk, mesa_dir, mesa_cols, WR_cond=wr_cond, no_WR=no_WR)

                if matched_chunk_df is not None and not matched_chunk_df.empty:
                    output_path = os.path.join(output_dir, f"matched_chunk_{i+1}.parquet")
                    matched_chunk_df.to_parquet(output_path)
                    print(f"Chunk {i+1} successfully processed and saved to {output_path}")

            print("\nAll chunks processed!")
            df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

    # --- SSE LOGIC ---
    elif mode == "SSE":
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

        elif num_stars < 1e7:
            os.makedirs(output_dir, exist_ok=True)
            df_sse = match_to_legacy_sims_fast(q_df, sse_dir, WR_cond=wr_cond, no_WR=no_WR, no_compactOBJ=no_compactOBJ)

            if df_sse is not None and not df_sse.empty:
                output_path = os.path.join(output_dir, "matched_chunk.parquet")
                df_sse.to_parquet(output_path)
                print(f"Processed and saved to {output_path}")
                df_stars = dd.read_parquet(output_path)

        else:
            os.makedirs(output_dir, exist_ok=True)
            num_chunks = int(np.ceil(num_stars / chunk_size))
            print(f"Processing {num_stars} stars in {num_chunks} chunks of size {chunk_size} using SSE tracks...")

            for i in range(num_chunks):
                start_index = i * chunk_size
                end_index = (i + 1) * chunk_size
                df_chunk = q_df.iloc[start_index:end_index]

                print(f"\n--- Processing Chunk {i+1}/{num_chunks} ---")
                matched_chunk_df = match_to_legacy_sims_fast(df_chunk, sse_dir, WR_cond=wr_cond, no_WR=no_WR, no_compactOBJ=no_compactOBJ)

                if matched_chunk_df is not None and not matched_chunk_df.empty:
                    output_path = os.path.join(output_dir, f"matched_chunk_{i+1}.parquet")
                    matched_chunk_df.to_parquet(output_path)
                    print(f"Chunk {i+1} successfully processed and saved to {output_path}")

            print("\nAll chunks processed!")
            df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

    # --- BSE LOGIC ---
    elif mode == "BSE":
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

        elif num_stars < 1e7:
            os.makedirs(output_dir, exist_ok=True)
            df_bse = match_to_legacy_sims_BSE_fast(
                q_df,
                bse_dir,
                WR_cond=wr_cond,
                q_dist_type=q_dist_type,
                p_dist_type=p_dist_type,
                q_grid=bse_grid_q,
                p_grid=bse_grid_logp,
                no_WR=no_WR,
                no_compactOBJ=no_compactOBJ
            )

            if df_bse is not None and not df_bse.empty:
                output_path = os.path.join(output_dir, "matched_chunk.parquet")
                df_bse.to_parquet(output_path)
                print(f"Processed and saved to {output_path}")
                df_stars = dd.read_parquet(output_path)

        else:
            os.makedirs(output_dir, exist_ok=True)
            num_chunks = int(np.ceil(num_stars / chunk_size))
            print(f"Processing {num_stars} stars in {num_chunks} chunks of size {chunk_size} using BSE tracks...")

            for i in range(num_chunks):
                start_index = i * chunk_size
                end_index = (i + 1) * chunk_size
                df_chunk = q_df.iloc[start_index:end_index]

                print(f"\n--- Processing Chunk {i+1}/{num_chunks} ---")
                matched_chunk_df = match_to_legacy_sims_BSE_fast(
                    df_chunk,
                    bse_dir,
                    WR_cond=wr_cond,
                    q_dist_type=q_dist_type,
                    p_dist_type=p_dist_type,
                    q_grid=bse_grid_q,
                    p_grid=bse_grid_logp,
                    no_WR=no_WR,
                    no_compactOBJ=no_compactOBJ
                )

                if matched_chunk_df is not None and not matched_chunk_df.empty:
                    output_path = os.path.join(output_dir, f"matched_chunk_{i+1}.parquet")
                    matched_chunk_df.to_parquet(output_path)
                    print(f"Chunk {i+1} successfully processed and saved to {output_path}")

            print("\nAll chunks processed!")
            df_stars = dd.read_parquet(os.path.join(output_dir, '*.parquet'))
            print(f"Dataset has {len(df_stars)} rows and {len(df_stars.columns)} columns.")

    return df_stars
