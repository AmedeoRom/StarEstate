import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import random
import pickle
from scipy import integrate, optimize
from scipy.stats import norm, truncnorm
from scipy.interpolate import interp1d
from functools import lru_cache

#==================================================================
# OPTIMIZATION: INVERSE CDF SAMPLER SETUP
#==================================================================
def setup_inverse_cdf_sampler(cdf_function, domain_min, domain_max, args=(), resolution=100000):
    """
    Pre-calculates the inverse CDF to create a fast interpolation function for sampling.
    This is a performance optimization to avoid repeated numerical integration.
    """
    u_values = np.linspace(0.0000001, 0.9999999, resolution)
    def inverse_cdf_slow(u):
        # Find the value x such that CDF(x) = u
        return optimize.root_scalar(lambda x: cdf_function(x, *args) - u, bracket=[domain_min, domain_max]).root
    x_values = [inverse_cdf_slow(u) for u in u_values]
    # Create an interpolation function from the pre-calculated values
    sampler = interp1d(u_values, x_values, bounds_error=False, fill_value=(x_values[0], x_values[-1]))
    return sampler

#==================================================================
# SAMPLER SETUP FUNCTIONS FOR EACH PARAMETER
#==================================================================
def setup_high_alpha_r_sampler(R_min=1.5, R_max=20.0, resolution=10000):
    """Sets up the fast sampler for the high-alpha (thick disk) radial distribution."""
    print("Setting up high-alpha (thick disk) radial sampler...")
    # Exponential disk model for the thick disk component
    Rdistr = lambda R: np.exp(-R / (1/0.43)) * R / (1/0.43)**2
    normalization_constant, _ = integrate.quad(Rdistr, R_min, R_max)
    Rdistr_normalized = lambda R: Rdistr(R) / normalization_constant
    CDF_R = lambda R: integrate.quad(Rdistr_normalized, R_min, R)[0]
    return setup_inverse_cdf_sampler(CDF_R, R_min, R_max, resolution=resolution)

def setup_age_sampler(age_min_Gyr, age_max_Gyr, resolution=10000):
    """Sets up the fast sampler for a specific age distribution range based on a star formation history."""
    # The star formation history (SFH) is a function of LOOKBACK time (which is equivalent to AGE)
    SFH = lambda t: np.exp(-(12.0 - t) / 6.8)

    t_min, t_max = age_min_Gyr, age_max_Gyr

    norm_const, _ = integrate.quad(SFH, t_min, t_max)
    if norm_const <= 0: return lambda u: np.full_like(u, t_min) # Handle cases with no star formation

    sfh_normalized = lambda t: SFH(t) / norm_const
    CDF_t = lambda t: integrate.quad(sfh_normalized, t_min, t)[0]

    return setup_inverse_cdf_sampler(CDF_t, t_min, t_max, resolution=resolution)

def setup_v_sampler(z_min=0, z_max=10, resolution=10000):
    """Sets up the sampler for vertical distance from the Galactic plane."""
    Vdistr = lambda z: (1 / 0.95) * np.exp(-np.abs(z) / 0.95)
    V_norm_const, _ = integrate.quad(Vdistr, z_min, z_max)
    Vdistr_normalized = lambda z: Vdistr(z) / V_norm_const
    CDF_V = lambda z: integrate.quad(Vdistr_normalized, z_min, z)[0]
    return setup_inverse_cdf_sampler(CDF_V, z_min, z_max, resolution=resolution)

def setup_imf_sampler(Mzams_min, Mzams_max, resolution=10000):
    """Sets up the sampler for the Initial Mass Function (IMF) using a broken power-law."""
    alpha1, alpha2, alpha3 = -1.3, -2.2, -2.7
    m_break1, m_break2 = 0.5, 1.0
    def integral_power_law(m, p):
        m = np.maximum(m, 1e-9)
        return (m**(p + 1)) / (p + 1)
    A1_unnorm = 1.0
    A2_unnorm = A1_unnorm * (m_break1**alpha1) / (m_break1**alpha2)
    A3_unnorm = A2_unnorm * (m_break2**alpha2) / (m_break2**alpha3)
    integral1, integral2, integral3 = 0, 0, 0
    if Mzams_min < m_break1:
        upper = min(m_break1, Mzams_max)
        integral1 = A1_unnorm * (integral_power_law(upper, alpha1) - integral_power_law(Mzams_min, alpha1))
    if Mzams_min < m_break2 and Mzams_max >= m_break1:
        lower = max(m_break1, Mzams_min)
        upper = min(m_break2, Mzams_max)
        integral2 = A2_unnorm * (integral_power_law(upper, alpha2) - integral_power_law(lower, alpha2))
    if Mzams_max >= m_break2:
        lower = max(m_break2, Mzams_min)
        integral3 = A3_unnorm * (integral_power_law(Mzams_max, alpha3) - integral_power_law(lower, alpha3))
    total_norm = integral1 + integral2 + integral3
    A1, A2, A3 = A1_unnorm/total_norm, A2_unnorm/total_norm, A3_unnorm/total_norm
    c1 = 0
    if Mzams_min < m_break1:
        c1 = A1 * (integral_power_law(min(m_break1, Mzams_max), alpha1) - integral_power_law(Mzams_min, alpha1))
    c2 = 0
    if Mzams_min < m_break2 and Mzams_max >= m_break1:
        c2 = A2 * (integral_power_law(min(m_break2, Mzams_max), alpha2) - integral_power_law(max(m_break1, Mzams_min), alpha2))
    def cdf_m(m):
        if m < Mzams_min: return 0.0
        if m < m_break1:
            return A1 * (integral_power_law(m, alpha1) - integral_power_law(Mzams_min, alpha1))
        elif m < m_break2:
            return c1 + A2 * (integral_power_law(m, alpha2) - integral_power_law(max(m_break1, Mzams_min), alpha2))
        else:
            return c1 + c2 + A3 * (integral_power_law(m, alpha3) - integral_power_law(max(m_break2, Mzams_min), alpha3))
    return setup_inverse_cdf_sampler(cdf_m, Mzams_min, Mzams_max, resolution=resolution)

def setup_age_dependent_r_sampler(R_min=1.5, R_max=20.0, age_bins=np.arange(0, 14, 0.5), resolution=10000):
    """Sets up binned radial samplers for the low-alpha (thin disk) population, where the scale length evolves with age."""
    print(f"Setting up binned low-alpha (thin disk) radial samplers for age bins: {age_bins}")
    # Model for how the radial scale length of the thin disk changes with stellar age
    Rdistr_lowZ = lambda R, age: np.exp(-R / (4*(1-0.3*age/8))) * R / (4*(1-0.3*age/8))**2 if (4*(1-0.3*age/8)) > 0 else 0
    sampler_dict = {}
    for i in range(len(age_bins)):
        age_mid_bin = age_bins[i] + (age_bins[1]-age_bins[0])/2 if i+1 < len(age_bins) else age_bins[i]
        norm_const, _ = integrate.quad(Rdistr_lowZ, R_min, R_max, args=(age_mid_bin,))
        if norm_const <= 0: continue
        Rdistr_normalized = lambda R: Rdistr_lowZ(R, age_mid_bin) / norm_const
        CDF_R = lambda R: integrate.quad(Rdistr_normalized, R_min, R)[0]
        sampler_dict[age_bins[i]] = setup_inverse_cdf_sampler(CDF_R, R_min, R_max, resolution=resolution)
    return sampler_dict, age_bins

#==================================================================
# SAMPLER CONTAINER CLASS
#==================================================================
class SamplerContainer:
    """A container class to hold all the pre-calculated sampler objects for efficient reuse."""
    def __init__(self, Mzams_min, Mzams_max, resolution=10000):
        print(f"Setting up all samplers with resolution={resolution}... (This may take a moment)")
        # Create samplers for each component based on their age (lookback time) ranges
        self.age_sampler_thin = setup_age_sampler(age_min_Gyr=0.0, age_max_Gyr=8.0, resolution=resolution)
        self.age_sampler_thick = setup_age_sampler(age_min_Gyr=8.0, age_max_Gyr=12.0, resolution=resolution)
        self.imf_sampler = setup_imf_sampler(Mzams_min, Mzams_max, resolution=resolution)
        self.v_sampler = setup_v_sampler(resolution=resolution)
        self.r_sampler_low_alpha_dict, self.r_age_bins = setup_age_dependent_r_sampler(resolution=resolution)
        self.r_sampler_high_alpha = setup_high_alpha_r_sampler(resolution=resolution)
        print("Samplers created successfully.")

def save_samplers(sampler_obj, filepath):
    """Saves the sampler container object to a file using pickle."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'wb') as f: pickle.dump(sampler_obj, f)
    print(f"Samplers saved to {filepath}")

def load_samplers(filepath):
    """Loads a sampler container object from a file."""
    with open(filepath, 'rb') as f: sampler_obj = pickle.load(f)
    print(f"Samplers loaded from {filepath}")
    return sampler_obj

#==================================================================
# FAST DRAWING FUNCTIONS (called during generation)
#==================================================================
def draw_samples_R_binned(ages, sampler_dict, age_bins):
    """Draws radial samples for stars based on their age, using the pre-binned samplers."""
    bin_indices = np.digitize(ages, age_bins) - 1
    bin_indices = np.maximum(0, bin_indices)
    corresponding_bins = age_bins[bin_indices]
    random_uniforms = np.random.rand(len(ages))
    radii = [sampler_dict[bin_val](u) for bin_val, u in zip(corresponding_bins, random_uniforms)]
    return np.array(radii)

def draw_spiral_arm_positions(radii, ages, arm_params_list, arm_membership_prob):
    """
    IMPROVED: Assigns stars to spiral arms based on their radius and age.
    This function now incorporates the key finding from the research: the arm-to-interarm
    density contrast is a strong function of stellar age. Younger stars are more
    tightly confined to the arms (lower scatter), while older stars are more dispersed.
    """
    num_stars = len(radii)
    arm_names = [arm['name'] for arm in arm_membership_prob]
    weights = [arm['probability'] for arm in arm_membership_prob]
    chosen_arm_names = random.choices(arm_names, weights=weights, k=num_stars)

    final_angles_deg = []

    # --- NEW: Age-dependent angular scatter ---
    # Young stars (<100 Myr) are "dynamically cold" -> low scatter
    # Intermediate stars are "warmer" -> medium scatter
    # Old stars (>1 Gyr) are "dynamically hot" -> high scatter, less confined to arms
    min_scatter_deg = 2.0  # For the youngest stars, tightly confined
    max_scatter_deg = 15.0 # For the oldest stars, very dispersed

    # Simple linear mapping from age to scatter. More complex functions could be used.
    angular_scatters_deg = min_scatter_deg + (ages / 2.0) * (max_scatter_deg - min_scatter_deg)

    for i in range(num_stars):
        R, arm_name, scatter_deg = radii[i], chosen_arm_names[i], angular_scatters_deg[i]

        # Find the correct segment of the arm for the star's given radius
        chosen_segment = arm_params_list[arm_name]

        # Calculate the mean angle of the arm at that radius using the logarithmic spiral equation
        # R = R_ref * exp((phi - phi_ref) * tan(p))
        # which rearranges to: phi = phi_ref + ln(R / R_ref) / tan(p)
        p_rad = math.radians(chosen_segment['pitch_angle_deg'])
        R_ref = chosen_segment['ref_radius_kpc']
        phi_ref_rad = math.radians(chosen_segment['ref_angle_deg'])

        try:
            # The term 1/tan(p) defines how tightly wound the arm is.
            mean_theta_rad = phi_ref_rad + (math.log(R / R_ref)) / math.tan(p_rad)
        except (ValueError, ZeroDivisionError):
            # Fallback for invalid math operations (e.g., log(negative number))
            mean_theta_rad = random.uniform(0, 2 * math.pi)

        # Apply the age-dependent scatter to the mean angle
        scatter_rad = math.radians(scatter_deg)
        final_theta_rad = random.normalvariate(mean_theta_rad, scatter_rad)
        final_angles_deg.append(math.degrees(final_theta_rad))

    return chosen_arm_names, np.array(final_angles_deg)


def draw_samples_FeH(Rs, ages, tm=12.0, feh_scatter_std=0.1, FeH_min=-3.0, FeH_max=0.5):
    """Draws metallicity ([Fe/H]) samples based on radius and age."""
    t = np.clip(ages, 0, tm)
    Fm, NablaFeH, RnowFeH, gammaFeH = -1.0, 0.075, 8.7, 0.3
    ft = (1 - t / tm)**gammaFeH
    mean_feh = Fm + NablaFeH * Rs - (Fm + NablaFeH * RnowFeH) * ft
    a = (FeH_min - mean_feh) / feh_scatter_std
    b = (FeH_max - mean_feh) / feh_scatter_std
    return truncnorm.rvs(a, b, loc=mean_feh, scale=feh_scatter_std)

#==================================================================
# MAIN GENERATOR FUNCTION
#==================================================================
def generate_star_population(num_stars, samplers, disk_type="all", galaxy_shape="elliptical",
                             arm_params_list=None, arm_membership_prob=None, seed=None,
                             FeH_min=-3.0, FeH_max=0.5):
    print(f"Generating {num_stars} stars for the '{disk_type}' disk model...")

    if disk_type == "thin":
        ages = samplers.age_sampler_thin(np.random.rand(num_stars))
        population_type = np.full(num_stars, "thin_disk")
    elif disk_type == "thick":
        ages = samplers.age_sampler_thick(np.random.rand(num_stars))
        population_type = np.full(num_stars, "thick_disk")
    elif disk_type == "all":
        # Approximate ratio of thin to thick disk stars
        num_thin = int(num_stars * 0.85)
        num_thick = num_stars - num_thin
        ages_thin = samplers.age_sampler_thin(np.random.rand(num_thin))
        ages_thick = samplers.age_sampler_thick(np.random.rand(num_thick))
        ages = np.concatenate([ages_thin, ages_thick])
        population_type = np.concatenate([np.full(num_thin, "thin_disk"), np.full(num_thick, "thick_disk")])
        # Shuffle the combined arrays to mix the populations
        p = np.random.permutation(num_stars)
        ages, population_type = ages[p], population_type[p]
    else:
        raise ValueError("disk_type must be 'thin', 'thick', or 'all'")

    masses = samplers.imf_sampler(np.random.rand(num_stars))
    vertical_distances = samplers.v_sampler(np.random.rand(num_stars))

    radii = np.zeros(num_stars)
    thin_mask = (population_type == "thin_disk")
    thick_mask = (population_type == "thick_disk")

    if np.any(thin_mask):
        radii[thin_mask] = draw_samples_R_binned(ages[thin_mask], samplers.r_sampler_low_alpha_dict, samplers.r_age_bins)
    if np.any(thick_mask):
        radii[thick_mask] = samplers.r_sampler_high_alpha(np.random.rand(np.sum(thick_mask)))

    metallicities = draw_samples_FeH(radii, ages, FeH_min=FeH_min, FeH_max=FeH_max)

    if galaxy_shape == "spiral":
        # Pass ages to the spiral arm drawing function
        arm_positions, alpha_angles = draw_spiral_arm_positions(radii, ages, arm_params_list, arm_membership_prob)
    else:
        arm_positions = ["none"] * num_stars
        alpha_angles = np.random.uniform(0, 360, size=num_stars)

    print("Generation complete.")

    data_dict = {
        "Mzams": masses, "Age [Gyr]": ages, "Fe/H": metallicities,
        "Radial Distance [kpc]": radii, "Vertical Distance [kpc]": vertical_distances,
        "Arm Position": arm_positions, "alpha_angle_deg": alpha_angles,
        "population_type": population_type, "random_seed": np.full(num_stars, seed),
        "disk_type": disk_type
    }
    return pd.DataFrame(data_dict)



#==================================================================
# DATA QUANTIZATION FUNCTION
#==================================================================
def quantize_population(df, mass_grid, feh_grid, in_place=True):

    """
    Snaps the mass and metallicity of a stellar population to the closest
    values in provided grids using a fast, vectorized NumPy approach.

    Args:
        df (pd.DataFrame): The input DataFrame with 'Mzams' and 'Fe/H' columns.
        mass_grid (np.ndarray): Sorted grid of mass values.
        feh_grid (np.ndarray): Sorted grid of metallicity values.
        in_place (bool): If True, modifies the original DataFrame to save memory.
                         If False, works on a copy. Defaults to False.
    """
    print("Quantizing mass and metallicity using vectorized NumPy...")

    """
    Snaps the mass and metallicity of a stellar population to the closest
    values in provided grids.
    """
    quantized_df = df if in_place else df.copy()

        # --- Vectorized quantization logic ---
    def find_nearest(grid, values):
        """Finds the nearest grid value for each value in an array."""
        # Ensure grid is sorted (np.searchsorted requires this)
        grid = np.sort(grid)

        # Find the insertion indices for all values at once
        # 'left' means if a value is equal to a grid point, the index of that grid point is returned
        indices = np.searchsorted(grid, values, side='left')

        # Handle edge cases: values smaller than the first grid point
        indices = np.maximum(indices, 1)
        # Handle edge cases: values larger than the last grid point
        indices = np.minimum(indices, len(grid) - 1)

        # Find the two nearest grid points for each value
        left_neighbors = grid[indices - 1]
        right_neighbors = grid[indices]

        # Return the closer of the two neighbors for each value
        return np.where(
            np.abs(values - left_neighbors) <= np.abs(values - right_neighbors),
            left_neighbors,
            right_neighbors
        )

    quantized_df['Mzams'] = find_nearest(mass_grid, quantized_df['Mzams'].values)
    quantized_df['Fe/H'] = find_nearest(feh_grid, quantized_df['Fe/H'].values)


    print("Quantization complete.")
    return quantized_df

#==================================================================
# PLOTTING FUNCTION
#==================================================================
def plot_distributions(df, disk_type="thin"):
    print("\nPlotting distributions...")
    plot_columns = [col for col in df.columns if pd.api.types.is_numeric_dtype(df[col]) and col != 'random_seed']
    plt.figure(figsize=(12, 5 * len(plot_columns)))
    for i, column in enumerate(plot_columns):
        plt.subplot(len(plot_columns), 1, i + 1)
        bins = 100 if column == "Mzams" else 50
        plt.hist(df[column], bins=bins, density=True, alpha=0.7, label=f'Sampled {column}')
        plt.title(f'Distribution of {column} ({disk_type.capitalize()} Disk)')
        if column == 'Mzams': plt.xlabel("$M_\\mathrm{ZAMS}$ [$M_\\odot$]"); plt.yscale('log')
        else: plt.xlabel(column)
        plt.ylabel('Probability Density'); plt.legend(); plt.grid(True, linestyle='--', alpha=0.6)
    plt.subplots_adjust(hspace=0.5)
    if 'Arm Position' in df.columns and df['Arm Position'].iloc[0] != 'none':
        plt.figure(figsize=(10, 6)); df['Arm Position'].value_counts().plot(kind='bar', alpha=0.8)
        plt.title('Star Counts per Spiral Arm'); plt.ylabel('Number of Stars'); plt.xticks(rotation=45, ha='right'); plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout(); plt.show()

#==================================================================
# SCRIPT EXECUTION
#==================================================================
if __name__ == "__main__":
    # --- Workflow Configuration ---
    MODE = 'generate_on_the_fly'  # Options: 'generate_on_the_fly', 'create_and_save_samplers', 'load_samplers_and_generate'
    SAMPLER_FILE = 'output/samplers/thin_disk_samplers.pkl'
    SEED = 42  # Set to None for a random run

    # --- Simulation Configuration ---
    NUM_STARS = 10000
    DISK_TYPE = "thin"
    GALAXY_SHAPE = "elliptical"
    MZAMS_MIN, MZAMS_MAX = 0.08, 150.0
    T_MAX_AGE = 13.6
    SAMPLER_RESOLUTION = 10000 # Controls accuracy vs. setup time for samplers

    # --- Spiral Arm Configuration ---
    SFR_PROBABILITY_BY_ARM = [
    {'name': 'Scutum-Centaurus', 'min_probability_pct': 35, 'max_probability_pct': 45},
    {'name': 'Perseus', 'min_probability_pct': 25, 'max_probability_pct': 35},
    {'name': 'Local (Orion) Arm', 'min_probability_pct': 10, 'max_probability_pct': 15},
    {'name': 'Sagittarius', 'min_probability_pct': 5, 'max_probability_pct': 10},
    {'name': 'Norma', 'min_probability_pct': 5, 'max_probability_pct': 10}]
    ARM_PARAMS = [
    {'name': 'Scutum-Centaurus', 'a_kpc': 3.5, 'b': 1/math.tan(math.radians(12)), 'theta0_rad': math.radians(20)},
    {'name': 'Perseus', 'a_kpc': 4.9, 'b': 1/math.tan(math.radians(12)), 'theta0_rad': math.radians(150)},
     {'name': 'Local (Orion) Arm', 'a_kpc': 4.0, 'b': 1/math.tan(math.radians(10)), 'theta0_rad': math.radians(100)},
      {'name': 'Sagittarius', 'a_kpc': 4.0, 'b': 1/math.tan(math.radians(12)), 'theta0_rad': math.radians(280)},
      {'name': 'Norma', 'a_kpc': 3.0, 'b': 1/math.tan(math.radians(12)), 'theta0_rad': math.radians(200)}]
    ANGULAR_SCATTER_DEG = 5.0

    # --- Seed the random number generators for reproducibility ---
    if SEED is not None:
        np.random.seed(SEED)
        random.seed(SEED)

    # --- Main Workflow Logic ---
    samplers = None
    if MODE == 'create_and_save_samplers':
        samplers = SamplerContainer(DISK_TYPE, MZAMS_MIN, MZAMS_MAX, resolution=SAMPLER_RESOLUTION)
        save_samplers(samplers, SAMPLER_FILE)

    else: # For 'generate_on_the_fly' or 'load_samplers_and_generate'
        if MODE == 'load_samplers_and_generate':
            if not os.path.exists(SAMPLER_FILE):
                raise FileNotFoundError(f"Sampler file not found: {SAMPLER_FILE}. Please run in 'create_and_save_samplers' mode first.")
            samplers = load_samplers(SAMPLER_FILE)
        else: # generate_on_the_fly
            samplers = SamplerContainer(DISK_TYPE, MZAMS_MIN, MZAMS_MAX, resolution=SAMPLER_RESOLUTION)

        population_df = generate_star_population(
            num_stars=NUM_STARS,
            samplers=samplers,
            disk_type=DISK_TYPE,
            galaxy_shape=GALAXY_SHAPE,
            arm_params_list=ARM_PARAMS,
            sfr_probs_list=SFR_PROBABILITY_BY_ARM,
            angular_scatter_deg=ANGULAR_SCATTER_DEG,
            t_max_age=T_MAX_AGE,
            seed=SEED
        )
        print("\nGenerated Population Head:")
        print(population_df.head())
        print("\nGenerated Population Stats:")
        print(population_df.describe())

        plot_distributions(population_df, disk_type=DISK_TYPE)


    # --- Quantization Configuration ---
    QUANTIZE_DATA = True

    # 1. Generate the stellar population
    population_df = generate_star_population(
        num_stars=NUM_STARS,
        disk_type=DISK_TYPE,
        Mzams_min=MZAMS_MIN,
        Mzams_max=MZAMS_MAX,
        galaxy_shape=GALAXY_SHAPE,
        arm_params_list=ARM_PARAMS,
        sfr_probs_list=SFR_PROBABILITY_BY_ARM,
        angular_scatter_deg=ANGULAR_SCATTER_DEG
    )

    # 2. Optionally save the original data
    if SAVE_TO_CSV:
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)
        filepath = os.path.join(OUTPUT_DIR, f"star_population_{DISK_TYPE}_{NUM_STARS}.csv")
        population_df.to_csv(filepath, index=False)
        print(f"Original data saved to {filepath}")

    # 3. Optionally quantize the data
    quantized_df = None
    if QUANTIZE_DATA:
        quantized_df = quantize_population(
            population_df,
            mass_grid=stellar_mass_grid,
            feh_grid=metallicity_grid
        )
        if SAVE_TO_CSV:
            quantized_filepath = os.path.join(OUTPUT_DIR, f"star_population_quantized_{DISK_TYPE}_{NUM_STARS}.csv")
            quantized_df.to_csv(quantized_filepath, index=False)
            print(f"Quantized data saved to {quantized_filepath}")

    # 4. Optionally perform MESA matching
    if PERFORM_MESA_MATCHING and quantized_df is not None:
        final_matched_df = match_to_mesa_sims(
            quantized_df,
            base_path=MESA_BASE_PATH,
            mesa_cols_to_extract=MESA_COLS_TO_EXTRACT
        )
        if final_matched_df is not None:
            print("\nFinal matched data with MESA parameters:")
            print(final_matched_df.head())
            if SAVE_TO_CSV:
                final_filepath = os.path.join(OUTPUT_DIR, f"final_matched_data_{DISK_TYPE}_{NUM_STARS}.csv")
                final_matched_df.to_csv(final_filepath, index=False)
                print(f"Final matched data saved to {final_filepath}")

    # 5. Plot the original distributions
    plot_distributions(population_df, disk_type=DISK_TYPE)
