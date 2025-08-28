# StarEstate
A population synthesis tool to generate entire populations of single stars for the Milky Way or elliptical galaxies.

IMF is taken from a Kroupa's broken power law.
Lookback time, metallicity and galactic position distribution taken from Wagg et al. (2022), with the following changes:

- The metallicity distribution was renormalized to have a maximum of [Fe/H] = 0.5, i.e. roughly 3 times Solar metallicity
- Depending on the chosen galaxy shape:
  - **"elliptical"** --> nothing changes in the galactic distribution
  - **"spiral"** --> every star is given a probability to be on a specific Milky Way arm. The code also give each arm a specific angle distribution within a set angular scatter (**ANGULAR_SCATTER_DEG**) to position stars within them (Reid et al. 2019)
 
## Samplers and efficient calculations

Doing the **Inverse Transform Sampling** for the distrbutions of every population (MODE = 'generate_on_the_fly') requires extensive computational resources, since the code needs to numerically solve integrals for each star.

The solution around this problem is the use of samplers that are saved in a dedicated directory and that already contain all the numerical solutions for the given initial distribution, and then loading it each time a new population needs to be drawn. This reduces the computational time of a population ~10<sup>6</sup> stars in a spiral galaxy from hours to few minutes. In the directory, two samplers for a spiral galaxy, one for an IMF between 0.11 and 220 Msun, and another between 50 and 300 Msun, is delivered as an example.

The recommended route would be to run the first cell in Section 1 with 

> MODE = 'create_and_save_samplers'

To set the desired IMF, position, age, and metallicity samplers. Then every other time with

> MODE = 'load_samplers_and_generate'

To load the samplers to generate the population.

## Generate the population

With the generate_star_population function one will be able to draw a number of stars equal to the variable **NUM_STARS** from the given sampled distributions, while defining whether the thick, thin or both parts of the galactic disk are included and optionally saving the population as a CSV file. 
Here below the generated plots for a population of 10<sup>6</sup> stars (the different alpha angles representing different galactic arms)

With the use of pre-generated samplers, **it took less than a minute to generate all the stars**

<img width="761" height="1720" alt="image" src="https://github.com/user-attachments/assets/14fbf11e-b4ce-4e1d-afeb-45df449852d9" />

From then on the file can be loaded within Section 1.1

### Galaxy visualization

Without yet combining the drawn stars with MESA models, one could visualize in Section 1.2 the distribution of stars. Here below:

**The disk of a spiral galaxy**, with the yellow cross representing the galactic centre
<img width="1760" height="899" alt="image" src="https://github.com/user-attachments/assets/05c772d4-7a85-40b0-b3d3-7f5b24fad2d3" />


## Quantization

The masses and metallicities of stars drawn by the given distributions may not match the discrete sample of MESA tracks that one owns. 
For this reason, they must be binned to the closest match in the MESA sample in Section 1.3, where the available stellar masses and metallicities are listed.

$\textcolor{red}{Metallicities\ are\ expressed\ as\ [Fe/H]}$ and therefore are converted, following Bertelli et al. 1994a, from fraction of solar metallicity onto [Fe/H]

Finally, the quantized dataframe is saved in its respective folder and ready to use. Here below an example with only 3 MESA metallicities

<img width="1194" height="1167" alt="image" src="https://github.com/user-attachments/assets/70c5dfaa-76c1-4b5f-bb83-10ad3738934e" />

## MESA processing

In Section 1.3, the code combines MESA tracks with the drawn population, for which at each age it gives the simulated stellar parameters. Any object past its total lifetime is automatically discarded.
The code also give each star a stellar type following the Morgan-Keenan (MK) spectral classification, which uses temperature for spectral type (O, B, A, F, G, K, M) and luminosity/evolutionary state to distinguish between main-sequence dwarfs, giants, and supergiants. Additionally, the code determines whether a star is a YSG, RSG, HG, naked He, or Wolf-Rayet (WR) star. For the WR star phase, a series of different methods can be choosen, based on one or more conditions like surface-H abundance, closeness to the free-electron scattering Eddington limit, or the wind efficiency factor (see more in the upcoming paper Romagnolo+ 2025b).

**Currently the code takes metallicity folders that are named in the format "1Zsun", "0.5Zsun", etc. and then searches the different mass folders (e.g. "0.15", "20", "300"), where LOGS directories are, within them**.

If I will have time I will make it more customizable. **Just remember to follow this structure.**

The final outcome will be a dataframe with the entire population of stars, ready to analyze. 
See e.g. the synthetic population from the upcoming Romagnolo+ (2025b) below for different MESA models and WR conditions

<img width="1169" height="880" alt="image" src="https://github.com/user-attachments/assets/f72db2e4-7f65-4b90-8656-d762a287533b" />



