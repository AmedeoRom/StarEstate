# StarEstate
A population synthesis tool to generate entire populations of single stars for the Milky Way or elliptical galaxies.

IMF is taken from a Kroupa's broken power law.
Lookback time, metallicity and galactic position distribution taken from Wagg et al. (2022), with the following changes:

- The metallicity distribution was renormalized to have a maximum of [Fe/H] = 0.5, i.e. roughly 3 times Solar metallicity
- Depending on the chosen galaxy shape:
  - **"elliptical"** --> nothing changes in the galactic distribution
  - **"spiral"** --> every star is given a probability to be on a specific Milky Way arm. The code also give each arm a specific angle distribution within a set angular scatter (**ANGULAR_SCATTER_DEG**) to position stars within them
 
## Samplers and efficient calculations

Doing the **Inverse Transform Sampling** for the distrbutions of every population (MODE = 'generate_on_the_fly') requires extensive computational resources, since the code needs to numerically solve integrals for each star.

The solution around this problem is the use of samplers that are saved in a dedicated directory and that already contain all the numerical solutions for the given initial distribution, and then loading it each time a new population needs to be drawn. This reduces the computational time of a population ~10<sup>6</sup> stars in a spiral galaxy from hours to few minutes. In the directory, two samplers for a spiral galaxy, one for an IMF between 0.11 and 220 Msun, and another between 50 and 300 Msun, is delivered as an example.

The recommended route would be to run the first cell in Section 1 with 

> MODE = 'create_and_save_samplers'

To set the desired IMF, position, age, and metallicity samplers. Then every other time with

> MODE = 'load_samplers_and_generate'

To load the samplers to generate the population.

## Generate the population

With the generate_star_population function one will be able to draw a number of stars equal to the variable **NUM_STARS** from the given sampled distributions, while optionally saving the population as a CSV file. 
Here below the generated plots for a population of 10<sup>6</sup> stars (the different alpha angles representing different galactic arms)

<img width="761" height="1720" alt="image" src="https://github.com/user-attachments/assets/14fbf11e-b4ce-4e1d-afeb-45df449852d9" />

From then on the file can be loaded within Section 1.1

### Galaxy visualization

Without yet combining the drawn stars with MESA models, one could visualize in Section 1.2 the distribution of stars. Here below:

**the disk of a spiral galaxy**
<img width="1760" height="899" alt="image" src="https://github.com/user-attachments/assets/05c772d4-7a85-40b0-b3d3-7f5b24fad2d3" />


## Quantization

The masses and metallicities of stars drawn by the given distributions may not match the discrete sample of MESA tracks that one owns. 
For this reason, they must be binned to the closest match in the MESA sample in Section 1.4, where the available stellar masses and metallicities are listed.

<font color="red;">Metallicities are expressed as [Fe/H]</font>
