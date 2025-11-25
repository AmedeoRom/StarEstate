# StarEstate
A population synthesis tool to generate entire populations of single stars for the Milky Way or elliptical galaxies.

IMF is taken from a Kroupa's broken power law.
Lookback time, metallicity and galactic position distribution taken from [Wagg et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...937..118W/abstract), with the following changes:

- The metallicity distribution was renormalized to have a maximum of [Fe/H] = 0.5, i.e. roughly 3 times Solar metallicity
- Depending on the chosen galaxy shape:
  - **"elliptical"** --> nothing changes in the galactic distribution. Uniform distribution of radial angles.
  - **"spiral"** --> every star is given a probability to be on a specific Milky Way arm. Deteails on spiral arms modelling at the end


## Samplers and efficient calculations

Doing the **Inverse Transform Sampling** for the distrbutions of every population (MODE = 'generate_on_the_fly') requires extensive computational resources, since the code needs to numerically solve integrals for each star.

The solution around this problem is the use of samplers that are saved in a dedicated directory and that already contain all the numerical solutions for the given initial distribution, and then loading it each time a new population needs to be drawn. This reduces the computational time of a population ~10<sup>6</sup> stars in a spiral galaxy from hours to few minutes. In the directory, a samplers for a spiral galaxy for an IMF between 0.11 and 100 $M_\odot$ and a resolution of 25,000 bins (across 4 different independent dimensions this represents $25,000^4$ $\approx$ 3.9 x $10^{17}$ combinations) is delivered as an example.

**If samplers for a specific mass range and resolution are not saved yet**, the recommended route would be to run the first cell in Section 1 with 

> MODE = 'create_and_save_samplers'

To set the desired IMF, position, age, and metallicity samplers.
Once the samplers are saved, then every other time the code can be run with

> MODE = 'load_samplers_and_generate'

To load the samplers to generate the population.

## Generate the population

### Spiral Galaxy Modelling

The code gives a membership probability for each star to be within one of the major Milky Way arms: Scutum-Centaurus, Perseus, Sagittarius-Carina, Norma-Outer, Local (Orion Spur). This is fully a guess work from the literature and model variations are invited. For each arm, as shown in the table below, various parameters from the literature are adopted to generate a distribution of angles as a function of age and radial positions.

In particular, the code defines introduces the concept of dynamical temperature, or velocity dispertion, which is why different stellar populations trace the spiral arms with varying fidelity (e.g. [Mackereth+ 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489..176M/abstract)):

- **Young Tracers (Age < 100 Myr)**: Populations with very low velocity dispersion are considered dynamically "cold." This includes molecular gas clouds, HII regions, and the massive, young O and B-type stars recently formed from them. These components have orbits that are nearly circular and are easily perturbed by the gravitational potential of the spiral arms. As a result, they are tightly confined to the arms.
- **Intermediate-Age Tracers (Age ~100 Myr – 1 Gyr)**: Stars in this age range, such as Classical Cepheids and young open clusters, are dynamically "warmer." Over their lifetimes, they have had time to drift from their birth sites within the arms. Their velocity dispersions have increased due to gravitational scattering off molecular clouds and the spiral arms themselves. Consequently, while they still show a clear concentration in the arms, their distribution is more dispersed than that of the youngest stars.
- **Old Tracers (Age > 1 Gyr)**: The general population of older disk stars. They are dynamically "hot," with large velocity dispersions. Their orbits are less affected by the relatively weak perturbation of the spiral arms. As a result, they trace a much smoother, lower-contrast pattern that reflects only the deepest parts of the gravitational potential well—the dominant two-armed stellar pattern.

| Arm Name | Pitch Angle (p) [deg] | Reference Radius [kpc] | Reference Angle β<sub>kink</sub> [deg] | Radial Range | Arm Membership Probability
|----------|----------|----------| ----------| ----------| ----------|
| Scutum-Centaurus | 2.0<sup>1</sup> |	3.14<sup>1, 2</sup> |	25<sup>1</sup>	| 3.0 - 16.0<sup>1, 3</sup> | 0.30 |
| Sagittarius-Carina | 13.1<sup>4</sup>	| 4.93<sup>4</sup> |	-45<sup>4</sup>	| 4.0 - 16.0<sup>1, 3</sup> | 0.25 |
| Perseus | 9.5<sup>5</sup>	| 9.94<sup>6</sup>	| 150<sup>1</sup>	| 6.0 - 18.0<sup>1</sup> | 0.20 |
| Sagittarius-Carina | 13.0<sup>4</sup> |	4.0<sup>7</sup>	| -100<sup>1</sup>	| 3.5 - 20<sup>1, 6, 8</sup> | 0.15 |
| Local (Orion Spur) | 10.1<sup>1</sup>	| 8.15<sup>9</sup>	| 0<sup>10</sup>	| 6.0 - 9.0<sup>11</sup> | 0.10 |

1. [Reid + (2019)](https://iopscience.iop.org/article/10.3847/1538-4357/ab4a11)
2. Starting point consistent with models where arms begin a few kpc from the center
3. It starts near the end of the central bar and extends far into the disk
4. [Vallée (2017)](https://arxiv.org/pdf/1711.05228)
5. Average of two values around the curvature kink in [Reid+ (2019)](https://iopscience.iop.org/article/10.3847/1538-4357/ab4a11)
6. [Bobylev & Bajkova (2014)](https://academic.oup.com/mnras/article/437/2/1549/1105964?login=false)
7. [Reid et al. (2019)](https://iopscience.iop.org/article/10.3847/1538-4357/ab4a11) place at a "kink" radius of 4.46 kpc
8. Combination of the inner Norma arm with the very distant Outer arm, which extends beyond the currently mapped regions of the Galaxy
9. Sun Distance
10. Galactocentric coordinate system
11. [Xu+ (2016)](https://www.science.org/doi/10.1126/sciadv.1600878)

### Code use and features

With the generate_star_population function one will be able to draw a number of stars equal to the variable **NUM_STARS** from the given sampled distributions, while defining whether the thick, thin or both parts of the galactic disk are included and optionally saving the population as a CSV file.  

With a good laptop, the samplers with a resolution of 25,000 bins were calculated and saved in less than an hour.
While using them, **it took 30 seconds to generate 10<sup>6</sup> stars, 4 minutes to generate 10<sup>7</sup>, and 40 minutes to generate 10<sup>8</sup>**

Here below the generated plots for a population of 10<sup>6</sup> stars.

<img width="784" height="1720" alt="Distributions" src="https://github.com/user-attachments/assets/4ddc1e8c-c610-4779-b1cc-078f27797ee2" />


From then on the file can be loaded within Section 1.1

### Galaxy morphology visualization

Without yet combining the drawn stars with evolutionary models, one could visualize in Section 1.2 the distribution of stars with a dynamical size adjustment for the scatter point sizes depending on the size of the stellar sample.
*Making a scatter plot of many data points can considerably slow the code.* For $10^8$ stars a good laptop can take nearly an hour.

> **Resolution comparison: [Click here to compare 1M vs 100M Star plots in Results.md](Results.md#1-galaxy-positions)**


## Quantization

The masses and metallicities of stars drawn by the given distributions likely don't match the discrete sample of stellar tracks that one owns. 
For this reason, they must be binned to the closest match in the star sample in Section 1.3, where the available stellar masses and metallicities are listed.

$\textcolor{red}{Metallicities\ are\ expressed\ as\ [Fe/H]}$ and therefore are converted, following Bertelli et al. 1994a, from fraction of solar metallicity to [Fe/H]

Finally, the quantized dataframe is saved in its respective folder and ready to use. Here below an example with only 3 initial metallicities

<img width="1651" height="1654" alt="quantization" src="https://github.com/user-attachments/assets/cb6e186c-d384-42c2-9d5e-7738af63a944" />


## Stellar tracks processing

<!--Example stellar tracks from [Romagnolo+ (2024)](https://ui.adsabs.harvard.edu/abs/2024ApJ...964L..23R/abstract) that can be used for a tutorial are available here **(link not yet available; still under review in Zenodo)**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17669656.svg)](https://doi.org/10.5281/zenodo.17669656).
Please download the `MESA_tracks.zip` file, extract it, and place the contents inside a `MESA Simulations` directory.-->
MESA tracks for this tutorial will be available soon after the submission of the code paper.
SSE tracks from the [Romagnolo+ (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..706R/abstract) StarTrack models can be found in the SSE folder. 
Both sets of tracks were taken until the simulated star becomes a white dwarf, exceeds 100 times the Hubble time, or ends its core-helium burning phase.

In Section 1.3, the code combines stellar tracks with the drawn population, for which at each age it gives the simulated stellar parameters. Any object past its total lifetime is automatically discarded.
The code also give each star a stellar type following the Morgan-Keenan (MK) spectral classification, which uses temperature for spectral type (O, B, A, F, G, K, M) and luminosity/evolutionary state to distinguish between main-sequence dwarfs, giants, and supergiants. Additionally, the code determines whether a star is a YSG, RSG, HG, naked He, or Wolf-Rayet (WR) star. For the WR star phase, a series of different methods can be choosen, based on one or more conditions like surface-H abundance, closeness to the free-electron scattering Eddington limit, or the wind efficiency factor (see more in the upcoming paper Romagnolo+ 2025b).

**Currently the code takes metallicity folders that are named in the format "1Zsun", "0.5Zsun", etc. and then searches the different mass folders (e.g. "0.15", "20", "300"), where LOGS directories are, within them. Remember to follow this structure.**

The final outcome will be a dataframe with the entire population of stars, ready to be analyzed. 
Integrated in the file, an Hertzprung-Russel diagram plotter is included to show the metallicity and mass distribution of the retrieved stars.
<!-- **Resolution and code comparison: [View the SSE vs. MESA benchmarks in Results.md](Results.md#2-hertzsprung-russell-hr-diagrams)**-->
> **Resolution comparison: [View the SSE benchmarks in Results.md](Results.md#2-hertzsprung-russell-hr-diagrams)**

