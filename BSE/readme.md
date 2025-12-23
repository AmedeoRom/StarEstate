# StarEstate BSE Wrapper

This repository contains a Python-driven wrapper for the **StarTrack** rapid binary evolution code. This setup is designed to generate deterministic evolutionary tracks for a specific **grid** of binary parameters.

## Prerequisites

To run this code, you need a standard Linux/Unix environment with a **GCC Compiler** (The Python script invokes `gcc` to compile the C code on the fly) and Standard C libraries (specifically `libm` for math operations).

## How to Run the Code

1.  **Prepare the Environment**: Ensure all files (`run_startrack.py`, `binary.c`, `singl.c`, `sinbin.h`) are in the same directory.
2.  **Execute the Script**: Run the main automation script.

    ```bash
    python3 run_startrack.py
    ```

### What happens when you run this?
The script performs the following automated steps:
1.  Creates a build directory (`StarTrack_Build`) and an output directory (`Simulations`).
2.  Iterates through a predefined grid of parameters (Metallicity, Mass, Mass Ratio, Period).
3.  **Patches** the C source code with specific values for the current iteration.
4.  **Compiles** a temporary executable (`startrack.out`).
5.  **Runs** the simulation for that specific binary system.
6.  **Saves** the evolution data into a structured directory tree.

---

## Changes to the StarTrack Setup

This setup fundamentally changes how StarTrack is controlled and executed.

| Feature | Standard StarTrack | StarEstate Grid Wrapper |
| :--- | :--- | :--- |
| **Initialization** | **Stochastic**: Uses Monte Carlo sampling to pick Initial Mass Function (IMF), eccentricity, and separations based on probability distributions. | **Deterministic**: Specific initial values are forced into the code by the Python wrapper. Random generation is bypassed. |
| **Configuration** | Parameters set manually in `sinbin.h` or input files before compilation. | Parameters are injected programmatically into `binary.c` and `sinbin.h` **during runtime** for every single system. |
| **Execution** | Runs thousands of systems in a single execution loop inside C (faster). | Runs **one** system per execution. The Python script handles the looping and repeatedly re-compiles the C code (probably there is a better way of doing this). |
| **Output** | Large data files containing the entire population. | Individual folders for each system configuration containing specific evolution tracks. |

---

## How Binary Evolution Tracks are Created

The binary evolution tracks are generated via Hot-Patching in `run_startrack.py`.

### 1. The Parameter Grid
The script defines a 4-dimensional grid of initial conditions:
* **Metallicity ($Z$):** 0.1 $Z_{\odot}$, 0.3 $Z_{\odot}$, 1.0 $Z_{\odot}$.
* **Primary Mass ($M_a$):** 0.11 $M_{\odot}$ to 100 $M_{\odot}$.
* **Mass Ratio ($q$):** 0.1 to 1.0 (where $M_b = q \times M_a$).
* **Orbital Period ($P$):** A list of specific periods (in days).

### 2. The Injection Process
For every combination of parameters, the Python script reads the raw `binary.c` and `sinbin.h` files and performs string replacement using Regular Expressions (`re` module):

* **In `sinbin.h`**: It updates the metallicity variable `ZZ`.
* **In `binary.c`**: It locates the `main()` function and overwrites the initial variable assignments:
    * `Mzamsa = [Current Grid Value];`
    * `Mzamsb = [Calculated Secondary Mass];`
    * `per0 = [Current Grid Period];`

### 3. Compilation and Simulation
The wrapper runs the following command for every single grid point:
```bash
gcc -O3 -w -o startrack.out binary.c singl.c -lm
./startrack.out
