import os
import shutil
import subprocess
import re
import sys

# --- Configuration ---
# 1. Metallicities to test
METALLICITY_MAP = {
    0.1 * 0.0142: "0.1Zsun",
    0.3 * 0.0142: "0.3Zsun",
    0.0142: "1Zsun"
}

# 2. Primary Masses (Mzamsa)
MASSES = [
    0.11, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    20, 30, 40, 50, 60, 70, 80, 90, 100
]

# 3. Mass Ratios (q = Mzamsb / Mzamsa)
# Range: [0.1 ... 1.0]
Q_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# 4. Initial Orbital Periods (per0) in days
# User requested folder depth for per0. Define your grid here.
PERIODS = [0.15,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5]

# Minimum allowed mass for a star (Solar masses)
MIN_STAR_MASS = 0.11

SOURCE_FILES = ["binary.c", "singl.c", "sinbin.h"]
OUTPUT_ROOT = "Simulations"
BUILD_DIR = "StarTrack_Build"

# --- Helper Functions ---

def read_file(path):
    with open(path, 'r') as f:
        return f.read()

def write_file(path, content):
    with open(path, 'w') as f:
        f.write(content)

def calculate_log_teff(L, R):
    """
    Calculates log10(Teff) based on L and R in solar units.
    Teff = 5778 * (L / R^2)^0.25
    log(Teff) = log(5778) + 0.25*log(L) - 0.5*log(R)
    log10(5778) approx 3.7617
    """
    if L <= 0 or R <= 0: return 0.0
    return 3.7617 + 0.25 * math.log10(L) - 0.5 * math.log10(R)

def process_output_line(line):
    """
    Parses a raw line from StarTrack, calculates logs, and formats for:
    t Ka Kb Ma Mb Ra Rb logLa logLb logTeffa logTeffb
    """
    parts = line.split()
    # We assume StarTrack output (evolution.dat) has at least:
    # t K1 K2 M1 M2 R1 R2 L1 L2 (Indices 0-8)
    # Adjust indices below if your binary.c output order differs!
    try:
        if len(parts) < 9: return None

        t  = float(parts[0])
        k1 = int(float(parts[1])) # cast float->int just in case
        k2 = int(float(parts[2]))
        m1 = float(parts[3])
        m2 = float(parts[4])
        r1 = float(parts[5])
        r2 = float(parts[6])
        l1 = float(parts[7])
        l2 = float(parts[8])

        # Derived Values
        logL1 = math.log10(l1) if l1 > 0 else -99.0
        logL2 = math.log10(l2) if l2 > 0 else -99.0
        logT1 = calculate_log_teff(l1, r1)
        logT2 = calculate_log_teff(l2, r2)

        # Formatting aligned columns
        return (f"{t:10.4e} {k1:2d} {k2:2d} {m1:8.4f} {m2:8.4f} "
                f"{r1:8.4f} {r2:8.4f} {logL1:8.4f} {logL2:8.4f} "
                f"{logT1:8.4f} {logT2:8.4f}")
    except ValueError:
        return None

def patch_sinbin(content, z_val, m_primary, q_val):
    """
    Rewrites sinbin.h for Mass and Metallicity only.
    Period is now handled in binary.c.
    """
    m_secondary = m_primary * q_val

    # 1. Set Binary Mode
    content = re.sub(r'#define\s+BINARY\s+\d+', '#define BINARY 1', content)

    # 2. Set only 1 system to be tested
    content = re.sub(r'#define\s+num_tested\s+\d+', '#define num_tested 1', content)

    # 3. Set Metallicity
    content = re.sub(r'#define\s+ZZ\s+[0-9\.eE+-]+', f'#define ZZ {z_val:.6f}', content)

    # 4. Set Primary Mass (Mmina = Mmaxa = m_primary)
    content = re.sub(r'#define\s+Mmina\s+[0-9\.eE+-]+', f'#define Mmina {m_primary:.4f}', content)
    content = re.sub(r'#define\s+Mmaxa\s+[0-9\.eE+-]+', f'#define Mmaxa {m_primary + 0.001:.4f}', content)

    # 5. Set Secondary Mass (Mminb = Mmaxb = m_secondary)
    content = re.sub(r'#define\s+Mminb\s+[0-9\.eE+-]+', f'#define Mminb {m_secondary:.4f}', content)
    content = re.sub(r'#define\s+Mmaxb\s+[0-9\.eE+-]+', f'#define Mmaxb {m_secondary + 0.001:.4f}', content)

    content = re.sub(r'#define\s+num_tested\s+\d+', '#define num_tested 1', content)
    content = re.sub(r'#define\s+BINOUT\s+\d+', '#define BINOUT 0', content)

    # Increase Hubble time by 10x
    match = re.search(r'#define\s+hub_val\s+([\d\.eE+-]+)', content)
    if match:
        val = float(match.group(1)) * 10.0
        content = re.sub(r'#define\s+hub_val\s+[\d\.eE+-]+', f'#define hub_val {val:.1f}', content)

    return content



def patch_binary_c(content, period_log_val, ma, mb):
    """
    Patches binary.c:
    1. Sets per0 = pow(10, val) to convert logP to Days.
    2. HARDCODES Mzamsa and Mzamsb by intercepting the IMF generation block.
    3. Injects ZAMS print.
    4. Injects Evolution logging with stopping condition.
    """

    # --- 1. Inject Data Logging & Header at Time Increment ---

    # User instruction: replace "per0=get_SS(0.15,5.5,-0.55);" with "per0 = {value};"
    # We look for a line assigning to per0 using get_SS or get_B
    per0_pattern = r'per0\s*=\s*get_SS\s*\([^;]+;'

    if re.search(per0_pattern, content):
        content = re.sub(per0_pattern, f'per0 = {period_log_val};', content, count=1)
    else:
        # Fallback: Look for just assignment if function name differs
        # This regex is broader, use with caution or rely on specific function name if known.
        # Given the user provided the specific line, we stick to that first.
        pass

    # Finding 't += dt;' is safer than the loop start because it's guaranteed
    # to be in the execution flow (after variable declarations).

    # --- Fix Masses (Mzamsa, Mzamsb) ---
    # We find the "else if(IMF==1 && SS==1)" block and preempt it with our own "else if(1)"
    # This effectively comments out the original random mass generation and inserts ours.
    imf_block_pattern = r'(else\s+if\s*\(\s*IMF\s*==\s*1\s*&&\s*SS\s*==\s*1\s*\)\s*\{)'

    # We insert a block that always runs (else if (1)), sets our masses, and then
    # puts the original check in an "else if" so the syntax remains valid but original code is skipped.
    mass_injection = f'''
    else if (1) {{
        Mzamsa = {ma};
        Mzamsb = {mb};
        /* q0, qmin calc if needed by other parts of code, though usually just Mzamsb is enough */
        q0 = Mzamsb / Mzamsa;
    }}
    \\1
    '''
    # Note: \\1 puts the original "else if(IMF==1...)" back, so it becomes the *next* else if branch.
    if re.search(imf_block_pattern, content):
        content = re.sub(imf_block_pattern, mass_injection, content, count=1)
    else:
        print("Warning: Could not find IMF generation block to patch Masses.")


    time_inc_pattern = r'(t\s*\+=\s*dt\s*;)'

    def logging_and_clamp_lambda(m):
        # We use a raw string r'''...''' and strictly formatted C code.
        # we use \\n to ensure a literal backslash+n is written to the C file
        # which the C compiler then interprets as a newline character.

        c_code_block = r'''
        /* --- INJECTED OUTPUT --- */
        {
            static int _header_done = 0;
            if(!_header_done) {
                 setvbuf(stdout, NULL, _IONBF, 0);
                 printf("t Ka Kb Ma Mb Ra Rb a logLa logLb logTeffa logTeffb\n");
                 _header_done = 1;
            }
        }

        if(1) {
            double _lgL1 = (La > 1.0e-10) ? log10(La) : -99.0;
            double _lgL2 = (Lb > 1.0e-10) ? log10(Lb) : -99.0;
            double _lgT1 = (La > 0 && Ra > 0) ? 3.7617 + 0.25*_lgL1 - 0.5*log10(Ra) : 0.0;
            double _lgT2 = (Lb > 0 && Rb > 0) ? 3.7617 + 0.25*_lgL2 - 0.5*log10(Rb) : 0.0;

            /* Print formatted line */
            printf("%12.5e %2d %2d %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f\n",
                   t, Ka, Kb, Ma, Mb, Ra, Rb, a, _lgL1, _lgL2, _lgT1, _lgT2);

            /* --- STOPPING CONDITION --- */
            /* Stop if both stars are evolved beyond normal phases (Type > 9) */
            if (Ka > 9 && Kb > 9) {
                 stop1 = 1;
            }
        }
        /* ----------------------- */
        '''

        # Clamp logic + Original t+=dt
        # We perform the clamping, then the original increment.
        clamp_logic = r' if(dt < 0.1) dt = 0.1; if(dt > 500.0) dt = 500.0; '

        return c_code_block + clamp_logic + m.group(1)

    if re.search(time_inc_pattern, content):
        content = re.sub(time_inc_pattern, logging_and_clamp_lambda, content)
    else:
        print("Warning: Could not find 't += dt' to inject logging. Output may be missing.")

    # --- 2. Fix Mzamsb Generator (Use Linear/Uniform for Grid) ---
    # Change get_M(Mminb, Mmaxb) -> get_M_linear(Mminb, Mmaxb)
    # This ensures the secondary mass is chosen exactly as specified, not weighted by IMF
    content = re.sub(r'(Mzamsb\s*=\s*)get_M(\s*\(\s*Mminb\s*,\s*Mmaxb\s*\)\s*;)',
                     r'\1get_M_linear\2', content)

    return content

def main():
    if os.path.exists(BUILD_DIR):
        shutil.rmtree(BUILD_DIR)
    os.makedirs(BUILD_DIR)

    # 1. Read Raw Source Files
    try:
        raw_sinbin = read_file("sinbin.h")
        raw_binary = read_file("binary.c")
    except FileNotFoundError as e:
        print(f"Error: Missing source file: {e}")
        return

    # Copy singl.c unchanged
    if os.path.exists("singl.c"):
        shutil.copy("singl.c", os.path.join(BUILD_DIR, "singl.c"))

    print(f"--- Starting StarTrack Grid ---")
    print(f"Structure: {OUTPUT_ROOT}/Metallicity/Mzamsa/q/per0/")

    # 2. Parameter Loops
    for z_val, z_name in METALLICITY_MAP.items():
        for ma in MASSES:
            for q in Q_VALUES:
                mb = ma * q

                # Check: When Mzams(b) is < mminb, ignore the simulation
                if mb < MIN_STAR_MASS:
                    continue

                for per0 in PERIODS:

                    dest_dir = os.path.join(OUTPUT_ROOT, z_name, f"{ma}", f"q_{q}", f"per_{per0}")
                    os.makedirs(dest_dir, exist_ok=True)
                    final_output = os.path.join(dest_dir, "evolution.dat")

                    if os.path.exists(final_output):
                        print(f"[Skip] Exists: {dest_dir}")
                        continue

                    print(f"[Run] Z={z_name} Ma={ma} q={q} (Mb={mb:.2f}) log(P/days)={per0}")

                    # A. Patch Binary.c with SPECIFIC PERIOD
                    # This must be done inside the loop now, as per0 changes per iteration
                    mod_binary = mod_binary = patch_binary_c(raw_binary, per0, ma, mb)
                    write_file(os.path.join(BUILD_DIR, "binary.c"), mod_binary)

                    # B. Patch Header (Masses & Z only)
                    mod_sinbin = patch_sinbin(raw_sinbin, z_val, ma, q)
                    write_file(os.path.join(BUILD_DIR, "sinbin.h"), mod_sinbin)

                    # C. Compile
                    compile_cmd = ["gcc", "-O3", "-w", "-o", "startrack.out", "binary.c", "singl.c", "-lm"]
                    try:
                        subprocess.run(compile_cmd, cwd=BUILD_DIR, capture_output=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"  -> Compile Error: {e.stderr.decode()}")
                        return

                    # D. Run
                    try:
                        with open(final_output, "w") as outfile:
                            subprocess.run(
                                ["./startrack.out"],
                                cwd=BUILD_DIR,
                                stdout=outfile,
                                stderr=subprocess.DEVNULL,
                                timeout=60
                            )
                    except subprocess.TimeoutExpired:
                        print("  -> Timeout!")
                    except Exception as e:
                        print(f"  -> Run Error: {e}")

    print("\nAll simulations completed.")

if __name__ == "__main__":
    main()
