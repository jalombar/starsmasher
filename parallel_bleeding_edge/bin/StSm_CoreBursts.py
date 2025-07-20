#!/usr/bin/env python3
"""
SPH BURST LUMINOSITY ANALYZER
================================================================================
Purpose:
Analyzes SPH simulations of degenerate core collisions to detect and characterize
burst luminosity events when gas envelopes interact.

Key Features:
- Detects envelope interactions between sink particles
- Calculates burst luminosities in both erg/s and L_sun
- Tracks interacting masses and relative velocities
- Outputs detailed ASCII data with conversion factors
- Generates customizable log/linear plots

REQUIREMENTS:
1. INPUT FILES:
   - Must be SPLASH-generated ASCII files named out*.sph.ascii
   - Must contain these columns (via --convert=7,1,2,3,4,5,6,8,14):
     7: Particle mass
     1: x position
     2: y position
     3: z position
     4: vx velocity
     5: vy velocity
     6: vz velocity
     8: Density
     14: Particle type (itype)

2. HOW TO GENERATE INPUT FILES:
   Run this command before using this script:
   $ splash to ascii -f starsmasher --convert=7,1,2,3,4,5,6,8,14 out*.sph
   Or use the provided StSm_TrajectoryAscii.sh script

PHYSICS DESCRIPTION:
Bursts occur when gas envelopes (defined by 6 densest particles) of two sink
particles overlap. Luminosity is calculated as:
L = η × (0.5 × M_int × v_rel²) / Δt
where:
- η = 0.1 (radiative efficiency)
- M_int = interacting gas mass
- v_rel = relative velocity
- Δt = 2 × min(r1,r2)/v_rel (duration)

UNITS:
All calculations done in code units (R_sun, M_sun, G=1) with physical unit
conversions provided in output file header.
"""

import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import rc

# =============================================
# PHYSICAL CONSTANTS
# =============================================
R_sun = 6.9599e10          # cm
M_sun = 1.9891e33          # g
G = 6.6739e-8              # cm³ g⁻¹ s⁻²
L_sun_erg = 3.828e33       # erg/s
sec_per_day = 86400
days_per_year = 365.25

# Derived units
t_sun = np.sqrt(R_sun**3 / (G * M_sun))  # 1593.6 seconds
days_per_code_time = t_sun / sec_per_day  # 0.01844 days

# =============================================
# GLOBAL STYLE PARAMETERS
# =============================================
GLOBAL_FONT_SIZE = 35
GLOBAL_LINE_WIDTH = 3
GLOBAL_TICK_SIZE = 30
OUTPUT_FILE = "burst_analysis.dat"

rc('text', usetex=True)
rc('font', family='serif', size=GLOBAL_FONT_SIZE)
rc('axes', titlesize=GLOBAL_FONT_SIZE, labelsize=GLOBAL_FONT_SIZE)
rc('xtick', labelsize=GLOBAL_TICK_SIZE)
rc('ytick', labelsize=GLOBAL_TICK_SIZE)

def get_user_choice(prompt, options):
    """Get validated user input"""
    while True:
        print(prompt)
        for i, opt in enumerate(options, 1):
            print(f"{i}: {opt}")
        choice = input(f"Select [1-{len(options)}]: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(options):
            return options[int(choice)-1]
        print(f"Invalid choice. Please enter 1-{len(options)}")

def get_simulation_time(filename):
    """Extract simulation time from 4th line of file"""
    with open(filename, 'r') as f:
        try:
            for _ in range(3): next(f)  # Skip 3 lines
            time_line = next(f).strip()
            return float(time_line.split('#')[-1].strip().split()[0])
        except Exception as e:
            print(f"Warning: Could not parse time in {filename}: {str(e)}")
            return None

def read_sph_data(filename):
    """Read SPH data with robust time extraction"""
    try:
        time_code = get_simulation_time(filename)
        if time_code is None:
            return None, None
        
        data = np.loadtxt(filename)
        return {
            'mass': data[:, 0],    # Column 7
            'x': data[:, 1],       # Column 1
            'y': data[:, 2],       # Column 2  
            'z': data[:, 3],       # Column 3
            'vx': data[:, 4],      # Column 4
            'vy': data[:, 5],      # Column 5
            'vz': data[:, 6],      # Column 6
            'density': data[:, 7], # Column 8
            'itype': data[:, 8].astype(int)  # Column 14
        }, time_code
    except Exception as e:
        print(f"Error reading {filename}: {str(e)}")
        return None, None

def calculate_burst(data, time_code):
    """Calculate burst luminosity for a single snapshot"""
    sink_mask = data['itype'] == 2
    gas_mask = data['itype'] == 1
    
    if np.sum(sink_mask) != 2:
        return None
    
    # Convert to physical units (cm and cm/s)
    sink_pos = np.column_stack((data['x'][sink_mask], data['y'][sink_mask])) * R_sun
    sink_vel = np.column_stack((data['vx'][sink_mask], data['vy'][sink_mask])) * R_sun / t_sun
    gas_pos = np.column_stack((data['x'][gas_mask], data['y'][gas_mask])) * R_sun
    gas_dens = data['density'][gas_mask]
    gas_mass = data['mass'][gas_mask] * M_sun
    
    # Process both sinks
    results = []
    for i in range(2):
        distances = np.linalg.norm(gas_pos - sink_pos[i], axis=1)
        densest_idx = np.argsort(gas_dens)[-6:]
        sorted_dist = np.sort(distances[densest_idx])
        r_3rd = sorted_dist[2]  # 3rd farthest particle
        
        within_r3rd = distances <= r_3rd
        rho = np.sum(gas_mass[within_r3rd]) / (4/3 * np.pi * r_3rd**3)
        
        results.append({
            'position': sink_pos[i],
            'velocity': sink_vel[i],
            'radius': r_3rd,
            'density': rho
        })
    
    # Calculate interaction parameters
    d_vec = results[0]['position'] - results[1]['position']
    d = np.linalg.norm(d_vec)
    v_rel = np.abs(np.dot(results[0]['velocity'] - results[1]['velocity'], d_vec/d))
    
    if d >= results[0]['radius'] + results[1]['radius']:
        return None
    
    # Burst physics
    delta = results[0]['radius'] + results[1]['radius'] - d
    V_int = np.pi * delta**2 * (results[0]['radius'] + results[1]['radius'])/2
    M_int = 0.5 * (results[0]['density'] + results[1]['density']) * V_int
    delta_t = 2 * min(results[0]['radius'], results[1]['radius']) / v_rel
    L_erg = 0.1 * (0.5 * M_int * v_rel**2) / delta_t  # 10% efficiency
    
    return {
        'time_code': time_code,
        'time_days': time_code * days_per_code_time,
        'time_years': time_code * days_per_code_time / days_per_year,
        'L_erg': L_erg,
        'L_sun': L_erg / L_sun_erg,
        'M_int': M_int / M_sun,
        'v_rel': v_rel / 1e5,  # km/s
        'duration': delta_t
    }

def write_ascii_output(bursts):
    """Write comprehensive burst data to ASCII file"""
    with open(OUTPUT_FILE, 'w') as f:
        # Detailed header
        f.write("# SPH BURST LUMINOSITY ANALYSIS\n")
        f.write("# ========================================================================\n")
        f.write("# PHYSICAL INTERPRETATION:\n")
        f.write("# Bursts occur when gas envelopes (defined by 6 densest particles) of two\n")
        f.write("# sink particles overlap. Luminosity is calculated as:\n")
        f.write("# L = η × (0.5 × M_int × v_rel²) / Δt\n")
        f.write("# where η=0.1 (radiative efficiency), M_int=interacting mass,\n")
        f.write("# v_rel=relative velocity, Δt=2×min(r1,r2)/v_rel (duration)\n#\n")
        f.write("# CONVERSION FACTORS (code → physical):\n")
        f.write(f"# 1 code_time = {t_sun:.6e} seconds\n")
        f.write(f"# 1 code_time = {days_per_code_time:.6e} days\n")
        f.write(f"# 1 code_time = {days_per_code_time/days_per_year:.6e} years\n")
        f.write(f"# 1 code_length = {R_sun:.6e} cm (1 R_sun)\n")
        f.write(f"# 1 code_mass = {M_sun:.6e} g (1 M_sun)\n")
        f.write(f"# 1 code_luminosity = {L_sun_erg:.6e} erg/s (1 L_sun)\n#\n")
        f.write("# COLUMN DESCRIPTIONS:\n")
        f.write("# 1: Time (code units) - Direct from simulation\n")
        f.write("# 2: Time (days) - code_time × {:.6e}\n".format(days_per_code_time))
        f.write("# 3: Time (years) - days / 365.25\n")
        f.write("# 4: Luminosity (erg/s) - 10% of kinetic energy over duration\n")
        f.write("# 5: Luminosity (L_sun) - erg/s / {:.6e}\n".format(L_sun_erg))
        f.write("# 6: Interacting mass (M_sun) - Mass in collision volume\n")
        f.write("# 7: Relative velocity (km/s) - Projected approach velocity\n")
        f.write("# 8: Duration (seconds) - 2×min_radius / v_rel\n#\n")
        f.write("# DATA FORMAT:\n")
        f.write("# All values space-delimited, scientific notation with 6 decimals\n")
        f.write("# Time[code] Time[days] Time[years] L[erg/s] L[Lsun] M_int[Msun] v_rel[km/s] Duration[s]\n")
        
        # Write data
        for burst in bursts:
            f.write(f"{burst['time_code']:.6e} ")
            f.write(f"{burst['time_days']:.6e} ")
            f.write(f"{burst['time_years']:.6e} ")
            f.write(f"{burst['L_erg']:.6e} ")
            f.write(f"{burst['L_sun']:.6e} ")
            f.write(f"{burst['M_int']:.6e} ")
            f.write(f"{burst['v_rel']:.6e} ")
            f.write(f"{burst['duration']:.6e}\n")
    
    print(f"\nBurst data written to {OUTPUT_FILE}")
    print(f"Total bursts detected: {len(bursts)}")

def load_ascii_output():
    """Load existing burst data from ASCII file"""
    if not os.path.exists(OUTPUT_FILE):
        return None
    
    bursts = []
    with open(OUTPUT_FILE, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            try:
                values = list(map(float, line.strip().split()))
                bursts.append({
                    'time_code': values[0],
                    'time_days': values[1],
                    'time_years': values[2],
                    'L_erg': values[3],
                    'L_sun': values[4],
                    'M_int': values[5],
                    'v_rel': values[6],
                    'duration': values[7]
                })
            except:
                continue
    
    return bursts if bursts else None

def analyze_simulation(sph_files):
    """Process all snapshots to detect bursts"""
    bursts = []
    print("\nProcessing snapshots...")
    
    for i, filename in enumerate(sorted(sph_files), 1):
        data, time_code = read_sph_data(filename)
        if data is None or time_code is None:
            continue
            
        burst = calculate_burst(data, time_code)
        if burst:
            bursts.append(burst)
            if i % 50 == 0 or i == len(sph_files):
                print(f"Processed {i}/{len(sph_files)} files, {len(bursts)} bursts found")
    
    if bursts:
        write_ascii_output(bursts)
    else:
        print("No bursts detected in any snapshot")
    
    return bursts

def plot_results(bursts, time_unit='years', xscale='linear', yscale='log'):
    """Create customizable luminosity plot"""
    if not bursts:
        print("No bursts to plot")
        return
    
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Set time units
    if time_unit == 'code':
        times = [b['time_code'] for b in bursts]
        xlabel = 'Time (code units)'
    elif time_unit == 'days':
        times = [b['time_days'] for b in bursts]
        xlabel = 'Time (days)'
    else:  # years
        times = [b['time_years'] for b in bursts]
        xlabel = 'Time (years)'
    
    luminosities = [b['L_erg'] for b in bursts]
    
    # Plotting
    ax.plot(times, luminosities, 'k-', linewidth=GLOBAL_LINE_WIDTH)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel, fontsize=GLOBAL_FONT_SIZE)
    ax.set_ylabel(r'Luminosity (erg s$^{-1}$)', fontsize=GLOBAL_FONT_SIZE)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=GLOBAL_TICK_SIZE)
    
    plt.tight_layout()
    plt.show()

def main():
    # Print full header documentation
    print(__doc__)
    
    # Check for input files
    sph_files = glob.glob('out*.sph.ascii')
    if not sph_files:
        print("\nERROR: No out*.sph.ascii files found")
        print("Please generate them first using:")
        print("  splash to ascii -f starsmasher --convert=7,1,2,3,4,5,6,8,14 out*.sph")
        return
    
    # Check for existing output
    existing_data = None
    if os.path.exists(OUTPUT_FILE):
        print(f"\nFound existing output file: {OUTPUT_FILE}")
        choice = get_user_choice(
            "What would you like to do?",
            ["Use existing data", "Recalculate from snapshots"]
        )
        if choice == "Use existing data":
            existing_data = load_ascii_output()
            if not existing_data:
                print("Error: Could not read existing file, recalculating...")
                existing_data = None
    
    # Process data if needed
    bursts = existing_data if existing_data else analyze_simulation(sph_files)
    
    if bursts:
        # Get user preferences
        print("\nSelect plot options:")
        time_unit = get_user_choice(
            "Time units for plotting:",
            ["years", "days", "code units"]
        )
        xscale = get_user_choice(
            "X-axis scale:",
            ["linear", "log"]
        )
        yscale = get_user_choice(
            "Y-axis scale:",
            ["linear", "log"]
        )
        
        plot_results(bursts, time_unit, xscale, yscale)

if __name__ == "__main__":
    main()
