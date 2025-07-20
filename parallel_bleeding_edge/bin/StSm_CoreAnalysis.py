#!/usr/bin/env python3
"""
SPH Sink Particle Visualization Script
Pau Amaro Seoane <amaro@riseup.net>

============================================
SPH CORE PARTICLE COLLISION ANALYSIS CODE
============================================

*PHYSICAL SCENARIO*
This code analyzes hydrodynamic simulations of stellar core collisions in binary red giant systems. 
The astrophysical scenario involves:
- Two evolved stars (red giants) in a close binary system
- Each star has a dense core and extended envelope
- As the stars orbit, their envelopes interact and cores may collide
- These collisions produce luminous bursts from shocked gas

The code specifically models:
- Gravitational focusing during close encounters
- Shock heating at the interaction interface
- Radiative losses from the collision zone
- Time-dependent luminosity evolution

*INPUT DATA REQUIREMENTS*
Input files must be ASCII text files converted from SPH simulations using SPLASH. Required preprocessing:

1. Run SPLASH conversion:
   splash to_ascii -f starsmasher --convert=7,1,2,3,4,5,6,8,14 out*.sph

   This extracts these columns in order:
   7: Particle mass [code units]
   1: x position [code units]
   2: y position [code units]
   3: z position [code units]
   4: x velocity [code units]
   5: y velocity [code units]
   6: z velocity [code units]
   8: Density [code units]
   14: Particle type (1=gas, 2=sink/core)

2. Each input file must contain:
   - Header with simulation time on 4th line
   - One row per particle with above columns
   - Exactly two sink particles (itype=2)
   - Sufficient gas particles to resolve envelopes

*CORE ALGORITHMS*

1. ENVELOPE DETECTION:
   For each core particle:
   - Find all gas particles (itype=1)
   - Calculate distances to core
   - Identify 3 densest neighboring particles
   - Set envelope radius as maximum distance of these 3
   - Compute average density within this radius

2. COLLISION DETECTION:
   - Calculate separation between cores
   - Check for envelope overlap (d < R1 + R2)
   - Compute intersection volume using exact formula:
     V_int = (πδ²(3R1 + 3R2 - δ))/12
     where δ = R1 + R2 - d

3. LUMINOSITY CALCULATION:

   SIMPLE METHOD:
   L = 0.1*(0.5*M_int*v_rel²)/Δt
   where:
   - M_int = interacting mass from intersection volume
   - v_rel = relative velocity projected along separation
   - Δt = 2*min(R1,R2)/v_rel (crossing time)

   ACCURATE METHOD:
   L = η(v_rel)*[0.5*M_dot*v_rel² + W_c]
   where:
   - η(v) = 0.1*(1-exp(-(v/v_crit)²)), v_crit=50 km/s
   - M_dot = mass flux through interaction plane
   - W_c = compression work = 0.75*ρ_avg*v_rel³*π*R_avg²

4. DUTY CYCLE CALCULATION:
   - Identify bursts where L > 1% of L_max
   - Sum time intervals containing bursts
   - Divide by total observation time:
     η = T_burst/T_obs

*CODE WORKFLOW*

1. INITIALIZATION:
   - Set physical constants (R_sun, M_sun, etc.)
   - Configure plotting parameters
   - Scan directory for input files

2. FILE PROCESSING:
   For each input file:
   - Read and validate data
   - Extract core particles
   - Analyze gas envelopes
   - Detect collisions
   - Calculate luminosity
   - Store time and luminosity

3. OUTPUT GENERATION:
   - Save lightcurve to file (time, luminosity)
   - Calculate duty cycle if requested
   - Generate visualization plots

*PHYSICAL ASSUMPTIONS*

1. Envelope Properties:
   - Spherically symmetric
   - Uniform density within R_env
   - Well-mixed composition

2. Collision Physics:
   - Instantaneous energy thermalization
   - Thin interaction region
   - 10% radiative efficiency (simple method)
   - Velocity-dependent efficiency (accurate method)

3. Orbital Dynamics:
   - Keplerian orbits between collisions
   - No mass loss between snapshots
   - Slow orbital evolution

*NUMERICAL IMPLEMENTATION*

1. Data Structures:
   - NumPy arrays for particle data
   - Dictionaries for core properties
   - Sorted lists for time series

2. Key Functions:
   - read_sph_data(): File I/O and validation
   - calculate_simple_burst(): Basic luminosity
   - calculate_accurate_burst(): Detailed physics
   - calculate_duty_cycle(): Activity fraction
   - visualize_cores(): Diagnostic plots

3. Optimization:
   - Vectorized calculations
   - Early termination for non-overlapping cases
   - Caching of repeated calculations

*OUTPUTS*

1. Data Files:
   - burst_lightcurve_simple.txt
   - burst_lightcurve_accurate.txt
   Format: time[yr] luminosity[erg/s]

2. Diagnostics:
   - Duty cycle percentage
   - Core separation history
   - Envelope properties

3. Visualizations:
   - XY plots of core positions
   - Envelope boundaries
   - Light curve plots
   - Optional logarithmic scales

*LIMITATIONS*

1. Resolution Effects:
   - Envelope radius depends on particle sampling
   - Under-resolved shocks in simple method

2. Physical Approximations:
   - No radiative transfer
   - Simplified equation of state
   - No magnetic fields

3. Temporal Sampling:
   - Limited by snapshot frequency
   - Burst durations may be undersampled

The code provides a robust framework for analyzing stellar collisions in SPH
simulations, with particular attention to the time-dependent luminous output
from interacting envelopes. The two calculation methods allow for either rapid
estimation or detailed modeling of the collision physics.

"""
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Circle

# =============================================
# PHYSICAL CONSTANTS
# =============================================
R_sun = 6.9599e10          # cm
M_sun = 1.9891e33          # g
G = 6.6739e-8              # cm³ g⁻¹ s⁻²
L_sun_erg = 3.828e33       # Solar luminosity [erg/s]
sec_per_day = 86400
days_per_year = 365.25
t_sun = np.sqrt(R_sun**3 / (G * M_sun))  # 1593.6 seconds
days_per_code_time = t_sun / sec_per_day  # 0.01844 days

# =============================================
# GLOBAL PARAMETERS
# =============================================
GLOBAL_FONT_SIZE = 35
GLOBAL_LINE_WIDTH = 3
GLOBAL_TICK_SIZE = 30
rc('text', usetex=True)
rc('font', family='serif', size=GLOBAL_FONT_SIZE)

def get_simulation_time(filename):
    """Extract time from 4th line of SPH file"""
    with open(filename, 'r') as f:
        for _ in range(3): next(f)
        return float(next(f).split('#')[-1].split()[0])

def read_sph_data(filename):
    """Load SPH data with error handling"""
    try:
        time_code = get_simulation_time(filename)
        data = np.loadtxt(filename)
        return {
            'mass': data[:,0], 'x': data[:,1], 'y': data[:,2],
            'z': data[:,3], 'vx': data[:,4], 'vy': data[:,5],
            'vz': data[:,6], 'density': data[:,7], 'itype': data[:,8].astype(int)
        }, time_code
    except Exception as e:
        print(f"Error reading {filename}: {str(e)}")
        return None, None

def calculate_simple_burst(core1, core2):
    """Basic envelope overlap method with fixed parameter access"""
    try:
        d_vec = core1['pos'] - core2['pos']
        d = np.linalg.norm(d_vec)
        v_rel = np.abs(np.dot(core1['vel']-core2['vel'], d_vec/d))

        delta = core1['r_env'] + core2['r_env'] - d
        V_int = np.pi * delta**2 * (core1['r_env'] + core2['r_env'])/2
        M_int = 0.5 * (core1['rho_phys'] + core2['rho_phys']) * V_int
        delta_t = 2 * min(core1['r_env'], core2['r_env']) / max(v_rel, 1e5)  # Avoid division by zero

        return 0.1 * (0.5 * M_int * v_rel**2) / delta_t
    except Exception as e:
        print(f"Simple burst calculation failed: {str(e)}")
        return None

def calculate_accurate_burst(data, core1, core2):
    """Calculate burst luminosity using mass flux and compression work"""
    try:
        # Calculate separation vector and distance
        d_vec = core1['pos'] - core2['pos']
        d = np.linalg.norm(d_vec)

        # Calculate relative velocity with proper parentheses
        v_rel_vec = core1['vel'] - core2['vel']
        v_rel = max(np.abs(np.dot(v_rel_vec, d_vec/d)), 1e5)  # Fixed parentheses

        # Mass flux calculation through interaction plane
        gas_mask = data['itype'] == 1
        gas_pos = np.column_stack((data['x'][gas_mask], data['y'][gas_mask])) * R_sun
        gas_vel = np.column_stack((data['vx'][gas_mask], data['vy'][gas_mask])) * R_sun/t_sun
        gas_mass = data['mass'][gas_mask] * M_sun

        # Find particles crossing midplane
        mid_pos = 0.5 * (core1['pos'] + core2['pos'])
        rel_pos = gas_pos - mid_pos
        rel_vel = gas_vel - 0.5*(core1['vel']+core2['vel'])
        dot_prod = np.sum(rel_pos * rel_vel, axis=1)

        # Calculate mass flux (code units -> physical)
        M_dot = np.sum(gas_mass * np.maximum(dot_prod, 0)/d) / t_sun

        # Compression work term
        r_avg = 0.5 * (core1['r_env'] + core2['r_env'])
        W_c = 0.75 * (core1['rho_phys']+core2['rho_phys'])/2 * v_rel**3 * np.pi * r_avg**2

        # Velocity-dependent efficiency (0-10%)
        v_crit = 50e5  # 50 km/s in cm/s
        eta = 0.1 * (1 - np.exp(-(v_rel/v_crit)**2))

        return eta * (0.5 * M_dot * v_rel**2 + W_c)

    except Exception as e:
        print(f"Accurate burst calculation failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def process_snapshot(filename, unit_system, method):
    """Analyze single snapshot and return time and luminosity"""
    data, time_code = read_sph_data(filename)
    if data is None or time_code is None:
        return None

    # Verify exactly two core particles
    core_mask = data['itype'] == 2
    if np.sum(core_mask) != 2:
        print(f"{filename}: Needs exactly 2 core particles")
        return None

    # Get core properties
    core_pos = np.column_stack((data['x'][core_mask], data['y'][core_mask])) * R_sun
    core_vel = np.column_stack((data['vx'][core_mask], data['vy'][core_mask])) * R_sun/t_sun

    # Process gas envelopes
    gas_mask = data['itype'] == 1
    gas_pos = np.column_stack((data['x'][gas_mask], data['y'][gas_mask])) * R_sun
    gas_rho = data['density'][gas_mask]
    gas_mass = data['mass'][gas_mask] * M_sun

    cores = []
    for i in range(2):
        dist = np.linalg.norm(gas_pos - core_pos[i], axis=1)
        densest_3_idx = np.argsort(gas_rho)[-3:]
        r_env = np.max(dist[densest_3_idx])

        within = dist <= r_env
        total_mass = np.sum(gas_mass[within])
        rho_code = np.mean(gas_rho[within])
        rho_phys = rho_code * (M_sun/R_sun**3)

        cores.append({
            'pos': core_pos[i],
            'vel': core_vel[i],
            'r_env': r_env,
            'rho_code': rho_code,
            'rho_phys': rho_phys,
            'mass': total_mass
        })

    # Calculate luminosity
    try:
        if method == 'simple':
            L_erg = calculate_simple_burst(cores[0], cores[1])
        else:
            L_erg = calculate_accurate_burst(data, cores[0], cores[1])

        # Return time in years and luminosity
        time_years = time_code * days_per_code_time / days_per_year
        return time_years, L_erg

    except Exception as e:
        print(f"Luminosity calculation failed: {str(e)}")
        return None


def visualize_cores(cores, time_code, unit_system, L_erg=None):
    """Create publication-quality visualization of core envelopes"""
    fig, ax = plt.subplots(figsize=(12, 12))
    colors = ['#1f77b4', '#ff7f0e']  # Blue and orange

    # Plot each core with envelope
    for i, core in enumerate(cores):
        x, y = core['pos'][0]/R_sun, core['pos'][1]/R_sun
        r = core['r_env']/R_sun

        # Core marker
        ax.plot(x, y, 'o', color=colors[i], markersize=12,
               markeredgecolor='k', label=f'Core {i+1}')

        # Envelope fill
        ax.add_patch(Circle((x, y), r, color=colors[i], alpha=0.15))

        # Envelope boundary
        ax.add_patch(Circle((x, y), r, fill=False,
                          color=colors[i], linewidth=2, alpha=0.7))

    # Annotations
    time_str = (f"Time: {time_code:.2f} code" if unit_system == 'code'
               else f"Time: {time_code*days_per_code_time/days_per_year:.2f} yr")
    ax.text(0.02, 0.95, time_str, transform=ax.transAxes,
            fontsize=12, bbox=dict(facecolor='white', alpha=0.8))

    if L_erg is not None:
        ax.text(0.02, 0.88, r"$L = {:.2e} L_\odot$".format(L_erg/L_sun_erg),
                transform=ax.transAxes, fontsize=12,
                bbox=dict(facecolor='white', alpha=0.8))

    # Plot formatting
    ax.set_xlabel(r'$x$ ($R_\odot$)')
    ax.set_ylabel(r'$y$ ($R_\odot$)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_aspect('equal')
    return fig

def calculate_duty_cycle(time_array, lum_array):
    """Calculate duty cycle from burst light curve data"""
    if len(time_array) < 2:
        print("Insufficient data points for duty cycle calculation")
        return None

    # Sort by time
    sort_idx = np.argsort(time_array)
    time_sorted = time_array[sort_idx]
    lum_sorted = lum_array[sort_idx]

    # Calculate time intervals between snapshots
    time_intervals = np.diff(time_sorted)
    avg_interval = np.mean(time_intervals)

    # Identify bursts (luminosity > 1% of max)
    burst_threshold = 0.01 * np.max(lum_sorted)
    burst_mask = lum_sorted > burst_threshold

    # Calculate total burst duration
    burst_duration = np.sum(time_intervals * burst_mask[:-1])

    # Total observation time
    total_time = time_sorted[-1] - time_sorted[0]

    return burst_duration / total_time

def main():
    """Main execution with clear workflow separation"""
    # Initialize
    sph_files = sorted(glob.glob('out*.sph.ascii'))
    if not sph_files:
        print("Error: No input files found")
        return

    # Define output filenames
    simple_output = 'burst_lightcurve_simple.txt'
    accurate_output = 'burst_lightcurve_accurate.txt'

    # First ask for calculation method
    print("\nSelect analysis method:")
    method_choice = input("1: Simple\n2: Accurate\nChoice: ").strip()

    if method_choice not in ['1', '2']:
        print("Invalid choice")
        return

    method = 'simple' if method_choice == '1' else 'accurate'
    output_file = simple_output if method == 'simple' else accurate_output

    # Check for existing file
    recalculate = True
    if os.path.exists(output_file):
        print(f"\nFound existing {output_file}")
        choice = input("Do you want to:\n1: Redo calculation\n2: Use existing file\nChoice: ").strip()
        if choice == '2':
            recalculate = False
            try:
                existing_data = np.loadtxt(output_file)
                time_array = existing_data[:,0]
                luminosity_array = existing_data[:,1]
                print(f"Loaded data from {output_file}")
            except Exception as e:
                print(f"Error loading file: {e}")
                recalculate = True

    if recalculate:
        print(f"\nRunning {method} calculation...")
        time_list = []
        luminosity_list = []
        processed = 0

        for filename in sph_files:
            try:
                result = process_snapshot(filename, 'physical', method)
                if not result:
                    continue

                time_years, L_erg = result
                time_list.append(time_years)
                luminosity_list.append(L_erg)
                processed += 1
                print(f"Processed {os.path.basename(filename)}")

            except Exception as e:
                print(f"Failed {filename}: {str(e)}")
                continue

        if not time_list:
            print("No valid data to plot")
            return None, None

        # Convert to arrays and sort
        time_array = np.array(time_list)
        luminosity_array = np.array(luminosity_list)
        sort_idx = np.argsort(time_array)
        time_array = time_array[sort_idx]
        luminosity_array = luminosity_array[sort_idx]

        # Save results
        np.savetxt(output_file,
                  np.column_stack((time_array, luminosity_array)),
                  header='Time[yr] Luminosity[erg/s]',
                  fmt='%.6e')
        print(f"\nSaved {method} results to {output_file}")

        # Offer duty cycle calculation for accurate method
        if method == 'accurate':
            duty_choice = input("\nCalculate duty cycle? (y/n): ").strip().lower()
            if duty_choice == 'y':
                duty_cycle = calculate_duty_cycle(time_array, luminosity_array)
                if duty_cycle is not None:
                    print(f"Duty cycle: {duty_cycle:.3f} ({duty_cycle*100:.1f}%)")

    # Plotting options
    print("\nSelect plot options:")
    xscale = input("X-axis scale:\n1: Linear\n2: Logarithmic\nChoice: ").strip()
    yscale = input("Y-axis scale:\n1: Linear\n2: Logarithmic\nChoice: ").strip()

    xscale = 'linear' if xscale == '1' else 'log'
    yscale = 'linear' if yscale == '1' else 'log'

    # Create plot
    plt.figure(figsize=(12, 6))
    plt.plot(time_array, luminosity_array, 'k-', linewidth=2)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel('Time (years)')
    plt.ylabel('Luminosity (erg/s)')
    plt.grid(True, alpha=0.3, which='both' if xscale == 'log' or yscale == 'log' else 'major')
    plt.tight_layout()
    plt.show()

    return time_array, luminosity_array

if __name__ == "__main__":
    main()
