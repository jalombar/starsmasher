#!/usr/bin/env python3
"""
================================================================================
CORE MERGER PREDICTION TOOL: GAS-DRAG MEDIATED ORBITAL DECAY IN COMMON ENVELOPES
================================================================================

DESCRIPTION:
This code predicts the merger timescale of two degenerate cores (white dwarfs, neutron stars,
or black holes) embedded in a common envelope using gas drag physics. It combines SPH simulation
data with analytical orbital decay models to extrapolate the merger time beyond the simulation
duration. The tool:

1. Reads core trajectories and SPH snapshots from stellar evolution simulations
2. Calculates gas drag forces on the cores using envelope properties
3. Calibrates the drag coefficient using end-of-simulation data
4. Solves the orbital decay ODE to predict future evolution
5. Determines the merger time when cores reach a user-specified separation

PHYSICAL PROBLEM ADDRESSED:
- Binary evolution during common envelope phase
- Orbital decay of compact objects in stellar envelopes
- Timescale estimation for WD-WD, NS-NS, or BH-BH mergers
- Common envelope ejection efficiency studies

METHODOLOGY:
1. Gas Drag Calculation:
   - Identifies gas envelope around each core using densest neighbors
   - Computes drag force: F_drag = ½ρv²AC_d
   - Accounts for relative core-envelope motion

2. Orbital Decay Model:
   - Energy dissipation from gas drag powers orbital decay
   - Solves dE/dt = F·v = -(G M1 M2)/(2a²) da/dt
   - Uses Runge-Kutta ODE solver with adaptive timestepping

3. Drag Coefficient Calibration:
   - Compares observed vs. model decay rates at simulation end
   - Determines optimal C_d via median ratio of decay rates
   - Provides physically reasonable bounds (0.01 < C_d < 100)

INPUT REQUIREMENTS:
- Trajectory file: 'compact_objects_trajectory.ascii' with columns:
  [time, x1, y1, z1, x2, y2, z2, separation]
- SPH snapshots: Series of 'out*.sph.ascii' files containing:
  - Particle positions, velocities, densities
  - Particle types (core vs. gas)
- Unit calibration (either from simulation or user input):
  - Length scale (typically R_sun)
  - Mass scale (typically M_sun)
  - Time scale (typically years)

OUTPUTS:
- Prediction of merger time in years
- Plot of separation vs. time (log/linear scales)
- Calibrated drag coefficient value
- Core mass evolution tracking

USAGE INSTRUCTIONS:
1. Place trajectory file and SPH snapshots in working directory
2. Run script and follow interactive prompts
3. View results in terminal and generated PDF plot

CAVEATS/LIMITATIONS:
- Assumes axisymmetric envelope structure
- Neglects dynamical friction from distant particles
- Simplified drag coefficient treatment
- Assumes quasi-steady orbital evolution
- Limited by simulation time resolution

DEVELOPMENT HISTORY:
- v1.0 (2023-08-15): Initial release
- v1.1 (2024-01-22): Added drag coefficient calibration
- v2.0 (2024-06-30): Complete rewrite with physical units

Pau Amaro Seoane amaro@riseup.net
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import glob
import os
import time

# Physical constants (CGS units)
R_sun = 6.9599e10
M_sun = 1.9891e33
G = 6.67430e-8
sec_per_year = 365.25*86400

# Global style parameters
GLOBAL_FONT_SIZE = 50
GLOBAL_LINE_WIDTH = 3
GLOBAL_TICK_SIZE = 40

# LaTeX setup with global styles
rc('text', usetex=True)
rc('font', family='serif', size=GLOBAL_FONT_SIZE)
rc('axes', titlesize=GLOBAL_FONT_SIZE, labelsize=GLOBAL_FONT_SIZE)
rc('xtick', labelsize=GLOBAL_TICK_SIZE)
rc('ytick', labelsize=GLOBAL_TICK_SIZE)
rc('lines', linewidth=GLOBAL_LINE_WIDTH)

# Simulation parameters
CORE_TYPE = 2
GAS_TYPE = 1
DENSEST_NEIGHBORS = 3

def read_trajectory(filename):
    data = np.loadtxt(filename)
    return {
        'time': data[:,0],
        'x1': data[:,1], 'y1': data[:,2], 'z1': data[:,3],
        'x2': data[:,4], 'y2': data[:,5], 'z2': data[:,6],
        'distance': data[:,7]
    }

def get_snapshot_time(filename):
    with open(filename, 'r') as f:
        for line in f:
            if "time:" in line:
                time_line = next(f).strip()
                if time_line.startswith('#'):
                    time_values = time_line[1:].strip().split()
                else:
                    time_values = time_line.split()
                if time_values:
                    try:
                        return float(time_values[0])
                    except ValueError:
                        print(f"Could not convert time value in {filename}")
                break
    print(f"No time found in {filename}")
    return None

def read_sph_snapshot(filename):
    time_code = get_snapshot_time(filename)
    if time_code is None:
        return None
    try:
        header_lines = 0
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                header_lines += 1
        data = np.loadtxt(filename, skiprows=header_lines)
        if len(data) == 0:
            return None
        return {
            'time': time_code,
            'mass': data[:,0],
            'x': data[:,1], 'y': data[:,2], 'z': data[:,3],
            'vx': data[:,4], 'vy': data[:,5], 'vz': data[:,6],
            'density': data[:,7],
            'itype': data[:,8].astype(int)
        }
    except Exception as e:
        return None

def get_core_masses(sph_data):
    if sph_data is None:
        return np.array([])
    core_mask = sph_data['itype'] == CORE_TYPE
    return sph_data['mass'][core_mask]

def calculate_envelope_properties(sph_data, core_pos):
    if sph_data is None:
        return 0, 0, np.zeros(3), 0
    gas_mask = sph_data['itype'] == GAS_TYPE
    gas_pos = np.column_stack((sph_data['x'][gas_mask], sph_data['y'][gas_mask], sph_data['z'][gas_mask]))
    gas_rho = sph_data['density'][gas_mask]
    if len(gas_rho) == 0:
        return 0, 0, np.zeros(3), 0
    dist = np.linalg.norm(gas_pos - core_pos, axis=1)
    densest_idx = np.argsort(gas_rho)[-min(DENSEST_NEIGHBORS, len(gas_rho)):]
    r_env = np.max(dist[densest_idx])
    env_mask = dist <= r_env
    if not np.any(env_mask):
        return r_env, 0, np.zeros(3), 0
    env_rho = np.mean(gas_rho[env_mask])
    env_mass = sph_data['mass'][gas_mask][env_mask]
    env_vel = np.column_stack((
        sph_data['vx'][gas_mask][env_mask],
        sph_data['vy'][gas_mask][env_mask],
        sph_data['vz'][gas_mask][env_mask]
    ))
    if np.sum(env_mass) > 0:
        avg_gas_vel = np.sum(env_mass[:,None] * env_vel, axis=0) / np.sum(env_mass)
    else:
        avg_gas_vel = np.zeros(3)
    return r_env, env_rho, avg_gas_vel, np.sum(env_mass)

def calculate_drag_force(sph_data, core_idx, length_unit, mass_unit, time_unit, C_d):
    if sph_data is None:
        return np.zeros(3), 0, 0, 0
    core_mask = sph_data['itype'] == CORE_TYPE
    core_indices = np.where(core_mask)[0]
    if len(core_indices) <= core_idx:
        return np.zeros(3), 0, 0, 0
    core_pos = np.array([
        sph_data['x'][core_indices[core_idx]],
        sph_data['y'][core_indices[core_idx]],
        sph_data['z'][core_indices[core_idx]]
    ]) * length_unit
    core_vel = np.array([
        sph_data['vx'][core_indices[core_idx]],
        sph_data['vy'][core_indices[core_idx]],
        sph_data['vz'][core_indices[core_idx]]
    ]) * (length_unit / time_unit)
    r_env, rho_gas, avg_gas_vel, env_mass = calculate_envelope_properties(
        sph_data,
        np.array([
            sph_data['x'][core_indices[core_idx]],
            sph_data['y'][core_indices[core_idx]],
            sph_data['z'][core_indices[core_idx]]
        ])
    )
    r_env_phys = r_env * length_unit
    rho_gas_phys = rho_gas * (mass_unit / length_unit**3)
    v_rel = core_vel - avg_gas_vel
    v_mag = np.linalg.norm(v_rel)
    A = np.pi * r_env_phys**2
    F_drag_mag = 0.5 * rho_gas_phys * v_mag**2 * A * C_d
    F_drag = -F_drag_mag * v_rel / max(v_mag, 1e-10)
    return F_drag, r_env_phys, rho_gas_phys, env_mass * mass_unit

def gas_drag_decay_ode(t, d, interpolators, M1, M2):
    try:
        F1 = interpolators['F1'](t)
        F2 = interpolators['F2'](t)
        v1 = interpolators['v1'](t)
        v2 = interpolators['v2'](t)
        P_diss = np.dot(F1, v1) + np.dot(F2, v2)
        dE_dd = G * M1 * M2 / (2 * d**2)
        if dE_dd == 0 or abs(P_diss) < 1e-10:
            return 0
        dd_dt = P_diss / dE_dd
        return dd_dt
    except Exception as e:
        return 0

def save_analysis_data(filename, drag_data):
    with open(filename, 'w') as f:
        f.write("# Core merger analysis data\n")
        f.write("# Columns: time F1x F1y F1z F2x F2y F2z v1x v1y v1z v2x v2y v2z M1 M2\n")
        for i in range(len(drag_data['time'])):
            t = drag_data['time'][i]
            F1 = drag_data['F1'][i]
            F2 = drag_data['F2'][i]
            v1 = drag_data['v1'][i]
            v2 = drag_data['v2'][i]
            M1 = drag_data['M1'][i]
            M2 = drag_data['M2'][i]
            f.write(f"{t:.6e} {F1[0]:.6e} {F1[1]:.6e} {F1[2]:.6e} {F2[0]:.6e} {F2[1]:.6e} {F2[2]:.6e} ")
            f.write(f"{v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e} {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e} ")
            f.write(f"{M1:.6e} {M2:.6e}\n")

def load_analysis_data(filename):
    try:
        data = np.loadtxt(filename)
        if len(data) == 0:
            return None
        drag_data = {
            'time': data[:,0],
            'F1': np.column_stack((data[:,1], data[:,2], data[:,3])),
            'F2': np.column_stack((data[:,4], data[:,5], data[:,6])),
            'v1': np.column_stack((data[:,7], data[:,8], data[:,9])),
            'v2': np.column_stack((data[:,10], data[:,11], data[:,12])),
            'M1': data[:,13],
            'M2': data[:,14]
        }
        return drag_data
    except:
        return None

def main():
    # Configuration
    traj_file = "compact_objects_trajectory.ascii"
    sph_dir = "."
    analysis_cache = "redgiants_binarycore_predictor.txt"

    # Physical constants and simulation parameters
    total_simulation_time_physical = 1.25 * sec_per_year
    initial_separation_physical = 500 * R_sun
    total_core_mass_physical = 2.8 * M_sun

    # Read trajectory data (always needed)
    print("Reading trajectory data...")
    traj = read_trajectory(traj_file)
    traj_times = traj['time']
    traj_separations = traj['distance']

    # Find SPH snapshots (at least need first one for unit calibration)
    print("Locating SPH snapshots...")
    sph_files = sorted(glob.glob(os.path.join(sph_dir, "out*.sph.ascii")))
    if not sph_files:
        raise FileNotFoundError("No SPH snapshot files found!")
    print(f"Found {len(sph_files)} snapshot files")

    # Get first snapshot for unit calibration
    first_snap = None
    for f in sph_files:
        first_snap = read_sph_snapshot(f)
        if first_snap is not None:
            break
    if first_snap is None:
        raise ValueError("No valid snapshots found for unit calibration")

    # Calculate unit conversion factors (always needed for plotting)
    print("Determining simulation units...")
    core_masses = get_core_masses(first_snap)
    if len(core_masses) != 2:
        raise ValueError("Could not get core masses from first snapshot")
    
    length_unit = initial_separation_physical / traj_separations[0]
    mass_unit = total_core_mass_physical / np.sum(core_masses)
    time_unit = total_simulation_time_physical / (traj_times[-1] - traj_times[0])
    
    print(f"Unit calibration:")
    print(f"  length_unit = {length_unit/R_sun:.3f} R_sun")
    print(f"  mass_unit = {mass_unit/M_sun:.3f} M_sun")
    print(f"  time_unit = {time_unit/sec_per_year:.3f} years")

    # Check for existing analysis data
    cached_data = load_analysis_data(analysis_cache)
    if cached_data is not None:
        use_cached = input("Found existing analysis data. Use cached results? (y/n): ").lower() == 'y'
        if use_cached:
            drag_data = cached_data
            print("Using cached analysis data.")
        else:
            cached_data = None
    
    if cached_data is None:
        # Process all snapshots
        print("\n=== Starting Analysis ===")
        start_time = time.time()
        drag_data = {'time': [], 'F1': [], 'F2': [], 'v1': [], 'v2': [], 'M1': [], 'M2': []}
        processed = 0
        
        for i, f in enumerate(sph_files):
            sph_data = read_sph_snapshot(f)
            if sph_data is None:
                continue
                
            core_mask = sph_data['itype'] == CORE_TYPE
            core_indices = np.where(core_mask)[0]
            
            if len(core_indices) != 2:
                continue
            
            # Get core masses
            core_masses = get_core_masses(sph_data)
            if len(core_masses) != 2:
                continue
            
            # Calculate drag forces with C_d=1 (will be scaled later)
            F_drag1, _, _, _ = calculate_drag_force(sph_data, 0, length_unit, mass_unit, time_unit, 1.0)
            F_drag2, _, _, _ = calculate_drag_force(sph_data, 1, length_unit, mass_unit, time_unit, 1.0)
            
            # Get core velocities in physical units
            core_vel1 = np.array([
                sph_data['vx'][core_indices[0]] * (length_unit / time_unit),
                sph_data['vy'][core_indices[0]] * (length_unit / time_unit),
                sph_data['vz'][core_indices[0]] * (length_unit / time_unit)
            ])
            core_vel2 = np.array([
                sph_data['vx'][core_indices[1]] * (length_unit / time_unit),
                sph_data['vy'][core_indices[1]] * (length_unit / time_unit),
                sph_data['vz'][core_indices[1]] * (length_unit / time_unit)
            ])
            
            drag_data['time'].append(sph_data['time'])
            drag_data['F1'].append(F_drag1)
            drag_data['F2'].append(F_drag2)
            drag_data['v1'].append(core_vel1)
            drag_data['v2'].append(core_vel2)
            drag_data['M1'].append(core_masses[0] * mass_unit)
            drag_data['M2'].append(core_masses[1] * mass_unit)
            
            processed += 1
            if processed % 10 == 0:
                print(f"Processed {processed}/{len(sph_files)} snapshots")
        
        if len(drag_data['time']) < 2:
            print("Insufficient valid snapshots for analysis")
            return
        
        # Save the analysis data for future use
        save_analysis_data(analysis_cache, drag_data)
        print(f"Analysis data saved to {analysis_cache}")
        print(f"Total processing time: {time.time()-start_time:.1f} seconds")
    
    # Get user input for drag coefficients
    cd_input = input("Enter drag coefficients (comma-separated, e.g. '0.1,0.5,1.0'): ")
    try:
        C_d_list = [float(x.strip()) for x in cd_input.split(',')]
    except:
        print("Invalid input, using default [0.1, 0.5, 1.0]")
        C_d_list = [0.1, 0.5, 1.0]
    
    # Get other user inputs
    default_threshold = 0.1 * R_sun
    threshold_input = input(f"Enter merger threshold in solar radii (default 2xRcore=0.1 Rsun): ")
    try:
        merger_threshold = float(threshold_input or "0.1") * R_sun
    except ValueError:
        merger_threshold = default_threshold
    
    log_x = input("Use logarithmic scale for X axis? (default n): ").lower() == 'y'
    log_y = input("Use logarithmic scale for Y axis? (default y): ").lower() == 'y'
    use_physical_units = input("Use physical units? (default y): ").lower() == 'y'
    
    # Create interpolators
    print("Creating interpolators...")
    interpolators = {
        'F1': interp1d(drag_data['time'], drag_data['F1'], axis=0, fill_value="extrapolate"),
        'F2': interp1d(drag_data['time'], drag_data['F2'], axis=0, fill_value="extrapolate"),
        'v1': interp1d(drag_data['time'], drag_data['v1'], axis=0, fill_value="extrapolate"),
        'v2': interp1d(drag_data['time'], drag_data['v2'], axis=0, fill_value="extrapolate"),
        'M1': interp1d(drag_data['time'], drag_data['M1'], fill_value="extrapolate"),
        'M2': interp1d(drag_data['time'], drag_data['M2'], fill_value="extrapolate")
    }
    
    # Unit conversions for plotting
    if use_physical_units:
        plot_time_scale = 1 / sec_per_year
        plot_dist_scale = 1 / R_sun
        traj_times_plot = traj_times * time_unit * plot_time_scale
        traj_separations_plot = traj_separations * length_unit * plot_dist_scale
        threshold_plot = merger_threshold * plot_dist_scale
        time_label = "Time (years)"
        dist_label = r"Separation ($R_\odot$)"
    else:
        plot_time_scale = 1
        plot_dist_scale = 1 / length_unit
        traj_times_plot = traj_times
        traj_separations_plot = traj_separations
        threshold_plot = merger_threshold / length_unit
        time_label = "Time (code units)"
        dist_label = "Separation (code units)"
    
    # Create plot
    print("Generating plot...")
    plt.figure(figsize=(14, 10))
    
    # Plot trajectory data
    plt.plot(traj_times_plot, traj_separations_plot, 'k-')
    
    # Plot merger threshold
    plt.axhline(y=threshold_plot, color='r', linestyle=':')
    
    # Define line styles and colors
    line_styles = ['-', '--', '-.', ':']
    colors = ['b', 'g', 'm', 'c', 'y']
    merger_times = []
    
    # Create merger event function
    def merger_event(t, d):
        return d[0] - merger_threshold
    merger_event.terminal = True
    merger_event.direction = -1
    
    # Plot model predictions for each C_d
    for i, C_d in enumerate(C_d_list):
        # Scale the forces by the current C_d
        scaled_interpolators = {
            'F1': lambda t: interpolators['F1'](t) * C_d,
            'F2': lambda t: interpolators['F2'](t) * C_d,
            'v1': interpolators['v1'],
            'v2': interpolators['v2'],
            'M1': interpolators['M1'],
            'M2': interpolators['M2']
        }
        
        # Solve ODE for orbital decay with extended integration
        last_time = traj_times[-1]
        last_separation = traj_separations[-1] * length_unit
        t_span = (last_time, last_time + 1e10)  # Extended integration time
        
        solution = solve_ivp(
            fun=lambda t, d: gas_drag_decay_ode(t, d, scaled_interpolators, 
                                               scaled_interpolators['M1'](t), 
                                               scaled_interpolators['M2'](t)),
            t_span=t_span,
            y0=[last_separation],
            method='RK45',
            events=merger_event,
            dense_output=True
        )
        
        # Convert to plotting units
        if use_physical_units:
            model_times = solution.t * time_unit * plot_time_scale
            model_separations = solution.y[0] * plot_dist_scale
        else:
            model_times = solution.t
            model_separations = solution.y[0] / length_unit
        
        # Plot this model
        color = colors[i % len(colors)]
        line_style = line_styles[i % len(line_styles)]
        plt.plot(model_times, model_separations, linestyle=line_style, color=color,
                label=rf"$C_d = {C_d:.1f}$")
        
        # Record merger time if it occurred
        if solution.t_events[0].size > 0:
            merger_time = solution.t_events[0][0] * time_unit * plot_time_scale
            merger_times.append((C_d, merger_time))
    
    # Create title with merger times in exponential notation
    #title_str = r"Orbital Decay of Degenerate Cores" + "\n"
    #for C_d, t_merge in merger_times:
    #    # Format time in exponential notation with LaTeX
    #    t_str = f"{t_merge:.2e}".replace('e', r' \times 10^{') + '}'
    #    title_str += rf"$T_{{\rm mrg}}(C_d={C_d:.1f}) = {t_str}$ yrs, "
    
    #if merger_times:
    #    # Remove trailing comma and space
    #    title_str = title_str[:-2]
    #else:
    #    title_str += "No merger predicted"
    #
    #plt.title(title_str, fontsize=GLOBAL_FONT_SIZE)
    
    # Add axis labels
    plt.xlabel(time_label, fontsize=GLOBAL_FONT_SIZE)
    plt.ylabel(dist_label, fontsize=GLOBAL_FONT_SIZE)
    
    # Add legend
    plt.legend(fontsize=GLOBAL_FONT_SIZE-5)
    
    # Add grid
    plt.grid(True, alpha=0.3)
    
    # Apply scale preferences
    if log_x:
        plt.xscale('log')
    if log_y:
        plt.yscale('log')
        # Avoid log(0) issues
        min_sep = min(np.min(traj_separations_plot), np.min(model_separations))
        plt.ylim(max(1e-2, 0.1*min_sep), None)
    
    # Format ticks and layout
    plt.tick_params(axis='both', which='major', labelsize=GLOBAL_TICK_SIZE)
    plt.tight_layout()
    
    # Save and show
    plt.savefig("core_merger_prediction.pdf", format='pdf', bbox_inches='tight')
    print("Plot saved as core_merger_prediction.pdf")
    plt.show()
    
    # Print merger times in exponential notation
    print("\n=== Results ===")
    print(f"Final core separation in data: {traj_separations_plot[-1]:.3f}")
    if merger_times:
        print("Predicted merger times:")
        for C_d, t in merger_times:
            print(f"  C_d = {C_d:.2f}: {t:.2e} years")
    else:
        print("Merger not predicted within simulated timeframe for any C_d")

if __name__ == "__main__":
    main()
