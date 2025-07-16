#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import glob
import sys
import os
from matplotlib.colors import LinearSegmentedColormap

# =============================================
# PHYSICAL UNIT CONVERSION FACTORS (from log0.sph)
# =============================================
CODE_UNIT_MASS = 1.9891e33    # g (1 Msun)
CODE_UNIT_LENGTH = 6.9599e10  # cm (1 Rsun)
PC_PER_CM = 1/3.086e18        # Convert cm to parsecs
RSUN_PER_CM = 1/6.9599e10     # Convert cm to solar radii

# Core radius parameters
CORE_RADIUS_RSUN = 0.05       # Default core radius in solar radii
CORE_RADIUS_PC = CORE_RADIUS_RSUN * CODE_UNIT_LENGTH * RSUN_PER_CM * PC_PER_CM

# =============================================
# GLOBAL STYLE PARAMETERS
# =============================================
GLOBAL_FONT_SIZE = 35
GLOBAL_LINE_WIDTH = 3
GLOBAL_TICK_SIZE = 30

# LaTeX setup with global styles
rc('text', usetex=True)
rc('font', family='serif', size=GLOBAL_FONT_SIZE)
rc('axes', titlesize=GLOBAL_FONT_SIZE, labelsize=GLOBAL_FONT_SIZE)
rc('xtick', labelsize=GLOBAL_TICK_SIZE)
rc('ytick', labelsize=GLOBAL_TICK_SIZE)

def get_unit_preferences():
    """Get user preferences for unit systems"""
    print("\nAvailable unit systems:")
    print("1: Normalized units (R/R0, T/Ttot)")
    print("2: Physical units (parsecs, solar masses)")
    print("3: Code units (simulation native)")
    print("4: Physical units (solar radii, solar masses)")
    
    while True:
        choice = input("Select unit system [1/2/3/4]: ").strip()
        if choice in ['1', '2', '3', '4']:
            return int(choice)
        print("Invalid choice. Please enter 1, 2, 3, or 4")

def create_fading_colormap(base_color, n_segments=100):
    """Create a colormap that fades from base_color to transparent"""
    return LinearSegmentedColormap.from_list(
        'fading', 
        [(*base_color, 1.0), (*base_color, 0.1)],
        N=n_segments
    )

def read_sph_files(file_pattern):
    """Read all SPH files and extract positions of itype=2 particles"""
    files = sorted(glob.glob(file_pattern))
    if not files:
        print("\nERROR: No ASCII files found matching pattern 'out*.sph.ascii'")
        print("You can create these files by running:")
        print("  StSm_TrajectoryAscii.sh")
        print("which executes:")
        print("  `which splash` to ascii -f starsmasher --convert=7,1,2,3,4,5,6,14 out*.sph")
        sys.exit(1)

    particle1 = {'x': [], 'y': [], 'z': []}
    particle2 = {'x': [], 'y': [], 'z': []}
    times = []

    for fname in files:
        with open(fname, 'r') as f:
            # Get time from file
            time = 0.0
            for line in f:
                if line.startswith('# time:'):
                    time_line = next(f).strip()
                    if time_line.startswith('#'):
                        time_values = time_line[1:].strip().split()
                    else:
                        time_values = time_line.split()
                    if time_values:
                        try:
                            time = float(time_values[0])
                        except ValueError:
                            pass
                    break

            current_particles = {'x': [], 'y': [], 'z': []}
            for line in f:
                if not line.startswith('#') and line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 8 and int(parts[-1]) == 2:
                        current_particles['x'].append(float(parts[1]))
                        current_particles['y'].append(float(parts[2]))
                        current_particles['z'].append(float(parts[3]))

            if len(current_particles['x']) >= 2:
                times.append(time)
                particle1['x'].append(current_particles['x'][0])
                particle1['y'].append(current_particles['y'][0])
                particle1['z'].append(current_particles['z'][0])
                particle2['x'].append(current_particles['x'][1])
                particle2['y'].append(current_particles['y'][1])
                particle2['z'].append(current_particles['z'][1])

    return (np.array(times),
            np.array(particle1['x']), np.array(particle1['y']), np.array(particle1['z']),
            np.array(particle2['x']), np.array(particle2['y']), np.array(particle2['z']))

def calculate_distances(x1, y1, z1, x2, y2, z2):
    """Calculate distances between corresponding points"""
    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

def write_trajectory_file(times, x1, y1, z1, x2, y2, z2, distances, filename):
    """Write trajectory data to file with distance column"""
    with open(filename, 'w') as f:
        f.write('# Trajectory of compact objects with separation distance\n')
        f.write('# time x1 y1 z1 x2 y2 z2 distance\n')
        for i in range(len(times)):
            f.write(f'{times[i]:.8E} {x1[i]:.8E} {y1[i]:.8E} {z1[i]:.8E} '
                    f'{x2[i]:.8E} {y2[i]:.8E} {z2[i]:.8E} '
                    f'{distances[i]:.8E}\n')

def read_trajectory_file(filename):
    """Read existing trajectory file with distance column"""
    times, x1, y1, z1, x2, y2, z2, distances = [], [], [], [], [], [], [], []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 8:
                    times.append(float(parts[0]))
                    x1.append(float(parts[1]))
                    y1.append(float(parts[2]))
                    z1.append(float(parts[3]))
                    x2.append(float(parts[4]))
                    y2.append(float(parts[5]))
                    z2.append(float(parts[6]))
                    distances.append(float(parts[7]))
    return (np.array(times),
            np.array(x1), np.array(y1), np.array(z1),
            np.array(x2), np.array(y2), np.array(z2),
            np.array(distances))

def plot_2d_trajectory(x, y, xlabel, ylabel, d0, times, color, linestyle, unit_choice):
    """Plot 2D trajectory with fading color"""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Apply unit conversion
    if unit_choice == 2:  # Physical units (pc)
        x_plot = x * CODE_UNIT_LENGTH * PC_PER_CM
        y_plot = y * CODE_UNIT_LENGTH * PC_PER_CM
        d0_plot = d0 * CODE_UNIT_LENGTH * PC_PER_CM
    elif unit_choice == 4:  # Physical units (Rsun)
        x_plot = x * CODE_UNIT_LENGTH * RSUN_PER_CM
        y_plot = y * CODE_UNIT_LENGTH * RSUN_PER_CM
        d0_plot = d0 * CODE_UNIT_LENGTH * RSUN_PER_CM
    else:  # Normalized or code units
        x_plot = x
        y_plot = y
        d0_plot = d0

    cmap = create_fading_colormap(color)
    norm = plt.Normalize(times.min(), times.max())

    for i in range(len(x_plot)-1):
        ax.plot(x_plot[i:i+2]/d0_plot, y_plot[i:i+2]/d0_plot,
               linestyle=linestyle,
               color=cmap(norm(times[i])),
               linewidth=GLOBAL_LINE_WIDTH)

    ax.set_xlabel(xlabel, fontsize=GLOBAL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=GLOBAL_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=GLOBAL_TICK_SIZE)
    plt.tight_layout()
    plt.show()

def plot_both_trajectories(x1, y1, x2, y2, xlabel, ylabel, d0, times, unit_choice):
    """Plot both trajectories with fading colors and different line styles"""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Apply unit conversion
    if unit_choice == 2:  # Physical units (pc)
        x1_plot = x1 * CODE_UNIT_LENGTH * PC_PER_CM
        y1_plot = y1 * CODE_UNIT_LENGTH * PC_PER_CM
        x2_plot = x2 * CODE_UNIT_LENGTH * PC_PER_CM
        y2_plot = y2 * CODE_UNIT_LENGTH * PC_PER_CM
        d0_plot = d0 * CODE_UNIT_LENGTH * PC_PER_CM
    elif unit_choice == 4:  # Physical units (Rsun)
        x1_plot = x1 * CODE_UNIT_LENGTH * RSUN_PER_CM
        y1_plot = y1 * CODE_UNIT_LENGTH * RSUN_PER_CM
        x2_plot = x2 * CODE_UNIT_LENGTH * RSUN_PER_CM
        y2_plot = y2 * CODE_UNIT_LENGTH * RSUN_PER_CM
        d0_plot = d0 * CODE_UNIT_LENGTH * RSUN_PER_CM
    else:  # Normalized or code units
        x1_plot = x1
        y1_plot = y1
        x2_plot = x2
        y2_plot = y2
        d0_plot = d0

    cmap1 = create_fading_colormap((0.2, 0.4, 0.8))  # Blue-ish
    cmap2 = create_fading_colormap((0.8, 0.4, 0.2))  # Red-ish
    norm = plt.Normalize(times.min(), times.max())

    for i in range(len(x1_plot)-1):
        ax.plot(x1_plot[i:i+2]/d0_plot, y1_plot[i:i+2]/d0_plot,
               linestyle='-',
               color=cmap1(norm(times[i])),
               linewidth=GLOBAL_LINE_WIDTH,
               label='Core 1' if i == 0 else "")

    for i in range(len(x2_plot)-1):
        ax.plot(x2_plot[i:i+2]/d0_plot, y2_plot[i:i+2]/d0_plot,
               linestyle='--',
               color=cmap2(norm(times[i])),
               linewidth=GLOBAL_LINE_WIDTH,
               label='Core 2' if i == 0 else "")

    ax.set_xlabel(xlabel, fontsize=GLOBAL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=GLOBAL_FONT_SIZE)
    ax.legend(fontsize=GLOBAL_FONT_SIZE-4)
    ax.tick_params(axis='both', which='major', labelsize=GLOBAL_TICK_SIZE)
    plt.tight_layout()
    plt.show()

def plot_distance_vs_time(times, distances, d0, t_tot, unit_choice):
    """Plot distance between particles with user-selected options"""
    # Get core radii sum from user
    core_radius = None
    if unit_choice in [2, 4]:  # Only for physical units
        use_core_line = get_user_choice(
            "Do you want to plot a line for the sum of core radii? [y/n]: ", ['y', 'n'])
        if use_core_line == 'y':
            if unit_choice == 2:  # pc
                default_radius = 2 * CORE_RADIUS_PC
                prompt = f"Enter sum of core radii (pc) [default=2xRcore={default_radius:.2e}pc]: "
            else:  # Rsun
                default_radius = 2 * CORE_RADIUS_RSUN
                prompt = f"Enter sum of core radii (Rsun) [default=2xRcore={default_radius:.2f}]: "
            
            try:
                input_str = input(prompt)
                core_radius = float(input_str) if input_str else default_radius
            except ValueError:
                print(f"Invalid input. Using default value: {default_radius}")
                core_radius = default_radius

    if unit_choice == 1:  # Normalized units
        x = times / t_tot
        y = distances / d0
        xlabel = r'$T/T_{\rm tot}$'
        ylabel = r'$R/d_0$'
    elif unit_choice == 2:  # Physical units (pc)
        x = times  # No time conversion available
        y = distances * CODE_UNIT_LENGTH * PC_PER_CM
        xlabel = 'Time (code units)'
        ylabel = 'Distance (pc)'
        print(f"\nPhysical units conversion:")
        print(f"  1 code length = {CODE_UNIT_LENGTH*PC_PER_CM:.3e} pc")
        print(f"  1 code mass = {CODE_UNIT_MASS/1.9891e33:.3f} Msun")
    elif unit_choice == 4:  # Physical units (Rsun)
        x = times
        y = distances * CODE_UNIT_LENGTH * RSUN_PER_CM
        xlabel = 'Time (code units)'
        ylabel = r'Distance ($R_\odot$)'
        print(f"\nPhysical units conversion:")
        print(f"  1 code length = {CODE_UNIT_LENGTH*RSUN_PER_CM:.3f} Rsun")
        print(f"  1 code mass = {CODE_UNIT_MASS/1.9891e33:.3f} Msun")
    else:  # Code units
        x = times
        y = distances
        xlabel = 'Time (code units)'
        ylabel = 'Distance (code units)'

    log_choice = get_user_choice(
        "Log scale options:\n"
        "  (1) Linear x and y\n"
        "  (2) Log x\n"
        "  (3) Log y\n"
        "  (4) Log x and y\n"
        "Your choice [1/2/3/4]: ", ['1', '2', '3', '4'])

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.plot(x, y, '-', color='black', linewidth=GLOBAL_LINE_WIDTH)
    
    # Plot core radius line if requested
    if core_radius is not None:
        if unit_choice == 2:  # pc
            label = r'$2 \times R_{\rm core} = %.2e\,{\rm pc}$' % core_radius
        else:  # Rsun
            label = r'$2 \times R_{\rm core} = %.2f\,R_\odot$' % core_radius
            
        ax.axhline(y=core_radius, color='gray', linestyle='--', 
                  alpha=0.5, linewidth=GLOBAL_LINE_WIDTH-1,
                  label=label)
        ax.legend(fontsize=GLOBAL_FONT_SIZE-4)

    ax.set_xlabel(xlabel, fontsize=GLOBAL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=GLOBAL_FONT_SIZE)

    if log_choice in ['2', '4']:
        ax.set_xscale('log')
    if log_choice in ['3', '4']:
        ax.set_yscale('log')

    ax.tick_params(axis='both', which='major', labelsize=GLOBAL_TICK_SIZE)
    plt.tight_layout()
    plt.show()

def get_user_choice(prompt, valid_options):
    """Get validated user input"""
    while True:
        choice = input(prompt).lower()
        if choice in valid_options:
            return choice
        print(f"Invalid option. Please choose from: {', '.join(valid_options)}")

def main():
    try:
        print("=== Core Merger Analysis ===")
        print(f"Using conversion factors from log0.sph:")
        print(f"  1 code length = {CODE_UNIT_LENGTH:.3e} cm")
        print(f"  1 code mass = {CODE_UNIT_MASS:.3e} g (1 Msun)")
        
        unit_choice = get_unit_preferences()
        output_file = 'compact_objects_trajectory.ascii'

        if os.path.exists(output_file):
            print(f"\nFound existing trajectory file: {output_file}")
            choice = get_user_choice(
                "Do you want to:\n"
                "  (r)egenerate from SPH files\n"
                "  (p)lot from existing file\n"
                "  (q)uit\n"
                "Your choice [r/p/q]: ", ['r', 'p', 'q'])

            if choice == 'q':
                return
            elif choice == 'p':
                times, x1, y1, z1, x2, y2, z2, distances = read_trajectory_file(output_file)
                if len(times) == 0:
                    print("ERROR: No valid data in trajectory file")
                    return
            else:  # 'r'
                times, x1, y1, z1, x2, y2, z2 = read_sph_files('out*.sph.ascii')
                distances = calculate_distances(x1, y1, z1, x2, y2, z2)
                write_trajectory_file(times, x1, y1, z1, x2, y2, z2, distances, output_file)
                print(f"Trajectory data written to {output_file}")
        else:
            times, x1, y1, z1, x2, y2, z2 = read_sph_files('out*.sph.ascii')
            distances = calculate_distances(x1, y1, z1, x2, y2, z2)
            write_trajectory_file(times, x1, y1, z1, x2, y2, z2, distances, output_file)
            print(f"Trajectory data written to {output_file}")

        if len(times) == 0:
            print("ERROR: No valid trajectory data")
            return

        d0 = distances[0] if len(distances) > 0 else 1.0
        t_tot = times[-1] if len(times) > 0 else 1.0
        print(f"\nInitial separation: {d0:.3f} code units")
        print(f"Total time: {t_tot:.3f} code units")

        plot_choice = get_user_choice(
            "Select plot type:\n"
            "  (1) Particle trajectory\n"
            "  (2) Distance vs time\n"
            "Your choice [1/2]: ", ['1', '2'])

        if plot_choice == '1':
            particle_choice = get_user_choice(
                "Select which particle to plot (1/2/both) and plane (xy/xz/yz) [e.g., 'bothxy' or '2xy']: ",
                ['1xy', '1xz', '1yz', '2xy', '2xz', '2yz', 'bothxy', 'bothxz', 'bothyz'])

            particle = particle_choice.replace('xy', '').replace('xz', '').replace('yz', '')
            plane = particle_choice[-2:]

            if particle == 'both':
                if plane == 'xy':
                    plot_both_trajectories(x1, y1, x2, y2, r'$x/d_0$', r'$y/d_0$', d0, times, unit_choice)
                elif plane == 'xz':
                    plot_both_trajectories(x1, z1, x2, z2, r'$x/d_0$', r'$z/d_0$', d0, times, unit_choice)
                else:
                    plot_both_trajectories(y1, z1, y2, z2, r'$y/d_0$', r'$z/d_0$', d0, times, unit_choice)
            else:
                if particle == '1':
                    x, y, z = x1, y1, z1
                    color = (0.2, 0.4, 0.8)
                    linestyle = '-'
                else:
                    x, y, z = x2, y2, z2
                    color = (0.8, 0.4, 0.2)
                    linestyle = '--'

                if plane == 'xy':
                    plot_2d_trajectory(x, y, r'$x/d_0$', r'$y/d_0$', d0, times, color, linestyle, unit_choice)
                elif plane == 'xz':
                    plot_2d_trajectory(x, z, r'$x/d_0$', r'$z/d_0$', d0, times, color, linestyle, unit_choice)
                else:
                    plot_2d_trajectory(y, z, r'$y/d_0$', r'$z/d_0$', d0, times, color, linestyle, unit_choice)
        else:
            plot_distance_vs_time(times, distances, d0, t_tot, unit_choice)

    except Exception as e:
        print(f"ERROR: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
