#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator, LogLocator

def read_history_data(filename='history.data'):
    """Read history.data and extract columns with proper indexing"""
    if not os.path.exists(filename):
        print("\nERROR: File 'history.data' not found!")
        print("You must be in the LOGS folder of your MESA simulation")
        print("Current working directory:", os.getcwd())
        raise FileNotFoundError(f"'{filename}' not found")

    star_age, center_T, center_rho = [], [], []
    data_started = False
    header_line = ""

    with open(filename, 'r') as f:
        for line in f:
            # Find the header line containing column names
            if 'log_cntr_T' in line and 'log_cntr_Rho' in line:
                header_line = line
                break

        if not header_line:
            raise ValueError("Could not find required columns in header")

        # Now find the start of numerical data
        for line in f:
            if line.strip() and not line.strip().startswith('#'):
                try:
                    # First numerical line indicates data start
                    cols = line.split()
                    float(cols[0])
                    data_started = True
                    break
                except (ValueError, IndexError):
                    continue

        if not data_started:
            raise ValueError("Could not find numerical data in file")

        # Process all numerical data lines
        while data_started:
            try:
                cols = line.split()
                if len(cols) < 45:  # Need at least 45 columns
                    line = next(f)
                    continue

                age = float(cols[2])       # Column 3: star_age (1-based)
                log_T = float(cols[43])    # Column 44: log_cntr_T (1-based)
                log_rho = float(cols[42])  # Column 43: log_cntr_Rho (1-based)

                # Convert log values to linear scale
                star_age.append(age)
                center_T.append(10**log_T)
                center_rho.append(10**log_rho)

                line = next(f)
            except (ValueError, IndexError, StopIteration):
                break

    return np.array(star_age), np.array(center_T), np.array(center_rho)

def save_plot_data(time, temp, rho, filename='MESA_Central_Temp_rho_Time.txt'):
    """Save the data used for plotting to a file"""
    with open(filename, 'w') as f:
        f.write("# Time (yr)    Central_T (K)    Central_Rho (g/cm^3)\n")
        f.write("# Values calculated as 10^log_cntr_T and 10^log_cntr_Rho\n")
        for t, T, ρ in zip(time, temp, rho):
            f.write(f"{t:.6e}    {T:.6e}    {ρ:.6e}\n")
    print(f"\nData written to {filename}")
    print("First 5 rows:")
    print(f"{time[0]:.6e}    {temp[0]:.6e}    {rho[0]:.6e}")
    print(f"{time[1]:.6e}    {temp[1]:.6e}    {rho[1]:.6e}")
    print(f"{time[2]:.6e}    {temp[2]:.6e}    {rho[2]:.6e}")
    print(f"{time[3]:.6e}    {temp[3]:.6e}    {rho[3]:.6e}")
    print(f"{time[4]:.6e}    {temp[4]:.6e}    {rho[4]:.6e}")

def plot_central_properties():
    # Set up LaTeX rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('axes', titlesize=30)
    plt.rc('axes', labelsize=30)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)

    try:
        star_age, center_T, center_rho = read_history_data()

        # Save the data to file
        save_plot_data(star_age, center_T, center_rho)

        # User options
        xlog = input("Log scale for x-axis (time)? (y/n): ").strip().lower() == 'y'
        y1log = input("Log scale for left y-axis (temperature)? (y/n): ").strip().lower() == 'y'
        y2log = input("Log scale for right y-axis (density)? (y/n): ").strip().lower() == 'y'

        fig, ax1 = plt.subplots(figsize=(14, 8))
        ax2 = ax1.twinx()

        # Plot central temperature on left axis (solid black)
        ax1.plot(star_age, center_T, 'k-', linewidth=2.5)
        ax1.set_xlabel(r'Time (yr)', fontsize=30)
        ax1.set_ylabel(r'Central T (K)', fontsize=30)

        # Plot central density on right axis (dashed black)
        ax2.plot(star_age, center_rho, 'k--', linewidth=3.5, alpha=0.3)
        ax2.set_ylabel(r'Central density (g cm$^{-3}$)', fontsize=30)

        # Apply log scales
        if xlog:
            ax1.set_xscale('log')
            ax1.xaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))
        else:
            ax1.xaxis.set_minor_locator(AutoMinorLocator(5))

        if y1log:
            ax1.set_yscale('log')
            ax1.yaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))

        if y2log:
            ax2.set_yscale('log')
            ax2.yaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))

        # Configure ticks and grid
        for ax in [ax1, ax2]:
            ax.tick_params(which='both', width=1.5)
            ax.tick_params(which='major', length=8)
            ax.tick_params(which='minor', length=4)
            ax.grid(True, which='major', alpha=0.3)
            ax.grid(True, which='minor', alpha=0.1)
            ax.minorticks_on()

        plt.tight_layout()
        plt.show()

    except FileNotFoundError:
        return
    except Exception as e:
        print(f"An error occurred: {e}")

plot_central_properties()
