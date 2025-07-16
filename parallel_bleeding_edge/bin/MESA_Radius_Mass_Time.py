#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator, LogLocator

def read_history_data(filename='history.data'):
    """Read history.data and extract columns"""
    if not os.path.exists(filename):
        print("\nERROR: File 'history.data' not found!")
        print("You must be in the LOGS folder of your MESA simulation")
        print("Current working directory:", os.getcwd())
        raise FileNotFoundError(f"'{filename}' not found")

    star_age, mass, log_R = [], [], []

    with open(filename, 'r') as f:
        for line in f:
            if 'star_age' in line and 'star_mass' in line and 'log_R' in line:
                break

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                cols = line.split()
                star_age.append(float(cols[2]))    # Column 3: star_age
                mass.append(float(cols[4]))       # Column 5: star_mass
                log_R.append(float(cols[38]))     # Column 39: log_R
            except (ValueError, IndexError):
                continue

    return np.array(star_age), np.array(mass), np.array(log_R)

def plot_stellar_evolution():
    # Set up LaTeX rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('axes', titlesize=30)
    plt.rc('axes', labelsize=30)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)

    try:
        star_age, mass, log_R = read_history_data()
        R = 10**log_R  # Convert to linear solar radii

        # User options
        xlog = input("Log scale for x-axis (time)? (y/n): ").strip().lower() == 'y'
        y1log = input("Log scale for left y-axis (radius)? (y/n): ").strip().lower() == 'y'
        y2log = input("Log scale for right y-axis (mass)? (y/n): ").strip().lower() == 'y'

        fig, ax1 = plt.subplots(figsize=(14, 8))
        ax2 = ax1.twinx()

        # Plot radius on left axis (solid black line)
        ax1.plot(star_age, R, 'k-', linewidth=2.5)
        ax1.set_xlabel(r'Time (yr)', fontsize=30)
        ax1.set_ylabel(r'Radius ($R_{\odot}$)', fontsize=30)

        # Plot mass on right axis (dashed black line with alpha)
        ax2.plot(star_age, mass, 'k--', linewidth=3.5, alpha=0.3)
        ax2.set_ylabel(r'Mass ($M_{\odot}$)', fontsize=30)

        # Apply log scales
        if xlog:
            ax1.set_xscale('log')
            ax1.xaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))
        else:
            ax1.xaxis.set_minor_locator(AutoMinorLocator(5))

        # Configure left y-axis (radius)
        if y1log:
            ax1.set_yscale('log')
            ax1.yaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))
        else:
            ax1.yaxis.set_minor_locator(AutoMinorLocator(5))

        # Configure right y-axis (mass)
        if y2log:
            ax2.set_yscale('log')
            ax2.yaxis.set_minor_locator(LogLocator(subs=np.linspace(0.1, 0.9, 9)))
        else:
            ax2.yaxis.set_minor_locator(AutoMinorLocator(5))

        # Configure ticks and grid (all black)
        for ax in [ax1, ax2]:
            ax.tick_params(axis='y', colors='k')
            ax.tick_params(which='both', width=1.5, colors='k')
            ax.tick_params(which='major', length=8)
            ax.tick_params(which='minor', length=4)
            ax.grid(True, which='major', alpha=0.4)
            ax.grid(True, which='minor', alpha=0.2)
            ax.minorticks_on()

        plt.tight_layout()
        plt.show()

    except FileNotFoundError:
        return
    except Exception as e:
        print(f"An error occurred: {e}")

plot_stellar_evolution()
