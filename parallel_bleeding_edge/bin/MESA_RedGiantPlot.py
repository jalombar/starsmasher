#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator, LogLocator

# Global style parameters
GLOBAL_FONT_SIZE = 26
GLOBAL_LINE_WIDTH = 4
GLOBAL_TICK_SIZE = 20

def read_history_data(filename='history.data'):
    """Read MESA history.data file and extract evolution data."""
    if not os.path.exists(filename):
        print("\nERROR: File 'history.data' not found!")
        print("You must be in the LOGS folder of your MESA simulation")
        print("Current working directory:", os.getcwd())
        raise FileNotFoundError(f"'{filename}' not found")

    star_age, log_L, log_Teff, log_R, log_g, he_core_mass, log_LH = [], [], [], [], [], [], []

    with open(filename, 'r') as f:
        for line in f:
            if 'star_age' in line and 'log_L' in line and 'log_Teff' in line:
                break

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                cols = line.split()
                star_age.append(float(cols[2]))
                log_L.append(float(cols[37]))
                log_Teff.append(float(cols[36]))
                log_R.append(float(cols[38]))
                log_g.append(float(cols[39]))
                he_core_mass.append(float(cols[48]))
                log_LH.append(float(cols[26]))
            except (ValueError, IndexError):
                continue

    return (
        np.array(star_age), 
        np.array(log_L), 
        np.array(log_Teff), 
        np.array(log_R), 
        np.array(log_g), 
        np.array(he_core_mass), 
        np.array(log_LH)
    )

def plot_evolution(data):
    """Create clean evolution plots with consistent styling."""
    star_age, log_L, log_Teff, log_R, log_g, he_core_mass, log_LH = data
    
    # Apply global style settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=GLOBAL_FONT_SIZE)
    plt.rc('axes', titlesize=GLOBAL_FONT_SIZE)
    plt.rc('axes', labelsize=GLOBAL_FONT_SIZE)
    plt.rc('xtick', labelsize=GLOBAL_TICK_SIZE)
    plt.rc('ytick', labelsize=GLOBAL_TICK_SIZE)
    plt.rc('lines', linewidth=GLOBAL_LINE_WIDTH)

    R = 10**log_R
    Teff = 10**log_Teff

    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    plt.subplots_adjust(wspace=0.35, hspace=0.4)

    # Panel 1: HR Diagram
    ax1 = axes[0, 0]
    ax1.plot(log_Teff, log_L, 'k-')
    ax1.set_xlabel(r'$\log\,T_{\rm eff}$ (K)')
    ax1.set_ylabel(r'$\log\,L/L_{\odot}$')
    ax1.invert_xaxis()
    ax1.grid(True, alpha=0.3)
    ax1.set_title('HR diagram')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())

    # Panel 2: Radius Evolution
    ax2 = axes[0, 1]
    ax2.plot(star_age, R, 'k-')
    ax2.set_xlabel(r'Time (yr)')
    ax2.set_ylabel(r'Radius ($R_{\odot}$)')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Radius evolution')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(LogLocator(subs='auto'))

    # Panel 3: He Core Mass
    ax3 = axes[1, 0]
    ax3.plot(star_age, he_core_mass, 'k-')
    ax3.set_xlabel(r'Time (yr)')
    ax3.set_ylabel(r'He core mass ($M_{\odot}$)')
    ax3.grid(True, alpha=0.3)
    ax3.set_title('He core growth')
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())

    # Panel 4: Nuclear Luminosity
    ax4 = axes[1, 1]
    ax4.plot(star_age, log_LH, 'k-')
    ax4.set_xlabel(r'Time (yr)')
    ax4.set_ylabel(r'$\log\,L_{\rm H}/L_{\odot}$')
    ax4.grid(True, alpha=0.3)
    ax4.set_title('Nuclear luminosity')
    ax4.xaxis.set_minor_locator(AutoMinorLocator())
    ax4.yaxis.set_minor_locator(AutoMinorLocator())

    # Panel 5: Surface Gravity
    ax5 = axes[2, 0]
    ax5.plot(star_age, log_g, 'k-')
    ax5.set_xlabel(r'Time (yr)')
    ax5.set_ylabel(r'$\log\,g$ (cm/s$^2$)')
    ax5.grid(True, alpha=0.3)
    ax5.set_title('Surface gravity')
    ax5.xaxis.set_minor_locator(AutoMinorLocator())
    ax5.yaxis.set_minor_locator(AutoMinorLocator())

    # Panel 6: Temperature
    ax6 = axes[2, 1]
    ax6.plot(star_age, Teff, 'k-')
    ax6.set_xlabel(r'Time (yr)')
    ax6.set_ylabel(r'$T_{\rm eff}$ (K)')
    ax6.set_yscale('log')
    ax6.grid(True, alpha=0.3)
    ax6.set_title('Effective temperature')
    ax6.xaxis.set_minor_locator(AutoMinorLocator())
    ax6.yaxis.set_minor_locator(LogLocator(subs='auto'))

    plt.tight_layout()
    plt.show()

def main():
    try:
        data = read_history_data()
        plot_evolution(data)

    except Exception as e:
        print(f"\nERROR: {e}")

if __name__ == "__main__":
    main()
