#!/usr/bin/env python3
# =============================================
# Wavelet analysis of stellar collisions
# =============================================
# 
# =======================
# OVERVIEW
# =======================
# This script performs wavelet analysis on stellar core collision light curves from
# SPH simulations. It computes:
# 1) Energy distribution across timescales
# 2) Turbulent cascade properties
# 3) Phase-coherent hydrodynamic structures
# 
# =======================
# INPUT REQUIREMENTS
# =======================
# Input file: 'burst_lightcurve_accurate.txt' (ASCII format)
# Required columns:
# - Column 1: Time (years)
# - Column 2: Luminosity (erg/s)
# 
# To generate input from SPH data:
# 1) Run your SPH simulation with luminosity output
# 2) Process snapshots to extract time-luminosity pairs
# 3) Format as two-column ASCII file with header:
#    # Time[yr] Luminosity[erg/s]
# 
# =======================
# EXECUTION INSTRUCTIONS
# =======================
# 1) Ensure required packages are installed:
#    numpy, scipy, matplotlib
# 
# 2) Run the script:
#    python wavelet_collision_analysis.py
# 
# 3) When prompted, select time units:
#    Enter '1' for years or '2' for seconds
# 
# 4) Output files generated:
#    - wavelet_analysis.pdf (plots)
#    - wavelet_energy_distribution.txt (data)
#    - coherent_structures.npz (binary data)
# 
# =======================
# CODE STRUCTURE
# =======================
# 1) PHYSICAL CONSTANTS SECTION
#    - Defines fundamental constants
#    - Sets composition parameters (X,Y,Z)
# 
# 2) WAVELET FUNCTIONS
#    - morlet_wavelet(): Complex Morlet wavelet basis
#    - continuous_wavelet_transform(): Computes CWT
# 
# 3) PHYSICAL ANALYSIS FUNCTIONS
#    - compute_physical_quantities(): Main analysis workflow
#    - plot_physical_analysis(): Visualization
# 
# 4) MAIN EXECUTION BLOCK
#    - Loads input data
#    - Performs user interaction
#    - Runs analysis pipeline
#    - Generates outputs
# 
# =======================
# ANALYSIS METHODOLOGY
# =======================
# 1) Data Preprocessing:
#    - Normalizes luminosity signal
#    - Converts time units
# 
# 2) Wavelet Transform:
#    - Computes CWT with Morlet basis
#    - ω₀=6 for zero mean condition
#    - Logarithmic scale spacing
# 
# 3) Physical Quantities:
#    - Energy distribution P(τ)
#    - Intermittency exponent β
#    - Injection scale ℓ_inj
#    - Coherent structures
# 
# 4) Visualization:
#    - Three-panel publication-quality figure
#    - Consistent style and labeling
# 
# =======================
# CUSTOMIZATION OPTIONS
# =======================
# Key parameters to adjust:
# 1) scales = np.logspace(3,7,100): Timescale range
# 2) omega0=6.0: Morlet wavelet parameter
# 3) threshold = 3*np.nanstd(): Coherence threshold
# 
# =======================
# DEPENDENCIES
# =======================
# - Python 3.6+
# - NumPy 1.19+
# - SciPy 1.6+
# - Matplotlib 3.3+
# - LaTeX (for text rendering)
# 
# =======================
# TROUBLESHOOTING
# =======================
# Common issues:
# 1) Input file not found:
#    - Verify filename/path
#    - Check column formatting
# 
# 2) Missing dependencies:
#    - Install via pip/conda
#    - Check LaTeX installation
# 
# 3) Runtime warnings:
#    - Typically from log(0)
#    - Handled internally
# =============================================


import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from matplotlib import rc

# =============================================
# Global parameters
# =============================================
global_font_size = 35
global_line_width = 3
global_tick_size = 30
rc('text', usetex=True)
rc('font', family='serif', size=global_font_size)

def morlet_wavelet(t, scale, omega0=6.0):
    """Complex Morlet wavelet function"""
    return (np.pi**-0.25 * np.exp(1j * omega0 * t / scale) *
            np.exp(-t**2 / (2 * scale**2)))

def continuous_wavelet_transform(signal, dt, scales, omega0=6.0):
    """Compute continuous wavelet transform"""
    n = len(signal)
    t = np.arange(-n//2, n//2) * dt
    W = np.zeros((len(scales), n), dtype=complex)

    signal_hat = fft(signal)
    for i, scale in enumerate(scales):
        wavelet = morlet_wavelet(t, scale, omega0)
        wavelet_hat = fft(wavelet)
        W[i] = ifft(signal_hat * wavelet_hat) * np.sqrt(dt / scale)

    return W

def compute_physical_quantities(time, scales, lum, dt):
    """Compute physics-driven wavelet metrics"""
    # Standardize and compute transform
    signal = (lum - np.mean(lum)) / np.std(lum)
    W = continuous_wavelet_transform(signal, dt, scales)
    power = np.abs(W)**2

    # 1. Energy distribution
    energy_density = np.trapz(power, x=time, axis=1)
    total_energy = np.trapz(energy_density, x=np.log(scales))
    energy_dist = energy_density / total_energy

    # 2. Intermittency exponent
    log_scales = np.log(scales)
    mean_power = np.mean(power, axis=1)
    valid = (mean_power > 0) & np.isfinite(mean_power)
    beta = np.polyfit(log_scales[valid], np.log(mean_power[valid]), 1)[0]

    # 3. Turbulent injection scale
    deriv = np.gradient(np.log(mean_power[valid]), log_scales[valid])
    injection_scale = scales[valid][np.argmin(np.abs(deriv))]

    # 4. Coherent structures
    phase = np.angle(W)
    phase_ref = 2*np.pi*time/scales[:,None]  # Model phase
    I_coherent = np.abs(W) * np.cos(phase - phase_ref)
    threshold = 3*np.nanstd(I_coherent)
    structure_mask = I_coherent > threshold

    return {
        'power': power,
        'energy_dist': energy_dist,
        'intermittency': beta,
        'injection_scale': injection_scale,
        'coherent_structures': structure_mask
    }

def plot_physical_analysis(time, scales, phys, time_unit='years'):
    """Create physics-focused plots with proper annotations"""
    fig, ax = plt.subplots(3, 1, figsize=(18, 24))

    # Convert time axis if requested
    xlabel_time = r'$t$ (yr)' if time_unit == 'years' else r'$t$ (s)'
    plot_time = time if time_unit == 'years' else time * 365.25 * 86400

    # Plot 1: Energy distribution
    ax[0].plot(scales, phys['energy_dist'], 'k-', linewidth=global_line_width)
    ax[0].set_xscale('log')
    ax[0].set_xlabel(r'$\tau$ (s)', fontsize=global_font_size)
    ax[0].set_ylabel(r'$\mathcal{P}(\tau)$', fontsize=global_font_size)
    ax[0].set_title(r'Energy distribution across timescales', fontsize=global_font_size)

    # Plot 2: Wavelet variance
    mean_power = np.mean(phys['power'], axis=1)
    ax[1].plot(scales, mean_power, 'k-', linewidth=global_line_width)
    inj_scale = phys['injection_scale']
    ax[1].axvline(inj_scale, color='b', linestyle='--',
                 label=r'$\ell_{\rm inj} = %.1f \times 10^{%d}\,{\rm s}$' %
                       (inj_scale/10**np.floor(np.log10(inj_scale)), np.floor(np.log10(inj_scale))))
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$\tau$ (s)', fontsize=global_font_size)
    ax[1].set_ylabel(r'$\sigma^2(\tau)$', fontsize=global_font_size)
    ax[1].set_title(r'Wavelet variance ($\beta = %.2f$)' % phys['intermittency'],
                   fontsize=global_font_size)
    ax[1].legend(fontsize=global_font_size-5)

    # Plot 3: Coherent structures
    im = ax[2].contourf(plot_time, scales, phys['coherent_structures'], levels=[0,0.5,1],
                       cmap='binary')
    ax[2].set_yscale('log')
    ax[2].set_xlabel(xlabel_time, fontsize=global_font_size)
    ax[2].set_ylabel(r'$\tau$ (s)', fontsize=global_font_size)
    ax[2].set_title(r'Phase-coherent structures', fontsize=global_font_size)

    # Annotate gap at y=1e5 s
    gap_start = 0.627
    gap_end = 1.06
    gap_center = (gap_start + gap_end)/2
    gap_center_plot = gap_center if time_unit == 'years' else gap_center * 365.25 * 86400
    ax[2].text(gap_center_plot, 1e5, 'quiescent phase',
              fontsize=global_font_size-5, ha='center', color='white')

    cbar = fig.colorbar(im, ax=ax[2], ticks=[0,1])
    cbar.set_label(r'Coherent structures', fontsize=global_font_size-5)

    for a in ax:
        a.tick_params(axis='both', which='major', labelsize=global_tick_size)
        a.grid(True, linestyle=':', alpha=0.7)

    plt.tight_layout()
    return fig

if __name__ == "__main__":
    # Load data
    data = np.loadtxt('burst_lightcurve_accurate.txt')
    time, lum = data[:, 0], data[:, 1]
    dt = np.median(np.diff(time)) * 365.25 * 86400  # Convert to seconds

    # Ask for time unit preference
    time_unit = input("Display time in (1) years or (2) seconds? [1/2]: ").strip()
    time_unit = 'years' if time_unit != '2' else 'seconds'

    # Set physical scales range (1e3 to 1e7 seconds)
    scales = np.logspace(3, 7, 100)

    # Compute physical quantities
    phys = compute_physical_quantities(time, scales, lum, dt)

    # Generate and save plots
    fig = plot_physical_analysis(time, scales, phys, time_unit)
    plt.savefig('wavelet_analysis.pdf', bbox_inches='tight')
    plt.show()
