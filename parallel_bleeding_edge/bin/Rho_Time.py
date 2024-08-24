#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Constants for conversion:
# time, converted to seconds, has a conversion factor = sqrt(Rsun^3/(G*Msun)))
# rho, converted to g/cm^3, has a conversion factor = Msun/Rsun^3
time_conversion_factor = 1593.0  # seconds
density_conversion_factor = 5.91  # g/cm^3
seconds_in_year = 3.154e7  # seconds per year
seconds_in_day = 86400  # seconds per day

# Extract and convert the time data
time_data_years = []
time_data_days = []
with open('energy0.sph', 'r') as energy_file:
    for line in energy_file:
        time_in_seconds = float(line.split()[0]) * time_conversion_factor
        time_in_years = time_in_seconds / seconds_in_year
        time_in_days = time_in_seconds / seconds_in_day
        time_data_years.append(time_in_years)
        time_data_days.append(time_in_days)

# Extract and convert the density data
density_data = []
with open('log0.sph', 'r') as log_file:
    for line in log_file:
        if "rhomax= " in line:
            density = float(line.split()[1]) * density_conversion_factor
            density_data.append(density)

# Ensure both columns have the same number of rows
min_length = min(len(time_data_years), len(density_data))
time_data_years = time_data_years[:min_length]
time_data_days = time_data_days[:min_length]
density_data = density_data[:min_length]

# Save the combined data to a file
with open('DENS_TIME.txt', 'w') as out_file:
    out_file.write("#1: Time (yrs) #2: Density (g/cm^3)\n")
    for t, d in zip(time_data_years, density_data):
        out_file.write(f"{t}\t{d}\n")

# Plot the data using Matplotlib with LaTeX rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax1 = plt.subplots(figsize=(10, 6))

# Main plot
ax1.plot(time_data_years, density_data, 'bo-', alpha=0.7)
ax1.set_xlabel(r'Time (yrs)', fontsize=24)
ax1.set_ylabel(r'Density $\rho$ (g/cm$^3$)', fontsize=24)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

# Upper x-axis with time in days
ax2 = ax1.twiny()
ax2.plot(time_data_days, density_data, 'bo-', alpha=0.0)  # Use alpha=0 to avoid showing the line again
ax2.set_xlabel(r'Time (days)', fontsize=24)
ax2.tick_params(axis='x', labelsize=20)

# Title and grid
#ax1.grid(True)
plt.tight_layout()

# Display the plot
plt.show()

