# StarSmasher post-pocessing code: Charles Gibson's sph-to-mesa

This repository contains post-processing codes developed in Python to analyze
output data from the StarSmasher simulation software. The primary focus of
these codes is to construct 1-D stellar models based on the 3-D data produced
by StarSmasher. The generated output files can be utilized to evolve the
results of Smoothed Particle Hydrodynamics (SPH) collisions within the MESA
(Modules for Experiments in Stellar Astrophysics) framework.

The code `sph-to-mesa` provided in the `tools2` directory is a reproduction of
the program originally developed by Charles Gibson. This work serves as a
foundation for the current repository, contributing to the ongoing development
and refinement of post-processing tools for stellar simulations.

We have decided to include it here to have everything centralised, but do not
get it wrong; the author is Charles Gibson and he deserves all of the credit
for the work. This is his original repository:

<a href="https://github.com/charlie-gibson/sph-to-mesa">https://github.com/charlie-gibson/sph-to-mesa</a>

## Features and Limitations

### Particle Management and SPH Data

Currently, the code does not support "jump runs" in the SPH simulation. This
limitation necessitates the presence of all particles throughout the entire
simulation. A potential workaround is the `cc` variable in the StarSmasher
output, which stores parent-star data for the simulation, thereby aiding in
post-processing.

### Spherical Symmetry in MESA

MESA operates under the assumption of spherical symmetry, while StarSmasher and
similar SPH models accommodate 3-D asymmetries. Consequently, it is crucial
that the SPH-derived star approximates a spherical shape. Although the provided
code can generate files that MESA may read successfully, the subsequent stellar
evolution modeled by MESA will assume spherical symmetry, even if the
underlying data represent non-spherical distributions. This assumption could
lead to inaccuracies in the stellar evolution, particularly when dealing with
complex structures, such as extended tails or non-uniform stellar shapes.

### Equation of State (EOS) and Metallicity Considerations

To ensure compatibility with MESA's definitions of thermodynamic values, it is
recommended to use MESA EOS tables. The current implementation primarily
supports a metallicity of $Z=0.02$. While this is sufficient for many cases,
extending this functionality to other metallicities requires the creation of
additional EOS table files, a process that, while straightforward, can be
labor-intensive. This becomes particularly important for older giant stars
(e.g., Helium-burning stars), where metallicities can vary significantly
between the core and outer layers. Ideally, MESA EOS tables would reflect this
variability across different metallicity values.

It is worth noting that the analytic EOS method can provide a reasonable
approximation for thermodynamic values used by MESA. However, this method may
lack precision in cases involving extreme models, such as very high-density
scenarios where degeneracy pressures are significant, or very
low-density/temperature particles where recombination and ionization effects
become critical.


