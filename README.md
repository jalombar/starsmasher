# StarSmasher - a Smoothed Particle Hydrodynamics code.

Please contact Jamie Lombardi (jamie.lombardi@allegheny.edu) with any questions.

## Overview
StarSmasher (originally StarCrash) is a hydrodynamics code developed originally by Rasio (1991), updated and maintained as described in Lombardi et al. (1999), Faber & Rasio (2000), and Lombardi et al. (2006). 
_SPH_ is a Lagrangian particle method that approximates a continuous fluid as discrete nodes, each carrying various parameters such as mass, position, velocity, pressure, and temperature. 
In an SPH simulation the resolution scales with the particle density, and StarSmasher is able to handle both equal-mass and equal number-density particle models.

StarSmasher solves for hydro forces by calculating the pressure for each particle as a function of the particle's properties - density, internal energy, and internal properties (e.g. temperature and mean molecular weight).
The code now implements variational equations of motion and libraries to calculate the gravitational forces between particles using direct summation on NVIDIA graphics cards as described in Gaburov et al. (2010b). 
Using a direct summation instead of a tree-based algorithm for gravity increases the accuracy of the gravity calculations at the cost of speed (Gaburov et al. 2010a). 
The code uses a cubic spline (Monaghan & Lattanzio 1985) for the smoothing kernel and an artificial viscosity prescription coupled with a Balsara Switch (Balsara 1995) to prevent unphysical interparticle penetration. 
The code also implements an artificial relaxation force to the equations of motion to add a drag term to the calculated accelerations during relaxation integrations.

Although this Starsmasher code is not well documented, there is documentation for Starcrash, which is a previous version of this [code](http://ciera.northwestern.edu/StarCrash/).
Although Starcrash is not the same as Starsmasher, the variable names, input files, output files, and parallelization strategy is very similar.
The Starcrash documentation is available as the file usersmanual.pdf in the misc subdirectory and may be somewhat helpful to peruse.
You will be most interested in pages 31 - 36, which talks about input and output.
Remember, however, that this Starsmasher SPH code is different in many ways.
The algorithms of Starsmasher are mostly described in [Evghenii et al. (2010)](http://adsabs.harvard.edu/abs/2010MNRAS.402..105G).
The AV scheme is described in Ponce et al. (2011).

## Roadmap
- [ ] Integrate equations of state to handle planet-planet collisions
- [ ] Link to data analysis libraries
- [ ] Improve documentation
- [ ] Update README
  - [ ] Add in references to papers and algorithms

## Installation
[Installing StarSmasher](https://github.com/jalombar/starsmasher/blob/master/documentation/installation.md)

## Example calculation
[Star-star fly-by calculation](https://github.com/jalombar/starsmasher/blob/master/documentation/walkthroughs/star_star_flyby.md)

## Publications using StarSmasher
The first iteration of StarSmasher appeared in 1991, used to study the collisions of giant stars with compact objects.
Since then, StarSmasher has been used to perform a variety of hydrodynamical calculations from common-envelope evolution and collisions in  galactic clusters involving binary stars at various stages and compact objects, and more recently, planet-planet collisions in Kepler Multis.

[List of publications](https://github.com/jalombar/starsmasher/blob/master/documentation/publications.md)

## DATA VISUALIZATION
I recommend Price's SPLASH software (Price 2007).
It will need to be customized to read in our output file type.

## Citations
Balsara, D. S. 1995, Journal of Computational Physics, 121, 357

Gaburov, E., Bédorf, J., & Portegies Zwart, S. 2010, Procedia Computer Science, volume 1, p. 1119-1127

Gaburov, E., Lombardi, J. C., Jr., & Portegies Zwart, S. 2010, MNRAS, 402, 105

Faber, J. A., & Rasio, F. A. 2000, Phys. Rev. D, 62, 064012

Lombardi, J. C., Sills, A., Rasio, F. A., & Shapiro, S. L. 1999, Journal of Computational Physics, 152, 687

Lombardi, J. C., Jr., Proulx, Z. F., Dooley, K. L., et al. 2008, ApJ, 640, 44

Monaghan, J. J., & Lattanzio, J. C. 1985, A&A, 149, 135

Price, D. J. 2007, PASA, 24, 159

Rasio, F. A. 1991, Ph.D. Thesis

