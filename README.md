# StarSmasher - a Smoothed-Particle Hydrodynamics code.

Please contact Jamie Lombardi (jamie.lombardi@allegheny.edu) with any questions.

## Overview
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
[Installing StarSmasher](https://github.com/jalombar/starsmasher/blob/bleeding-edge/INSTALLATION.md)

## Example calculation
[Star-star fly-by calculation](https://github.com/jalombar/starsmasher/blob/bleeding-edge/WALKTHROUGH.md)

## DATA VISUALIZATION

I recommend Price's SPLASH software.
It will need to be customized to read in our output file type.
