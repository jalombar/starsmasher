# StarSmasher - a Smoothed-Particle Hydrodynamics code.

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
[Installing StarSmasher](https://github.com/jalombar/starsmasher/blob/bleeding-edge/documentation/installation.md)

## Example calculation
[Star-star fly-by calculation](https://github.com/jalombar/starsmasher/blob/bleeding-edge/documentation/walkthroughs/star_star_flyby.md)

## Papers using StarSmasher
### Algorithms
[An Efficient Radiative Cooling Approximation for Use in Hydrodynamic Simulations](http://adsabs.harvard.edu/abs/2015MNRAS.447...25L)

[Smoothed Particle Hydrodynamics Calculations of Stellar Interactions](http://adsabs.harvard.edu/abs/1999JCoAM.109..213R)

[Test of Spurious Transport in Smoothed Particle Hydrodynamics](http://adsabs.harvard.edu/abs/1999JCoPh.152..687L)

[The Importance of Realistic Starting Models for Hydrodynamics Simulations of Stellar Collisions](http://adsabs.harvard.edu/abs/1997ApJ...484L..51S)

### Planet collisions
[Hydrodynamics of Collisions Between Sub-Neptunes](https://arxiv.org/abs/1707.01628)

[Dynamics and Collisional Evolution of Closely Packed Planetary Systems](http://adsabs.harvard.edu/abs/2017MNRAS.470.4145H)

### Binary Stability (Especially Common-Envelope Evolution)
[The Influence of Neutrinos on r-process Nucleosynthesis in the Ejecta of Black Hole-Neutron Star Mergers](http://adsabs.harvard.edu/abs/2017MNRAS.464.3907R)

[Stability and Coalescence of Massive Twin Binaries](http://adsabs.harvard.edu/abs/2015ApJ...806..135H)

[Recombination Energy in Double White Dwarf Formation](http://adsabs.harvard.edu/abs/2015MNRAS.450L..39N)

[Testing Common Envelopes on Double White Dwarf Binaries](http://adsabs.harvard.edu/abs/2015ASPC..493..469N)

[Implications of the Delayed 2013 Outburst of ESO 243-49 HLX-1](http://adsabs.harvard.edu/abs/2014ApJ...793..105G)

[V1309 Sco--Understanding a Merger](http://adsabs.harvard.edu/abs/2014ApJ...786...39N)

[Identification of the Long-Sought Common-Envelope Events](http://adsabs.harvard.edu/abs/2013Sci...339..433I)

[Twin Binaries: Studies of Stability, Mass Transfer, and Coalescence](http://adsabs.harvard.edu/abs/2011ApJ...737...49L)

[Binary Neutron Star Mergers](http://adsabs.harvard.edu/abs/2012LRR....15....8F)

[Dynamical Evolution of Black Hole-Neutron Star Binaries in General Relativity: Simulations of Tidal Disruption](http://adsabs.harvard.edu/abs/2006PhRvD..73b4012F)

[On the Formation and Evolution of Common Envelope Systems](http://adsabs.harvard.edu/abs/1996ApJ...471..366R)

[Hydrodynamics of Binary Coalescence. 2: Polytropes with Gamma=5/3](http://adsabs.harvard.edu/abs/1995ApJ...438..887R)

[Hydrodynamics of Binary Coalescence. 1: Polytropes with Stiff Equations of State](http://adsabs.harvard.edu/abs/1994ApJ...432..242R)

[Hydrodynamical Evolution of Coalescing Binary Neutron Stars](http://adsabs.harvard.edu/abs/1992ApJ...401..226R)


### Tides and Collisions
[Micro-Tidal Disruption Events by Stellar Compact Objects and the Production of Ultra-Long GRBs](http://adsabs.harvard.edu/abs/2016ApJ...823..113P)

[Formation of Black Hole X-Ray Binaries in Globular Clusters](http://iopscience.iop.org/article/10.1088/0004-637X/717/2/948/pdf)

[Tidal Breakup of Binary Stars at the Galactic Center. II. Hydrodynamic Simulations](http://adsabs.harvard.edu/abs/2011ApJ...731..128A)

[On the Origin of the Metallicity Dependence in Dynamically formed Extragalactic Low-mass X-Ray Binaries](http://adsabs.harvard.edu/abs/2010ApJ...717..948I)

[On the Onset of Runaway Stellar Collisions in Dense Star Clusters - II. Hydrodynamics of Three-Body Interactions](http://adsabs.harvard.edu/abs/2010MNRAS.402..105G)

[Mixing in Massive Stellar Mergers](http://adsabs.harvard.edu/abs/2008MNRAS.383L...5G)

[Stellar Collisions and Ultracompact X-Ray Binary Formation](http://adsabs.harvard.edu/abs/2006ApJ...640..441L)

[Formation of Ultracompact X-Ray Binaries in Dense Star Clusters](http://adsabs.harvard.edu/abs/2005ApJ...621L.109I)

[Tidal Interactions and Disruptions of Giant Planets on Highly Eccentric Orbits](http://adsabs.harvard.edu/abs/2005Icar..175..248F)

[Modeling Collision Products of Triple-Star Mergers](http://adsabs.harvard.edu/abs/2003MNRAS.345..762L)

[Stellar Collisions and the Interior Structure of Blue Stragglers](http://adsabs.harvard.edu/abs/2002ApJ...568..939L)

[Evolution of Stellar Collision Products in Globular Clusters. II. Off-Axis Collisions](http://adsabs.harvard.edu/abs/2001ApJ...548..323S)

[Collisions of Main-Sequence Stars and the Formation of Blue Stragglers in Globular Clusters](http://adsabs.harvard.edu/abs/1997AAS...191.9801L)

[Evolution of Stellar Collision Products in Globular Clusters. I. Head-On Collisions](http://adsabs.harvard.edu/abs/1997ApJ...487..290S)

[Collisions of Main-Sequence Stars and the Formation of Blue Stragglers in Globular Clusters](http://adsabs.harvard.edu/abs/1996ApJ...468..797L)

[On Blue Straggler Formation by Direct Collisions of Main Sequence Stars](http://adsabs.harvard.edu/abs/1995ApJ...445L.117L)

[Collisions of Giant Stars with Compact Objects](http://adsabs.harvard.edu/abs/1991ApJ...377..559R)

### Other Dynamical Events
[Accretion Disks around Kicked Black Holes: Post-kick Dynamics](http://adsabs.harvard.edu/abs/2012ApJ...745...71P)

### Theses
[Smoothed Particle Hydrodynamic Simulations of Stellar Collisions](http://adsabs.harvard.edu/abs/1998PhDT........16L)

[Hydrodynamical Calculations of Stellar Interactions](http://adsabs.harvard.edu/abs/1991PhDT........11R)


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

