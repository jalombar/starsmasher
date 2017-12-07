# DRAFT

## Introduction

Hello there, welcome to our public release of *StarSmasher*- a Smoothed Particle Hydrodynamics (SPH) based, parallel code. *StarSmasher* has been developed and maintained by James Lombardi with help from many others. *StarSmasher* is a code that evolved from a code created and maintained by Joshua Faber, James Lombardi, and Fred Rasio called *StarCrash*, and much of this guide is directly from the guide for *StarCrash*. My name is Sarah Seitanakis, and I will help you get to know our code-- so that you, too, can become a StarSmasher. 

The following sections of this guide will include installation instructions, a discription of portions of the code (such as some background on SPH, parallelization techniques, and physical effects), and we will wrap it up with an overview of how to properly use the code- and even add your own routines. 
  
While this code is being freely distributed, we ask that you please give appropriate credit to those who have spent a great deal of time developing it. If you are going to publish work done with this code, please see Sec. 8 for the proper references to cite, as well as the license agreement we are using, the rather ubiquitous Gnu General Public License (GPL).



## Installation Instructions

The installation instructions are layed out in [this file](./installation.md). I believe in you. :thumbsup:



## SPH: Smoothed Particle Hydrodynamics

#### History & The Basics

The basic Smoothed Particle Hydrodynamics (SPH) method was created by Lucy (1977) and
Gingold & Monaghan (1977) in order to study fission in rotating stars. It has since been used to
study, among other astrophysical topics, large scale structure in the universe, galaxy formation,
star formation, supernovae, solar system formation, tidal disruption of stars by massive black holes,
and stellar collisions; see Rasio & Shapiro (1992), Monaghan (1992), and Rasio & Lombardi (1999)
for a more complete list of references. Our particular code has been used primarily in the study of
stellar collisions and mergers, including merging compact object binaries (Rasio & Shapiro 1992,
1994, 1995), collisions involving main sequence stars (Rasio & Shapiro 1991; Lai et al. 1993a; Sills
et al. 1997, 2001), blue-straggler formation (Lombardi et al. 1995, 1996, 2002), and most recently
post-Newtonian (PN) and relativistic studies of binary neutron star (NS) systems (Faber & Rasio
2000; Faber et al. 2001; Faber & Rasio 2002). 
Because of its Lagrangian nature, SPH presents some clear advantages over more traditional
grid-based methods for calculations of stellar interactions. Most importantly, fluid advection, even
for stars with a sharply defined surface such as NS, is accomplished without difficulty in SPH,
since the particles simply follow their trajectories in the flow. In contrast, to track accurately the
orbital motion of two stars across a large 3D grid can be quite tricky, and the stellar surfaces
then require a special treatment (to avoid “bleeding”). SPH is also very computationally efficient,
since it concentrates the numerical elements (particles) where the fluid is at all times, not wasting
any resources on emty regions of space. For this reason, with given computational resources, SPH
provides higher averaged spatial resolution than grid-based calculations, although Godunov-type
schemes such as PPM typically provide better resolution of shock fronts (this is certainly not a
decisive advantage for binary coalescence calculations, where no strong shocks ever develop). SPH
also makes it easy to track the hydrodynamic ejection of matter to large distances from the central
dense regions. Sophisticated nested-grid algorithms are necessary to accomplish the same with
grid-based methods.

#### SPH Density: Kernels, Sums, and all that Jazz

A good description of the particular SPH algorithm implemented in our code can be found in
Rasio & Shapiro (1992). We are going to summarize the most important features for you here.


[\alpha]

