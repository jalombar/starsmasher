# How to collide a compact object with a star

Before trying these simulations, you may want to practice with StarSmasher following the other tutorials. If you want to collide a compact object (a black hole, neutron star, or white dwarf modelled by a point mass) with a star, then look here!

First of all, you need a different code from the one present in the folder parallel_bleeding_edge. This modified code in src's folder is optimized for this kind of simulation and is much faster than the classical one when we are talking about a BH/NS collision with a MS star.

But now let's move on to the simulations. To collide a black hole or neutron star with a star, you first need to relax your star (transforming it from a 1D text file to a 3D SPH model). To do that, please, follow the tutorial "Creating a Mesa's star", BUT using this version of the code with the required files "sph.init" & "sph.input" located in the folders "relaxation_files" & "collision_files". Once that the relaxation ends, a file called "sph.passivelyAdvected" will be generated. This file keeps track of information needed for updating particle smoothing lengths. You'll need to copy that file to the directory where you run the collision; otherwise, the collision won't start.

Rename the relaxed star (for example, the last out*.sph file) "sph.start1u".  If there is no sph.start2u file, then a compact object of mass mbh (which can be set in sph.input) will be used as the second object.  That is, when Starsmahser notices that in your folder you don't have two stars (sph.start1u & sph.start2u), it will create a black hole / neutron star, a point mass, in the simulation. To decide the mass of your black hole / neutron star, you must write in sph.input " mbh=10 ", or whatever the mass you want (in solar masses). 

When you have relaxed your star, as explained in the tutorial "Creating a Mesa's star", you can start your simulation, then type in your terminal (that must be opened inside your folder):

```
mpirun -np N test_gpu_sph
```

Where N is the number of processes you want to use, usually less than or equal to the number of CPU cores on your machine.

When your star touches the black hole, the simulation will instantly become slower, as simulating the close interaction between a star and a pointmass black hole / neutron star is computationally expensive.  You'll also want to use the GPU version of the code to simulate this kind of collision.

# Visualizing data

In StarSmasher, your compact object is treated like a point mass with gravitational softening. In SPLASH, the software used to visualize StarSmasher's snapshots, you cannot see your point mass BH/NS in the default mode. Let's change this. First, let's visualize all your snapshots as we explained in the "visualizing_data's tutorial":

```
splash -f starsmasher out*.sph
```

Now, before visualizing the data in the usual way, type:

```
o
1
yes
yes
yes
```

Now check that the BH/NS is visible by, for example, entering the following.

```
2
1
8
0
/xw
```

The point mass is now shown.
