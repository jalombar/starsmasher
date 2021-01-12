# How to make collide a Black Hole (BH) / Neutron star (NS) with a Main Sequence Star

Before making these simulations, we suggest you to make practice with StarSmasher following the other tutorials. If you want to make collide a Black Hole/Neutron Star with a Main Sequence star or viceversa, you are welcome!

First of all, you need a different code from the one present in the folder parallel_bleeding_edge. This modified code in src's folder is optimized for this kinds of simulation and is much faster than the classical one when we are talking about BH/NS collision with a MS star. The classical code will go very slow when the encounter happen, this will go faster instead.

It's also been reported that compiling this code is not easier. This one is been tested in Keeneland and Ubuntu 20.04. We are sure that the compilation will fail in Ubuntu 19 or lower version.

But now lets pass to the simulations. To make collide a black hole / neutron star to a Main sequence star or "viceversa", you first need to relax your star (transforming it from a 1D text file to a 3D SPH Star). To do that, please, follow the tutorial "Creating a Mesa's star", BUT using this version of the code with the files "sph-init" & "sph-input" related in the folders "relaxation_files" & "collision_files". If you won't use these files, the simulation won't start. Once that the relaxation ended, it will be generated a file called "sph.passivelyAdvected". This file keeps track of information needed for updating particle smoothing lengths. You'll need to copy that file to the directory where you run the collision. If you won't do, the collision won't start.

Once you have your star, now you have to decide if it gets hit by the black hole/the neutron star (1) or if it will collide with the black hole/ neutron star (2).

If you want (1), you have to rename the relaxed star "sph.start1u" If you want (2), you have to rename the relaxed star "sph.start2u"

When starsmahser notice that in your folder you don't have two stars (sph.start1u & sph.start2u), it will create a black hole / neutron star, a point mass, in the simulation. To decide the mass of your black hole / neutron star, you must write in sph.input " mbh=10 ", or whatever the mass you want (in solar masses). 

When you have relaxed your star, as explained in the tutorial "Creating a Mesa's star", you can start your simulation, then type in your terminal (that must be opened inside your folder):

```
mpirun -np N test_gpu_sph
```

Where N are the numbers of cores (for section) of your machine.

You will likely need a supercomputer to simulate this kind of collisions, even if your star has 5k particles and you builded the gpu version of the program. That's because simulating the interaction between a star and a pointmass black hole / neutron star will require a very high computing capability. If you are simulating in a laptop, you will see that, when your star will touch the black hole, the simulation will instantly become slower (unless you have a supercomputer or a GPU cloud).



# Visualizing data

In StarSmasher, your black hole or Neutron star is treated like a point mass. In SPLASH, the software used to visualize StarSmasher's snapshots, as default, you cannot see your point mass BH/NS. To do that, first, let's visualize all your snapshots as we explained in the "visualizing_data's tutorial":

```
jsplash out*.sph
```

Now, before visualizing the grid as in the classical way, type:

```
o
yes
yes
yes
```

Now the BH/NS is visible. check with:

```
2
1
8
0
/xw
```

The point mass is now evidenced with a white circumference.
