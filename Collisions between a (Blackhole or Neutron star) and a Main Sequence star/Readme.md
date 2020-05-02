To make collide a black hole / neutron star to a Main sequence star or "viceversa", you first need to relax your star (transforming it from a 1D text file to a 3D SPH Star). To do that, please, follow the tutorial "Creating a Mesa's star".

Once you have your star, now you have to decide if it gets hit by the black hole/the neutron star (1) or if it will collide with the black hole/ neutron star.

If you want (1), you have to rename the relaxed star "sph.start1u"
If you want (2), you have to rename the relaxed star "sph.start2u"

When starsmahser notice that in your folder you don't have two stars (sph.start1u & sph.start2u), it will create a black hole / neutron star, a point mass, in the simulation. To decide the mass of your black hole / neutron star, you must writhe in sph.input  " mbh=10 ", or whatever the mass you want. In this tutorial there is an example with 10. The unit of N is in solar mass. 
In this simulatio sph.imput is setted with a bimpact (the periastron separation) of 0. It means hat the black hole will it the center of the star or viceversa.

When you have your star, as explained in the tutorial "Creating a Mesa's star", you can start your simulation, then type in your terminal (that must be opened inside your folder):

mpirun -np N test_gpu_sph


Where N are the numbers of cores (for section) of your machine.

You will likely need a supercomputer to simulate this kind of collisions, even if your star has 10k particles and you builded the gpu version of the program. That's because simulating the interaction between a star and a pointmass black hole / neutron star will require  a very high computing capability. If you are simulating in a laptop, you will se that when your star will touch the black hole, the simulation will instantly become slower (unless you have a supercomputer).
