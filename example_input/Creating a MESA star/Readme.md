# How to create a star in Starsmasher and use it to make a collision


To create a Star in Starsmasher and use it, at the moment the only way to do that is using MESA stellar evolution code.
Mesa is a code that is able to create stars having mainly the mass and the metallicty of the Stars.

To get Mesa go to the following link              http://mesa.sourceforge.net/

and be sure to follow all the tutorials to use it!

At the end of the Mesa's simulation we have to look at the folder called LOGS and then at the files called "ProfileN.data". There will be many, each number "1, 2, 3" and so on" is a "snapshot" of your star. You have to look obviously at the profile you are interested in (that differ in radius, age, helium concentration and so on!).

Once you have your file, you have to edit your sph.input file. Be sure that, if your profile file is called, for example, "start1". Inside the sph.input file there is written 

profilefile='start1'

You can change the name inside '*' as you want. Be just sure that the name of your file and the name written inside sph.input are the same!

Be also sure that there is written  "stellarevolutioncodetype=1". In this way Starsmasher will recognize that this is a MESA file.

Be sure that the nrelax=1. In this way Starsmasher will convert your mesa stars in SPH particles. The number of SPH particles is setted in the 4th line of sph.input file. James Lombardi suggest a number of 50000 particles for solar mass. In this way the simulation will be realistic.

To make the star, after that in your folder there are all the files needed (like in this example folder), type in your terminal:

```
mpirun -np N test_gpu_sph
```
Where N are the numbers of cores (cpu) of your machine.

Obviously be sure that everything (number of cpus, graphic cards and so on) is setted for your machine, or mpirun will go in error.

At the end of the simulation, you will have many files called out0000.sph, out0001.sph and so on. The last file should be a "out0300.sph"; it obviously depend when you decide to finish your simulation (specified in tf=N inside sph.input).

To use your star in a simulation, create a new folder copying the "collision" folder that is in "example_imput" folder. In that folder you have to delete the files called "sph.start1u and sph.start2u" and then paste the last file of your relaxation run. That file contains you star.

If you want that your star is the target of the collision, rename it "sph.start1u". If you want that is the impactor, call it "sph.start2u.

At this point set the collision parameters in the sph.input for the that collision that you want.

then start the simulation:

```
mpirun -np N test_gpu_sph
```

Where N is the number of cores (cpu) of your machine.


The simulation will create many out*.sph files. you can see them following the tutorial inside the "visualizing data" folder of the Starsmasher repository!

In this folder there are also 2 already done MESA stars!    Start1 is a SUN (1 solar mass and 1 solar radius), Start2 is a red dwarf (0.5 solar mass and 0.5 solar radius). You can use them to practice with Starsmasher! 

There is also another file in this folder called "tips and tricks". In this one there is a list of tricks that you can add inside your sph.imput before relaxing the star, for example, to make it rotate!

For any problems please contact Francesco Radica at "francyrad97@gmail.com".

