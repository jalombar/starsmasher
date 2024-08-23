Compile 
========

Compile the source and place the executable in the `bin/` directory


```$ cd src/SPHgrav_lib2
   $ make
   $ cd ..
   $ make
   $ mv ./*_sph ../bin/
```

Create 3D models from MESA 1D models
=====================================

Get into MESA/ and use your 1D profile file created with MESA
as input for the sph.input (last rows). This will create an SPH
profile from the 1D MESA file, say profile9.data. It needs to relax.

```
   $ cp MESA_profiles/LOGS/profile9.data ./MESA/
   $ sed -i.bak "s|^profilefile=.*|profilefile='profile9.data'|" sph.input
   $ ./run.sh
```

Then wait for it to finish.

Take the last snapshot as your input 3D file for the
stellar collision. Rename it

```
   $ mv out0300.sph sph.start1u
```    

Repeat for the second profile, corresponding to the second star.

```
   $ mv out0300.sph sph.start2u
```

Brace for the collision
========================

Create a collision directory,

```
   $ mkdir Collision
```

Copy all necessary collision files into that directory,

```
   $ cp coll_files/* Collision
```


Move the 3D MESA SPH initial profiles into the Collision directory

```
   $ mv MESA/sph.start* Collision/
```


Get into the directory and start the simulation,

```
   $ cd Collision/
   $ ./run.sh
```
