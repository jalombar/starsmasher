# Directory `tools2/` 

This directory contains more complex tools than the main `tools/` one,
where we have configuration files and light-weight tools. In `tools2/`
we will add codes which go beyond the small scripts of `bin/`.

# Charles Gibson's sph-to-mesa

These codes developed in Python analyse output data from StarSmasher
collisions. The primary focus of these codes is to construct 1-D stellar models
based on the 3-D data produced by StarSmasher. The generated output files can
be utilized to evolve the results of Smoothed Particle Hydrodynamics (SPH)
collisions within the MESA (Modules for Experiments in Stellar Astrophysics)
framework.

The code `sph-to-mesa` provided in the `tools2/` directory is a reproduction of
the program originally developed by Charles Gibson.  I have adapted it to make
it work with the rolling version of `StarSmasher`.

This is his original repository of the code (not compatible with the rolling version):

<a href="https://github.com/charlie-gibson/sph-to-mesa">https://github.com/charlie-gibson/sph-to-mesa</a>

## A step by step guide

I am using `sph-to-mesa` to produce a file which I can further evolve with `MESA`. In order to do that,
it is a good idea to follow these steps.

### Adapt your scripts

First let `splot.py` where your `sph-to-mesa` files are. Edit `splot.py` and add the path
to line 25:

```
sys.path.insert(0, f'/path/to/sph-to-mesa/python_splot/') # location of the main folder
```

Then move the script `splot.py` to any folder which is in your path. For
instance, in `$HOME/bin`.

### Create `sph.composition`

In the folder in which you have created the 3D models for the collision simulation,
run `splot.py` for the last file and choose option 2. 

For instance,

```
splot.py 
Which option would you like to do?
    Write ascii output                         [0] 
    Generate MESA files - no collision         [1] 
    Generate composition.sph                   [2] 
    Generate MESA files - collision            [3] 
    Determine components                       [4] 
    Generate MESA files - existing bestfit     [5] 
    pa plot                                    [6] 
    v plot                                     [7] 
    Force Hydrostatic Equilibrium              [8] 
    Calculate Eccentricity                     [9] 
: 2
Starting out file (int): 10
Ending out file (int, must be >= the start value): 10
Step size between files: 1
```

I have chosen file 10 because it is the last one I have, and I ask splot to
create the composition file just for this one last file because the components
will not have varied.

Once you run it, you will get some screen output plus a matplotlib plot,
giving you the fraction of elements as a function of the radius.

You will get a similar output to this

```
------------------- sph.input ------------------
Ignoring invalid line 0:  &input

Ignoring invalid line 11: ! bimpact=1.d30, ! the impact parameter for collisions

Ignoring invalid line 25: ! the courant numbers cn1, cn2, cn3, and cn4 are for sph particles:

Ignoring invalid line 25: ! dt_sph=1/(1/dt1 + 1/dt2 + 1/dt3 + 1/dt4)

Ignoring invalid line 29: ! the courant numbers cn5, cn6, and cn7 are for a particle i that is a compact object (co):

Ignoring invalid line 29: ! dt_co=1/(1/dt5 + 1/dt6)

Ignoring invalid line 32: ! the final timestep dt is the minimum of dt_sph and dt_co for all particles i

Ignoring invalid line 36:  &end

tf              100.0
dtout           10.0
n               20000
nnopt           23
nav             3
alpha           1.0
beta            2.0
ngr             3
nrelax          1
trelax          1.0
sep0            20000.0
vinf2           1.d30
equalmass       0.0
treloff         50.0
tresplintmuoff  0.0
nitpot          1
tscanon         0.0
sepfinal        1.d30
nintvar         2
ngravprocs      -2
qthreads        0.0
gflag           1.0
mbh             5.0
runit           6.9599d10
munit           1.9891d33
cn1             .6d0
cn2             0.06d0
cn3             0.06d0
cn4             1.d30
cn5             0.0005d0
cn6             0.01d0
cn7             4.d0
computeexclusivemode 0.0
ppn             16.0
profilefile     profile8.data
stellarevolutioncodetype 0.0



--------------------------Output file 10--------------------------
Looking for out0010.sph
DATA READ IN
VELOCITIES UPDATED
Output data analysis completed.
```

And you will have got your `sph.composition` for your first star.

Now you should repeat the same steps for the second star. Since the name of the
file is going to be the same one, it is a good idea to have two separated
folders where you create and relax the 3D model for each individual star, so
that the outputs are not overwritten.

Finally, you concatenate the two `sph.composition` files into a single one,

```
cat 1st_star/sph.composition 2nd_star/sph.composition > sph.composition
```



## Note: Features and limitations

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


