HDF5 snapshots
==============

The snapshots are ``out*.sph`` files, a binary format of the code's own, which
:doc:`SPLASH <visualisation>` reads.  The code can write them in HDF5 as well.
An HDF5 file is read from Python in a couple of lines, with ``h5py`` or with the
reader that comes with the code, and a column can be read without reading the
rest: there is no format to decode by hand and no SPLASH needed to look at a
run.

The ``out*.sph`` files stay the default, and ``restartrad.sph`` is always
written in the native format.

Switching it on
---------------

Two things, one when the code is built and one when it runs.

**At build time**, ``HDF5=1``::

   $ make HDF5=1            # the GPU version
   $ make cpu HDF5=1        # the CPU version

This needs the HDF5 library with its headers (``libhdf5-dev`` on Debian and
Ubuntu) and ``mpicc``.  ``make config HDF5=1`` shows which HDF5 the build found.
If it found the wrong one, or none, say where it is with ``HDF5ROOT``, the
directory that holds ``include`` and ``lib``, or give ``HDF5_CFLAGS`` and
``HDF5_LIBS`` yourself::

   $ make HDF5=1 HDF5ROOT=/opt/hdf5
   $ make HDF5=1 HDF5_CFLAGS='-I/opt/hdf5/include' HDF5_LIBS='-L/opt/hdf5/lib -lhdf5'

The writer uses only the C interface of HDF5 and only its serial form, so any
HDF5 from 1.10 on will do, whether it was built with MPI or not.  Without
``HDF5=1`` the code is built exactly as before, with nothing new to link.

**At run time**, ``nhdf5`` in ``sph.input``:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - ``nhdf5``
     - Effect
   * - 0
     - Only ``out*.sph``, as always.  This is the default.
   * - 1
     - ``out*.sph`` and ``out*.h5``, with the same number and at the same moment.
   * - 2
     - Only ``out*.h5``.

A run that asks for ``nhdf5`` above 0 with a build made without ``HDF5=1`` stops
at the start with ``nhdf5>0 needs a build with make HDF5=1``, rather than going
on and writing nothing.

The name is ``nhdf5`` and not ``enable_hdf5`` because the integer settings of
the code are all called ``n...`` (``ncooling``, ``nkernel``, ``nintvar``), and
because this one has three values.

.. warning::

   With ``nhdf5=2`` there are no ``out*.sph`` files.  Those are also what a
   collision starts from (``startfile1`` and the others read that format), what
   SPLASH reads, and what the tools in ``splot_routines`` read.  Use ``nhdf5=2``
   to analyse a run, and ``nhdf5=1`` for a run whose snapshots you may want to
   start another calculation from.

What is written
---------------

``out0001.h5`` goes with ``out0001.sph``: the same number, with as many digits
as the native file has, written at the same point of the run and from the same
arrays, after the same MPI gathers, so that the two hold the same numbers
exactly.  Rank 0 writes it.  The file is written under the name
``out0001.h5.tmp`` and renamed when it is complete, so that a run that is
killed in the middle of a write never leaves a half-written ``out0001.h5``.

``log*.sph`` gets a line for every snapshot written, in either format, with the
wall time of the write::

   duout: wall time writing out0001.sph (native):   <seconds> s
   duout: wall time writing out0001.h5 (hdf5):   <seconds> s

The file has two groups.  This is version 1 of the layout.

``/header``
~~~~~~~~~~~

Attributes.  The first 25 are the entries of the header of ``out*.sph``, with the
same names and the same types: ``int32`` for the integers and ``float64`` for the
``real*8``:

``ntot nnopt hco hfloor sep0 tf dtout nout nit t nav alpha beta tjumpahead ngr
nrelax trelax dt omega2 ncooling erad lastnumstars epoteat ekineat einteat``

Then come these, which the native header does not have:

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Attribute
     - Type
     - Meaning
   * - ``format``
     - string
     - ``StarSmasher HDF5``
   * - ``format_version``
     - int32
     - 1
   * - ``munit``, ``runit``, ``gravconst``
     - float64
     - The three numbers that give the code units their size in cgs: grams in the
       unit of mass, centimetres in the unit of length, and :math:`G` in cgs.  See
       :ref:`code units <code-units>`.
   * - ``tunit``
     - float64
     - The unit of time in seconds, :math:`\sqrt{R_\mathrm{unit}^3/(G M_\mathrm{unit})}`.
   * - ``nintvar``
     - int32
     - What ``u`` holds: the specific internal energy (2), the entropic variable
       :math:`A` (1) or :math:`\ln A` (3).
   * - ``neos``, ``nkernel``, ``gflag``, ``nusegpus``, ``nselfgravity``, ``nblock``, ``dynhco``
     - int32
     - The settings of the same name in ``sph.input``, as the run has them.
   * - ``gam``, ``reat``, ``boxx``, ``boxy``, ``boxz``
     - float64
     - Likewise.
   * - ``hdf5_library_version``
     - string
     - The version of HDF5 that wrote the file.

``/particles``
~~~~~~~~~~~~~~

One dataset of length ``ntot`` for every column of ``out*.sph``, named like the
Fortran variable, in the order of the native file.  Each has two attributes,
``description`` and ``units``; the units are written in terms of ``munit``,
``runit`` and ``tunit`` from ``/header``.

.. list-table::
   :header-rows: 1
   :widths: 16 12 40 32

   * - Dataset
     - Type
     - Meaning
     - Units
   * - ``x``, ``y``, ``z``
     - float64
     - position
     - ``runit``
   * - ``am``
     - float64
     - mass
     - ``munit``
   * - ``hp``
     - float64
     - smoothing length :math:`h` (the kernel reaches :math:`2h`)
     - ``runit``
   * - ``rho``
     - float64
     - density
     - ``munit/runit^3``
   * - ``vx``, ``vy``, ``vz``
     - float64
     - velocity
     - ``runit/tunit``
   * - ``vxdot``, ``vydot``, ``vzdot``
     - float64
     - acceleration
     - ``runit/tunit^2``
   * - ``u``
     - float64
     - internal energy variable: see ``nintvar`` above
     - ``(runit/tunit)^2`` for ``nintvar=2``
   * - ``udot``
     - float64
     - time derivative of ``u``
     - ``(runit/tunit)^2/tunit`` for ``nintvar=2``
   * - ``grpot``
     - float64
     - gravitational potential at the particle
     - ``(runit/tunit)^2``
   * - ``meanmolecular``
     - float64
     - mean molecular weight times the atomic mass unit
     - g
   * - ``cc``
     - int32
     - integer label of the component (star) the setup routines gave the particle
     -
   * - ``divv``
     - float64
     - divergence of the velocity
     - ``1/tunit``
   * - ``ueq``, ``tthermal``
     - float64
     - equilibrium specific internal energy and thermal timescale of the
       cooling.  Only when ``ncooling`` is not 0, as in ``out*.sph``
     - ``(runit/tunit)^2``, ``tunit``

A snapshot takes 140 bytes per particle (148 in ``out*.sph``, which frames every
record with two four-byte markers), 156 with ``ncooling``.

Nothing is compressed, on purpose: see the measurements below.

How long it takes
~~~~~~~~~~~~~~~~~

Measured on one compute node of the UPSCALE cluster (gfortran 14.2, Open MPI
5.0.7, HDF5 1.14.5, four MPI ranks), writing to its network file system (NFS).
Every value is the median of three runs.  The time to write is the one the code
logs for each snapshot, taken as the median over the snapshots of a run.  The
times to read are the median of five reads in Python, after one more to warm
up, with the file in the cache of the node, so they leave the disk out.  Sizes
are in MB of :math:`10^6` bytes.

.. list-table::
   :header-rows: 1
   :widths: 34 14 14 19 19

   * - Particles, format
     - Write [ms]
     - Size [MB]
     - Read all [ms]
     - Read x, y, z, rho [ms]
   * - 99 954, ``out*.sph``
     - 61
     - 14.8
     - 7.4
     - 3.1
   * - 99 954, ``out*.h5``
     - 25
     - 14.0
     - 5.0
     - 2.7
   * - 99 954, ``out*.h5``, deflate 1
     - 271
     - 10.6
     - 42
     - 10.5
   * - 99 954, ``out*.h5``, deflate 4
     - 299
     - 10.5
     - 44
     - 10.9
   * - 499 860, ``out*.sph``
     - 423
     - 74.0
     - 51
     - 21.8
   * - 499 860, ``out*.h5``
     - 238
     - 70.0
     - 21.5
     - 6.4
   * - 499 860, ``out*.h5``, deflate 1
     - 948
     - 34.7
     - 165
     - 49
   * - 499 860, ``out*.h5``, deflate 4
     - 1099
     - 34.0
     - 175
     - 46
   * - 1 999 440, ``out*.sph``
     - 1076
     - 295.9
     - 278
     - 120
   * - 1 999 440, ``out*.h5``
     - 906
     - 279.9
     - 78
     - 12.9
   * - 1 999 440, ``out*.h5``, deflate 1
     - 5015
     - 200.0
     - 615
     - 101
   * - 1 999 440, ``out*.h5``, deflate 4
     - 5603
     - 198.4
     - 613
     - 100

What the numbers do and do not say:

* The 99 954 particle star was run with gravity.  The two larger ones were run
  *without* gravity (``ngr=0``) and as pure hydrodynamics, because gravity by
  direct summation on CPUs is too slow for two million particles.  In those
  snapshots ``grpot`` and ``udot`` are all zero, which compress far better than
  real data, so their deflate rows flatter the compression.  Only the first
  star says how well real data compress.
* Writing an ``out*.h5`` takes about as long as writing the ``out*.sph``, or
  less, on that file system.  Times to write on NFS vary from run to run (for the
  two million particle star, 1.07 to 2.10 s for ``out*.sph`` and 0.36 to 1.02 s
  for ``out*.h5``, over nine snapshots), so read it as "not slower", not as a
  factor.
* The clear difference is in reading: the whole file 1.5 to 3.6 times faster,
  and a few columns up to nine times faster (13 ms against 120 ms for
  :math:`x, y, z, \rho` of two million particles), since only an HDF5 file lets
  you read a column without reading the rest.  The reader of ``out*.sph`` has to
  read the whole file in any case.
* Compression was measured and taken out of the code.  Deflate, with the shuffle
  filter, at level 1 or 4, brought the file of real data down to about 75% of its
  size, and made writing it about eleven times slower and reading it eight to
  nine times slower.  A quarter of the space is not worth that, so the datasets
  are not compressed and there is no option for it.
* One machine, one file system, a warm cache, one run length.  The ``log*.sph``
  of any run has the time of each write, so the writing can be checked on
  another file system.

Reading it
----------

From Python, with the reader in ``parallel_bleeding_edge/bin`` (it needs numpy,
and h5py for HDF5 files; put the directory on ``PYTHONPATH`` or copy the file next
to your script)::

   import StSm_snapshot as ss

   header, p = ss.read('out0003.h5')            # a dict of numpy arrays
   print(header['t'], p['am'].sum())
   header, p = ss.read('out0003.h5', columns=['x', 'y', 'z', 'rho'])  # only these

``ss.read`` takes ``out*.sph`` as well, in both layouts, and returns the same
structure; for a native file it checks the length marker of every record and the
``ntot`` that ends the file.  The other functions::

   cu = ss.code_units(header)            # what a code unit is worth in cgs
   rho_cgs = p['rho'] * cu['density']    # g/cm^3 (the units are in an HDF5 header;
                                         # for out*.sph pass munit= and runit=)
   df = ss.to_dataframe(p)               # a pandas DataFrame
   ss.compare('out0003.sph', 'out0003.h5')   # bit for bit: (18, [], []) if equal

``to_dataframe`` gives the mass as ``m`` and the smoothing length as ``h``, the
names that `sarracen <https://sarracen.readthedocs.io>`_, which renders SPH data
in Python, takes.  (The kernel of the code is a Wendland one by default, and may
not be one that your version of sarracen has: check that before trusting a
rendering to the last digit.)  ``compare`` returns the number of columns that
both files hold, the names of the columns that differ in any bit, and the
entries of the header that differ.

``StSm_snapshot_plots.py`` is an example to copy::

   $ python3 StSm_snapshot_plots.py out0003.h5 plots/star

It draws the density of every particle against its distance from the centre of
mass, and the column density on the x-y plane, as two PNG files.

Without the reader, ``h5py`` does it directly::

   import h5py
   with h5py.File('out0003.h5', 'r') as f:
       print(f['header'].attrs['t'], f['particles/rho'][:].max())
       print(f['particles/rho'].attrs['units'])

If the HDF5 tools are installed, ``h5dump -H out0003.h5`` shows the structure,
and ``h5ls -r`` lists it.

.. note::

   HDF5 locks the files it opens.  Some network file systems do not support
   that, and reading then fails with a message about locking.  The writer of the
   code switches the locking off for its own file; to read the file in such a
   place, set ``HDF5_USE_FILE_LOCKING=FALSE`` in the environment first.
