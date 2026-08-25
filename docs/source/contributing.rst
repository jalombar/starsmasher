Contributing
============

Issues and pull requests are welcome at `github.com/jalombar/starsmasher
<https://github.com/jalombar/starsmasher>`_.

Before opening a pull request, run the test suite::

   python3 tests/run_tests.py

Working on the Fortran
----------------------

The code is a mixture of fixed-form (``.f``) and free-form (``.f90``) source,
which has a few consequences worth knowing before they cost you an afternoon.

**Comment markers in fixed-form source belong in column 1.**  A ``c`` in column 1
starts a comment; a ``c`` after leading whitespace is the first letter of a
variable.  This bites tooling more often than people: a script that treats
anything matching ``\s*c`` as a comment will silently discard lines beginning
``cn1 =`` or ``computeexclusivemode =``.

**Do not put fixed-form comments in** ``starsmasher.h``.  The header is included
by ``advance.f90``, which is free form, so a ``c`` in column 1 of the header is a
syntax error there even though it would be a comment in every ``.f`` file that
includes it.  Use ``!``, which means the same thing in both.

**Restoring a file from a copy can defeat** ``make``.  ``cp -p`` preserves the
original modification time, so the object file still looks newer than the source
and is not rebuilt: the old code keeps running while you debug the new one.
``touch`` the file after any such restore.  The makefiles do list
``starsmasher.h`` as a dependency of every object, so editing the header does
trigger the rebuilds you would expect.

Changing the input namelist
---------------------------

:doc:`reference/sph_input` is generated from ``init.f``.  After adding a setting
or changing a default, regenerate it::

   python3 docs/generate_reference.py

and commit the result.  The documentation workflow fails if you do not, so this
is caught before review rather than after.

Give the new variable a trailing comment in the default-initialisation block of
``init.f``.  That comment becomes its description on the reference page, so it is
worth writing for a reader rather than for yourself.

Building the documentation
--------------------------

::

   python3 -m venv env
   env/bin/pip install -r docs/requirements.txt
   env/bin/sphinx-build -b html docs/source docs/build

The workflow builds with ``-W``, so a Sphinx warning fails the build.  If you
quote code containing ``|`` or ``*`` in a table cell, escape it.

Writing documentation
---------------------

Anything version- or machine-specific -- paths, terminal output, file names,
timings, particle counts -- should be produced by running it and pasted from real
output, not carried over from a previous version or written from memory.  Several
pages here quote figures from specific runs for exactly this reason.

Where a number depends on the model, say so and give the range rather than
implying a constant.
