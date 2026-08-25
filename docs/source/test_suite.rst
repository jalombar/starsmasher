Test suite
==========

::

   python3 tests/run_tests.py

The suite runs short calculations and checks them against things that must be
true, rather than against stored output.  A golden-file test tells you that
something changed; these tell you whether the answer is right.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Test
     - Checks
   * - ``virial``
     - A relaxed :math:`n=1.5` polytrope satisfies
       :math:`W = -\frac{3}{5-n}\frac{GM^2}{R}`, and :math:`U = -W/2`.
   * - ``energy``
     - Total energy is conserved through the run.
   * - ``quiet``
     - A relaxed star stays put: kinetic energy remains negligible.
   * - ``ranks``
     - The answer does not depend on how many MPI ranks are used.
   * - ``cpu_gpu``
     - The CPU and GPU builds agree.
   * - ``guard_ngravprocs``
     - Asking for more gravity processes than there are GPUs is clamped rather
       than fatal.

Options::

   python3 tests/run_tests.py --list            # what exists
   python3 tests/run_tests.py virial energy     # named tests only
   python3 tests/run_tests.py --np 2            # fewer ranks
   python3 tests/run_tests.py --exe /path/to/..._sph   # skip the build
   python3 tests/run_tests.py --with-cpu        # also build the CPU version

Each build happens in a copy of the source under the suite's own scratch
directory, so a test run leaves the repository untouched.

Why these tests
---------------

The last three exist because each corresponds to a real defect.

Results once depended on the number of MPI ranks: the CPU gravity kernel
accumulates into particles a rank does not own, and several code paths combined
that with a gather over each rank's own slice, which silently discarded the rest.
The gravitational potential of a stellar model came out 30 per cent wrong on four
ranks and correct on one.

The CPU and GPU builds once disagreed for the same reason, since only the CPU
path was affected.

Asking for more gravity processes than there were GPUs produced a segmentation
fault inside the CUDA library, with a stack trace that gave no indication the
cause was a setting in ``sph.input``.

None of these would have been caught by comparing against saved output, because
there was no known-good output to compare against.  All three are caught by
asking whether the answer is physically right, and whether it is the same answer
under configurations that ought not to matter.

Adding a test
-------------

Write a function taking a context and returning a one-line description of what
it verified, raise ``Fail`` with the numbers when it does not hold, and add it to
the ``TESTS`` dictionary in ``tests/run_tests.py``.  Prefer an invariant that has
to be true over a value that happened to be produced.
