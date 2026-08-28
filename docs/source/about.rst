About StarSmasher
=================

StarSmasher is a smoothed particle hydrodynamics code for collisions and close
interactions between stars and planets.  It descends from StarCrash, a parallel SPH code written by
Joshua Faber, James Lombardi, and Frederic Rasio, which descends the SPH code written for Rasio's PhD thesis. StarSmasher is currently maintained by, primarily, James
Lombardi at Allegheny College.

Papers
------

The algorithms, and in particular the equations of motion described in
:doc:`reference/equations_of_motion`, are set out in

* Gaburov, E., Lombardi, J. C., Jr., & Portegies Zwart, S. 2010,
  *On the onset of runaway stellar collisions in dense star clusters -- II.
  Hydrodynamics of three-body interactions*, MNRAS, 402, 105
  (`ADS <https://ui.adsabs.harvard.edu/abs/2010MNRAS.402..105G/abstract>`_,
  `arXiv:0904.0997 <https://arxiv.org/abs/0904.0997>`_)

The artificial viscosity follows Balsara (1995) and Ponce et al. (2011).

A fuller list, including work that has used the code, is kept in the repository
in `papers.md
<https://github.com/jalombar/starsmasher/blob/master/documentation/papers.md>`_
and `publications.md
<https://github.com/jalombar/starsmasher/blob/master/documentation/publications.md>`_.

Citing
------

If you publish work done with StarSmasher, please cite the relevant papers above.

Licence
-------

`GNU General Public License v3
<https://github.com/jalombar/starsmasher/blob/master/LICENSE>`_.

Contact
-------

Questions and feedback are welcome as an issue at `github.com/jalombar/starsmasher
<https://github.com/jalombar/starsmasher>`_, which is the best place for
anything others might benefit from seeing.

.. raw:: html

   <p>By email:
     <span id="ss-contact">jamie&#46;lombardi <em>at</em> allegheny <em>dot</em> edu</span>
   </p>
   <script>
   // Assembled at display time so the address is never present as a single
   // string in the served HTML.  Without scripting the reader still sees a
   // human-readable form above.
   (function () {
     var el = document.getElementById('ss-contact');
     if (!el) { return; }
     var at = String.fromCharCode(64), dot = String.fromCharCode(46);
     var who = ['jamie', 'lombardi'].join(dot);
     var where = ['allegheny', 'edu'].join(dot);
     var a = document.createElement('a');
     a.href = 'mai' + 'lto:' + who + at + where;
     a.textContent = who + at + where;
     el.parentNode.replaceChild(a, el);
   })();
   </script>
