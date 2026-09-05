This builds a detached binary in a corotating frame and relaxes it, so that both
stars settle into the tidally distorted, synchronously rotating shapes they would
have at this separation.  Whatever rotation the stars were given during their own
relaxations does not matter, because the corotating frame sets it.

It needs two relaxed stars, under the names `startfile1` and `startfile2` give.
The two that ship with the collision example will do:

    cp ../collision/sph.start1u .
    cp ../collision/sph.start2u .

Those are a 8.0 Msun star of radius 3.17 and a 0.3 Msun star of radius 0.69, in
code units.  At `sep0=8` the primary fills about 60 per cent of its Roche lobe, so
the binary is comfortably detached.  Closing to `sep0=5` leaves it filling about
98 per cent, which is close enough to contact that the outer layers will start to
stream, so move `sep0` in gradually rather than in one step if that is where you
are heading.

Leave `sph.start2u` out and the second body becomes a single point mass of mass
`mbh` instead, which is what you want when the companion is small compared with
the smoothing lengths of the gas it moves through.

`sepfinal` is set equal to `sep0` because this run is not meant to scan.  Leaving
`sepfinal` at its default stops the run during setup, which is the intended guard
against scanning somewhere by accident.  To walk the separation inward instead
and trace a quasi-equilibrium sequence, set `sepfinal` smaller than `sep0` and
give `tscanon`.  The tutorial "A corotating scan" in the documentation covers
that.

The drag switches off at `treloff`, after which the binary is released and orbits
in the inertial frame.  Set `treloff` beyond `tf` if you would rather it stayed
corotating for the whole run.
