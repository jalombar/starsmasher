c     The limits of the step of one gas particle (see dtlimit in tstep.f): from hpi, vxi, vyi, vzi, uiji (uijmax), axi, ayi, azi, ui and udoti it makes dtivel, dtiacc, dtiu,
c     dtiacc4 and blim, the step they give together.
      dtivel=cn1*hpi/sqrt(max(uiji,(vxi**2+vyi**2+vzi**2)/16))
      dtiacc=cn2*hpi**0.5d0/(axi**2+ayi**2+azi**2)**0.25d0
      dtiu=cn3*ui/dabs(udoti)
      dtiacc4=cn4*sqrt(uiji)/(axi**2+ayi**2+azi**2)**0.5d0
      blim=(1.d0/dtivel+1.d0/dtiacc+1.d0/dtiu+1.d0/dtiacc4)**(-1.d0)
