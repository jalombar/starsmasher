      subroutine gettemperature(q,r,x3)
c     Temperature x3 for an equation of state with both ideal gas and
c     radiation pressure: u_gas + u_rad = u, with u_gas proportional to T and
c     u_rad proportional to T^4, gives the quartic
c
c        x^4 + q x + r = 0,   q = qconst*rho/mu > 0,   r = -u*rho/arad < 0.
c
c     f(x) = x^4+qx+r increases and is convex for x>0, so it has exactly one
c     positive root, and Newton's method started above that root comes down
c     onto it monotonically.  Both x^4 and qx are at most -r at the root, so
c     x0 = min((-r)^(1/4), -r/q) lies above it, and by no more than a factor
c     of about 2 (at the root one of the two terms is at least -r/2), so
c     Newton converges to machine precision within 6 or so iterations.
c
c     This replaces the closed-form solution (Stillwell's reduction to a
c     resolvent cubic, by Scott Fleming and James Lombardi, 2002-2003), which
c     loses accuracy as radiation comes to dominate: the error in T reaches
c     ~1e-4 when the gas holds 1e-6 of the internal energy, ~2% at 1e-7 and
c     ~100% below 1e-9, and roundoff can make the square root at its last
c     step negative, returning NaN.  That NaN, from a hot, nearly isolated
c     particle orbiting close to a black hole, stopped a run by making the
c     timestep NaN.

      implicit none
      real*8 q,r,x3,dx
      integer it

c     no internal energy (or a negative one, which has no temperature); a
c     NaN goes on through, so that it shows
      if(r.ge.0.d0) then
         x3=0.d0
         return
      endif

      x3=(-r)**0.25d0
      if(q.gt.0.d0) x3=min(x3,-r/q)
      do it=1,50
         dx=(x3**4+q*x3+r)/(4.d0*x3**3+q)
         x3=x3-dx
         if(abs(dx).le.1.d-15*x3) return
      enddo

      end


      subroutine getlna_from_u(ucgs,rhocgs,mucgs,alna)
c     Buoyancy of Gaburov, Lombardi & Portegies Zwart (2008, MNRAS 383, L5),
c     their eq. 2, from the specific internal energy.  The specific entropy of
c     a monatomic ideal gas plus radiation is
c
c        s - s0 = (3k/2mu) [ ln(k T/(rho^(2/3) mu)) + (8/3)(1-beta)/beta ]
c               = (3k/2mu) ln A,
c
c     so ln A is the true specific entropy up to a factor and an additive
c     constant that depends only on composition.  s0 never enters the equation
c     of state, and it is fixed for a particle because sph does not mix
c     microscopically, so integrating ln A integrates the true entropy without
c     needing s0 -- which is just as well, since s0 needs the individual
c     species and the code carries only the mean molecular weight.
c
c     We store ln A rather than A because A = (kT/(rho^(2/3) mu)) exp((8/3)
c     (1-beta)/beta) overflows double precision once beta drops below about
c     0.004, which is exactly the radiation-dominated regime this variable is
c     meant to serve.
      include 'starsmasher.h'
      real*8 ucgs,rhocgs,mucgs,alna,temperature,pgas,prad,betaone

      call gettemperature(qconst*rhocgs/mucgs,-ucgs*rhocgs/arad,
     $     temperature)
      pgas=rhocgs*boltz*temperature/mucgs
      prad=arad*temperature**4/3.d0
      betaone=pgas/(pgas+prad)
      alna=log(boltz*temperature/(rhocgs**(2.d0/3.d0)*mucgs))
     $     +8.d0/3.d0*(1.d0-betaone)/betaone

      return
      end


      subroutine getT_from_lna(alna,rhocgs,mucgs,tguess,temperature)
c     Invert the buoyancy for the temperature.  Writing y=ln T and using
c     (1-beta)/beta = P_rad/P_gas = a mu T^3/(3 rho k), eq. 2 becomes
c
c        y + c exp(3y) = D,   c = 8 a mu/(9 rho k),
c        D = ln A - ln(k/(rho^(2/3) mu)).
c
c     The internal energy gives a quartic that gettemperature solves in closed
c     form because gas and radiation each contribute a power of T.  Entropy
c     contributes ln T from the gas and T^3 from the radiation, so there is no
c     polynomial form; it is a Lambert W problem, and the W argument overflows
c     for a radiation-dominated envelope.  Solving on the log scale instead
c     avoids that, and both terms rise with y, so the root is unique and a
c     bracketed Newton converges from anywhere.
c
c     tguess>0 is used as the starting point.  The previous step's temperature
c     is a very good one and makes this one or two iterations; pass a
c     non-positive value when there is nothing to go on.
      include 'starsmasher.h'
      real*8 alna,rhocgs,mucgs,tguess,temperature
      real*8 clna,dlna,ylna,ylolna,yhilna,flna,fplna,e3ylna,ynewlna
      integer itlna

      clna=8.d0*arad*mucgs/(9.d0*rhocgs*boltz)
      dlna=alna-log(boltz/(rhocgs**(2.d0/3.d0)*mucgs))

clna     ylna=D always overshoots, since flna(D)=clna*exp(3D)>0.  Walk down from there for
clna     a lower bound; flna tends to ylna-D as ylna falls, so this always terminates.
      yhilna=dlna
      ylolna=dlna-1.d0
      do itlna=1,200
         e3ylna=exp(min(3.d0*ylolna,7.d2))
         if(ylolna+clna*e3ylna-dlna.lt.0.d0) goto 10
         ylolna=dlna-2.d0**itlna
      enddo
      stop 'getT_from_lna: could not bracket the root'
 10   continue

      ylna=0.5d0*(ylolna+yhilna)
      if(tguess.gt.0.d0) then
         ylna=log(tguess)
         if(ylna.le.ylolna .or. ylna.ge.yhilna) ylna=0.5d0*(ylolna+yhilna)
      endif

      do itlna=1,100
         e3ylna=exp(min(3.d0*ylna,7.d2))
         flna=ylna+clna*e3ylna-dlna
         if(flna.gt.0.d0) then
            yhilna=ylna
         else
            ylolna=ylna
         endif
         fplna=1.d0+3.d0*clna*e3ylna
         ynewlna=ylna-flna/fplna
clna     keep Newton inside the bracket; fall back on bisection when itlna leaves
         if(ynewlna.le.ylolna .or. ynewlna.ge.yhilna) ynewlna=0.5d0*(ylolna+yhilna)
         if(abs(ynewlna-ylna).lt.1.d-14*abs(ynewlna)) then
            ylna=ynewlna
            goto 20
         endif
         ylna=ynewlna
      enddo
 20   continue

      temperature=exp(ylna)

      return
      end


      real*8 function uofstored(i)
c     The code-unit specific internal energy of particle i, whatever variable
c     u(i) currently holds.  Callers that want only the energy -- the energy
c     budget in output, changetf's bookkeeping -- can use this rather than
c     repeating the conversion inline, which is how the nintvar=1 version of it
c     came to be written out in five separate places.
      include 'starsmasher.h'
      integer i
      real*8 rhocgsi,tempi

      if(nintvar.eq.1) then
         uofstored=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
      else if(nintvar.eq.3) then
         rhocgsi=rho(i)*munit/runit**3.d0
         call getT_from_lna(u(i),rhocgsi,meanmolecular(i),-1.d0,tempi)
         uofstored=(1.5d0*boltz*tempi/meanmolecular(i)
     $        +arad*tempi**4/rhocgsi)/(gravconst*munit/runit)
      else
         uofstored=u(i)
      endif

      return
      end
