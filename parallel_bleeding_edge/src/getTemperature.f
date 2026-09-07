      subroutine gettemperature(q,r,x3)
c     subroutine to solve 4th order equations to determine the temperature x3
c     for an equation of state with both ideal gas and radiation pressure.
c     written by scott fleming 10/04/02 and james lombardi 2002-2003

c     the fourth order equation comes from u_gas+ u_rad = u, with
c     u_gas proportional to t and u_rad proportional to t^4

c     in general, we can transform a 4th order equation to x^4+px^2+qx+r=0
c     (see pages 57-58 of stillwell's "mathematics and its history" text)
c     but we fortunately don't even have an x^2 term (that is, p=0).

c     follow stillwell, we can transform this into a cubic equation:
c     first solve for y by using the fact that b^2-4ac=0
c     equation is then:  y^3=ry+q^2/8
c     using the solution of cubic equations found in stillwell page 55:

      implicit none
      real*8 q,r,k,b,piece1,piece2
      real*8 y1,y2,yy,aa,b2,c2,x3,kh
c      real*8 piece2old,x3old

      k = 0.125d0*q**2
      kh=0.5d0*k
      if(kh**2-(r/3.d0)**3.le.0.d0)then
         write(69,*) k,r,kh**2-(r/3.d0)**3
         stop 'bad input: imaginary results?'
      endif

      piece1 = kh+(kh**2-(r/3.d0)**3)**0.5d0
c      piece2old = kh-(kh**2-(r/3.d0)**3)**0.5d0
      piece2 = (r/3.d0)**3.d0/piece1

c      write(69,*)piece2old,piece2

      y1 = piece1**(1.d0/3.d0)
      
c     fortran can't handle cube roots of neg. #'s
      y2 = -dabs(piece2)**(1.d0/3.d0)
      yy = y1+y2

c     equation to solve: (x^2+p+y)^2=ax^2+bx+c
c     now take square root of both sides with:

      aa = 2.d0*yy
      b = -q
c      c = -r+y**2

c     re-writing ax^2+bx+c as a square then solving the equation we
c     obtain 2 results:
c     x^2 + (-(a^(1/2)))x + (-b/(2(a)^(1/2))+p+y) = 0 (1)
c     or
c     x^2 + (a^(1/2))x + (b/(2(a)^(1/2))+p+y) = 0     (2)

c     our solution we're interested in:
      b2 = aa**0.5d0
      c2 = 0.5d0*b/b2 + yy

c     therefore, we once again have x^2+bx+c=0, and our answer we want is:
c      x3old = 0.5d0*(-b2 + (b2**2-4.d0*c2)**0.5d0)
      x3 = -2.d0*c2/(b2 + (b2**2-4.d0*c2)**0.5d0)

c      write(69,*) 'temperature components', x3old, x3

      if(piece1.lt.0.d0) write(69,*)
     $     'piece 1 lt 0',k,r,piece1,piece2
c      if(piece1.eq.-piece2)then
      if(b2.eq.0.d0)then
c         write(69,*)
c     $        'piece 1 eq -piece 2 (rad pressure dominates)',
c     $        k,r,piece1,piece2,b2,c2,x3,(-r)**0.25d0
         x3=(-r -q*( -r-q*(-r)**0.25d0 )**0.25d0)**0.25d0
c         write(69,*) x3
c         write(69,*) piece1,piece2
      endif
      if(piece2.ge.0.d0) x3=-(r+(r/q)**4)/q
      
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
