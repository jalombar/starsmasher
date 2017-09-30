      subroutine tabulinit
c     initialize look-up tables:
      include 'starsmasher.h'
      integer i
      real*8 u2,dig,ig,w,dw,dwdh,dphidh,dphi332hdh,dphic42hdh
      real*8 dw332h,dw332hdh,w332h,dwc42h,dwc42hdh,wc42h

      if(myrank.eq.0) write(69,*)'gflag=',gflag

      if(ncooling.ne.0)then
c         if(myrank.eq.0) write(69,*)'prepare to use opacity tables'
c         call prepare_opacity
c         if(myrank.eq.0) write(69,*)'about to read sph.opacity'
         call readinkappatable
      endif
      if(neos.eq.2)call readineostable ! use tabulated equation of state
      if (myrank.eq.0) write(69,*)'nkernel is',nkernel
c     compute normalization constant:
      ctab=dble(ntab-1)/4.d0
c     compute tabulations of w(u) and dw(u)/du:
      do i=1,ntab
         u2 = 4.d0*dble(i-1)/dble(ntab-1)

         if(u2.lt.0 .or. u2.gt.4d0) then
            write(69,*)'u=r/h is out of range: u^2=',u2
            stop
         endif

         gtab(i)    = ig(sqrt(u2))
         dgtab(i)   = dig(sqrt(u2))
         dgdhtab(i) = dig(sqrt(u2))*u2

c         gtab(i)    = w(sqrt(u2))
c         dgtab(i)   = dwtab(sqrt(u2))
c         dgdhtab(i) = u2*dwtab(sqrt(u2))
         
         if (nkernel.eq.0) then
c     cubic spline
            wtab(i)=w(sqrt(u2))
            dwtab(i)=dw(sqrt(u2))
            dwdhtab(i)=dwdh(sqrt(u2))
            dphidhtab(i)=dphidh(sqrt(u2))
         elseif (nkernel.eq.1) then
c     Wendland 3,3 kernel
            wtab(i)=w332h(sqrt(u2))
            dwtab(i)=dw332h(sqrt(u2))
            dwdhtab(i)=dw332hdh(sqrt(u2))
            dphidhtab(i)=dphi332hdh(sqrt(u2))
         elseif (nkernel.eq.2) then
c     Wendland C4 kernel
            wtab(i)=wc42h(sqrt(u2))
            dwtab(i)=dwc42h(sqrt(u2))
            dwdhtab(i)=dwc42hdh(sqrt(u2))
            dphidhtab(i)=dphic42hdh(sqrt(u2))
         else
            write(69,*) 'nkernel must equal 0, 1, or 2'
            stop
         endif
      enddo

      return
      end
***********************************************************************
      real*8 function iw(u)
      implicit none
      real*8 u
      if (u < 0)  then
         iw = 0
      else if (u .lt. 1) then
         iw = u**3/3.d0 - 0.3d0*u**5 + 0.125d0*u**6
      else if (u .lt. 2) then
         iw = -1.d0/60 + 2.d0/3*u**3 - 0.75d0*u**4
     $        + 0.3d0*u**5 - u**6/24.d0
      else
         iw = 1/4.d0
      end if
      iw = iw * 4
      return
      end
***********************************************************************
      real*8 function diw(u)
      implicit none
      real*8 u, pi_const, w
      parameter(pi_const = 3.1415926535897932384626d0)
      diw = 4*pi_const * u*u * w(u)
      return
      end
***********************************************************************
      real*8 function ig(u)
c     if gflag=0 in the sph.input file, then the function ig is the g
c     defined by equation (a2) and shown in figure a1 of gaburov et
c     al. (2010).  the default value of gflag is 1, in which case the g
c     function is the same as in figure a1 except that g=1 for u<1.  we
c     find that gflag=1 works better in cases where there are massive
c     compact objects or core particles.
      
c     you can change g if you want to, but there is no compelling
c     reason to do so.  the g function is used in determining smoothing
c     lengths (see equation a1).  g should be a smoothly varying
c     function, but to understand g pretend that we were to let it be
c     the step function with value 1 for u<2 and value 0 for u>2.  then
c     n_i on the left hand side of eq. (a1) would just be the number of
c     neighbors for particle i.  so g helps count the neighboring
c     particles, although for some smoothly varying g function not all
c     neighbors are counted with a full value 1 but rather as the local
c     value of g.  for our g functions the actual number of
c     neighbors is about 1.9*n_i, where n_i is the nnopt value set in
c     sph.input or assigned by default in init.f.

      implicit none
      integer gflag
      common/flag/gflag
      real*8 alpha_n, u, iw
      parameter(alpha_n = 4.d0)

      if (u < 1) then
         if(gflag.eq.0)then
            ig = iw(alpha_n * u)
         else
            ig = 1.d0
         endif
      else if (u < 2) then
         ig = iw(alpha_n * (2 - u))
      else
         ig = 0
      end if

      return
      end
***********************************************************************
      real*8 function dig(u)
      implicit none
      integer gflag
      common/flag/gflag
      real*8 alpha_n, u, diw
      parameter(alpha_n = 4.d0)

      if (u .eq. 0) then
         dig = 0
      else if (u < 1) then
         if(gflag.eq.0)then
            dig = +alpha_n/u * diw(alpha_n * u)
         else
            dig = 0.d0
         endif
      else if (u < 2) then
         dig = -alpha_n/u * diw(alpha_n * (2 - u))
      else
        dig = 0
      end if
      end
***********************************************************************
      real*8 function w(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W is the usual cubic spline kernel with compact support 2h
      implicit none
      real*8 u
      if (u < 0) then
         w = 1
      else if (u.lt.1.d0) then
         w=1.d0+u**2.d0*(-1.5d0+0.75d0*u)
      else if (u.lt.2.d0) then
         w=0.25d0*(2.d0-u)**3
      else
         w=0.d0
      endif
      w=w/3.1415926535897932384626d0

      return
      end
************************************************************************
      real*8 function dw(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dr*h^5/r, where h=smoothing
c     length and W is the usual cubic spline kernel w/ compact support 2h
      implicit none
      real*8 u
      if (u.lt.1.d0) then
         dw=-3.d0+2.25d0*u
      else if (u.lt.2.d0) then
         dw=-0.75d0*(2.d0-u)**2/u
      else
         dw=0.d0
      endif
      dw=dw/3.1415926535897932384626d0
      
      return
      end
************************************************************************
      real*8 function dwdh(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dh*h^4, where h=smoothing
c     length and W is the usual cubic spline kernel w/ compact support 2h
      implicit none
      real*8 u,nu,sigma

      nu=3.d0
      sigma=1.d0/3.1415926535897932384626d0
      
      if (u.lt.1.d0) then
         dwdh=sigma*(-nu+3.d0*u**2*(2.d0+nu)/2.d0
     $        -3.d0*u**3*(3.d0+nu)/4.d0)
      else if (u.lt.2.d0) then
c         dwdh=sigma*(-2.d0*nu+3.d0*u*(1.d0+nu)-3.d0*u**2*(2.d0+nu)/2.d0
c     $        +u**3*(3.d0+nu)/4.d0)
         dwdh=sigma*0.25d0*(2.d0-u)**2.d0*((3.d0+nu)*u-2.d0*nu)
      else
         dwdh=0.d0
      endif
      
      return
      end
************************************************************************
      real*8 function dphidh(u)
c     Cubic spline kernel
c     The dphidh(u) function is related to d(phi)/dh, where phi is the
c     specific gravitational potential given by eq. (A1) of Hernquist &
c     Katz (1989,HK) or introduced in equation (A17) of Gaburov et
c     al. (2010).  In HK's equation (A1), epsilon=smoothing length h.
c     Our function g equals HK's function f but with an overall minus
c     sign in the front.  Specifically, if you differentiate their eq.
c     (A1) with respect to h, you get dphidh multiplied by -1/h^2.  In
c     other words, dphidh is -h^2 times the derivative of HK's eq. (A1)
c     with respect to h.  The dphidh values are tabulated in the
c     array dphidhtab.  This array is used in advance.f90 to get the
c     psi_i defined by eq. (A17) of Gaburov et al. (2010).  Specifically,
c     psi_i is the bonet1_psi being calculated in advance.f90 with
c     the line
c              bonet1_psi = bonet1_psi + am(j)*dphidhtab(itab)
c     which is inside a loop over neighbors of i.  The bonet1_psi value
c     is divided later in advance.f90 by h2=h_i^2 to account for dphidh
c     being dimensionless:
c          bonet1_psi = bonet1_psi/h2
c     If we change phi (that is, we change
c     how gravity is softened) then we would have to change the code in
c     sphgrav_lib/grav_force_direct.cu.  The paper
c               http://adsabs.harvard.edu/abs/2001mnras.324..273d
c     says that changing phi
c     could help improve the gravitational softening. The kernel used for
c     softening gravity doesn't need to be the same one used for
c     calculating hydrodynamic forces, although that might be nice.

      implicit none
      real*8 u

      if (u.le.1.d0) then
         dphidh=-2*u**2+1.5d0*u**4-0.6d0*u**5+1.4d0
      else if (u.lt.2.d0) then
c         dphidh=-4*u**2+4*u**3-1.5d0*u**4+0.2d0*u**5+1.6d0
         dphidh=(0.1d0+0.2d0*u)*(2.d0-u)**4.d0
      else
         dphidh=0.d0
      endif
      
      return
      end
************************************************************************
      real*8 function dphi332hdh(u)
c     Wendland 3,3 kernel with compact support 2h
      implicit none
      real*8 u,q
      
      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dphi332hdh=(35*(1 - q)**9*(7 + 63*q + 237*q**2 + 453*q**3 + 384*q**4))/128d0
      
      return
      end
************************************************************************
      real*8 function dphic42hdh(u)
c     Wendland C4 kernel with compact support 2h
      implicit none
      real*8 u,q
      
      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dphic42hdh=(55*(1 - q)**7*(1 + 7*q + 19*q**2 + 21*q**3))/32d0
      
      return
      end
************************************************************************
cC Wendland 33 kernel
c      real*8 function w33(q)
c      implicit none
c      real*8 q,sigW
c      sigW =1365.d0/(64.d0*3.1415926535897932384626d0)
c
c      if (q.lt.0.d0) then
c          w33 = 0.d0
c      elseif (q.lt.1.d0) then
c          w33 = sigW*(1.d0-q)**8.d0*(32.d0*q**3.d0+25.d0*q**2.d0+8.d0*q+1.d0)
c      else
c          w33 = 0.d0
c      endif
c
c      return
c      end
c
c      real*8 function dw33(q)
c      implicit none
c      real*8 q,sigW
c      sigW =1365.d0/(64.d0*3.1415926535897932384626d0)
c
c      if ((q.gt.0.d0).and.(q.lt.1.d0)) then
c         dw33 = sigW*(-8.d0*(1.d0-q)**7.d0*(32.d0*q**3.d0+25.d0*q**2.d0+
c     $                 8.d0*q+1.d0)/q +
c     $                 (1.d0-q)**8.d0*(96.d0*q+50.d0+8.d0/q))
c      else
c         dw33 = 0.d0
c      endif
c
c      return
c      end
c
c 
c      real*8 function dw33dh(q)
c      implicit none
c      real*8 q,sigW
c      sigW =1365.d0/(64.d0*3.1415926535897932384626d0)
c
c      if (q.lt.1.d0) then
c         dw33dh =  sigW*(-3.d0*(1.d0-q)**8.d0*
c     &                        (32.d0*q**3.d0+25.d0*q**2.d0+8.d0*q+1.d0)
c     &                   -(1.d0-q)**8.d0*(96.d0*q**3.d0+50.d0*q**2.d0+8.d0*q)
c     &                   +8.d0*(1.d0-q)**7.d0*q*
c     &                   (32.d0*q**3.d0+25.d0*q**2.d0+8.d0*q+1.0d0))
c      else
c         dw33dh = 0.d0
c      endif
c
c      return
c      end

      real*8 function w332h(u)  ! NOTE... function of u=r/h
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W=Wendland 3,3 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =1365.d0/(512.d0*3.1415926535897932384626d0)) ! Note the 512
!     We're using 512 instead of 64 because the output is W*h^3 and not W*(2h)^3
      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      w332h = sigW*(1.d0-q)**8.d0*(32.d0*q**3.d0+25.d0*q**2.d0+8.d0*q+1.d0)

      return
      end
 
      real*8 function dw332h(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dr*h^5/r, where h=smoothing
c     length & W=Wendland 3,3 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =1365.d0/(512.d0*3.1415926535897932384626d0)) ! Note the 512
!     Note: the output is dW/dr*h^5/r = d(Wh^3)/dr*h^2/r = 0.25 * d(Wh^3)/dr*(2h)^2/r
!     so below we have -11/2 instead of -22

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dw332h = -11.d0/2.d0*sigW*(1.d0-q)**7.d0*(1+7*q+16*q**2)

      return
      end
 
      real*8 function dw332hdh(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dh*h^4, where h=smoothing
c     length & W=Wendland 3,3 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =1365.d0/(512.d0*3.1415926535897932384626d0)) ! Note the 512

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dw332hdh =  sigW*(1.d0-q)**7.d0*(-3-21*q-29*q**2+133*q**3+448*q**4)

      return
      end

************************************************************************
cC Wendland c4 kernel
c      real*8 function wc4(q)
c      implicit none
c      real*8 q,sigW
c      sigW =495.d0/(32.d0*3.1415926535897932384626d0)
c
c      if (q.lt.0.d0) then
c          wc4 = 0.d0
c      elseif (q.lt.1.d0) then
c          wc4 = sigW*(1.d0-q)**6.d0*(35.d0/3d0*q**2.d0+6.d0*q+1.d0)
c      else
c          wc4 = 0.d0
c      endif
c
c      return
c      end

      real*8 function wc42h(u)  ! NOTE... function of u=r/h
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W=Wendland C4 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =495.d0/(256.d0*3.1415926535897932384626d0)) ! Note the 256 instead of 32 (this is because the output is W*h^3 and not W*(2h)^3

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      wc42h = sigW*(1.d0-q)**6.d0*(35.d0/3d0*q**2.d0+6.d0*q+1.d0)
      return
      end
 
      real*8 function dwc42h(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dr*h^5/r, where h=smoothing
c     length & W=Wendland C4 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =495.d0/(256.d0*3.1415926535897932384626d0)) ! Note the 256

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dwc42h = -14.d0/3.d0*sigW*(1.d0-q)**5.d0*(1+5*q)

      return
      end
 
      real*8 function dwc42hdh(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity dW/dh*h^4, where h=smoothing
c     length & W=Wendland C4 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =495.d0/(256.d0*3.1415926535897932384626d0)) ! Note the 256

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      dwc42hdh =  sigW*(1.d0-q)**5.d0*(-9-45*q+5*q**2+385*q**3)/3d0

      return
      end
