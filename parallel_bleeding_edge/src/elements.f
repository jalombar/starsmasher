      subroutine elements(m1,m2,xrel,yrel,zrel,vxrel,vyrel,vzrel,
     $     a,ecc,lx,ly,lz)
      implicit none
      real*8 m1,m2,xrel,yrel,zrel,vxrel,vyrel,vzrel
      real*8 e,l,lx,ly,lz
      real*8 g
      parameter(g=1.d0)
      real*8 mu,semilatusrectum,ecc,k
      real*8 r
      real*8 w12
      real*8 ecc2,a
      logical verbose
      parameter(verbose=.false.)

      r=(xrel**2+yrel**2+zrel**2)**0.5d0
      if(verbose) write(69,*)'separation r=',r

      w12=-g*m1*m2/r
      mu=m1*m2/(m1+m2)
c      write(69,*)'reduced mass mu=',mu
      e=0.5d0*mu*(vxrel**2+vyrel**2+vzrel**2)+w12
      if(verbose)write(69,'(a,3g13.4)')'total orbital energy=',e
      k=g*m1*m2
c      write(69,*)'force constant k=',k
      a=-0.5d0*k/e

      lz=mu*(xrel*vyrel-yrel*vxrel)
      lx=mu*(yrel*vzrel-zrel*vyrel)
      ly=mu*(zrel*vxrel-xrel*vzrel)
      l=(lx**2+ly**2+lz**2)**0.5d0
      if(verbose)write(69,'(a,4g13.4)')'total angular momentum=',
     $     l,lx,ly,lz


      semilatusrectum=l**2/(mu*k)
      if(verbose)write(69,*)'semi-latus rectum alpha=',semilatusrectum

      ecc2=1.d0+2.d0*e*l**2.d0/(mu*k**2.d0)
      if(ecc2.lt.0.d0) then
         write(69,'(a,9g17.9)')'whoops, ecc squared=',ecc2,e,l,mu,k,m1,m2
         ecc2=0.d0
         stop
      endif
      ecc=ecc2**0.5d0
      if(verbose) then
         write(69,*)'eccentricity ecc=',ecc
         write(69,*)'apastron sep rmax=',semilatusrectum/(1.d0-ecc)
         write(69,*)'periastron sep rmin=',semilatusrectum/(1.d0+ecc)
      endif

      end
