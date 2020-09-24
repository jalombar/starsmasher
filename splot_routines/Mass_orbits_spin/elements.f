      subroutine elements(m1,m2,xrel,yrel,zrel,vxrel,vyrel,vzrel,
     $     a,ecc,lx,ly,lz)
      implicit none
      real*8 m1,m2,xrel,yrel,zrel,vxrel,vyrel,vzrel
      real*8 E,l,lx,ly,lz
      real*8 G
      parameter(G=1.d0)
      real*8 mu,semilatusrectum,ecc,k
      real*8 r
      real*8 W12
      real*8 ecc2,a
      logical verbose
      parameter(verbose=.false.)

      r=(xrel**2+yrel**2+zrel**2)**0.5d0
      if(verbose) print *,'separation r=',r

      W12=-G*m1*m2/r
      mu=m1*m2/(m1+m2)
c      print *,'reduced mass mu=',mu
      E=0.5d0*mu*(vxrel**2+vyrel**2+vzrel**2)+W12
      if(verbose)write(6,'(a,3g13.4)')'Total orbital energy=',E
      k=G*m1*m2
c      print *,'force constant k=',k
      a=-0.5d0*k/E

      lz=mu*(xrel*vyrel-yrel*vxrel)
      lx=mu*(yrel*vzrel-zrel*vyrel)
      ly=mu*(zrel*vxrel-xrel*vzrel)
      l=(lx**2+ly**2+lz**2)**0.5d0
      if(verbose)write(6,'(a,4g13.4)')'Total angular momentum=',
     $     l,lx,ly,lz


      semilatusrectum=l**2/(mu*k)
      if(verbose)print *,'semi-latus rectum alpha=',semilatusrectum

      ecc2=1.d0+2.d0*E*l**2.d0/(mu*k**2.d0)
      if(ecc2.lt.0.d0) then
         write(6,'(a,9g17.9)')'WHOOPS, ecc squared=',ecc2,E,l,mu,k,m1,m2
         ecc2=0.d0
         sToP
      endif
      ecc=ecc2**0.5d0
      if(verbose) then
         print *,'eccentricity ecc=',ecc
         print *,'apastron sep rmax=',semilatusrectum/(1.d0-ecc)
         print *,'periastron sep rmin=',semilatusrectum/(1.d0+ecc)
      endif

      end
