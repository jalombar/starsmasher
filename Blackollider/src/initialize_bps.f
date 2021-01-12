      subroutine bps

c to follow the logic in this routine,
c it will be helpful to look at section 11.8 of the thornton and
c mation text.  in their notation, the
c inertial frame is the x' frame and the corotating frame is the x frame.
c we want to find a rotation matrix that brings coordinates from the x
c frame to the x' frame.  that is x'=(matrix)*x, where * is just matrix
c multiplication and where x and x' are 3d position (or velocity) vectors.
c in the text, they give the inverse of the matrix that we want.  that is,
c they give a 3x3 lambda matrix such that x=(lambda)*x'.  you can get the
c inverse of lambda from lambda by letting theta -> -theta, phi -> -psi,
c and psi -> -phi in equation (11.99) as you can convince yourself from
c equations (11.91)-(11.98).
c 
c you can get the rotation angles for any particular binary by looking at
c figure 11-9.  the x_3 axis is in the direction of the binary's angular
c momentum.  the x'_3 axis is the z-axis in the inertial frame (the
c z'-axis).  the line of nodes is along the intersection of the orbital
c plane and the z'=0 plane.  the stars are constructed to be separated
c along the x-axis in the corotating frame.... that's the x_1 axis in
c figure 11-9.  once you know what each axis in figure 11-9 represents in
c our problem, then appropriate dot and cross products can be used to find
c the three angles in terms of things you can calculate, like the angular
c momentum of the binary, etc.
c 
c the good thing is that i'm confident in the bps code now.  there are
c enough stop statements in it that if something goes wrong, it will
c break.  and that doesn't seem to be happening fortunately!

      include 'starsmasher.h'
      real*8 vxcm,vycm,vzcm,xcm,ycm,zcm
      integer n3,i,nchk,corepts
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      real*8 costh1,sinth1,costh2,sinth2,cosps1,sinps1,ps2,cosph1,sinph1
      real*8 xold,yold,zold,vxold,vyold,vzold
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2,am3
      real*8 vxcm1,vycm1,vzcm1,vxcm2,vycm2,vzcm2
      real*8 xcm3,ycm3,zcm3
      real*8 vxcm3,vycm3,vzcm3
      real*8 lambda11,lambda12,lambda13
      real*8 lambda21,lambda22,lambda23
      real*8 lambda31,lambda32,lambda33
      character*7 dummy
      real*8 am1chk,am2chk,am3chk
      real*8 lx,ly,lz,ltotbin
      real*8 lonx,lony,lonhatx,lonhaty,lontot
      real*8 crossx,crossy,crossz,crosstot
      real*8 deltaxbin,deltaybin,deltazbin,deltarbin
      real*8 deltavxbin,deltavybin,deltavzbin
      real*8 mubin
      real*8 deltax1,deltay1,deltaz1,deltavx1,deltavy1,deltavz1
      real*8 deltax2,deltay2,deltaz2,deltavx2,deltavy2,deltavz2
      real*8 deltax3,deltay3,deltaz3,deltavx3,deltavy3,deltavz3
      real*8 egsol,solrad
      parameter(egsol=1.9891d+33,solrad=6.9599d10)
      real*8 xcombin,ycombin,zcombin
      real*8 vxcombin,vycombin,vzcombin
      real*8 amtot
      real*8 xnew1,ynew1,znew1
      real*8 xnew2,ynew2,znew2
      real*8 deltarbinnew
      integer nbinary,n2
      logical tweakv
      real*8 bincomx,bincomy,bincomz
      real*8 bincomvx,bincomvy,bincomvz
      real*8 deltarsphbin

      write(69,*)'bps: reading in input.bs'
      corepts=0

      open(30,file='input.bs')
      read(30,*) dummy
      if(dummy.ne.'binary:')then
         write(69,*) 'input.bs should have binary first'
         stop
      endif
      read(30,*) am1chk
      read(30,*) deltax1,deltay1,deltaz1
      read(30,*) deltavx1,deltavy1,deltavz1
      read(30,*) am2chk
      read(30,*) deltax2,deltay2,deltaz2
      read(30,*) deltavx2,deltavy2,deltavz2
      read(30,*) dummy
      if(dummy.ne.'single:')then
         write(69,*) 'input.bs should have single star second'
         stop
      endif
      read(30,*) am3chk
      read(30,*) deltax3,deltay3,deltaz3
      read(30,*) deltavx3,deltavy3,deltavz3
      close(30)

      deltax1=deltax1*solrad/runit
      deltay1=deltay1*solrad/runit
      deltaz1=deltaz1*solrad/runit
      deltax2=deltax2*solrad/runit
      deltay2=deltay2*solrad/runit
      deltaz2=deltaz2*solrad/runit
      deltax3=deltax3*solrad/runit
      deltay3=deltay3*solrad/runit
      deltaz3=deltaz3*solrad/runit
      deltavx1=deltavx1*(egsol/munit*runit/solrad)**0.5d0
      deltavy1=deltavy1*(egsol/munit*runit/solrad)**0.5d0
      deltavz1=deltavz1*(egsol/munit*runit/solrad)**0.5d0
      deltavx2=deltavx2*(egsol/munit*runit/solrad)**0.5d0
      deltavy2=deltavy2*(egsol/munit*runit/solrad)**0.5d0
      deltavz2=deltavz2*(egsol/munit*runit/solrad)**0.5d0
      deltavx3=deltavx3*(egsol/munit*runit/solrad)**0.5d0
      deltavy3=deltavy3*(egsol/munit*runit/solrad)**0.5d0
      deltavz3=deltavz3*(egsol/munit*runit/solrad)**0.5d0

c     position and velocity of center of mass of binary:
      xcombin=(am1chk*deltax1+am2chk*deltax2)/(am1chk+am2chk)
      ycombin=(am1chk*deltay1+am2chk*deltay2)/(am1chk+am2chk)
      zcombin=(am1chk*deltaz1+am2chk*deltaz2)/(am1chk+am2chk)
      vxcombin=(am1chk*deltavx1+am2chk*deltavx2)/(am1chk+am2chk)
      vycombin=(am1chk*deltavy1+am2chk*deltavy2)/(am1chk+am2chk)
      vzcombin=(am1chk*deltavz1+am2chk*deltavz2)/(am1chk+am2chk)

      write(69,'(a,3g14.6)') 'center of mass position of binary:',
     $     xcombin*runit/solrad,
     $     ycombin*runit/solrad,
     $     zcombin*runit/solrad
      write(69,'(a,3g14.6)') 'center of mass velocity of binary:',
     $     vxcombin*(egsol/munit*runit/solrad)**(-0.5d0),
     $     vycombin*(egsol/munit*runit/solrad)**(-0.5d0),
     $     vzcombin*(egsol/munit*runit/solrad)**(-0.5d0)

      write(69,*)'star 1 mass [msun]=',am1chk
      write(69,*)'star 2 mass [msun]=',am2chk
      write(69,*)'star 3 mass [msun]=',am3chk

      deltaxbin=deltax2-deltax1
      deltaybin=deltay2-deltay1
      deltazbin=deltaz2-deltaz1
      deltarbin=(deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0
      write(69,'(a,4g14.6)')'input.bs: separation hat of binary=',
     $     deltaxbin/deltarbin,deltaybin/deltarbin,deltazbin/deltarbin,
     $     deltarbin

      deltavxbin=deltavx2-deltavx1
      deltavybin=deltavy2-deltavy1
      deltavzbin=deltavz2-deltavz1

      mubin=am1chk*am2chk/(am1chk+am2chk)
c      lx=mubin*(deltaybin*deltavzbin-deltazbin*deltavybin)
c      ly=mubin*(deltazbin*deltavxbin-deltaxbin*deltavzbin)
c      lz=mubin*(deltaxbin*deltavybin-deltaybin*deltavxbin)
      lx=deltaybin*deltavzbin-deltazbin*deltavybin
      ly=deltazbin*deltavxbin-deltaxbin*deltavzbin
      lz=deltaxbin*deltavybin-deltaybin*deltavxbin
      ltotbin=(lx**2+ly**2+lz**2)**0.5d0
      write(69,'(a,4g14.6)')'input.bs: angular momentum hat of binary=',
     $     lx/ltotbin,ly/ltotbin,lz/ltotbin,ltotbin

      costh1=lz/ltotbin ! dot product of lhat and zhat is cos(theta1)
      sinth1=(1.d0-costh1**2)**0.5d0
      write(69,*)'theta1=',acos(costh1)
      write(69,*)'cos(theta1),sin(theta1)=',costh1,sinth1
      if(dabs(sin(acos(costh1))/sinth1-1.d0).gt.1.d-15 .or.
     $     dabs(costh1**2+sinth1**2-1.d0).gt.1.d-15) then
         write(69,*) 'problem setting theta1'
         stop
      endif

c     lon=line of nodes (its z component is zero)
c     the line of nodes is perpendicular to both the z direction and the
c     angular momentum vector.  so get lon by taking (z hat)x(angular mom),
c     that is, the cross product of zhat and the binary's angular momentum
      lonx=-ly
      lony= lx
      lontot=(lonx**2+lony**2)**0.5d0
      lonhatx=lonx/lontot
      lonhaty=lony/lontot
c     the angle psi1 is between the separation vector r2-r1 and the line of nodes.
c     so we can get cos(psi1) from the dot product between these vectors:
      cosps1=(lonhatx*deltaxbin+lonhaty*deltaybin)/deltarbin
c     and we can get sin(psi1) from the cross product between these vectors:
      crossx= lony*deltazbin
      crossy=-lonx*deltazbin
      crossz= lonx*deltaybin-lony*deltaxbin
      crosstot=(crossx**2+crossy**2+crossz**2)**0.5d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this was bad!  it forces sin(ps1) to be positive
c     sinps1=crosstot/(lontot*deltarbin)  
c     this was bad!  it forces sin(ps1) to be positive
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      sinps1=(crossx*lx+crossy*ly+crossz*lz)/(lontot*deltarbin*ltotbin)

      write(69,*)'cos(psi1),sin(psi1)=',cosps1,sinps1
      if(dabs(cosps1**2+sinps1**2-1.d0).gt.1.d-15) then
         write(69,*) 'problem setting psi1'
         stop
      endif

c     the angle phi1 is between the x axis and the line of nodes.
c     so we can get cos(phi1) from the dot product between these vectors:
      cosph1=lonhatx
c     and we can get sin(phi1) from the cross product between these vectors...
c     because lonz=0, the result is particularly simple:
      sinph1=lonhaty
      write(69,*)'cos(phi1),sin(phi1)=',cosph1,sinph1
      if(dabs(cosph1**2+sinph1**2-1.d0).gt.1.d-15) then
         write(69,*) 'problem setting phi1'
         stop
      endif

      costh2=1.d0
      sinth2=0.d0
      ps2=0.d0
c      write(69,*)'rotation angles for binary (in radians):'
c      write(69,*)'phase phi1=',ph1,'theta1=',th1,' psi1=',ps1
      write(69,*)'rotation angles for star 2 (in radians):'
      write(69,*)'theta2=',acos(costh2),' psi2=',ps2

c     determine coefficients of rotation matrix for the binary
      lambda11= cosph1*cosps1-costh1*sinps1*sinph1
      lambda21= sinph1*cosps1+costh1*sinps1*cosph1
      lambda31= sinth1*sinps1
      lambda12=-cosph1*sinps1-costh1*cosps1*sinph1
      lambda22=-sinph1*sinps1+costh1*cosps1*cosph1
      lambda32= sinth1*cosps1
      lambda13= sinph1*sinth1
      lambda23=-cosph1*sinth1
      lambda33= costh1

      write(69,*)'the rotation matrix has been determined.'
      write(69,*)'lets test where it would place the binary components:'

      xold=-deltarbin*am2chk/(am1chk+am2chk)
      yold=0.d0
      zold=0.d0
      xnew1= lambda11*xold  + lambda12*yold  + lambda13*zold + xcombin
      ynew1= lambda21*xold  + lambda22*yold  + lambda23*zold + ycombin
      znew1= lambda31*xold  + lambda32*yold  + lambda33*zold + zcombin
      write(69,'(a,3g14.6)')'        star 1 will be at (natural units)',
     $     xnew1*runit/solrad,ynew1*runit/solrad,znew1*runit/solrad

      xold=deltarbin*am1chk/(am1chk+am2chk)
      yold=0.d0
      zold=0.d0
      xnew2= lambda11*xold  + lambda12*yold  + lambda13*zold + xcombin
      ynew2= lambda21*xold  + lambda22*yold  + lambda23*zold + ycombin
      znew2= lambda31*xold  + lambda32*yold  + lambda33*zold + zcombin
      write(69,'(a,3g14.6)')'        star 2 will be at (natural units)',
     $     xnew2*runit/solrad,ynew2*runit/solrad,znew2*runit/solrad

      if(dabs(deltax1/xnew1-1.d0).gt.1.d-4
     &     .or.dabs(deltay1/ynew1-1.d0).gt.1.d-4
     &     .or.dabs(deltaz1/znew1-1.d0).gt.1.d-4
     &     .or.dabs(deltax2/xnew2-1.d0).gt.1.d-4
     &     .or.dabs(deltay2/ynew2-1.d0).gt.1.d-4
     &     .or.dabs(deltaz2/znew2-1.d0).gt.1.d-4) then
         write(69,*)'rotation matrix has failed the test. :( '
         write(69,*)'input.bs value           attempted value'
         write(69,*)'----------------------------------------'
         write(69,*)'x1',deltax1,xnew1
         write(69,*)'y1',deltay1,ynew1
         write(69,*)'z1',deltaz1,znew1
         write(69,*)'x2',deltax2,xnew2
         write(69,*)'y2',deltay2,ynew2
         write(69,*)'z2',deltaz2,znew2
         write(69,*)'this is a serious problem that should not be'
         write(69,*)'fixed simply by relaxing the inequality in the if'
         write(69,*)'statement above.  there has been a problem'
         write(69,*)'determining the rotation angles.'
         stop
      endif

      deltarbinnew=((xnew2-xnew1)**2+
     $     (ynew2-ynew1)**2+(znew2-znew1)**2)**0.5d0
      if(dabs(deltarbinnew/deltarbin-1.d0).gt.1.d-15) then
         write(69,*)'rotation matrix will not preserve binary separation'
         stop
      endif

      write (69,*) 'bps: reading start files ...'
      open(12,file=startfile1,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) nbinary,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold,omega2

      omeg=omega2**0.5d0
      write(69,*)'omega=',omeg

c     let us test how rotation will affect motion'
c      vxold=0.d0
c      vyold=-deltarbin*am2chk/(am1chk+am2chk)*omeg
c      vzold=0.d0
c      vxnew1= lambda11*vxold  + lambda12*vyold  + lambda13*vzold
c     $     + vxcombin
c      vynew1= lambda21*vxold  + lambda22*vyold  + lambda23*vzold
c     $     + vycombin
c      vznew1= lambda31*vxold  + lambda32*vyold  + lambda33*vzold
c     $     + vzcombin
c      write(69,'(a,3g14.6)')'star 1 will move with (natural units)',
c     $     vxnew1*(egsol/munit*runit/solrad)**(-0.5d0),
c     $     vynew1*(egsol/munit*runit/solrad)**(-0.5d0),
c     $     vznew1*(egsol/munit*runit/solrad)**(-0.5d0)

c      vxold=0.d0
c      vyold=deltarbin*am1chk/(am1chk+am2chk)*omeg
c      vzold=0.d0
c      vxnew2= lambda11*vxold  + lambda12*vyold  + lambda13*vzold
c     $     + vxcombin
c      vynew2= lambda21*vxold  + lambda22*vyold  + lambda23*vzold
c     $     + vycombin
c      vznew2= lambda31*vxold  + lambda32*vyold  + lambda33*vzold
c     $     + vzcombin
c      write(69,'(a,3g14.6)')'star 2 will move with (natural units)',
c     $     vxnew2*(egsol/munit*runit/solrad)**(-0.5d0),
c     $     vynew2*(egsol/munit*runit/solrad)**(-0.5d0),
c     $     vznew2*(egsol/munit*runit/solrad)**(-0.5d0)

c      if(dabs(deltax1/xnew1-1.d0).gt.1.d-15
c     &     .or.dabs(deltay1/ynew1-1.d0).gt.1.d-15
c     &     .or.dabs(deltaz1/znew1-1.d0).gt.1.d-15
c     &     .or.dabs(deltax2/xnew2-1.d0).gt.1.d-15
c     &     .or.dabs(deltay2/ynew2-1.d0).gt.1.d-15
c     &     .or.dabs(deltaz2/znew2-1.d0).gt.1.d-15) then
c         write(69,*)'rotation will not lay stars down super accurately'
c         stop
c      endif

      am1=0.d0
      am2=0.d0
      write(69,*)'nbinary=n1+n2=',nbinary
c      amass1=nbinary
      n1=0
      n2=0
      xcm1=0.d0
      ycm1=0.d0
      zcm1=0.d0
      xcm2=0.d0
      ycm2=0.d0
      zcm2=0.d0
      vxcm1=0.d0
      vycm1=0.d0
      vzcm1=0.d0
      vxcm2=0.d0
      vycm2=0.d0
      vzcm2=0.d0
      amtot=0.d0

      do i=1,nbinary
         read (12) xold,yold,zold,am(i),hp(i),rho(i),vxold,vyold,
     $        vzold,vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)

c     place velocities at same time as everything else:
c         vxold=vxold-vxdot(i)*0.5d0*dtold
c         vyold=vyold-vydot(i)*0.5d0*dtold
c         vzold=vzold-vzdot(i)*0.5d0*dtold
         vxold=-omeg*yold
         vyold= omeg*xold
         vzold= 0.d0

         x(i)= lambda11*xold  + lambda12*yold  + lambda13*zold + xcombin
         y(i)= lambda21*xold  + lambda22*yold  + lambda23*zold + ycombin
         z(i)= lambda31*xold  + lambda32*yold  + lambda33*zold + zcombin
         vx(i)=lambda11*vxold + lambda12*vyold + lambda13*vzold+vxcombin
         vy(i)=lambda21*vxold + lambda22*vyold + lambda23*vzold+vycombin
         vz(i)=lambda31*vxold + lambda32*vyold + lambda33*vzold+vzcombin
         
c         if(xold.le.0.d0) then
         if(cc(i).eq.cc(1)) then
            n1=n1+1
            am1=am1+am(i)
            xcm1=xcm1+am(i)*x(i)
            ycm1=ycm1+am(i)*y(i)
            zcm1=zcm1+am(i)*z(i)
            vxcm1=vxcm1+am(i)*vx(i)
            vycm1=vycm1+am(i)*vy(i)
            vzcm1=vzcm1+am(i)*vz(i)
         else
c         else if(xold.gt.0.d0) then
            n2=n2+1
            am2=am2+am(i)
            xcm2=xcm2+am(i)*x(i)
            ycm2=ycm2+am(i)*y(i)
            zcm2=zcm2+am(i)*z(i)
            vxcm2=vxcm2+am(i)*vx(i)
            vycm2=vycm2+am(i)*vy(i)
            vzcm2=vzcm2+am(i)*vz(i)
         endif

      enddo
      read(12) nchk
      close(12)

      if (nchk.ne.nbinary) stop 'bps: problem with file'
      if (n1+n2.ne.nbinary) then
         write(69,*)'n1+n2.ne.nbinary???'
         stop
      endif

      xcm1=xcm1/am1
      ycm1=ycm1/am1
      zcm1=zcm1/am1
      vxcm1=vxcm1/am1
      vycm1=vycm1/am1
      vzcm1=vzcm1/am1
      xcm2=xcm2/am2
      ycm2=ycm2/am2
      zcm2=zcm2/am2
      vxcm2=vxcm2/am2
      vycm2=vycm2/am2
      vzcm2=vzcm2/am2


      deltarsphbin=((xcm2-xcm1)**2+(ycm2-ycm1)**2+(zcm2-zcm1)**2)**0.5d0
      if(dabs(deltarsphbin/deltarbin-1.d0).gt.1.d-4) then
         write(69,*)'separation from sph.start1u=',deltarsphbin
         write(69,*)'separation from input.bs=',deltarbin
         write(69,*)'the binary in sph.start1u has a much'
         write(69,*)'different separation than input.bs suggests.'
         write(69,*)'this is assumed to be a contact binary, so...'
         write(69,*)'velocities *wont* be adjusted to match input.bs'
         write(69,*)'velocities will be from sph.start1u (e=0 orbit)'
         write(69,*)'orientation of binary still comes from input.bs'
         tweakv=.false.
      else
         tweakv=.true.
      endif

      if(tweakv) then
         if(dabs(deltavx1-vxcm1).gt.1.d-2 .or.
     $        dabs(deltavy1-vycm1).gt.1.d-2 .or.
     $        dabs(deltavz1-vzcm1).gt.1.d-2 .or.
     $        dabs(deltavx2-vxcm2).gt.1.d-2 .or.
     $        dabs(deltavy2-vycm2).gt.1.d-2 .or.
     $        dabs(deltavz2-vzcm2).gt.1.d-2) then
            write(69,*)'warning: changing v by a lot from circular value'
            write(69,*) deltavx1-vxcm1
            write(69,*) deltavy1-vycm1
            write(69,*) deltavz1-vzcm1
            write(69,*) deltavx2-vxcm2
            write(69,*) deltavy2-vycm2
            write(69,*) deltavz2-vzcm2
            stop
         endif
         write(69,*)'changing v by a little from circular value',
     $        'to match input.bs'
         
         do i=1,n1
            vx(i)=vx(i)+deltavx1-vxcm1
            vy(i)=vy(i)+deltavy1-vycm1
            vz(i)=vz(i)+deltavz1-vzcm1
         enddo
         do i=n1+1,n1+n2
            vx(i)=vx(i)+deltavx2-vxcm2
            vy(i)=vy(i)+deltavy2-vycm2
            vz(i)=vz(i)+deltavz2-vzcm2
         enddo
      endif

      amtot=am1+am2
      xcm=(xcm1*am1+xcm2*am2)/amtot
      ycm=(ycm1*am1+ycm2*am2)/amtot
      zcm=(zcm1*am1+zcm2*am2)/amtot
      vxcm=(vxcm1*am1+vxcm2*am2)/amtot
      vycm=(vycm1*am1+vycm2*am2)/amtot
      vzcm=(vzcm1*am1+vzcm2*am2)/amtot

c      write(69,*) 'begin natural units...'
c      write(69,'(a,3g14.6)') 'center of mass position of binary:',
c     $     xcm*runit/solrad,
c     $     ycm*runit/solrad,
c     $     zcm*runit/solrad
cc      write(69,'(a,3g14.6)') 'center of mass velocity of binary:',
cc     $     vxcm/amtot*(egsol/munit*runit/solrad)**(-0.5d0),
cc     $     vycm/amtot*(egsol/munit*runit/solrad)**(-0.5d0),
cc     $     vzcm/amtot*(egsol/munit*runit/solrad)**(-0.5d0)
c      write(69,*) '...end natural units'
      if(dabs(xcm/xcombin-1.d0).gt.1.d-12 .or.
     $     dabs(ycm/ycombin-1.d0).gt.1.d-12 .or.
     $     dabs(zcm/zcombin-1.d0).gt.1.d-12)then
         write(69,*)'center of mass position of binary is not correct'
         write(69,*) dabs(xcm/xcombin-1.d0)
         write(69,*) dabs(ycm/ycombin-1.d0)
         write(69,*) dabs(zcm/zcombin-1.d0)
         stop
      endif
      if(dabs(vxcm/vxcombin-1.d0).gt.1.d-11 .or.
     $     dabs(vycm/vycombin-1.d0).gt.1.d-11 .or.
     $     dabs(vzcm/vzcombin-1.d0).gt.1.d-11)then
         write(69,*)'center of mass velocity of binary is not correct'
         write(69,*) vxcm/vxcombin-1.d0
         write(69,*) vycm/vycombin-1.d0
         write(69,*) vzcm/vzcombin-1.d0
         stop
      endif

      open(12,file=startfile2,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n3,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tjumpahead,ngrold,
     $     nrelaxold,trelaxold,dtold
c      amass2=n3
      am3=0.d0
      ntot=nbinary+n3
      if (ntot.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      do i=nbinary+1,ntot
         read (12) xold,yold,zold,am(i),hp(i),rho(i),vxold,vyold,
     $        vzold,vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)

c     place velocities at same time as everything else:
         vxold=vxold-vxdot(i)*0.5d0*dtold
         vyold=vyold-vydot(i)*0.5d0*dtold
         vzold=vzold-vzdot(i)*0.5d0*dtold

         am3=am3+am(i)
         x(i)=cos(ps2)*xold+costh2*sin(ps2)*yold+
     &        sin(ps2)*sinth2*zold
         y(i)=-sin(ps2)*xold+costh2*cos(ps2)*yold+
     &        cos(ps2)*sinth2*zold
         z(i)=-sinth2*yold+costh2*zold
         vx(i)=cos(ps2)*vxold+costh2*sin(ps2)*vyold+
     &        sin(ps2)*sinth2*vzold
         vy(i)=-sin(ps2)*vxold+costh2*cos(ps2)*vyold+
     &        cos(ps2)*sinth2*vzold
         vz(i)=-sinth2*vyold+costh2*vzold
         if(hp(i).le.0.d0) then
            corepts=corepts+1
         endif
      enddo
      n=ntot-corepts
      read(12) nchk
      close(12)
      if (nchk.ne.n3) stop 'bps: problem with file'

      write (69,'(a,5i7)')'bps: n1,n2,n3,n,ntot=',n1,n2,n3,n,ntot

      if(dabs(am1*munit/egsol/am1chk-1.d0).gt.1.d-8 .or.
     $     dabs(am2*munit/egsol/am2chk-1.d0).gt.1.d-8 .or.
     $     dabs(am3*munit/egsol/am3chk-1.d0).gt.1.d-8) then
         write(69,*)'mass(es) in sph.start?u does not match with'
         write(69,*)'mass in input.bs'
         write(69,*)am1*munit/egsol,'should be',am1chk
         write(69,*)am2*munit/egsol,'should be',am2chk
         write(69,*)am3*munit/egsol,'should be',am3chk
         if(.not.tweakv .and.
     $        dabs((am1+am2)*munit/egsol/(am1chk+am2chk)-1.d0).gt.1.d-8)
     $        then
            stop
         else
            write(69,*)'but presumably this is ok in this case'
         endif
      endif

      write(69,*) 'in code units ...'
      write(69,'(3(a,g14.6))') '  binary starts with: x=',
     $     xcombin,
     $     ',y=',ycombin,',z=',zcombin
      write(69,'(3(a,g13.6))') '                     vx=',
     $     vxcombin,
     $     ',vy=',vycombin,
     $     ',vz=',vzcombin
      write(69,'(3(a,g14.6))') 'intruder starts with: x=',
     $     deltax3,
     $     ',y=',deltay3,',z=',deltaz3
      write(69,'(3(a,g13.6))') '                     vx=',
     $     deltavx3,
     $     ',vy=',deltavy3,
     $     ',vz=',deltavz3

      write(69,*) 'now in natural units (g=msun=rsun=1)...'
      write(69,'(3(a,g14.6))') '  binary starts with: x=',
     $     xcombin*runit/solrad,
     $     ',y=',ycombin*runit/solrad,',z=',zcombin*runit/solrad
      write(69,'(3(a,g13.6))') '                     vx=',
     $     vxcombin*(egsol/munit*runit/solrad)**(-0.5d0),
     $     ',vy=',vycombin*(egsol/munit*runit/solrad)**(-0.5d0),
     $     ',vz=',vzcombin*(egsol/munit*runit/solrad)**(-0.5d0)
      write(69,'(3(a,g14.6))') 'intruder starts with: x=',
     $     deltax3*runit/solrad,
     $     ',y=',deltay3*runit/solrad,',z=',deltaz3*runit/solrad
      write(69,'(3(a,g13.6))') '                     vx=',
     $     deltavx3*(egsol/munit*runit/solrad)**(-0.5d0),
     $     ',vy=',deltavy3*(egsol/munit*runit/solrad)**(-0.5d0),
     $     ',vz=',deltavz3*(egsol/munit*runit/solrad)**(-0.5d0)

      do i=nbinary+1,ntot
         x(i)=x(i)+deltax3
         y(i)=y(i)+deltay3
         z(i)=z(i)+deltaz3
         vx(i)=vx(i)+deltavx3
         vy(i)=vy(i)+deltavy3
         vz(i)=vz(i)+deltavz3
      enddo

      xcm1=0.d0
      ycm1=0.d0
      zcm1=0.d0
      vxcm1=0.d0
      vycm1=0.d0
      vzcm1=0.d0
      am1=0.d0
      xcm2=0.d0
      ycm2=0.d0
      zcm2=0.d0
      vxcm2=0.d0
      vycm2=0.d0
      vzcm2=0.d0
      am2=0.d0
      xcm3=0.d0
      ycm3=0.d0
      zcm3=0.d0
      vxcm3=0.d0
      vycm3=0.d0
      vzcm3=0.d0
      am3=0.d0
      do i=1,n1
         xcm1=xcm1+am(i)*x(i)
         ycm1=ycm1+am(i)*y(i)
         zcm1=zcm1+am(i)*z(i)
         vxcm1=vxcm1+am(i)*vx(i)
         vycm1=vycm1+am(i)*vy(i)
         vzcm1=vzcm1+am(i)*vz(i)
         am1=am1+am(i)
      enddo
      do i=n1+1,n1+n2
         xcm2=xcm2+am(i)*x(i)
         ycm2=ycm2+am(i)*y(i)
         zcm2=zcm2+am(i)*z(i)
         vxcm2=vxcm2+am(i)*vx(i)
         vycm2=vycm2+am(i)*vy(i)
         vzcm2=vzcm2+am(i)*vz(i)
         am2=am2+am(i)
      enddo
      do i=n1+n2+1,ntot
         xcm3=xcm3+am(i)*x(i)
         ycm3=ycm3+am(i)*y(i)
         zcm3=zcm3+am(i)*z(i)
         vxcm3=vxcm3+am(i)*vx(i)
         vycm3=vycm3+am(i)*vy(i)
         vzcm3=vzcm3+am(i)*vz(i)
         am3=am3+am(i)
      enddo

      xcm1=xcm1*runit/solrad/am1
      ycm1=ycm1*runit/solrad/am1
      zcm1=zcm1*runit/solrad/am1
      xcm2=xcm2*runit/solrad/am2
      ycm2=ycm2*runit/solrad/am2
      zcm2=zcm2*runit/solrad/am2
      xcm3=xcm3*runit/solrad/am3
      ycm3=ycm3*runit/solrad/am3
      zcm3=zcm3*runit/solrad/am3
      vxcm1=vxcm1*(egsol/munit*runit/solrad)**(-0.5d0)/am1
      vycm1=vycm1*(egsol/munit*runit/solrad)**(-0.5d0)/am1
      vzcm1=vzcm1*(egsol/munit*runit/solrad)**(-0.5d0)/am1
      vxcm2=vxcm2*(egsol/munit*runit/solrad)**(-0.5d0)/am2
      vycm2=vycm2*(egsol/munit*runit/solrad)**(-0.5d0)/am2
      vzcm2=vzcm2*(egsol/munit*runit/solrad)**(-0.5d0)/am2
      vxcm3=vxcm3*(egsol/munit*runit/solrad)**(-0.5d0)/am3
      vycm3=vycm3*(egsol/munit*runit/solrad)**(-0.5d0)/am3
      vzcm3=vzcm3*(egsol/munit*runit/solrad)**(-0.5d0)/am3

      am1=am1*munit/egsol
      am2=am2*munit/egsol
      am3=am3*munit/egsol
      amtot=am1+am2+am3
      xcm=(xcm1*am1+xcm2*am2+xcm3*am3)/amtot
      ycm=(ycm1*am1+ycm2*am2+ycm3*am3)/amtot
      zcm=(zcm1*am1+zcm2*am2+zcm3*am3)/amtot
      vxcm=(vxcm1*am1+vxcm2*am2+vxcm3*am3)/amtot
      vycm=(vycm1*am1+vycm2*am2+vycm3*am3)/amtot
      vzcm=(vzcm1*am1+vzcm2*am2+vzcm3*am3)/amtot

      write(69,'(a,3g14.6)') 'binary:'
      write(69,*) am1
      write(69,*) xcm1,ycm1,zcm1
      write(69,*) vxcm1,vycm1,vzcm1
      write(69,*) am2
      write(69,*) xcm2,ycm2,zcm2
      write(69,*) vxcm2,vycm2,vzcm2
      write(69,*)'single:'
      write(69,*) am3
      write(69,*) xcm3,ycm3,zcm3
      write(69,*) vxcm3,vycm3,vzcm3

      bincomx=(am1*xcm1+am2*xcm2)/(am1+am2)
      bincomy=(am1*ycm1+am2*ycm2)/(am1+am2)
      bincomz=(am1*zcm1+am2*zcm2)/(am1+am2)
      bincomvx=(am1*vxcm1+am2*vxcm2)/(am1+am2)
      bincomvy=(am1*vycm1+am2*vycm2)/(am1+am2)
      bincomvz=(am1*vzcm1+am2*vzcm2)/(am1+am2)

      deltaxbin=xcm2-xcm1
      deltaybin=ycm2-ycm1
      deltazbin=zcm2-zcm1
      deltarbin=(deltaxbin**2+deltaybin**2+deltazbin**2)**0.5d0
      write(69,'(a,4g14.6)')'actual separation hat of binary=',
     $     deltaxbin/deltarbin,deltaybin/deltarbin,deltazbin/deltarbin,
     $     deltarbin

      deltavxbin=vxcm2-vxcm1
      deltavybin=vycm2-vycm1
      deltavzbin=vzcm2-vzcm1

      write(69,'(a,3g14.6)') 'binary c.o.m. position:',
     $     bincomx,bincomy,bincomz
      write(69,'(a,3g14.6)') 'binary c.o.m. velocity:',
     $     bincomvx,bincomvy,bincomvz

      lx=deltaybin*deltavzbin-deltazbin*deltavybin
      ly=deltazbin*deltavxbin-deltaxbin*deltavzbin
      lz=deltaxbin*deltavybin-deltaybin*deltavxbin
      ltotbin=(lx**2+ly**2+lz**2)**0.5d0
      write(69,'(a,4g14.6)')'actual angular momentum hat of binary=',
     $     lx/ltotbin,ly/ltotbin,lz/ltotbin,ltotbin

      write(69,'(a,3g14.6)') 'overall c.o.m. position:',xcm,ycm,zcm
      write(69,'(a,3g14.6)') 'overall c.o.m. velocity:',vxcm,vycm,vzcm

      write(69,*)'...end natural units'

      if(xcm**2+ycm**2+zcm**2.gt.1.d-12) then
         write(69,*)'center of mass is too far from origin',
     $        xcm**2+ycm**2+zcm**2
         stop
      endif
      if(vxcm**2+vycm**2+vzcm**2.gt.1.d-15) then
         write(69,*)'center of mass velocity is too large',
     $        vxcm**2+vycm**2+vzcm**2
         stop
      endif

      call stride_setup
      
c     prepare leap-frog scheme for first iteration:
      call lfstart
      write (69,*) 'bps:          ... done'
      
      return
      end
