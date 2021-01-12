      subroutine triple
      include 'starsmasher.h'
      real*8 vxcm,vycm,vzcm,xcm,ycm,zcm
      integer n2,i,nchk,corepts,n3
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      real*8 deltax1,deltay1,deltaz1,deltavx1,deltavy1,deltavz1
      real*8 deltax2,deltay2,deltaz2,deltavx2,deltavy2,deltavz2
      real*8 deltax3,deltay3,deltaz3,deltavx3,deltavy3,deltavz3
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2,am3
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      character*7 dummy
      real*8 am1chk,am2chk,am3chk
      real*8 egsol,solrad
      parameter(egsol=1.9891d+33,solrad=6.9599d10)
      
      corepts=0
      write (69,*) 'triple: reading start files ...'

      open(12,file=startfile1,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n1,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold
      am1=0.d0
      write(69,*)'n1=',n1
      amass1=n1
      if (n1.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      do i=1,n1
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)
c     place velocities at same time as everything else:
         vx(i)=vx(i)-vxdot(i)*0.5d0*dtold
         vy(i)=vy(i)-vydot(i)*0.5d0*dtold
         vz(i)=vz(i)-vzdot(i)*0.5d0*dtold
         am1=am1+am(i)
         if(hp(i).le.0.d0) then
            write(69,*)'only star three can have a core point'
            stop
         endif
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n1) stop 'triple: problem with file'

      open(12,file=startfile2,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n2,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tjumpahead,ngrold,
     $     nrelaxold,trelaxold,dtold
      write(69,*)'n2=',n2
      amass2=n2
      am2=0.d0
      ntot=n1+n2
      if (ntot.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      do i=n1+1,ntot
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)
c     place velocities at same time as everything else:
         vx(i)=vx(i)-vxdot(i)*0.5d0*dtold
         vy(i)=vy(i)-vydot(i)*0.5d0*dtold
         vz(i)=vz(i)-vzdot(i)*0.5d0*dtold
         am2=am2+am(i)
         if(hp(i).le.0.d0) then
            write(69,*)'only star three can have a core point'
            stop
         endif
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n2) stop 'triple: problem with file'

      open(12,file='sph.start3u',form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n3,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tjumpahead,ngrold,
     $     nrelaxold,trelaxold,dtold
      write(69,*)'n3=',n3
      am3=0.d0
      ntot=n1+n2+n3
      if (ntot.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      do i=n1+n2+1,ntot
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)
c     place velocities at same time as everything else:
         vx(i)=vx(i)-vxdot(i)*0.5d0*dtold
         vy(i)=vy(i)-vydot(i)*0.5d0*dtold
         vz(i)=vz(i)-vzdot(i)*0.5d0*dtold
         am3=am3+am(i)
         if(hp(i).le.0.d0) then
            corepts=corepts+1
         endif
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n3) stop 'triple: problem with file'

      n=ntot-corepts
      write (69,*)'triple: n=',n,'ntot=',ntot

      open(30,file='input.3s')
      read(30,*) dummy
      read(30,*) am1chk
      read(30,*) deltax1,deltay1,deltaz1
      read(30,*) deltavx1,deltavy1,deltavz1
      read(30,*) am2chk
      read(30,*) deltax2,deltay2,deltaz2
      read(30,*) deltavx2,deltavy2,deltavz2
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

      write(69,*)'star 1 mass [msun]=',am1*munit/egsol,am1chk
      write(69,*)'star 2 mass [msun]=',am2*munit/egsol,am2chk
      write(69,*)'star 3 mass [msun]=',am3*munit/egsol,am3chk

      if(dabs(am1*munit/egsol/am1chk-1.d0).gt.1.d-8 .or.
     $     dabs(am2*munit/egsol/am2chk-1.d0).gt.1.d-8 .or.
     $     dabs(am3*munit/egsol/am3chk-1.d0).gt.1.d-8) then
         write(69,*)'mass(es) in sph.start?u does not match with'
         write(69,*)'mass in input.3s'
         stop
      endif

      write(69,'(3(a,g14.6))') 'star 1 starts with: x=',deltax1,
     $     ',y=',deltay1,',z=',deltaz1
      write(69,'(3(a,g13.6))') '                    vx=',deltavx1,
     $     ',vy=',deltavy1,',vz=',deltavz1
      write(69,'(3(a,g14.6))') 'star 2 starts with: x=',deltax2,
     $     ',y=',deltay2,',z=',deltaz2
      write(69,'(3(a,g13.6))') '                    vx=',deltavx2,
     $     ',vy=',deltavy2,',vz=',deltavz2
      write(69,'(3(a,g14.6))') 'star 3 starts with: x=',deltax3,
     $     ',y=',deltay3,',z=',deltaz3
      write(69,'(3(a,g13.6))') '                    vx=',deltavx3,
     $     ',vy=',deltavy3,',vz=',deltavz3
      xcm=(am1*deltax1+am2*deltax2+am3*deltax3)/(am1+am2+am3)
      ycm=(am1*deltay1+am2*deltay2+am3*deltay3)/(am1+am2+am3)
      zcm=(am1*deltaz1+am2*deltaz2+am3*deltaz3)/(am1+am2+am3)
      vxcm=(am1*deltavx1+am2*deltavx2+am3*deltavx3)/(am1+am2+am3)
      vycm=(am1*deltavy1+am2*deltavy2+am3*deltavy3)/(am1+am2+am3)
      vzcm=(am1*deltavz1+am2*deltavz2+am3*deltavz3)/(am1+am2+am3)
      write(69,'(a,3g13.6)')'center of mass position:',xcm,ycm,zcm
      write(69,'(a,3g13.6)')'center of mass velocity:',vxcm,vycm,vzcm

      write(69,*)'n.b.:subtracting off these center of mass values!'

      do i=1,n1
         x(i)=x(i)+deltax1-xcm
         y(i)=y(i)+deltay1-ycm
         z(i)=z(i)+deltaz1-zcm
         vx(i)=vx(i)+deltavx1-vxcm
         vy(i)=vy(i)+deltavy1-vycm
         vz(i)=vz(i)+deltavz1-vzcm
      enddo
      do i=n1+1,n1+n2
         x(i)=x(i)+deltax2-xcm
         y(i)=y(i)+deltay2-ycm
         z(i)=z(i)+deltaz2-zcm
         vx(i)=vx(i)+deltavx2-vxcm
         vy(i)=vy(i)+deltavy2-vycm
         vz(i)=vz(i)+deltavz2-vzcm
      enddo
      do i=n1+n2+1,ntot
         x(i)=x(i)+deltax3-xcm
         y(i)=y(i)+deltay3-ycm
         z(i)=z(i)+deltaz3-zcm
         vx(i)=vx(i)+deltavx3-vxcm
         vy(i)=vy(i)+deltavy3-vycm
         vz(i)=vz(i)+deltavz3-vzcm
      enddo
      call stride_setup
      
c     prepare leap-frog scheme for first iteration:
      call lfstart
      write (69,*) 'triple:          ... done'
      
      return
      end
