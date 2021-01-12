      subroutine smbh
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
      real*8 au
      parameter(au=14959787066000.d0)
      real*8 displacex, displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace
      
      corepts=0
      if(myrank.eq.0) write (69,*) 'smbh: reading start files ...'

      open(12,file=startfile1,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n1,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold
      am1=0.d0
      if(myrank.eq.0) write (69,*) 'n1=',n1
      amass1=n1
      if (n1.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      do i=1,n1
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i)
c     place velocities at same time as everything else:
         vx(i)=vx(i)-vxdot(i)*0.5d0*dtold
         vy(i)=vy(i)-vydot(i)*0.5d0*dtold
         vz(i)=vz(i)-vzdot(i)*0.5d0*dtold
         am1=am1+am(i)
         if(hp(i).le.0.d0) then
            if(myrank.eq.0) write (69,*) 'star 1 has a core point ...'
c            write(69,*)'only star three can have a core point'
c            stop
         endif
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n1) stop 'smbh: problem with file'

      open(12,file=startfile2,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n2,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tjumpahead,ngrold,
     $     nrelaxold,trelaxold,dtold
      if(myrank.eq.0) write(69,*)'n2=',n2
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
     $        grpot(i),meanmolecular(i),
     $        cc(i)
c     place velocities at same time as everything else:
         vx(i)=vx(i)-vxdot(i)*0.5d0*dtold
         vy(i)=vy(i)-vydot(i)*0.5d0*dtold
         vz(i)=vz(i)-vzdot(i)*0.5d0*dtold
         am2=am2+am(i)
         if(hp(i).le.0.d0) then
            if(myrank.eq.0) write (69,*) 'star 2 has a core point ...'
c            write(69,*)'only star three can have a core point'
c            stop
         endif
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n2) stop 'smbh: problem with file'

      open(30,file='input.3s')
      read(30,*) dummy
      read(30,*) am3chk
      read(30,*) deltax3,deltay3,deltaz3
      read(30,*) deltavx3,deltavy3,deltavz3
      read(30,*) am1chk
      read(30,*) deltax1,deltay1,deltaz1
      read(30,*) deltavx1,deltavy1,deltavz1
      read(30,*) am2chk
      read(30,*) deltax2,deltay2,deltaz2
      read(30,*) deltavx2,deltavy2,deltavz2
      close(30)

      am1chk=3.d0*am1chk
      am2chk=3.d0*am2chk
      am3chk=3.d0*am3chk

      deltax1=deltax1*0.01d0*au/runit
      deltay1=deltay1*0.01d0*au/runit
      deltaz1=deltaz1*0.01d0*au/runit
      deltax2=deltax2*0.01d0*au/runit
      deltay2=deltay2*0.01d0*au/runit
      deltaz2=deltaz2*0.01d0*au/runit
      deltax3=deltax3*0.01d0*au/runit
      deltay3=deltay3*0.01d0*au/runit
      deltaz3=deltaz3*0.01d0*au/runit

      deltavx1=deltavx1*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavy1=deltavy1*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavz1=deltavz1*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavx2=deltavx2*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavy2=deltavy2*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavz2=deltavz2*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavx3=deltavx3*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavy3=deltavy3*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0
      deltavz3=deltavz3*(3.d0*egsol/munit*runit/(0.01d0*au))**0.5d0

      n3=1
      if(myrank.eq.0) write(69,*)'n3=',n3
      am3=am3chk
      ntot=n1+n2+n3
      if (ntot.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif
      x(ntot)=deltax3
      y(ntot)=deltay3
      z(ntot)=deltaz3
      vx(ntot)=deltavx3
      vy(ntot)=deltavy3
      vz(ntot)=deltavz3
      am(ntot)=am3
      hp(ntot)=hco
      u(ntot)=0.d0
      meanmolecular(ntot)=0.d0
      cc(ntot)=1
      corepts=0
      n=ntot-corepts
      if(myrank.eq.0) then
         write (69,*)'smbh: n=',n,'ntot=',ntot
         
         write(69,*)'star 1 mass [msun]=',am1*munit/egsol,am1chk
         write(69,*)'star 2 mass [msun]=',am2*munit/egsol,am2chk
         write(69,*)'star 3 mass [msun]=',am3*munit/egsol,am3chk
      endif

      if(dabs(am1*munit/egsol/am1chk-1.d0).gt.1.d-5 .or.
     $     dabs(am2*munit/egsol/am2chk-1.d0).gt.1.d-5 .or.
     $     dabs(am3*munit/egsol/am3chk-1.d0).gt.1.d-5) then
         write(69,*)'mass(es) in sph.start?u does not match with'
         write(69,*)'mass in input.3s'
         stop
      endif

      if(myrank.eq.0) then
         write(69,'(3(a,g14.6))') 'star 1 starts with: x=',deltax1,
     $        ',y=',deltay1,',z=',deltaz1
         write(69,'(3(a,g13.6))') '                    vx=',deltavx1,
     $        ',vy=',deltavy1,',vz=',deltavz1
         write(69,'(3(a,g14.6))') 'star 2 starts with: x=',deltax2,
     $        ',y=',deltay2,',z=',deltaz2
         write(69,'(3(a,g13.6))') '                    vx=',deltavx2,
     $        ',vy=',deltavy2,',vz=',deltavz2
         write(69,'(3(a,g14.6))') 'star 3 starts with: x=',deltax3,
     $        ',y=',deltay3,',z=',deltaz3
         write(69,'(3(a,g13.6))') '                    vx=',deltavx3,
     $        ',vy=',deltavy3,',vz=',deltavz3
      endif
      xcm=(am1*deltax1+am2*deltax2+am3*deltax3)/(am1+am2+am3)
      ycm=(am1*deltay1+am2*deltay2+am3*deltay3)/(am1+am2+am3)
      zcm=(am1*deltaz1+am2*deltaz2+am3*deltaz3)/(am1+am2+am3)
      vxcm=(am1*deltavx1+am2*deltavx2+am3*deltavx3)/(am1+am2+am3)
      vycm=(am1*deltavy1+am2*deltavy2+am3*deltavy3)/(am1+am2+am3)
      vzcm=(am1*deltavz1+am2*deltavz2+am3*deltavz3)/(am1+am2+am3)
c     Find the center of mass position of just the two stars and not the black hole:
      ndisplace=1
      displacex=(am1*deltax1+am2*deltax2)/(am1+am2)
      displacey=(am1*deltay1+am2*deltay2)/(am1+am2)
      displacez=(am1*deltaz1+am2*deltaz2)/(am1+am2)

      if(myrank.eq.0) then
         write(69,'(a,3g13.6)')'center of mass position:',xcm,ycm,zcm
         write(69,'(a,3g13.6)')'center of mass velocity:',vxcm,vycm,vzcm
         write(69,*)'n.b.:subtracting off these center of mass values!'
         write(69,*) 'Shifting origin by displace{x,y,z}=',
     $        displacex,displacey,displacez
      endif

      do i=1,n1
         x(i)=x(i)+(deltax1-xcm-displacex)
         y(i)=y(i)+(deltay1-ycm-displacey)
         z(i)=z(i)+(deltaz1-zcm-displacez)
         vx(i)=vx(i)+deltavx1-vxcm
         vy(i)=vy(i)+deltavy1-vycm
         vz(i)=vz(i)+deltavz1-vzcm
      enddo
      do i=n1+1,n1+n2
         x(i)=x(i)+(deltax2-xcm-displacex)
         y(i)=y(i)+(deltay2-ycm-displacey)
         z(i)=z(i)+(deltaz2-zcm-displacez)
         vx(i)=vx(i)+deltavx2-vxcm
         vy(i)=vy(i)+deltavy2-vycm
         vz(i)=vz(i)+deltavz2-vzcm
      enddo
      do i=n1+n2+1,ntot
         x(i)=x(i)+(deltax3-xcm-displacex)
         y(i)=y(i)+(deltay3-ycm-displacey)
         z(i)=z(i)+(deltaz3-zcm-displacez)
         vx(i)=vx(i)+deltavx3-vxcm
         vy(i)=vy(i)+deltavy3-vycm
         vz(i)=vz(i)+deltavz3-vzcm
      enddo
      call stride_setup
      
c     prepare leap-frog scheme for first iteration:
      call lfstart
      if(myrank.eq.0) write (69,*) 'smbh:          ... done'
      
      return
      end
