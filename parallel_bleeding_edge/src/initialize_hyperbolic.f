      subroutine hyperbolic
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 k,rdot,rdotcheck,costheta,sintheta,ltot,thetadot,
     $     semilatusrectum,mu,e0check,sinthetacheck,
     $     eorb0check
      real*8 altotint,r,amu,ak,vxcm,vycm,xcm,ycm
      integer n2,i,nchk,corepts
      integer nnoptold,nnoptold2,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      logical twofiles
      real*8 deltax1,deltay1,deltavx1,deltavy1
      real*8 deltax2,deltay2,deltavx2,deltavy2
      real*8 vinf,semilatusrectumprime,theta
      integer idumb
      real*8 costh1,sinth1,costh2,sinth2,ps1,ps2
      real*8 xold,yold,zold,vxold,vyold,vzold,ran1
      real*8 xnew,ynew,znew,vxnew,vynew,vznew
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 eorb0,xx
      real*8 x1,y1,z1,vx1,vy1,vz1,amone
      real*8 x2,y2,z2,vx2,vy2,vz2,amtwo
      real*8 x3,y3,z3,vx3,vy3,vz3,amthree,am4
      integer icomp(nmax)
      common/compbettercom3/amone,x1,y1,z1,vx1,vy1,vz1,
     $     amtwo,x2,y2,z2,vx2,vy2,vz2,
     $     amthree,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp
      real*8 divv(nmax)
      common/commdivv/divv
      common/orbitalelements/e0,semimajoraxis,impactparameter,vinf2
      real*8 bbh_m1,bbh_m2,bbh_rp,bbh_semimajoraxis,bbh_vinf2,
     $     bbh_e0,bbh_trueanomaly,bbh_argperi,bbh_inclination,bbh_longitude
      common/bbhinfo/bbh_m1,bbh_m2,bbh_rp,bbh_semimajoraxis,bbh_vinf2,
     $     bbh_e0,bbh_trueanomaly,bbh_argperi,bbh_inclination,bbh_longitude
      real*8 xcmparent,ycmparent,zcmparent,
     $     vxcmparent,vycmparent,vzcmparent
      real*8 hprms,bbh_sep0
      logical passivelyAdvected

      nrelax=0
      
      call cpu_time(time1)

      corepts=0
      if(myrank.eq.0) write (69,*)
     $     'hyperbolic: reading startfile1 ',
     $     trim(startfile1)

      open(12,file=startfile1,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
      read(12) n1,nnoptold,hcoold,hfloorold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tjumpahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold,omega2

      if(nnoptold.ne.nnopt) then
         if(myrank.eq.0) then
            write(69,*) 'NOTE: Currently nnopt=',nnopt
            write(69,*) '      The NNOPT in startfile1=',nnoptold
            write(69,*) '      Changing NNOPT to be',nnoptold
         endif
         nnopt=nnoptold
      endif

      am1=0.d0
      if(myrank.eq.0) write(69,*)'n1=',n1
      amass1=n1

      inquire(file=startfile2,exist=twofiles)
c      if(twofiles) then
      if(.false.) then
         if(.true.)then
            seconds = dabs(vinf2)/7.d0+1.d-7
         else
            call cpu_time(time2)
            seconds = time2-time1
         endif
         if(myrank.eq.0) write (69,'(a,g16.9,a)')
     $        'hyperbolic: ',seconds,' will help set random number'
         idumb=-nint(10000*(100000*seconds-int(100000*seconds)))
         if(myrank.eq.0) write(69,*)'idumb=',idumb
         
         costh1=2.d0*(ran1(idumb)-0.5d0)
         sinth1=sqrt(1.d0-costh1**2.d0)
         ps1=2.d0*pi*(ran1(idumb)-0.5d0)
         if(myrank.eq.0) write(69,*)
     $        'rotation angles for star 1 (in radians):'
         if(myrank.eq.0) write(69,*)'theta1=',acos(costh1),' psi1=',ps1
      else
         costh1=1.d0
         sinth1=0.d0
         ps1=0d0
      endif

      xcmparent=0d0
      ycmparent=0d0
      zcmparent=0d0
      vxcmparent=0d0
      vycmparent=0d0
      vzcmparent=0d0
      do i=1,n1
         icomp(i)=2
         read (12) xold,yold,zold,am(i),hp(i),rho(i),vxold,vyold,
     $        vzold,vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i),divv(i)
         ueq(i)=0
         tthermal(i)=1d30

c     place velocities at same time as everything else:
         vxold=vxold-vxdot(i)*0.5d0*dtold
         vyold=vyold-vydot(i)*0.5d0*dtold
         vzold=vzold-vzdot(i)*0.5d0*dtold

         am1=am1+am(i)
         x(i)=cos(ps1)*xold+costh1*sin(ps1)*yold+
     &        sin(ps1)*sinth1*zold
         y(i)=-sin(ps1)*xold+costh1*cos(ps1)*yold+
     &        cos(ps1)*sinth1*zold
         z(i)=-sinth1*yold+costh1*zold
         vx(i)=cos(ps1)*vxold+costh1*sin(ps1)*vyold+
     &        sin(ps1)*sinth1*vzold
         vy(i)=-sin(ps1)*vxold+costh1*cos(ps1)*vyold+
     &        cos(ps1)*sinth1*vzold
         vz(i)=-sinth1*vyold+costh1*vzold
         xcmparent=xcmparent+am(i)*x(i)
         ycmparent=ycmparent+am(i)*y(i)
         zcmparent=zcmparent+am(i)*z(i)
         vxcmparent=vxcmparent+am(i)*vx(i)
         vycmparent=vycmparent+am(i)*vy(i)
         vzcmparent=vzcmparent+am(i)*vz(i)
      enddo
      read(12) nchk
      close(12)
      if (nchk.ne.n1) stop 'hyperbolic: problem with sph.start1u file'
      xcmparent=xcmparent/am1
      ycmparent=ycmparent/am1
      zcmparent=zcmparent/am1
      vxcmparent=vxcmparent/am1
      vycmparent=vycmparent/am1
      vzcmparent=vzcmparent/am1
      amtwo=am1
      if(myrank.eq.0) then
         write(69,*)'Center of mass position & velocity of star 1'
         write(69,*)'(to be subtracted off):'
         write(69,*) xcmparent,ycmparent,zcmparent
         write(69,*) vxcmparent,vycmparent,vzcmparent
      endif
      do i=1,n1
         x(i)=x(i)-xcmparent
         y(i)=y(i)-ycmparent
         z(i)=z(i)-zcmparent
         vx(i)=vx(i)-vxcmparent
         vy(i)=vy(i)-vycmparent
         vz(i)=vz(i)-vzcmparent
      enddo
      passivelyAdvected = .false.

      if(twofiles) then
         if(myrank.eq.0) write (69,*)
     $        'hyperbolic: reading startfile2 ',
     $        trim(startfile2)
         if(.false.) then
            costh2=2.d0*(ran1(idumb)-0.5d0)
            sinth2=sqrt(1.d0-costh2**2.d0)
            ps2=2.d0*pi*(ran1(idumb)-0.5d0)      
            if(myrank.eq.0) write(69,*)
     $           'rotation angles for star 2 (in radians):'
            if(myrank.eq.0) write(69,*)'theta2=',acos(costh2),' psi2=',ps2
         else
            costh2=1.d0
            sinth2=0.d0
            ps2=0.d0
         endif
         open(12,file=startfile2,form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
         read(12) n2,nnoptold2,hcoold,hfloorold,sep0old,
     $        tfold,dtoutold,noutold,nitold,told,navold,
     $        alphaold,betaold,tjumpahead,ngrold,
     $        nrelaxold,trelaxold,dtold,omega2
         if(nnoptold2.ne.nnoptold) then
            write(69,*) 'ERROR: nnopt from star 1=',nnoptold
            write(69,*) '       nnopt from star 2=',nnoptold2
            stop
         endif
         if(nnoptold2.ne.nnopt) then
            write(69,*) 'ERROR: currently nnopt=',nnopt
            write(69,*) '       the NNOPT in startfile2=',nnoptold
            stop
         endif

         amass2=n2
         am2=0.d0
         ntot=n1+n2
         if (ntot.gt.nmax) then
            if(myrank.eq.0) write(69,*)'must increase nmax...'
            stop
         endif
         xcmparent=0d0
         ycmparent=0d0
         zcmparent=0d0
         vxcmparent=0d0
         vycmparent=0d0
         vzcmparent=0d0
         do i=n1+1,ntot
            icomp(i)=3
            read (12) xold,yold,zold,am(i),hp(i),rho(i),vxold,vyold,
     $           vzold,vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $           grpot(i),meanmolecular(i),
     $           cc(i),divv(i)
            ueq(i)=0
            tthermal(i)=1d30
            
c     place velocities at same time as everything else:
            vxold=vxold-vxdot(i)*0.5d0*dtold
            vyold=vyold-vydot(i)*0.5d0*dtold
            vzold=vzold-vzdot(i)*0.5d0*dtold
            
            am2=am2+am(i)
            x(i)=cos(ps2)*xold+costh2*sin(ps2)*yold+
     &           sin(ps2)*sinth2*zold
            y(i)=-sin(ps2)*xold+costh2*cos(ps2)*yold+
     &           cos(ps2)*sinth2*zold
            z(i)=-sinth2*yold+costh2*zold
            vx(i)=cos(ps2)*vxold+costh2*sin(ps2)*vyold+
     &           sin(ps2)*sinth2*vzold
            vy(i)=-sin(ps2)*vxold+costh2*cos(ps2)*vyold+
     &           cos(ps2)*sinth2*vzold
            vz(i)=-sinth2*vyold+costh2*vzold
            if(hp(i).le.0.d0) then
               corepts=corepts+1
            endif
            xcmparent=xcmparent+am(i)*x(i)
            ycmparent=ycmparent+am(i)*y(i)
            zcmparent=zcmparent+am(i)*z(i)
            vxcmparent=vxcmparent+am(i)*vx(i)
            vycmparent=vycmparent+am(i)*vy(i)
            vzcmparent=vzcmparent+am(i)*vz(i)
         enddo
         n=ntot-corepts
         read(12) nchk
         close(12)
         if (nchk.ne.n2) stop
     $        'hyperbolic: problem with sph.start2u file'
         xcmparent=xcmparent/am2
         ycmparent=ycmparent/am2
         zcmparent=zcmparent/am2
         vxcmparent=vxcmparent/am2
         vycmparent=vycmparent/am2
         vzcmparent=vzcmparent/am2
         if(myrank.eq.0) then
            write(69,*)'Center of mass position & velocity of star 2'
            write(69,*)'(to be subtracted off):'
            write(69,*) xcmparent,ycmparent,zcmparent
            write(69,*) vxcmparent,vycmparent,vzcmparent
         endif
         do i=n1+1,ntot
            x(i)=x(i)-xcmparent
            y(i)=y(i)-ycmparent
            z(i)=z(i)-zcmparent
            vx(i)=vx(i)-vxcmparent
            vy(i)=vy(i)-vycmparent
            vz(i)=vz(i)-vzcmparent
         enddo

         if(passivelyAdvected) then
            if(myrank.eq.0) write(69,*)'Reading in aa and bb values...',ntot
            open(99,file=advectedfile,status='old')
            do i=1,ntot
               read(99,*) aa(i),bb(i),dd(i)
            enddo
            close(99)
         endif

      else
c     There is only one start file, so the second object will be a single or
c     binary compact object
         if(passivelyAdvected) then
            if(myrank.eq.0) write(69,*)'Reading in aa and bb values...',n1
            open(99,file=advectedfile,status='old')
            do i=1,n1
               read(99,*) aa(i),bb(i),dd(i)
            enddo
            close(99)
         endif
            
         hprms=0d0
         do i=1,n1
            hprms=hprms+hp(i)**2
         enddo
         hprms=sqrt(hprms/n1)
         if(myrank.eq.0) then
            write(69,*) 'rms smoothing length h in star=',hprms
         endif
         
         if(myrank.eq.0) write (69,*)
     $        'startfile2 ',
     $        trim(startfile2),
     $        ' does not exist: will use'
         if(bbh_m2.le.0.d0) then
c     Second object will be a single compact object of mass mbh.  A binary is
c     asked for by giving bbh_m2 a positive mass in sph.input.  The default is
c     negative, so reaching the binary branch means the user set it deliberately.
            if(myrank.eq.0) write (69,*)
     $           'compact object particle instead'
            
            n2=1
            ntot=n1+n2
            i=ntot
            am2=mbh
            am(i)=am2
            x(i)=0.d0
            y(i)=0.d0
            z(i)=0.d0
            vx(i)=0.d0
            vy(i)=0.d0
            vz(i)=0.d0
            gx(i)=0.d0
            gy(i)=0.d0
            gz(i)=0.d0
            vxdot(i)=0.d0
            vydot(i)=0.d0
            vzdot(i)=0.d0
            u(i)=0.d0
            udot(i)=0.d0
            if(hco.gt.0.d0) then
               hp(i)=hco
            else
               hp(i)=hprms
            endif
            if(passivelyAdvected) then
               aa(i)=hp(i)
               bb(i)=0
               dd(i)=0
            endif
            cc(i)=0
            if(myrank.eq.0) write(69,*)'compact_object_mass',am(i)
            if(myrank.eq.0) write(69,*)'compact_object_smoothing_length',
     $           hp(i)
         else
c     Second object will be a compact object binary such as a binary black hole.
c     bbh_m2 being positive is what selected this branch, so bbh_m1 must have
c     been set too.  Stopping here is better than proceeding with a primary mass
c     the user never chose.
            if(bbh_m1.le.0.d0) then
               if(myrank.eq.0) then
                  write(69,*)'hyperbolic: bbh_m2 is set, so the second'
                  write(69,*)'hyperbolic: object is a compact-object binary,'
                  write(69,*)'hyperbolic: but bbh_m1 has not been given a'
                  write(69,*)'hyperbolic: positive mass.  Set bbh_m1 in'
                  write(69,*)'hyperbolic: sph.input, or leave bbh_m2 unset'
                  write(69,*)'hyperbolic: for a single point mass of mass mbh.'
               endif
               stop
            endif
            if(myrank.eq.0) write (69,*)
     $           'binary compact object instead'
            n2=2
            ntot=n1+n2
            am(ntot-1)=bbh_m1
            am(ntot)=bbh_m2
            am2=bbh_m1+bbh_m2
            do i=ntot-1,ntot
               x(i)=0.d0
               y(i)=0.d0
               z(i)=0.d0
               vx(i)=0.d0
               vy(i)=0.d0
               vz(i)=0.d0
               gx(i)=0.d0
               gy(i)=0.d0
               gz(i)=0.d0
               vxdot(i)=0.d0
               vydot(i)=0.d0
               vzdot(i)=0.d0
               u(i)=0.d0
               udot(i)=0.d0
               if(hco.gt.0.d0) then
                  hp(i)=hco
               else if (u(1).eq.0) then
                  hp(i)=hp(1) ! If first point is core point, use same softening
               else
                  hp(i)=hprms
               endif
               if(passivelyAdvected) then
                  aa(i)=hp(i)
                  bb(i)=0
                  dd(i)=0
               endif
               cc(i)=0
               if(myrank.eq.0) write(69,*)'compact_object_mass',
     $              am(i),'i=',i
               if(myrank.eq.0) write(69,*)
     $              'compact_object_smoothing_length',hp(i),
     $              'i=',i
            enddo

            k=bbh_m1*bbh_m2
            mu=k/(bbh_m1+bbh_m2)

            if(bbh_rp.lt.0d0 .and. bbh_semimajoraxis.eq.0.d0)then
c     presumably bbh_e0 and bbh_vinf2 have been set in sph.input:
               eorb0=0.5d0*mu*bbh_vinf2
               bbh_semimajoraxis=-0.5d0*k/eorb0
               rp=bbh_semimajoraxis*(1.d0-bbh_e0)
            else if(bbh_vinf2.ge.1d30 .and. bbh_semimajoraxis.eq.0.d0)then
c     presumably bbh_e0 and bbh_rp have been set in sph.input:
               bbh_semimajoraxis=bbh_rp/(1.d0-bbh_e0)
               eorb0=-0.5d0*k/bbh_semimajoraxis
               bbh_vinf2=2.d0*eorb0/mu
            else if(bbh_vinf2.ge.1d30 .and. bbh_e0.lt.0.d0)then
c     presumably bbh_semimajoraxis and bbh_rp have been set in sph.input:
               eorb0=-0.5d0*k/bbh_semimajoraxis
               bbh_vinf2=2.d0*eorb0/mu
               bbh_e0=1.d0-bbh_rp/bbh_semimajoraxis
            else if(bbh_rp.lt.0.d0 .and. bbh_vinf2.ge.1.d30)then
c     presumably bbh_semimajoraxis and bbh_e0 have been set in sph.input:
               bbh_rp=bbh_semimajoraxis*(1.d0-bbh_e0)
               eorb0=-0.5d0*k/bbh_semimajoraxis
               bbh_vinf2=2.d0*eorb0/mu
            else
c     presumably bbh_rp and bbh_vinf2 have been set in sph.input:
               eorb0=0.5d0*mu*bbh_vinf2
               bbh_semimajoraxis=-0.5d0*k/eorb0
               bbh_e0=1.d0-bbh_rp/bbh_semimajoraxis
            endif

            if(myrank.eq.0) then
               write(69,*) 'hyperbolic: bbh_rp=',
     $              bbh_rp,'bbh_v_inf2=',bbh_vinf2,'bbh_semimajoraxis=',
     $              bbh_semimajoraxis
               write(69,*)'hyperbolic: bbh_e_orb=',eorb0
               write(69,*)'hyperbolic: n=',n,'ntot=',ntot
            endif
            
            semilatusrectum=bbh_rp*(1.d0+bbh_e0)
            
c     equation (8.41) of marion and thornton
c            costheta = (semilatusrectum/bbh_sep0-1.d0)/bbh_e0
            theta=bbh_trueanomaly ! true anomaly should be from -pi to +pi
            costheta=cos(theta)
            bbh_sep0=semilatusrectum/(1d0+bbh_e0*costheta)
            
c     equation (8.40) of marion and thornton
            ltot = sqrt(semilatusrectum*mu*k)
c     equation (8.10) of marion and thornton
            thetadot = ltot/mu/bbh_sep0**2
c     equation (8.40) of marion and thornton
            eorb0check = (bbh_e0**2-1.d0)*mu*k*k/(ltot*ltot)*0.5d0
            if(myrank.eq.0) write(69,*)'1st simple bbh check:',
     $           eorb0,eorb0check
            
            sintheta = sin(theta)
            sinthetacheck = sqrt(1.d0-costheta**2)
            if( abs(abs(sintheta)-sinthetacheck ).gt.1.d-13) then
               write(69,*)'2nd simple bbh check fails',sintheta,sinthetacheck
               stop
            endif
            if(myrank.eq.0) write(69,*)'bbh semilatusrectum=',
     $           semilatusrectum,' bbh mu=',mu,
     $           ' bbh_sep0=',bbh_sep0,' bbh_e0=',bbh_e0
            if(myrank.eq.0) write(69,*) 'bbh cos=',costheta,
     $           ' sin=',sintheta
c     the minus signs on the position and velocity component equations
c     have been chosen so that the separation vector r equals r2-r1,
c     *not* r1-r2 as in marion and thornton.  this was done so that the
c     code would give the same initial conditions as twostars.f when the
c     eccentricity is 1
            deltax1 = -bbh_m2/(bbh_m1+bbh_m2)*bbh_sep0*costheta
            deltay1 = -bbh_m2/(bbh_m1+bbh_m2)*bbh_sep0*sintheta
            deltax2 = bbh_m1/(bbh_m1+bbh_m2)*bbh_sep0*costheta
            deltay2 = bbh_m1/(bbh_m1+bbh_m2)*bbh_sep0*sintheta
c     differentiating eqn. 8.41
            rdot=semilatusrectum/
     $           (1.d0+bbh_e0*costheta)**2*bbh_e0*sintheta*thetadot
c     another expr. for rdot (using eqn. 8.14) and compare
            rdotcheck=sqrt((eorb0-0.5d0*ltot**2/mu/bbh_sep0**2+k/bbh_sep0)/
     $           0.5d0/mu)
            if(myrank.eq.0) write(69,*)
     $           'bbh rdot(from semilatusrectum): ',rdot,
     $           '  magnitude of bbh rdot(from e): ',rdotcheck
c            if(rdot.ne.rdot)rdot=rdotcheck
            deltavx1=bbh_m2/(bbh_m1+bbh_m2)*
     $           (bbh_sep0*sintheta*thetadot-rdot*costheta)
            deltavy1=bbh_m2/(bbh_m1+bbh_m2)*
     $           (-bbh_sep0*costheta*thetadot-rdot*sintheta)
            deltavx2=bbh_m1/(bbh_m1+bbh_m2)*
     $           (-bbh_sep0*sintheta*thetadot+rdot*costheta)
            deltavy2=bbh_m1/(bbh_m1+bbh_m2)*
     $           (bbh_sep0*costheta*thetadot+rdot*sintheta)
            if(myrank.eq.0) then 
               write(69,*) 'compact object 1 starts with:'
               write(69,*) '                     x=',deltax1,
     $              ', y=',deltay1
               write(69,*) '                    vx=',deltavx1,
     $              ',vy=',deltavy1
               write(69,*) 'compact object 2 starts with:'
               write(69,*) '                     x=',deltax2,
     $              ', y=',deltay2
               write(69,*) '                    vx=',deltavx2,
     $              ',vy=',deltavy2
            endif
            xcm=(bbh_m1*deltax1+bbh_m2*deltax2)/(bbh_m1+bbh_m2)
            ycm=(bbh_m1*deltay1+bbh_m2*deltay2)/(bbh_m1+bbh_m2)
            vxcm=(bbh_m1*deltavx1+bbh_m2*deltavx2)/(bbh_m1+bbh_m2)
            vycm=(bbh_m1*deltavy1+bbh_m2*deltavy2)/(bbh_m1+bbh_m2)
            if(myrank.eq.0) write(69,*)'bbh center of mass position:',
     $           xcm,ycm
            if(myrank.eq.0) write(69,*)'bbh center of mass velocity:',
     $           vxcm,vycm
            ak=bbh_m1*bbh_m2
            mu=bbh_m1*bbh_m2/(bbh_m1+bbh_m2)
            r = sqrt((deltax1-deltax2)**2 + (deltay1-deltay2)**2)
            eorb0check = 0.5d0*(bbh_m1*(deltavx1**2+deltavy1**2) +
     $           bbh_m2*(deltavx2**2+deltavy2**2)) - bbh_m1*bbh_m2/r
            altotint = bbh_m1*(deltax1*deltavy1 - deltay1*deltavx1) +
     $           bbh_m2*(deltax2*deltavy2 - deltay2*deltavx2)
            e0check = sqrt(1.d0 + 2.d0*eorb0check*altotint**2/mu/ak**2)
            semilatusrectumprime = r+bbh_e0*(deltax2-deltax1)
            if(abs(r/bbh_sep0-1.d0).gt.1.d-8) then
               write(69,*)'hyperbolic: bbh_sep0 problem',r,bbh_sep0
               stop
            endif
            if(abs(eorb0check-eorb0).gt.1.d-8) then
               write(69,*)'hyperbolic: bbh eorb0 problem',
     $              eorb0check,eorb0
               write(69,*)'Be sure to set two of bbh_semimajoraxis,',
     $              'bbh_e0,bbh_vinf2,bbh_rp in sph.input'
               stop
            endif
            if(abs(altotint/ltot-1.d0).gt.1.d-8) then
               write(69,*)'hyperbolic: bbh angular momentum problem',
     $              altotint,ltot
               stop
            endif
            if(abs(e0check-bbh_e0).gt.1.d-7) then
               write(69,*)'hyperbolic: bbh eccentricity problem',
     $              e0check,bbh_e0
               stop
            endif
            if(abs(semilatusrectumprime/semilatusrectum-1.d0).gt.1.d-8)then
               write(69,*)'hyperbolic: bbh semilatusrectum problem',
     $              semilatusrectumprime,semilatusrectum
               stop
            endif
            if(bbh_vinf2.ge.0.d0) then
               vinf=bbh_vinf2**0.5d0
               if(myrank.eq.0) write(69,*)'bbh v_infinity=',vinf,
     $              '(code units)=',
     $              vinf*(gravconst*munit/runit)**0.5d0/1.d5,'km/s'
               if(myrank.eq.0) write(69,*)'converted with velocity unit',
     $              (gravconst*munit/runit)**0.5d0/1.d5,'km/s'
            endif

            do i=ntot-1,ntot
C     Rotations as described here https://youtu.be/N_njPJYaXaU?t=1983
c     First rotate around z-axis by the argument of periapsis
               if(i.eq.ntot-1) then
                  xold=deltax1
                  yold=deltay1
                  vxold=deltavx1
                  vyold=deltavy1
               else
                  xold=deltax2
                  yold=deltay2
                  vxold=deltavx2
                  vyold=deltavy2
               endif
               zold=0
               vzold=0
               xnew=cos(bbh_argperi)*xold-sin(bbh_argperi)*yold
               ynew=sin(bbh_argperi)*xold+cos(bbh_argperi)*yold
               znew=zold
               vxnew=cos(bbh_argperi)*vxold-sin(bbh_argperi)*vyold
               vynew=sin(bbh_argperi)*vxold+cos(bbh_argperi)*vyold
               vznew=vzold
c     Next rotate around y-axis by the inclination angle
               xold=xnew
               yold=ynew
               zold=znew
               vxold=vxnew
               vyold=vynew
               vzold=vznew
               xnew=cos(bbh_inclination)*xold-sin(bbh_inclination)*zold
               ynew=yold
               znew=sin(bbh_inclination)*xold+cos(bbh_inclination)*zold
               vxnew=cos(bbh_inclination)*vxold-sin(bbh_inclination)*vzold
               vynew=vyold
               vznew=sin(bbh_inclination)*vxold+cos(bbh_inclination)*vzold
c     Finally rotate around z-axis by the longitude of ascending node
               xold=xnew
               yold=ynew
               zold=znew
               vxold=vxnew
               vyold=vynew
               vzold=vznew
               x(i)=cos(bbh_longitude)*xold-sin(bbh_longitude)*yold
               y(i)=sin(bbh_longitude)*xold+cos(bbh_longitude)*yold
               z(i)=zold
               vx(i)=cos(bbh_longitude)*vxold-sin(bbh_longitude)*vyold
               vy(i)=sin(bbh_longitude)*vxold+cos(bbh_longitude)*vyold
               vz(i)=vzold
            enddo
            
         endif
         corepts=corepts+n2
      endif
      amthree=am2


c     From energy conservation 0.5*vinf^2=0.5*vp^2-G*M/rp
c     -> vinf^2=vp^2-2*G*M/rp
c            =(vinf*b/rp)^2-2*G*M/rp because of angular momentum conservation vinf*b = vp*rp
c     -> rp^2=b^2-2*rp*G*M/vinf^2
c            =b^2-2*rp*xx*b where we have defined xx=G*M/(vinf^2*b)
c     -> rp^2*(1.d0+xx^2)=b^2-2*rp*xx*b+rp^2*xx^2
c     -> rp*sqrt(1.d0+xx**2)=b-rp*xx
c     -> b=rp*xx+rp*sqrt(1.d0+xx**2)
c     Solving for rp yields the desired result
      if(rp.lt.0d0 .and. impactparameter.ge.0) then
         xx=(am1+am2)/(vinf2*impactparameter)
         if(myrank.eq.0) write(69,*)'dimensionless parameter x=',xx
         rp=impactparameter/(xx+sqrt(1.d0+xx**2))
      endif

      k=am1*am2
      mu=am1*am2/(am1+am2)

      if(rp.lt.0.d0 .and. semimajoraxis.eq.0.d0)then
c     presumably e0 and vinf2 have been set in sph.input:
         eorb0=0.5d0*mu*vinf2
         semimajoraxis=-0.5d0*k/eorb0
         rp=semimajoraxis*(1.d0-e0)
      else if(vinf2.ge.1d30 .and. semimajoraxis.eq.0.d0)then
c     presumably e0 and (rp or impactparameter) have been set in sph.input:
         semimajoraxis=rp/(1.d0-e0)
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
      else if(vinf2.ge.1d30 .and. e0.lt.0.d0)then
c     presumably semimajoraxis and (rp or impactparameter) have been set in sph.input:
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
         e0=1.d0-rp/semimajoraxis
      else if(rp.lt.0.d0 .and. vinf2.ge.1.d30)then
c     presumably semimajoraxis and e0 have been set in sph.input:
         rp=semimajoraxis*(1.d0-e0)
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
      else
c     presumably (rp or impactparameter) and vinf2 have been set in sph.input:
         eorb0=0.5d0*mu*vinf2
         semimajoraxis=-0.5d0*k/eorb0
         e0=1.d0-rp/semimajoraxis
      endif

      if(myrank.eq.0) then
         write(69,*) 'hyperbolic: impactparameter (if set)=',
     $        impactparameter,'v_inf2=',vinf2,'semimajoraxis=',semimajoraxis
         write(69,*)'hyperbolic: e_orb=',eorb0
         write (69,*)'hyperbolic: n=',n,'ntot=',ntot,'rp=',rp
         open(30,file='m1m2rp.sph')
c     make file m1m2rp.sph with masses and radius in solar units
         write(30,*) n1,n2,rp
         write(30,*) am1*munit/1.9884098706980504d33,am2*munit/1.9884098706980504d33,
     $        rp*runit/6.957d10     
         close(30)
      endif

c      if(myrank.eq.0) write(69,*)k,mu,semimajoraxis,e0

      semilatusrectum=rp*(1.d0+e0)

c     equation (8.41) of marion and thornton
      costheta = (semilatusrectum/sep0-1.d0)/e0
      if(abs(costheta).gt.1.d0 .or. e0.eq.0) then
         sep0=semilatusrectum/(1-e0)
         costheta = -1d0
         if(myrank.eq.0) then
            if(abs(costheta).gt.1.d0) then
               write(69,*)'bad initial conditions in sph.input:'
               write(69,*)'cos(theta)=',cos(theta)
               write(69,*)'sep0 cannot be larger than',
     $              semilatusrectum/(1-e0)
            endif
            write(69,*)'suggestion: reset sep0<=',sep0
            sToP
         endif
      endif

c     equation (8.40) of marion and thornton
      ltot = sqrt(semilatusrectum*mu*k)
c     equation (8.10) of marion and thornton
      thetadot = ltot/mu/sep0**2
c     equation (8.40) of marion and thornton
      eorb0check = (e0**2-1.d0)*mu*k*k/(ltot*ltot)*0.5d0
      if(myrank.eq.0) write(69,*)'1st simple check:',eorb0,eorb0check

      theta=-acos(costheta)
      sinthetacheck = sin(theta)
      sintheta = -sqrt(1.d0-costheta**2)
      if(abs(sintheta-sinthetacheck).gt.1.d-13) then
         write(69,*)'2nd simple check fails',sintheta,sinthetacheck,
     $        abs(sintheta-sinthetacheck)
         stop
      endif
      if(myrank.eq.0) write(69,*)'semilatusrectum=',semilatusrectum,
     $     ' mu=',mu,' sep0=',sep0,' e0=',e0
      if(myrank.eq.0) write(69,*) 'cos=',costheta,' sin=',sintheta
c     the minus signs on the position and velocity component equations
c     have been chosen so that the separation vector r equals r2-r1,
c     *not* r1-r2 as in marion and thornton.  this was done so that the
c     code would give the same initial conditions as twostars.f when the
c     eccentricity is 1
      deltax1 = -am2/(am1+am2)*sep0*costheta
      deltay1 = -am2/(am1+am2)*sep0*sintheta
      deltax2 = am1/(am1+am2)*sep0*costheta
      deltay2 = am1/(am1+am2)*sep0*sintheta
c     differentiating eqn. 8.41
      rdot=semilatusrectum/
     $     (1.d0+e0*costheta)**2*e0*sintheta*thetadot
c     another expr. for rdot (using eqn. 8.14) and compare
      rdotcheck=-sqrt((eorb0-0.5d0*ltot**2/mu/sep0**2+k/sep0)/0.5d0/mu)
      if(myrank.eq.0) write(69,*) 'rdot(from semilatusrectum): ',rdot,
     $     '  rdot(from e): ',rdotcheck
      if(rp.eq.0)rdot=rdotcheck
      deltavx1=am2/(am1+am2)*(sep0*sintheta*thetadot-rdot*costheta)
      deltavy1=am2/(am1+am2)*(-sep0*costheta*thetadot-rdot*sintheta)
      deltavx2=am1/(am1+am2)*(-sep0*sintheta*thetadot+rdot*costheta)
      deltavy2=am1/(am1+am2)*(sep0*costheta*thetadot+rdot*sintheta)
      if(myrank.eq.0) then 
        write(69,*) 'star 1 starts with: x=',deltax1,',y=',deltay1
        write(69,*) '                    vx=',deltavx1,',vy=',deltavy1
        write(69,*) 'star 2 starts with: x=',deltax2,',y=',deltay2
        write(69,*) '                    vx=',deltavx2,',vy=',deltavy2
      endif
      xcm=(am1*deltax1+am2*deltax2)/(am1+am2)
      ycm=(am1*deltay1+am2*deltay2)/(am1+am2)
      vxcm=(am1*deltavx1+am2*deltavx2)/(am1+am2)
      vycm=(am1*deltavy1+am2*deltavy2)/(am1+am2)
      if(myrank.eq.0) write(69,*)'center of mass position:',xcm,ycm
      if(myrank.eq.0) write(69,*)'center of mass velocity:',vxcm,vycm
      ak=am1*am2
      amu=am1*am2/(am1+am2)
      r = sqrt((deltax1-deltax2)**2 + (deltay1-deltay2)**2)
      eorb0check = 0.5d0*(am1*(deltavx1**2+deltavy1**2) +
     $     am2*(deltavx2**2+deltavy2**2)) - am1*am2/r
      altotint = am1*(deltax1*deltavy1 - deltay1*deltavx1) +
     $     am2*(deltax2*deltavy2 - deltay2*deltavx2)
      e0check = sqrt(1.d0 + 2.d0*eorb0check*altotint**2/amu/ak**2)
      semilatusrectumprime = r+e0*(deltax2-deltax1)
      if(abs(r/sep0-1.d0).gt.1.d-8) then
         write(69,*)'hyperbolic: sep0 problem',r,sep0
         stop
      endif
      if(abs(eorb0check-eorb0).gt.1.d-8) then
         write(69,*)'hyperbolic: eorb0 problem',eorb0check,eorb0
         stop
      endif
      if(abs(altotint/ltot-1.d0).gt.1.d-8) then
         write(69,*)'hyperbolic: angular momentum problem',altotint,ltot
         stop
      endif
      if(abs(e0check-e0).gt.1.d-7) then
         write(69,*)'hyperbolic: eccentricity problem',e0check,e0
         stop
      endif
      if(abs(semilatusrectumprime/semilatusrectum-1.d0).gt.1.d-8)then
         write(69,*)'hyperbolic: semilatusrectum problem',
     $        semilatusrectumprime,semilatusrectum
         stop
      endif
c      vinf=(2.d0*eorb0/mu)**0.5d0
      if(vinf2.ge.0.d0) then
         vinf=vinf2**0.5d0
         if(myrank.eq.0) write(69,*)'v_infinity=',vinf,'(to units)=',
     $        vinf*(gravconst*munit/runit)**0.5d0/1.d5,'km/s'
         if(myrank.eq.0) write(69,*)'converted with',
     $        (gravconst*munit/runit)**0.5d0/1.d5
      endif

      do i=1,n1
         x(i)=x(i)+deltax1
         y(i)=y(i)+deltay1
         vx(i)=vx(i)+deltavx1
         vy(i)=vy(i)+deltavy1
      enddo
      do i=n1+1,ntot
         x(i)=x(i)+deltax2
         y(i)=y(i)+deltay2
         vx(i)=vx(i)+deltavx2
         vy(i)=vy(i)+deltavy2
      enddo
      call stride_setup
      
c     prepare leap-frog scheme for first iteration:
      call lfstart
      if(myrank.eq.0) write (69,*) 'hyperbolic:          ... done'
      
      return
      end
