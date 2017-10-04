      subroutine hyperbolic
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 k,rdot,rdotcheck,costheta,sintheta,ltot,thetadot,
     $     semilatusrectum,mu,e0check,sinthetacheck,
     $     eorb0check
      real*8 altotint,r,amu,ak,vxcm,vycm,xcm,ycm
      integer n2,i,nchk,corepts
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      logical twofiles
      real*8 deltax1,deltay1,deltavx1,deltavy1
      real*8 deltax2,deltay2,deltavx2,deltavy2
      real*8 vinf,semilatusrectumprime,theta
      integer idumb
      real*8 costh1,sinth1,costh2,sinth2,ps1,ps2
      real*8 xold,yold,zold,vxold,vyold,vzold,ran1
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
      common/orbitalelements/e0,semimajoraxis

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
         ps1=0.d0
      endif

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
      enddo
      read(12) nchk
      close(12)

      amtwo=am1

      if (nchk.ne.n1) stop 'hyperbolic: problem with sph.start1u file'
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
         read(12) n2,nnoptold,hcoold,hfloorold,sep0old,
     $        tfold,dtoutold,noutold,nitold,told,navold,
     $        alphaold,betaold,tjumpahead,ngrold,
     $        nrelaxold,trelaxold,dtold,omega2
         amass2=n2
         am2=0.d0
         ntot=n1+n2
         if (ntot.gt.nmax) then
            if(myrank.eq.0) write(69,*)'must increase nmax...'
            stop
         endif
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
         enddo
         n=ntot-corepts
         read(12) nchk
         close(12)
         if (nchk.ne.n2) stop
     $        'hyperbolic: problem with sph.start2u file'
      else
         if(myrank.eq.0) write (69,*)
     $        'startfile2 ',
     $        trim(startfile2),
     $        ' does not exist: will use black hole instead'
         
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
            hp(i)=hp(n1)
         endif
         cc(i)=0
         if(myrank.eq.0) write(69,*)'black_hole_mass',am(i)
         if(myrank.eq.0) write(69,*)'black_hole_smoothing_length',hp(i)
      endif
      amthree=am2

      xx=(am1+am2)/(vinf2*bimpact)
      if(myrank.eq.0) write(69,*)'dimensionless parameter x=',xx
      rp=bimpact/(xx+sqrt(1.d0+xx**2))


      if(myrank.eq.0) then
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
        write(69,*) 'setting rp=bimpact'
      endif

      k=am1*am2
      mu=am1*am2/(am1+am2)

      if(bimpact.lt.0.d0 .and. semimajoraxis.eq.0.d0)then
c     presumably e0 and vinf2 have been set in sph.input:
         eorb0=0.5d0*mu*vinf2
         semimajoraxis=-0.5d0*k/eorb0
         rp=semimajoraxis*(1.d0-e0)
      else if(vinf2.ge.1d30 .and. semimajoraxis.eq.0.d0)then
c     presumably e0 and bimpact have been set in sph.input:
         rp=bimpact
         semimajoraxis=rp/(1.d0-e0)
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
      else if(vinf2.ge.1d30 .and. e0.lt.0.d0)then
c     presumably semimajoraxis and bimpact have been set in sph.input:
         rp=bimpact
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
         e0=1.d0-rp/semimajoraxis
      else if(bimpact.lt.0.d0 .and. vinf2.ge.1.d30)then
c     presumably semimajoraxis and e0 have been set in sph.input:
         rp=semimajoraxis*(1.d0-e0)
         eorb0=-0.5d0*k/semimajoraxis
         vinf2=2.d0*eorb0/mu
      else
c     presumably bimpact and vinf2 have been set in sph.input:
         rp=bimpact
         eorb0=0.5d0*mu*vinf2
         semimajoraxis=-0.5d0*k/eorb0
         e0=1.d0-rp/semimajoraxis
      endif

      if(myrank.eq.0) then
         write(69,*) 'hyperbolic: bimpact=',
     $        bimpact,'v_inf2=',vinf2,'semimajoraxis=',semimajoraxis
         write(69,*)'hyperbolic: e_orb=',eorb0
         write (69,*)'hyperbolic: n=',n,'ntot=',ntot,'rp=',rp
         open(30,file='m1m2rp.sph')
c     make file m1m2rp.sph with masses and radius in solar units
         write(30,*) n1,n2,rp
         write(30,*) am1*munit/1.98843d33,am2*munit/1.98843d33,
     $        rp*runit/6.9599d10     
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
            write(69,*)'reset sep0=',sep0
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
      if(myrank.eq.0) write(69,*)'semilatusrectum=',semilatusrectum,' mu=',mu,
     $     ' sep0=',sep0,' e0=',e0
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
      if(rdot.ne.rdot)rdot=rdotcheck
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
