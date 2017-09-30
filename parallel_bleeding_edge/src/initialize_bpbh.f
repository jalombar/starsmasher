      subroutine bpbh
      include 'starsmasher.h'
      real*8 k,rdot,rdotcheck,costheta,sintheta,ltot,thetadot,
     $     semilatusrectum,mu,e0check,sinthetacheck,
     $     eorb0check
      real*8 altotint,r,amu,ak,vxcm,vycm,xcm,ycm
      integer n2,i,nchk,corepts
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      real*8 deltax1,deltay1,deltavx1,deltavy1
      real*8 deltax2,deltay2,deltavx2,deltavy2
      real*8 vinf,semilatusrectumprime,theta
      real*8 th1,costh2,sinth2,ps1,ps2,ph1
      real*8 xold,yold,zold,vxold,vyold,vzold
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 eorb0
      real*8 lambda11,lambda12,lambda13
      real*8 lambda21,lambda22,lambda23
      real*8 lambda31,lambda32,lambda33
      real*8 bhmass

      write(69,*)'bpbh: rp=',rp,'v_inf2=',vinf2
      corepts=0
      write (69,*) 'bpbh: reading start files ...'
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

      open(30,file='sph.bpbh')
      read(30,*) ph1,th1,ps1
      read(30,*) bhmass
      close(30)

      if(sin(th1).lt.0.d0) then
         write(69,*)'expecting theta1 to be between 0 and pi'
         stop
      endif

      costh2=1.d0
      sinth2=0.d0
      ps2=0.d0
      write(69,*)'rotation angles for binary (in radians):'
      write(69,*)'phase phi1=',ph1,'theta1=',th1,' psi1=',ps1
      write(69,*)'rotation angles for star 2 (in radians):'
      write(69,*)'theta2=',acos(costh2),' psi2=',ps2

c     determine coefficients of rotation matrix for the binary
      lambda11= cos(ps1)*cos(ph1)-cos(th1)*sin(ph1)*sin(ps1)
      lambda21=-sin(ps1)*cos(ph1)-cos(th1)*sin(ph1)*cos(ps1)
      lambda31= sin(th1)*sin(ph1)
      lambda12= cos(ps1)*sin(ph1)+cos(th1)*cos(ph1)*sin(ps1)
      lambda22=-sin(ps1)*sin(ph1)+cos(th1)*cos(ph1)*cos(ps1)
      lambda32=-sin(th1)*cos(ph1)
      lambda13= sin(ps1)*sin(th1)
      lambda23= cos(ps1)*sin(th1)
      lambda33= cos(th1)

      do i=1,n1
         read (12) xold,yold,zold,am(i),hp(i),rho(i),vxold,vyold,
     $        vzold,vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)

c     place velocities at same time as everything else:
         vxold=vxold-vxdot(i)*0.5d0*dtold
         vyold=vyold-vydot(i)*0.5d0*dtold
         vzold=vzold-vzdot(i)*0.5d0*dtold

         am1=am1+am(i)
         x(i)= lambda11*xold  + lambda12*yold  + lambda13*zold
         y(i)= lambda21*xold  + lambda22*yold  + lambda23*zold
         z(i)= lambda31*xold  + lambda32*yold  + lambda33*zold
         vx(i)=lambda11*vxold + lambda12*vyold + lambda13*vzold
         vy(i)=lambda21*vxold + lambda22*vyold + lambda23*vzold
         vz(i)=lambda31*vxold + lambda32*vyold + lambda33*vzold
      enddo
      read(12) nchk
      close(12)

      if (nchk.ne.n1) stop 'bpbh: problem with file'


      n2=1
      ntot=n1+n2
      if (ntot.gt.nmax) then
         write(69,*)'must increase nmax...'
         stop
      endif

      amass2=n2
      i=n1+1
      x(i)=0.d0
      y(i)=0.d0
      z(i)=0.d0
      vx(i)=0.d0
      vy(i)=0.d0
      vz(i)=0.d0
      am(i)=bhmass*1.98843d33/munit
      am2=am(i)
      hp(i)=-999.d0
      rho(i)=0.d0
      vxdot(i)=0.d0
      vydot(i)=0.d0
      vzdot(i)=0.d0
      u(i)=0.d0
      udot(i)=0.d0
      gx(i)=0.d0
      gy(i)=0.d0
      gz(i)=0.d0
      grpot(i)=0.d0
      meanmolecular(i)=0.d0
      cc(i)=0
      corepts=corepts+1
      n=ntot-corepts


      write (69,*)'bpbh: n=',n,'ntot=',ntot

      open(30,file='m1m2rp.sph')
c     make file m1m2rp.sph with masses and radius in solar units
      write(30,*) n1,n2,rp
      write(30,*) am1*munit/1.98843d33,am2*munit/1.98843d33,
     $     rp*runit/6.9599d10     
      close(30)

      k=am1*am2
      mu=am1*am2/(am1+am2)
c      vinf=(2.d0*eorb0/mu)**0.5d0
      eorb0=0.5d0*mu*vinf2
      write(69,*)'bpbh: e_orb=',eorb0
      semimajoraxis=-0.5d0*k/eorb0
      e0=1.d0-rp/semimajoraxis
      semilatusrectum=rp*(1.d0+e0)
c     equation (8.40) of marion and thornton
      ltot = sqrt(semilatusrectum*mu*k)
c     equation (8.10) of marion and thornton
      thetadot = ltot/mu/sep0**2
c     equation (8.40) of marion and thornton
      eorb0check = (e0**2-1.d0)*mu*k*k/(ltot*ltot)*0.5d0
      write(69,*)'1st simple check:',eorb0,eorb0check
c     equation (8.41) of marion and thornton
      costheta = (semilatusrectum/sep0-1.d0)/e0
      theta=-acos(costheta)
      sinthetacheck = sin(theta)
      sintheta = -sqrt(1.d0-costheta**2)
      if(abs(sintheta-sinthetacheck).gt.1.d-15) then
         write(69,*)'2nd simple check fails',sintheta,sinthetacheck
         stop
      endif
      write(69,*)'semilatusrectum=',semilatusrectum,' mu=',mu,
     $     ' sep0=',sep0,' e0=',e0
      write(69,*) 'cos=',costheta,' sin=',sintheta
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
      rdotcheck=semilatusrectum/
     $     (1.d0+e0*costheta)**2*e0*sintheta*thetadot
c     another expr. for rdot (using eqn. 8.14) and compare
      rdot=-sqrt((eorb0-0.5d0*ltot**2/mu/sep0**2+k/sep0)/0.5d0/mu)
      write(69,*) 'rdot(from semilatusrectum): ',rdotcheck,
     $     '  rdot(from e): ',rdot
      deltavx1=am2/(am1+am2)*(sep0*sintheta*thetadot-rdot*costheta)
      deltavy1=am2/(am1+am2)*(-sep0*costheta*thetadot-rdot*sintheta)
      deltavx2=am1/(am1+am2)*(-sep0*sintheta*thetadot+rdot*costheta)
      deltavy2=am1/(am1+am2)*(sep0*costheta*thetadot+rdot*sintheta)
      write(69,*) 'star 1 starts with: x=',deltax1,',y=',deltay1
      write(69,*) '                    vx=',deltavx1,',vy=',deltavy1
      write(69,*) 'star 2 starts with: x=',deltax2,',y=',deltay2
      write(69,*) '                    vx=',deltavx2,',vy=',deltavy2
      xcm=(am1*deltax1+am2*deltax2)/(am1+am2)
      ycm=(am1*deltay1+am2*deltay2)/(am1+am2)
      vxcm=(am1*deltavx1+am2*deltavx2)/(am1+am2)
      vycm=(am1*deltavy1+am2*deltavy2)/(am1+am2)
      write(69,*)'center of mass position:',xcm,ycm
      write(69,*)'center of mass velocity:',vxcm,vycm
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
         write(69,*)'bpbh: sep0 problem'
         stop
      endif
      if(abs(eorb0check-eorb0).gt.1.d-8) then
         write(69,*)'bpbh: eorb0 problem'
         stop
      endif
      if(abs(altotint/ltot-1.d0).gt.1.d-8) then
         write(69,*)'bpbh: angular momentum problem'
         stop
      endif
      if(abs(e0check/e0-1.d0).gt.1.d-8) then
         write(69,*)'bpbh: eccentricity problem'
         stop
      endif
      if(abs(semilatusrectumprime/semilatusrectum-1.d0).gt.1.d-8)then
         write(69,*)'bpbh: semilatusrectum problem'
         stop
      endif
c      vinf=(2.d0*eorb0/mu)**0.5d0
      if(vinf2.ge.0.d0) then
         vinf=vinf2**0.5d0
         write(69,*)'v_infinity=',vinf,'(to units)=',
     $        vinf*(gravconst*munit/runit)**0.5d0/1.d5,'km/s'
         write(69,*)'converted with',(gravconst*munit/runit)**0.5d0/1.d5
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
      write (69,*) 'bpbh:          ... done'
      
      return
      end
