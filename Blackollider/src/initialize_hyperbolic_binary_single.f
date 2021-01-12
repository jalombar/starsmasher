      subroutine hyperbolic_binary_single
c     hyperbolic collision with a neutron star
      include 'starsmasher.h'
      real*8 k,rdot,rdotcheck,costheta,sintheta,ltot,thetadot,
     $     semilatusrectum,mu,e0check,sinthetacheck,
     $     eorb0check,eorb0
      real*8 altotint,r,amu,ak,vxcm,vycm,xcm,ycm
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hcoold,hfloorold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,trelaxold,dtold
      real*8 deltax1,deltay1,deltavx1,deltavy1
      real*8 deltax2,deltay2,deltavx2,deltavy2
      integer i,nchk,corepts,ccmax
      real*8 vinf,semilatusrectumprime,theta
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2,xx
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      real*8 hpcore
      integer idumb
      real*8 costh2,sinth2,ps2
      real*8 xold,yold,zold,vxold,vyold,vzold,ran1,hmax,hmin

      call cpu_time(time1)

      hpcore=1d30

c         write(69,*)'hns: rp=',rp
      if(vinf2.ne.0) then
         write(69,*) 'hbs: bimpact=',bimpact,'v_inf2=',vinf2
      endif
         corepts=0
         write (69,*) 'hbs: reading start files ...'
         open(12,file='sph.startu',form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
         read(12) n1,nnoptold,hcoold,hfloorold,sep0old,
     $        tfold,dtoutold,noutold,nitold,told,navold,
     $        alphaold,betaold,tjumpahead,ngrold,
     $        nrelaxold,trelaxold,dtold
         am1=0.d0
         n=n1+2
         ntot=n
         if (ntot.gt.nmax) then
            write(69,*)'must increase nmax...'
            stop
         endif
         ccmax=0
         do i=1,n1
            read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $           vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $           gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $           cc(i)
            ccmax=max(ccmax,cc(i))
            am1=am1+am(i)
            if(u(i).eq.0.d0) then
               hpcore=hp(i)
               corepts=corepts+1
            endif
         enddo
         read(12) nchk
         close(12)

      if (nchk.ne.n1) stop 'hbs: problem with file'

      if(hpcore.eq.1d30) then
         hmin=1.d30
         hmax=0.d0
         do i=1,n1
            hmin=min(hmin,hp(i))
            hmax=max(hmax,hp(i))
         enddo
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
         hpcore=hmin
      endif

c     make the bh-wd binary:
         corepts=corepts+2
         am(ntot-1)=mbh         ! a bh mass
         am(ntot)=0.6d0         ! a wd mass
         am2=am(ntot-1)+am(ntot)
         x(ntot-1)=-sepfinal*am(ntot)/am2
         x(ntot)=+sepfinal*am(ntot-1)/am2
         vy(ntot-1)=-(am2/sepfinal)**0.5d0*am(ntot)/am2
         vy(ntot)=+(am2/sepfinal)**0.5d0*am(ntot-1)/am2
         do i=ntot-1,ntot
            y(i)=0.d0
            z(i)=0.d0
            vx(i)=0.d0
            vz(i)=0.d0
            vxdot(i)=0.d0       ! accelerations and udot
            vydot(i)=0.d0       ! will be set when lfstart
            vzdot(i)=0.d0       ! is called
            u(i)=0.d0           !
            udot(i)=0.d0        ! 
            hp(i)=2.d0*hpcore
            cc(i)=ccmax
         enddo

         if(vinf2.ne.0.d0)then
            xx=(am1+am2)/(vinf2*bimpact)
            write(69,*)'dimensionless parameter x=',xx
            rp=bimpact/(xx+sqrt(1.d0+xx**2))
         else
            rp=bimpact
            xx=1.d10
            bimpact=rp*(xx+sqrt(1.d0+xx**2))
            vinf2=(am1+am2)/(xx*bimpact)
            write(69,*) 'hbs: bimpact=',bimpact,'v_inf2=',vinf2
         endif
         write (69,*)'hbs: n=',n,'ntot=',ntot,'rp=',rp

         k=am1*am2
         mu=am1*am2/(am1+am2)
         eorb0=0.5d0*mu*vinf2
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
     $        ' sep0=',sep0,' e0=',e0
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
     $        (1.d0+e0*costheta)**2*e0*sintheta*thetadot
c     another expr. for rdot (using eqn. 8.14) and compare
         rdot=-sqrt((eorb0-0.5d0*ltot**2/mu/sep0**2+k/sep0)/0.5d0/mu)
         write(69,*) 'rdot(from semilatusrectum): ',rdotcheck,
     $        '  rdot(from e): ',rdot
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
     $        am2*(deltavx2**2+deltavy2**2)) - am1*am2/r
         altotint = am1*(deltax1*deltavy1 - deltay1*deltavx1) +
     $        am2*(deltax2*deltavy2 - deltay2*deltavx2)
         e0check = sqrt(1.d0 + 2.d0*eorb0check*altotint**2/amu/ak**2)
         semilatusrectumprime = r+e0*(deltax2-deltax1)
         if(abs(r/sep0-1.d0).gt.1.d-8) then
            write(69,*)'hbs: sep0 problem'
            stop
         endif
         if(abs(eorb0check-eorb0).gt.1.d-8) then
            write(69,*)'hbs: eorb0 problem'
            stop
         endif
         if(abs(altotint/ltot-1.d0).gt.1.d-8) then
            write(69,*)'hbs: angular momentum problem'
            stop
         endif
         if(abs(e0check/e0-1.d0).gt.1.d-8) then
            write(69,*)'hbs: eccentricity problem'
            stop
         endif
         if(abs(semilatusrectumprime/semilatusrectum-1.d0).gt.1.d-8)then
            write(69,*)'hbs: semilatusrectum problem'
            stop
         endif
         vinf=(2.d0*eorb0/mu)**0.5d0
         write(69,*)'v_infinity=',vinf,'(sph units)=',
     $        vinf*(gravconst*munit/runit)**0.5d0/1.d5,'km/s'

      call cpu_time(time2)
      if(.false.)then
         seconds = dabs(vinf2)/7.d0+1.d-7
      else
         seconds = time2-time1
      endif
      write (69,'(a,g16.9,a)')
     $     'hyperbolic: ',seconds,' will help set random number'
      idumb=-nint(10000*(100000*seconds-int(100000*seconds)))
      write(69,*)'idumb=',idumb
      costh2=2.d0*(ran1(idumb)-0.5d0)
      sinth2=sqrt(1.d0-costh2**2.d0)
      ps2=2.d0*pi*(ran1(idumb)-0.5d0)      
      write(69,*)'rotation angles for binary (in radians):'
      write(69,*)'theta2=',acos(costh2),' psi2=',ps2

      do i=n1+1,ntot
         xold=x(i)
         yold=y(i)
         zold=z(i)
         vxold=vx(i)
         vyold=vy(i)
         vzold=vz(i)
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
      enddo

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
      write(69,*)'hbs:          ... done'
      return
      end
