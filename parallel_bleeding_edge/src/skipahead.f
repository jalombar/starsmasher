      subroutine jumpahead
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 am1,am2
      integer icomp(nmax)
      real*8 x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4
      common/compbettercom3/am1,x1,y1,z1,vx1,vy1,vz1,
     $     am2,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp
      real*8 ecc,eorb
      real*8 lambda11,lambda12,lambda13
      real*8 lambda21,lambda22,lambda23
      real*8 lambda31,lambda32,lambda33
      real*8 lorb,lorbx,lorby,lorbz
      real*8 lonx,lony,lonhatx,lonhaty,lontot
      real*8 g
      parameter(g=1.d0)
      real*8 mu,semilatusrectum,k
      real*8 xcm,ycm,zcm,vxcm,vycm,vzcm
      real*8 ecc2
      real*8 r12,r,rdot,rdotcheck,costheta,sintheta,thetadot,
     $     ecccheck,sinthetacheck,
     $     eorbcheck
      real*8 ltotcheck
      integer i
      real*8 deltax1,deltay1,deltavx1,deltavy1
      real*8 deltax2,deltay2,deltavx2,deltavy2
      real*8 vinf,semilatusrectumprime,theta
      real*8 costh1,sinth1,cosps1,sinps1,cosph1,sinph1
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 divv(nmax)
      common/commdivv/divv
      real*8 dxold1,dyold1,dzold1
      real*8 dvxold1,dvyold1,dvzold1
      real*8 dxold2,dyold2,dzold2
      real*8 dvxold2,dvyold2,dvzold2
      real*8 dxnew1,dynew1,dznew1
      real*8 dvxnew1,dvynew1,dvznew1
      real*8 dxnew2,dynew2,dznew2
      real*8 dvxnew2,dvynew2,dvznew2
      real*8 thrownawaymass
      integer mygravlength,ierr
      integer comm_worker
      common/gravworkers/comm_worker
      real*8 eccx, eccy, eccz
      integer bhcomp
      real*8 displacex, displacey,displacez
      integer ndisplace
      common/displace/displacex,displacey,displacez,ndisplace

c     Initializing bhcomp
      bhcomp = -1

      do i=1,n
         vx(i)=vx(i)-0.5d0*dt*vxdot(i)
         vy(i)=vy(i)-0.5d0*dt*vydot(i)
         vz(i)=vz(i)-0.5d0*dt*vzdot(i)
         if(u(i).ne.0.d0) then
            u(i)=u(i)-0.5d0*dt*udot(i)
         endif
      enddo
      t=t-0.5d0*dt

      if(myrank.eq.0)write(69,*)
     $     '***analyze system right before jump ahead:***'
      call rho_and_h
      call gravforce
      if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
         if(nusegpus.eq.1)then
            call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
         else
            call get_gravity_using_cpus
         endif
         mygravlength=ngrav_upper-ngrav_lower+1
         if(myrank.ne.0)then
            call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
     $           grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $           comm_worker, ierr)
         else
            call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $           grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $           comm_worker, ierr)
         endif
      endif
      call enout(.false.)
      if(myrank.eq.0)write(69,*)
     $     '***done analyzing system right before jump ahead***'

      if(throwaway)then
         call compbest3(.true.)
      else
c     When throwaway is set equal to .false. (the default value) then...
c     We will be keeping *all* mass, and assuming it is in two separate components, with
c     the point masses making up one component.  This will need to be fixed if we ever
c     have, for example, a red giant star with a point mass core.
         am1=0d0
         x1=0d0
         y1=0d0
         z1=0d0
         vx1=0d0
         vy1=0d0
         vz1=0d0
         am2=0d0
         x2=0d0
         y2=0d0
         z2=0d0
         vx2=0d0
         vy2=0d0
         vz2=0d0
         am3=0d0
         x3=0d0
         y3=0d0
         z3=0d0
         vx3=0d0
         vy3=0d0
         vz3=0d0
         am4=0d0
         do i=1,n
            if(u(i).eq.0) then
               am1=am1+am(i)
               x1=x1+am(i)*x(i)
               y1=y1+am(i)*y(i)
               z1=z1+am(i)*z(i)
               vx1=vx1+am(i)*vx(i)
               vy1=vy1+am(i)*vy(i)
               vz1=vz1+am(i)*vz(i)
               icomp(i)=1
            else
               am2=am2+am(i)
               x2=x2+am(i)*x(i)
               y2=y2+am(i)*y(i)
               z2=z2+am(i)*z(i)
               vx2=vx2+am(i)*vx(i)
               vy2=vy2+am(i)*vy(i)
               vz2=vz2+am(i)*vz(i)
               icomp(i)=2
            endif
         enddo
         if(am1.gt.0) then
            x1=x1/am1
            y1=y1/am1
            z1=z1/am1
            vx1=vx1/am1
            vy1=vy1/am1
            vz1=vz1/am1
         endif
         if(am2.gt.0) then
            x2=x2/am2
            y2=y2/am2
            z2=z2/am2
            vx2=vx2/am2
            vy2=vy2/am2
            vz2=vz2/am2
         endif
      endif

      if(myrank.eq.0) then
         if(am1.gt.0)then
            write(69,'(a,3g17.9)')'position of star 1=',x1,y1,z1
            write(69,'(a,3g17.9)')'velocity of star 1=',vx1,vy1,vz1
         endif
         
         if(am2.gt.0)then
            write(69,'(a,3g17.9)')'position of star 2=',x2,y2,z2
            write(69,'(a,3g17.9)')'velocity of star 2=',vx2,vy2,vz2
         endif
         
         if(am3.gt.0)then
            write(69,'(a,3g17.9)')'position of star 3=',x3,y3,z3
            write(69,'(a,3g17.9)')'velocity of star 3=',vx3,vy3,vz3
         endif
      endif

      if (am3.gt.min(am1,am2))then
         if(am3.gt.am1 .and.am1.le.am2) then
            if(myrank.eq.0)write(69,*)'renaming star 3 to star 1'
            am1=am3
            x1=x3
            y1=y3
            z1=z3
            vx1=vx3
            vy1=vy3
            vz1=vz3
            am3=0.d0
            do i=1,ntot
               if(icomp(i).eq.3)then
                  icomp(i)=1
               elseif(icomp(i).eq.1)then
                  icomp(i)=3
               endif
            enddo
         else
            if(myrank.eq.0)write(69,*)'renaming star 3 to star 2'
            am2=am3
            x2=x3
            y2=y3
            z2=z3
            vx2=vx3
            vy2=vy3
            vz2=vz3
            am3=0.d0
            do i=1,ntot
               if(icomp(i).eq.3)then
                  icomp(i)=2
               elseif(icomp(i).eq.2)then
                  icomp(i)=3
               endif
            enddo
         endif
      endif

      r12=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

      if(myrank.eq.0)write (69,'(3(a,g17.9))')
     $     'eccentricity: r12=',r12,' am1=',am1,' am2=',am2,
     $     ' am3=',am3

c calculate center of mass position and velocity
      xcm=(am1*x1+am2*x2)/(am1+am2)
      ycm=(am1*y1+am2*y2)/(am1+am2)
      zcm=(am1*z1+am2*z2)/(am1+am2)
      vxcm=(am1*vx1+am2*vx2)/(am1+am2)
      vycm=(am1*vy1+am2*vy2)/(am1+am2)
      vzcm=(am1*vz1+am2*vz2)/(am1+am2)
      if(myrank.eq.0) then
         write(69,'(a,3g17.9)') 'center of mass position=',xcm,ycm,zcm
         write(69,'(a,3g17.9)') 'center of mass velocity=',vxcm,vycm,vzcm
      endif

      mu=am1*am2/(am1+am2)
      eorb=0.5d0*mu*((vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2)-am1*am2/r12
      lorbz=mu*((x1-x2)*(vy1-vy2)-(y1-y2)*(vx1-vx2))
      lorbx=mu*((y1-y2)*(vz1-vz2)-(z1-z2)*(vy1-vy2))
      lorby=mu*((z1-z2)*(vx1-vx2)-(x1-x2)*(vz1-vz2))
      lorb=(lorbx**2+lorby**2+lorbz**2)**0.5d0

      k=g*am1*am2
c      write(69,*)'force constant k=',k

c     let's get the components of the laplace-runge-lenz vector, defined
c     as p x l - mu k rhat, where p is momentum, l is
c     angular momentum, mu is reduced mass, k=g*am1*am2.
c     the l-r-l vector is conserved in the kepler problem and points in
c     the direction of periapse.
c      lrlx=mu*((vy1-vy2)*lorbz-(vz1-vz2)*lorby)-mu*k*(x1-x2)/r12
c      lrly=mu*((vz1-vz2)*lorbx-(vx1-vx2)*lorbz)-mu*k*(y1-y2)/r12
c      lrlz=mu*((vx1-vx2)*lorby-(vy1-vy2)*lorbx)-mu*k*(z1-z2)/r12

c     actually, let's get the eccentricity vector, which is proportional
c     to the lrl vector but its magnitude is the eccentricity
      eccx=((vy1-vy2)*lorbz-(vz1-vz2)*lorby)/k-(x1-x2)/r12
      eccy=((vz1-vz2)*lorbx-(vx1-vx2)*lorbz)/k-(y1-y2)/r12
      eccz=((vx1-vx2)*lorby-(vy1-vy2)*lorbx)/k-(z1-z2)/r12

      if(myrank.eq.0) then
         write(69,*)'reduced mass mu=',mu
         write(69,'(a,3g13.4)')'total orbital energy=',eorb
         write(69,'(a,4g13.4)')'total angular momentum=',lorb,
     $        lorbx,lorby,lorbz
         write(69,'(a,4g13.4)')'components of eccentricity vector=',
     $        eccx,eccy,eccz
      endif

      semilatusrectum=lorb**2/(mu*k)
c      write(69,*)'semi-latus rectum alpha=',semilatusrectum

      ecc2=1.d0+2.d0*eorb*lorb**2.d0/(mu*k**2.d0)

      if(abs(ecc2-(eccx**2+eccy**2+eccz**2)).gt.1.d-14)then
         write(69,*)'eccentricity (vector) problem',
     $        ecc2-(eccx**2+eccy**2+eccz**2)
         stop
      endif

      if(ecc2.lt.0.d0) then
         write(69,'(a,5g17.9)')'whoops, ecc squared=',ecc2,eorb,lorb,mu,k
         stop
      endif
      ecc=ecc2**0.5d0

      if(myrank.eq.0) then
         write(69,*)'apastron separation rmax=',semilatusrectum/(1.d0-ecc)
         
         write(69,'(33a12)')
     $        't','e','alpha','eorb','lorb','r12',
     $        'am1','am2','am3','x1','y1','z1','x2','y2','z2'
         write(69,'(33d12.4)') t,ecc,semilatusrectum,eorb,lorb,
     $        r12,am1,am2,am3,x1,y1,z1,x2,y2,z2
      endif
         
      if(ecc.ne.ecc .or. ecc.lt.0.d0) stop

      sep0=0.5d0*r12
c      sep0=r12
      if(myrank.eq.0) write(69,*)'sep0=',sep0

c     equation (8.10) of marion and thornton
      thetadot = lorb/mu/sep0**2
c     equation (8.40) of marion and thornton
      eorbcheck = (ecc2-1.d0)*mu*k*k/(lorb*lorb)*0.5d0
      if(myrank.eq.0) write(69,*)'1st simple check:',eorb,eorbcheck
c     equation (8.41) of marion and thornton
      costheta = (semilatusrectum/sep0-1.d0)/ecc
c     theta=0 corresponds to periapse, theta=+-pi corresponds to apapse.
c     arccos returns a non-negative number, so the minus sign in the
c     following equation will move the star to the pre-periapse rather
c     than post-periapse phase of the orbit:
      theta=-acos(costheta)
      sinthetacheck = sin(theta)
      sintheta = -sqrt(1.d0-costheta**2)
      if(abs(sintheta-sinthetacheck).gt.1.d-15) then
         write(69,*)'2nd simple check fails',sintheta,sinthetacheck
         stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if, as a test, you want *not* to move particles then use this:
c     
c     sintheta=-sintheta
c
c     there is another part below that also would need changed
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(myrank.eq.0) then
         write(69,*)'semilatusrectum=',semilatusrectum,' mu=',mu,
     $        ' sep0=',sep0,' ecc=',ecc
         write(69,*) 'cos=',costheta,' sin=',sintheta
      endif
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
     $     (1.d0+ecc*costheta)**2*ecc*sintheta*thetadot
c     another expr. for rdot (using eqn. 8.14) and compare
      rdot=-sqrt((eorb-0.5d0*lorb**2/mu/sep0**2+k/sep0)/0.5d0/mu)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if, as a test, you want *not* to move particles then use this:
c     
c     rdot=-rdot
c
c     there is another part above that also would need changed
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(myrank.eq.0) then
         write(69,*) 'rdot(from semilatusrectum): ',rdotcheck,
     $        '  rdot(from e): ',rdot
      endif
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

      r = sqrt((deltax1-deltax2)**2 + (deltay1-deltay2)**2)
      eorbcheck = 0.5d0*(am1*(deltavx1**2+deltavy1**2) +
     $     am2*(deltavx2**2+deltavy2**2)) - am1*am2/r
      ltotcheck = am1*(deltax1*deltavy1 - deltay1*deltavx1) +
     $     am2*(deltax2*deltavy2 - deltay2*deltavx2)
      ecccheck = sqrt(1.d0 + 2.d0*eorbcheck*ltotcheck**2/mu/k**2)
      semilatusrectumprime = r+ecc*(deltax2-deltax1)
      if(abs(r/sep0-1.d0).gt.1.d-8) then
         write(69,*)'jumpahead: sep0 problem'
         stop
      endif
      if(abs(eorbcheck-eorb).gt.1.d-8) then
         write(69,*)'jumpahead: eorb problem'
         stop
      endif
      if(abs(ltotcheck/lorb-1.d0).gt.1.d-8) then
         write(69,*)'jumpahead: angular momentum problem'
         stop
      endif
      if(abs(ecccheck/ecc-1.d0).gt.1.d-8) then
         write(69,*)'jumpahead: eccentricity problem'
         stop
      endif
      if(abs(semilatusrectumprime/semilatusrectum-1.d0).gt.1.d-8)then
         write(69,*)'jumpahead: semilatusrectum problem'
         stop
      endif

      if(eorb.ge.0.d0) then
         vinf=(2.d0*eorb/mu)**0.5d0
         if(myrank.eq.0) write(69,*)'v_infinity=',vinf,'(code units)=',
     $        vinf*(gravconst*munit/runit)**0.5d0/1.d5,'km/s'
      endif

      if(myrank.eq.0)write(69,'(a,4g14.6)')
     $     'jumpahead: angular momentum hat of binary=',
     $     lorbx/lorb,lorby/lorb,lorbz/lorb
      costh1=lorbz/lorb ! dot product of lhat and zhat is cos(theta1)
      sinth1=(1.d0-costh1**2)**0.5d0
      if(myrank.eq.0) then
         write(69,*)'theta1=',acos(costh1)
         write(69,*)'cos(theta1),sin(theta1)=',costh1,sinth1
         write(69,*)'cos(theta1)**2+sin(theta1)**2-1=',costh1**2+sinth1**2-1.d0
         write(69,*)'sin(acos(cos(theta1)/sin(theta1)-1=',sin(acos(costh1))/sinth1-1.d0
      endif
      if(dabs(sin(acos(costh1))/sinth1-1.d0).gt.1.d-6 .or.
     $     dabs(costh1**2+sinth1**2-1.d0).gt.1.d-6) then
         write(69,*) 'problem setting theta1'
         stop
      endif

c     lon=line of nodes (its z component is zero)
c     the line of nodes is perpendicular to both the z direction and the
c     angular momentum vector.  so get lon by taking (z hat)x(angular mom),
c     that is, the cross product of zhat and the binary's angular momentum
      lonx=-lorby
      lony= lorbx
      lontot=(lonx**2+lony**2)**0.5d0
      lonhatx=lonx/lontot
      lonhaty=lony/lontot

c     the angle phi1 is between the x axis and the line of nodes.
c     so we can get cos(phi1) from the dot product between these
c     vectors... because lonz=0, the result is particularly simple:
      cosph1=lonhatx
c     and we can get sin(phi1) from the cross product between these
c     vectors... because lonz=0, the result is particularly simple:
      sinph1=lonhaty
      if(myrank.eq.0) write(69,*)'cos(phi1),sin(phi1)=',cosph1,sinph1
      if(dabs(cosph1**2+sinph1**2-1.d0).gt.1.d-15) then
         write(69,*) 'problem setting phi1'
         stop
      endif

c     the angle psi1 is determined by making sure the eccentricity vector
c     is not going to be changed by jumpahead.  in other words, the lambda
c     matrix below when applied to the vector (-ecc,0,0) should give the
c     eccentricity vector (eccx,eccy,eccz) determined above.

      sinps1=-eccz/(ecc*sinth1)
      cosps1=-(eccx*cosph1+eccy*sinph1)/ecc
      if(dabs(cosps1**2+sinps1**2-1.d0).gt.1.d-4) then
         if(myrank.eq.0) then
            write(69,*)'poor cos(psi1),sin(psi1)=',cosps1,sinps1
            write(69,*) 'there would have been a problem setting psi1',
     $           cosps1**2+sinps1**2-1.d0
         endif
      endif
      sinps1=sign((1d0-cosps1**2)**0.5d0,sinps1)
      if(myrank.eq.0) write(69,*)'cos(psi1),sin(psi1)=',cosps1,sinps1

c     determine coefficients of rotation matrix for the binary
c     the following comes from (11.99) of marion and thornton;
c     note, however, that we need the inverse of lambda to go
c     from the body system to the fixed system.  thus, phi and
c     psi are swapped in (11.99), and all the angles pick up a
c     minus sign.
      lambda11= cosph1*cosps1-costh1*sinps1*sinph1
      lambda21= sinph1*cosps1+costh1*sinps1*cosph1
      lambda31= sinth1*sinps1
      lambda12=-cosph1*sinps1-costh1*cosps1*sinph1
      lambda22=-sinph1*sinps1+costh1*cosps1*cosph1
      lambda32= sinth1*cosps1
      lambda13= sinph1*sinth1
      lambda23=-cosph1*sinth1
      lambda33= costh1

      if(myrank.eq.0) then
         write(69,*)'the rotation matrix has been determined:'
         write(69,*)lambda11,lambda12,lambda13
         write(69,*)lambda21,lambda22,lambda23
         write(69,*)lambda31,lambda32,lambda33
         write(69,*)
      endif

      dxold1=deltax1
      dyold1=deltay1
      dzold1=0.d0
      dvxold1=deltavx1
      dvyold1=deltavy1
      dvzold1=0.d0
      dxold2=deltax2
      dyold2=deltay2
      dzold2=0.d0
      dvxold2=deltavx2
      dvyold2=deltavy2
      dvzold2=0.d0
      dxnew1= lambda11*dxold1  + lambda12*dyold1  + lambda13*dzold1
      dynew1= lambda21*dxold1  + lambda22*dyold1  + lambda23*dzold1
      dznew1= lambda31*dxold1  + lambda32*dyold1  + lambda33*dzold1
      dvxnew1=lambda11*dvxold1 + lambda12*dvyold1 + lambda13*dvzold1
      dvynew1=lambda21*dvxold1 + lambda22*dvyold1 + lambda23*dvzold1
      dvznew1=lambda31*dvxold1 + lambda32*dvyold1 + lambda33*dvzold1
      dxnew2= lambda11*dxold2  + lambda12*dyold2  + lambda13*dzold2
      dynew2= lambda21*dxold2  + lambda22*dyold2  + lambda23*dzold2
      dznew2= lambda31*dxold2  + lambda32*dyold2  + lambda33*dzold2
      dvxnew2=lambda11*dvxold2 + lambda12*dvyold2 + lambda13*dvzold2
      dvynew2=lambda21*dvxold2 + lambda22*dvyold2 + lambda23*dvzold2
      dvznew2=lambda31*dvxold2 + lambda32*dvyold2 + lambda33*dvzold2

      amass1=0
      amass2=0

      do i=1,n
         if(icomp(i).eq.1) then
            amass1=amass1+1
            if(ndisplace.eq.0) then
               x(i)=x(i)+dxnew1-x1+xcm
               y(i)=y(i)+dynew1-y1+ycm
               z(i)=z(i)+dznew1-z1+zcm
            else
c     The idea here is that if ndisplace=1 then the point particle should
c     be displaced in such a way that the SPH particles remain near the origin.
               x(i)=x(i)+dxnew1-x1 - (dxnew2-x2)
               y(i)=y(i)+dynew1-y1 - (dynew2-y2)
               z(i)=z(i)+dznew1-z1 - (dznew2-z2)
            endif
            vx(i)=vx(i)+dvxnew1-vx1+vxcm
            vy(i)=vy(i)+dvynew1-vy1+vycm
            vz(i)=vz(i)+dvznew1-vz1+vzcm
         elseif(icomp(i).eq.2) then
            amass2=amass2+1
            if(ndisplace.eq.0) then
c     The idea here is that if ndisplace=1 then the SPH particles should just
c     have their coordinates stay near zero.
               x(i)=x(i)+dxnew2-x2+xcm
               y(i)=y(i)+dynew2-y2+ycm
               z(i)=z(i)+dznew2-z2+zcm
            endif
            vx(i)=vx(i)+dvxnew2-vx2+vxcm
            vy(i)=vy(i)+dvynew2-vy2+vycm
            vz(i)=vz(i)+dvznew2-vz2+vzcm
         endif
      enddo
      if(ndisplace.ne.0) then
c     The idea is that displace{x,y,z} are the components of the displacement vector that
c     keep track of how far particles had been displaced:
c     (x-coordinate if ndisplace had been zero) = (x-coordinate in code if ndisplace is 1) + displacex
c     (y-coordinate if ndisplace had been zero) = (y-coordinate in code if ndisplace is 1) + displacey
c     (z-coordinate if ndisplace had been zero) = (z-coordinate in code if ndisplace is 1) + displacez
c     If ndisplace had been zero, the icomp(2) particles would have been shifted by 
c     (dxnew2-x2+xcm, dynew2-y2+ycm, dznew2-z2+zcm),
c     but since they have not been shifted at all when ndisplace=0, we should record the change that otherwise
c     would have happened:
         displacex=displacex+dxnew2-x2+xcm
         displacey=displacey+dynew2-y2+ycm
         displacez=displacez+dznew2-z2+zcm
      endif


      if(myrank.eq.0)write(69,*)
     $     '***analyze system right after jump ahead:***'
      call rho_and_h
      call gravforce
      if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
         if(nusegpus.eq.1)then
            call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
         else
            call get_gravity_using_cpus
         endif
         mygravlength=ngrav_upper-ngrav_lower+1
         if(myrank.ne.0)then
            call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
     $           grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $           comm_worker, ierr)
         else
            call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $           grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $           comm_worker, ierr) 
         endif
      endif
      call enout(.false.)
      if(myrank.eq.0) write(69,*)
     $     '***done analyzing system right after jump ahead***'

      if(throwaway)then
         if(myrank.eq.0) write(69,*)'now we throw away some particles:'
         ntot=0
         amass1=0
         amass2=0
         thrownawaymass=0.d0

         do i=1,n
            if(u(i).eq.0) bhcomp=icomp(i)
         enddo

         do i=1,n
            if(  (bhcomp.eq.2 .and. icomp(i).eq.1) .or.
     $           (bhcomp.eq.1 .and. icomp(i).eq.2) .or. u(i).eq.0)then
               if(icomp(i).eq.1) then
                  amass1=amass1+1
               elseif(icomp(i).eq.2) then
                  amass2=amass2+1
               endif
               ntot=ntot+1
               x(ntot)=x(i)
               y(ntot)=y(i)
               z(ntot)=z(i)
               vx(ntot)=vx(i)
               vy(ntot)=vy(i)
               vz(ntot)=vz(i)
               am(ntot)=am(i)
               hp(ntot)=hp(i)
               rho(ntot)=rho(i)
               vxdot(ntot)=vxdot(i)
               vydot(ntot)=vydot(i)
               vzdot(ntot)=vzdot(i)
               u(ntot)=u(i)
               udot(ntot)=udot(i)
               gx(ntot)=gx(i)
               gy(ntot)=gy(i)
               gz(ntot)=gz(i)
               grpot(ntot)=grpot(i)
               meanmolecular(ntot)=meanmolecular(i)
               cc(ntot)=cc(i)
               divv(ntot)=divv(i)
               icomp(ntot)=icomp(i)
            else
               thrownawaymass=thrownawaymass+am(i)
            endif
         enddo
         if(myrank.eq.0) then
            write(69,*)'mass thrownaway=',thrownawaymass
            write(69,*)'new ntot=',ntot
            if(ntot.lt.nnopt)
     $           write(69,*)'not enough particles to continue'        
         endif

         if(ntot.lt.nnopt) stop            

         n=ntot
         n_upper=ntot

         if(myrank.eq.0)write(69,*)
     $        '***analyze system right after throw away:***'
         call gravquant
         call rho_and_h
         call gravforce
         if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
            if(nusegpus.eq.1)then
               call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
            else
               call get_gravity_using_cpus
            endif
            mygravlength=ngrav_upper-ngrav_lower+1
            if(myrank.ne.0)then
               call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
     $              grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $              comm_worker, ierr)
            else
               call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
     $              grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
     $              comm_worker, ierr)
            endif
         endif
         call enout(.false.)
         if(myrank.eq.0)write(69,*)
     $        '***done analyzing system right after throw away***'

      endif

      if(myrank.eq.0) then
         write (69,*) 'jumpahead:          ... done'

         open(30,file='m1m2rp.sph')
c     make file m1m2rp.sph with masses and radius
         write(30,*) amass1,amass2,rp
         write(30,*) am1,am2,rp
         close(30)
      endif

      call lfstart
      t=t+0.5d0*dt

      return
      end
