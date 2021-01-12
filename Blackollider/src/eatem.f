      subroutine eatem
      include 'starsmasher.h'
      integer i
      real*8 divv(nmax)
      common/commdivv/divv
      real*8 thrownawaymass,ribh2,coeff
      real*8 thrownawaymx,thrownawaymy,thrownawaymz
      real*8 thrownawaypx,thrownawaypy,thrownawaypz
      real*8 thrownawayfx,thrownawayfy,thrownawayfz
      integer no,thrownawaynum
      real*8 x1,y1,z1,vx1,vy1,vz1,am1
      real*8 x2,y2,z2,vx2,vy2,vz2,am2
      real*8 x3,y3,z3,vx3,vy3,vz3,am3,am4
      integer icomp(nmax)
      common/compbettercom3/am1,x1,y1,z1,vx1,vy1,vz1,
     $     am2,x2,y2,z2,vx2,vy2,vz2,
     $     am3,x3,y3,z3,vx3,vy3,vz3,am4,
     $     icomp
      real*8 thermaltime(nmax)
      common/tt/ thermaltime

      no=ntot
      ntot=0
      thrownawaymass=0.d0
      thrownawayfx=0.d0
      thrownawayfy=0.d0
      thrownawayfz=0.d0
      thrownawaypx=0.d0
      thrownawaypy=0.d0
      thrownawaypz=0.d0
      thrownawaymx=0.d0
      thrownawaymy=0.d0
      thrownawaymz=0.d0
      thrownawaynum=0
      if(u(no).ne.0.d0)then
         write(69,*)'code assumes that last particle is a black hole...'
         write(69,*)'stopping'
         stop
      endif
      do i=1,no-1               ! no-1: do not include the black hole yet
         ribh2=(x(i)-x(no))**2+(y(i)-y(no))**2+(z(i)-z(no))**2
         if(ribh2.ge.reat**2) then
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
            if(ncooling.ne.0) then
               ueq(ntot)=ueq(i)
               tthermal(ntot)=tthermal(i)
               thermaltime(ntot)=thermaltime(i)
            endif
         else
            thrownawaymass=thrownawaymass+am(i)
            thrownawayfx=thrownawayfx+am(i)*vxdot(i)
            thrownawayfy=thrownawayfy+am(i)*vydot(i)
            thrownawayfz=thrownawayfz+am(i)*vzdot(i)
            thrownawaypx=thrownawaypx+am(i)*vx(i)
            thrownawaypy=thrownawaypy+am(i)*vy(i)
            thrownawaypz=thrownawaypz+am(i)*vz(i)
            thrownawaymx=thrownawaymx+am(i)*x(i)
            thrownawaymy=thrownawaymy+am(i)*y(i)
            thrownawaymz=thrownawaymz+am(i)*z(i)
            thrownawaynum=thrownawaynum+1
         endif
      enddo

c     always keep the black hole:
      i=no
      ntot=ntot+1
      coeff=1.d0/(1.d0+thrownawaymass/am(i))
      x(ntot)=coeff*x(i)+thrownawaymx/(am(i)+thrownawaymass)
      y(ntot)=coeff*y(i)+thrownawaymy/(am(i)+thrownawaymass)
      z(ntot)=coeff*z(i)+thrownawaymz/(am(i)+thrownawaymass)
      vx(ntot)=coeff*vx(i)+thrownawaypx/(am(i)+thrownawaymass)
      vy(ntot)=coeff*vy(i)+thrownawaypy/(am(i)+thrownawaymass)
      vz(ntot)=coeff*vz(i)+thrownawaypz/(am(i)+thrownawaymass)
      vxdot(ntot)=coeff*vxdot(i)+thrownawayfx/(am(i)+thrownawaymass)
      vydot(ntot)=coeff*vydot(i)+thrownawayfy/(am(i)+thrownawaymass)
      vzdot(ntot)=coeff*vzdot(i)+thrownawayfz/(am(i)+thrownawaymass)
      am(ntot)=am(i)+thrownawaymass
      hp(ntot)=hp(i)
      rho(ntot)=rho(i)
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
      if(ncooling.ne.0)then
         ueq(ntot)=0d0
         tthermal(ntot)=0d0
         thermaltime(ntot)=0d0
      endif

      n=ntot

      if(thrownawaynum.gt.0) then
         if(myrank.eq.0) then
            write(69,*)'now we throw away some particles:'
            write(69,*)'black hole is located at',x(ntot),y(ntot),z(ntot)
            write(69,*)'black hole has 2*h=',2*hp(ntot)
            write(69,*)'mass thrownaway=',thrownawaymass,thrownawaynum
            write(69,*)'new ntot=',ntot
            write(43,*) t, thrownawaymass, thrownawaynum
         endif
         call gravquant
      endif
      
      return
      end
