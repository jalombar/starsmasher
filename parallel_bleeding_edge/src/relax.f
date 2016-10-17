      subroutine relax
      include 'starsmasher.h'
      integer i
      real*8 trelaxuse,trelaxmin,trelaxmax

c     note: trelax and other parameters must be set-up previously.
c     nrelax=0 not a relaxation calculation                    
c     nrelax=1 for simple drag force
c     nrelax=2 for drag + centrifugal
c     nrelax=3 for drag + centrifugal + coriolis 
c     separation (set-up binary configuration)

c      write(69,*)'entering relax with nrleax=',nrelax,gonedynamic

      if (nrelax.eq.1) then
         trelaxuse=trelax
         trelaxmin=1.d30
         trelaxmax=-1.d30

         do i=1,ntot
            if(trelax.lt.0.d0) then
               trelaxuse=-trelax/rho(i)**0.5d0
               trelaxmin=min(trelaxmin,trelaxuse)
               trelaxmax=max(trelaxmax,trelaxuse)
            endif
            vxdot(i)=vxdot(i)-(vx(i)+omega_spin*y(i))/trelaxuse
            vydot(i)=vydot(i)-(vy(i)-omega_spin*x(i))/trelaxuse
            vzdot(i)=vzdot(i)-vz(i)/trelaxuse
         enddo
         if(trelax.lt.0.d0 .and. myrank.eq.0)
     $        write(69,*)'trelax min, max:',trelaxmin, trelaxmax
      else if (nrelax.eq.2) then
         if(.not. gonedynamic) call getomega2
c     add in centrifugal acceleration, since velocity is measured in 
c     corotating frame
         do i=1,ntot
            vxdot(i)=vxdot(i)-vx(i)/trelax+omega2*x(i)                     
            vydot(i)=vydot(i)-vy(i)/trelax+omega2*y(i)                     
            vzdot(i)=vzdot(i)-vz(i)/trelax
         enddo
      else if (nrelax.eq.3) then
         if(.not. gonedynamic) call getomega2
c     add in centrifugal and coriolis acceleration, since velocity
c     is measured in corotating frame
         do i=1,ntot
            vxdot(i)=vxdot(i)-vx(i)/trelax+omega2*x(i)+2.d0*omeg*vy(i)
            vydot(i)=vydot(i)-vy(i)/trelax+omega2*y(i)-2.d0*omeg*vx(i)
            vzdot(i)=vzdot(i)-vz(i)/trelax
         enddo
      endif

      return
      end

      subroutine getomega2
      include 'starsmasher.h'
      real*8 momentofinertia, numerator
      integer i

c     calculate omega so that the total net force (pressure+
c     gravity+centrifugal) on the centers of mass is zero)
      momentofinertia=0.d0
      numerator=0.d0
      do i=1,ntot
         momentofinertia=momentofinertia+am(i)*(x(i)**2+y(i)**2)
         numerator=numerator+am(i)*(x(i)*vxdot(i)+y(i)*vydot(i))
      enddo
      omega2=-numerator/momentofinertia
      omeg=sqrt(omega2)
      if(myrank.eq.0) write (69,'(2(a,g15.7))') 'omega2=',omega2
      
      end
