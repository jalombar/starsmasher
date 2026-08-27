      subroutine relax
      include 'starsmasher.h'
      integer i
      real*8 trelaxuse,trelaxmin,trelaxmax
      real*8 trelaxeff,amtotauto,r2maxauto
      logical autodone
      save autodone
      data autodone /.false./

c     note: trelax and other parameters must be set-up previously.
c     nrelax=0 not a relaxation calculation                    
c     nrelax=1 for simple drag force
c     nrelax=2 for drag + centrifugal
c     nrelax=3 for drag + centrifugal + coriolis 
c     separation (set-up binary configuration)

c      write(69,*)'entering relax with nrleax=',nrelax,gonedynamic

c     trelax=0 in sph.input asks for the drag schedule to be taken from the model
c     rather than set by hand.  It is derived here, on the first call, rather than
c     in init.f: initialize_parent runs a full force evaluation -- and so reaches
c     this routine -- before it returns to init, so anything set there would arrive
c     one step late.  Doing it here also works for every initializer without
c     touching any of them.
c
c     From the dynamical time t_dyn = sqrt(R^3/M) in code units (G=1), with R the
c     outermost particle radius:
c
c       trelax  = 4 t_dyn   The fundamental oscillation period, which is the drag
c                           timescale that damps it fastest, scales with t_dyn.
c                           Measured on two stars that bracket the intended use:
c                           a 5.4 Msun AGB giant gives 3.91 t_dyn and a 0.4 Msun
c                           M dwarf 4.43, so 4 is within 10% of both.  The spread
c                           is real -- it reflects central concentration and the
c                           run of Gamma_1 -- but a 10% error in trelax moves the
c                           relaxed model's residual motion by only 2-3%, which is
c                           why a single constant is good enough here.
c       treloff = 10 trelax The oscillation envelope decays with an e-folding time
c                           equal to trelax itself (measured 10558 and 10672 in two
c                           independent runs), so ten of them buys ten e-foldings.
c
c     For a star whose envelope is soft enough that Gamma_1 approaches 4/3 the
c     period lengthens sharply and 4 t_dyn will underestimate it.  If that is
c     suspected, measure the period directly: relax briefly with trelax=1.d30 and
c     take the dominant frequency of the kinetic energy, which oscillates at twice
c     the mode frequency.
      if(trelax.eq.0.d0 .and. .not.autodone) then
         amtotauto=0.d0
         r2maxauto=0.d0
         do i=1,ntot
            amtotauto=amtotauto+am(i)
            r2maxauto=max(r2maxauto,x(i)**2+y(i)**2+z(i)**2)
         enddo
         tdynauto=sqrt(sqrt(r2maxauto)**3/amtotauto)
         trelax0auto=4.d0*tdynauto
         treloff=10.d0*trelax0auto
         autodone=.true.
         if(myrank.eq.0) then
            write(69,*)'relax: trelax=0, so the drag schedule is derived:'
            write(69,*)'relax:   radius  =',sqrt(r2maxauto)
            write(69,*)'relax:   mass    =',amtotauto
            write(69,*)'relax:   t_dyn   =',tdynauto
            write(69,*)'relax:   trelax  =',trelax0auto
            write(69,*)'relax:   treloff =',treloff
            if(t.gt.0.d0) then
               write(69,*)'relax: NOTE resuming at t=',t,': the schedule was'
               write(69,*)'relax: re-derived from the current star, which has'
               write(69,*)'relax: relaxed and so is slightly smaller than the'
               write(69,*)'relax: one it started from.'
            endif
         endif
      endif

      trelaxeff=trelax
      if(trelax.eq.0.d0) trelaxeff=trelax0auto

      if (nrelax.eq.1) then
         trelaxuse=trelaxeff
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
            vxdot(i)=vxdot(i)-vx(i)/trelaxeff+omega2*x(i)                     
            vydot(i)=vydot(i)-vy(i)/trelaxeff+omega2*y(i)                     
            vzdot(i)=vzdot(i)-vz(i)/trelaxeff
         enddo
      else if (nrelax.eq.3) then
         if(.not. gonedynamic) call getomega2
c     add in centrifugal and coriolis acceleration, since velocity
c     is measured in corotating frame
         do i=1,ntot
            vxdot(i)=vxdot(i)-vx(i)/trelaxeff+omega2*x(i)+2.d0*omeg*vy(i)
            vydot(i)=vydot(i)-vy(i)/trelaxeff+omega2*y(i)-2.d0*omeg*vx(i)
            vzdot(i)=vzdot(i)-vz(i)/trelaxeff
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
