c
c     update timestep
c***********************************************************************
      subroutine tstep
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 dtiacc,dtivel,dtvelmin(2),dtaccmin(2),dtiacc4,dtacc4min(2)
      real*8 dtvelcomin(2),dtacccomin(2)
      real*8 mydt
      integer i,j
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      real*8 rij,vij,vdotij
      real*8 dtumin(2),dtiu
      real*8 mydtvel(2),mydtacc(2),mydtu(2),mydtacc4(2)
      real*8 mydtvelco(2),mydtaccco(2)
      integer idtvel,idtacc,idtu,idtacc4, ierr
      integer idtvelco,idtaccco
      
      if(nrelax.ne.1) then
         call vdotsm
      endif
      mydtvel(1)=1.d30          ! dt1 (on this myrank process)
      mydtacc(1)=1.d30          ! dt2 (on this myrank process)
      mydtu(1)=1.d30            ! dt3 (on this myrank process)
      mydtacc4(1)=1.d30         ! dt4 (on this myrank process)
      mydtvelco(1)=1.d30        ! dt5 (on this myrank process)
      mydtaccco(1)=1.d30        ! dt6 (on this myrank process)
      mydtvel(2)=0              ! particle index for dominate dt1 particle
      mydtacc(2)=0              ! particle index for dominate dt2 particle
      mydtu(2)=0                ! particle index for dominate dt3 particle
      mydtacc4(2)=0             ! particle index for dominate dt4 particle
      mydtvelco(2)=0            ! particle index for dominate dt5 particle
      mydtaccco(2)=0            ! particle index for dominate dt6 particle
      mydt=1.d30                ! overall minimum dt (on this myrank process)

      do i=n_lower,n_upper
         if(u(i).ne.0.d0) then
            dtivel=cn1*hp(i)/sqrt(uijmax(i))
            dtiacc=cn2*hp(i)**0.5d0/((vxdot(i)-vxdotsm(i))**2
     $           +(vydot(i)-vydotsm(i))**2
     $           +(vzdot(i)-vzdotsm(i))**2)**0.25d0
            if(ncooling.eq.0) then
               dtiu=cn3*u(i)/dabs(udot(i))
            else
c               dtiu=cn3*u(i)/dabs(udot(i) + (ueq(i)-u(i))*(1-exp(-dth/tthermal(i))/dth )
               dtiu=-cn3*u(i)/(udot(i) + (ueq(i)-u(i))/tthermal(i))
               if(dtiu.le.0d0) dtiu=1d30
            endif
            dtiacc4=cn4*sqrt(uijmax(i))/((vxdot(i)-vxdotsm(i))**2
     $           +(vydot(i)-vydotsm(i))**2
     $           +(vzdot(i)-vzdotsm(i))**2)**0.5d0
            if(dtivel.lt.mydtvel(1)) then
               mydtvel(1)=dtivel
               mydtvel(2)=i
            endif
            if(dtiacc.lt.mydtacc(1)) then
               mydtacc(1)=dtiacc
               mydtacc(2)=i
            endif
            if(dtiu.lt.mydtu(1)) then
               mydtu(1)=dtiu
               mydtu(2)=i
            endif
            if(dtiacc4.lt.mydtacc4(1)) then
               mydtacc4(1)=dtiacc4
               mydtacc4(2)=i
            endif
            mydt=min(mydt,(1.d0/dtivel+1.d0/dtiacc+1.d0/dtiu
     $            +1.d0/dtiacc4)**(-1.d0))
         else
            do j=1,ntot
               if(j.ne.i)then

c     Could consider changing hp(i) to hp(j) in the next line, but if we
c     do that then the smoothing lengths will need to be shared over all
c     processes in advance.f90.  Currently, smoothing lengths are shared
c     only to the gravity processes since they are the only processes
c     that need them.
                  rij=((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
     $                 +cn7*hp(i)**2)**0.5d0
                  vij= ((vx(i)-vx(j))**2
     $                 +(vy(i)-vy(j))**2
     $                 +(vz(i)-vz(j))**2)**0.5d0
                  vdotij= ((vxdot(i)-vxdot(j))**2
     $                 +(vydot(i)-vydot(j))**2
     $                 +(vzdot(i)-vzdot(j))**2)**0.5d0
                  dtivel=cn5*rij/vij
                  dtiacc=cn6*(rij/vdotij)**0.5d0
                  
                  if(dtivel.lt.mydtvelco(1)) then
                     mydtvelco(1)=dtivel
                     mydtvelco(2)=j
                  endif
                  if(dtiacc.lt.mydtaccco(1)) then
                     mydtaccco(1)=dtiacc
                     mydtaccco(2)=j
                  endif
                  mydt=min(mydt,(1.d0/dtivel+1.d0/dtiacc)**(-1.d0))
               endif
            enddo
         endif
      enddo


c     mpi sync here dt should be min for all processes
      call mpi_allreduce(mydt,dt,1,mpi_double_precision,mpi_min, 
     $      mpi_comm_world,ierr)

      call mpi_reduce(mydtvel,dtvelmin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mydtacc,dtaccmin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mydtu,dtumin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mydtacc4,dtacc4min,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mydtvelco,dtvelcomin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)
      call mpi_reduce(mydtaccco,dtacccomin,1,mpi_2double_precision,
     $     mpi_minloc,0,mpi_comm_world,ierr)

      if(myrank.eq.0) then
         idtvel=nint(dtvelmin(2))
         idtacc=nint(dtaccmin(2))
         idtu=nint(dtumin(2))
         idtacc4=nint(dtacc4min(2))
         idtvelco=nint(dtvelcomin(2))
         idtaccco=nint(dtacccomin(2))
         write(69,'(a4,9g10.3)')
     $    'dts=',dtvelmin(1),dtaccmin(1),dtumin(1),dtacc4min(1),
     $        dtvelcomin(1),dtacccomin(1),dt
         write(69,'(a4,9g10.3)')
     $    'indx',idtvel,idtacc,idtu,idtacc4,idtvelco,idtaccco
      endif

      return
      end
