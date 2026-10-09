c
c     update timestep
c***********************************************************************
      subroutine tstep
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 dtiacc,dtivel,dtvelmin(2),dtaccmin(2),dtiacc4,dtacc4min(2)
      real*8 dtvelcomin(2),dtacccomin(2)
      real*8 mydt
      integer i,j,irank
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      real*8 rij,vij,vdotij
      real*8 dtumin(2),dtiu
      real*8 mydtvel(2),mydtacc(2),mydtu(2),mydtacc4(2)
      real*8 mydtvelco(2),mydtaccco(2)
      integer idtvel,idtacc,idtu,idtacc4, ierr
      integer idtvelco,idtaccco
c     block steps: smallest value of each column, and its particle, over the
c     substeps since the last full synchronization
      real*8 dtsblk(7),dtsnow(7)
      integer indblk(6),indnow(6),k
      save dtsblk,indblk
      data dtsblk/7*1.d30/, indblk/6*0/
      
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

      if(nblock.eq.1) then
         do i=1,ntot
            dtpart(i)=1.d30
         enddo
      endif
      do i=n_lower,n_upper
         if(u(i).ne.0.d0) then
c            dtivel=cn1*hp(i)/sqrt(uijmax(i))
            dtivel=cn1*hp(i)/sqrt(max(uijmax(i),
     $           (vx(i)**2+vy(i)**2+vz(i)**2)/16))
            dtiacc=cn2*hp(i)**0.5d0/((vxdot(i)-vxdotsm(i))**2
     $           +(vydot(i)-vydotsm(i))**2
     $           +(vzdot(i)-vzdotsm(i))**2)**0.25d0
            if(ncooling.eq.0) then
               if(nintvar.eq.3) then
c     u(i) holds ln A here, so u(i)/|udot(i)| is not a relative change: the
c     value of a logarithm depends on the units its argument was measured in,
c     and shifting them would silently rescale the timestep.  d(lnA) is already
c     the fractional change in A, so limit that directly.
                  dtiu=cn3/dabs(udot(i))
               else
                  dtiu=cn3*u(i)/dabs(udot(i))
               endif
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
            if(nblock.eq.1) dtpart(i)=min(dtpart(i),(1.d0/dtivel
     $           +1.d0/dtiacc+1.d0/dtiu+1.d0/dtiacc4)**(-1.d0))
         endif
      enddo


c     point particles (u=0): the pairwise limits cn5, cn6 (and the softening
c     term cn7) against every other particle.  Each rank takes the pairs whose
c     partner j it owns, so the work is shared instead of falling on the rank
c     that owns the point particle, and the minima are combined below.  hp of
c     a point particle is broadcast from its owner first, since otherwise only
c     the gravity processes are sure to hold its current value.
      do i=1,ntot
         if(u(i).ne.0.d0) cycle
         do irank=0,nprocs-1
            if(i.gt.displs(irank+1) .and.
     $           i.le.displs(irank+1)+recvcounts(irank+1)) then
               call mpi_bcast(hp(i),1,mpi_double_precision,irank,
     $              mpi_comm_world,ierr)
               exit
            endif
         enddo
         do j=n_lower,n_upper
            if(j.ne.i)then
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
c     block steps: the pairwise limit protects both partners, since gas
c     particles have no acceleration criterion of their own (cn2=cn4=1e30)
               if(nblock.eq.1) then
                  dtpart(i)=min(dtpart(i),
     $                    (1.d0/dtivel+1.d0/dtiacc)**(-1.d0))
                  dtpart(j)=min(dtpart(j),
     $                    (1.d0/dtivel+1.d0/dtiacc)**(-1.d0))
               endif
            endif
         enddo
      enddo


c     mpi sync here dt should be min for all processes
      call mpi_allreduce(mydt,dt,1,mpi_double_precision,mpi_min, 
     $      mpi_comm_world,ierr)
      if(dtforce.gt.0.d0) dt=dtforce
      if(nblock.eq.1) call mpi_allreduce(mpi_in_place,dtpart,ntot,
     $     mpi_double_precision,mpi_min,mpi_comm_world,ierr)

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
         dtsnow=(/dtvelmin(1),dtaccmin(1),dtumin(1),dtacc4min(1),
     $        dtvelcomin(1),dtacccomin(1),dt/)
         indnow=(/idtvel,idtacc,idtu,idtacc4,idtvelco,idtaccco/)
c     With block steps tstep runs every substep, so one pair of lines per
c     call would flood the log.  Keep the smallest value of each column (and
c     the particle that set it) and write them once per full synchronization.
         if(nblock.eq.1) then
            do k=1,6
               if(dtsnow(k).lt.dtsblk(k)) then
                  dtsblk(k)=dtsnow(k)
                  indblk(k)=indnow(k)
               endif
            enddo
            dtsblk(7)=min(dtsblk(7),dtsnow(7))
            if(.not.blksync) return
            dtsnow=dtsblk
            indnow=indblk
            dtsblk=1.d30
            indblk=0
         endif
         write(69,'(a4,9g10.3)') 'dts=',dtsnow
         write(69,'(a4,9g10.3)') 'indx',indnow
      endif

      return
      end
