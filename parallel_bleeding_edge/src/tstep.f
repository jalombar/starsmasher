c
c     update timestep
c***********************************************************************
      subroutine tstep
      use omp_utils
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
c     the minima over the particles are found in blocks of fixed size and
c     combined in block order, which gives the first particle that reaches the
c     minimum whatever the number of threads (the minima themselves would be
c     the same anyway; the index goes to the log)
      integer ib,nb,nown,npoint
      real*8 bvl,bvi,bal,bai,bul,bui,b4l,b4i,bdt,bcvl,bcvi,bcal,bcai,bdtp
c     start-up ramp (see below)
      logical ramp,firstcall
      real*8 dtlast
      save ramp,firstcall,dtlast
      data ramp/.false./, firstcall/.true./, dtlast/0.d0/
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
!$omp parallel do schedule(static) private(i)
         do i=1,ntot
            dtpart(i)=1.d30
         enddo
      endif
      nown=n_upper-n_lower+1
      nb=nblk(nown,redblock_fine)
      call reserve_rpart(9,nb)
!$omp parallel do schedule(static) private(ib,i,dtivel,dtiacc,dtiu,dtiacc4,bvl,bvi,bal,bai,bul,bui,b4l,b4i,bdt)
      do ib=1,nb
      bvl=1.d30
      bvi=0.d0
      bal=1.d30
      bai=0.d0
      bul=1.d30
      bui=0.d0
      b4l=1.d30
      b4i=0.d0
      bdt=1.d30
      do i=n_lower+(ib-1)*redblock_fine,
     $     min(n_lower+ib*redblock_fine-1,n_upper)
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
            if(dtivel.lt.bvl) then
               bvl=dtivel
               bvi=i
            endif
            if(dtiacc.lt.bal) then
               bal=dtiacc
               bai=i
            endif
            if(dtiu.lt.bul) then
               bul=dtiu
               bui=i
            endif
            if(dtiacc4.lt.b4l) then
               b4l=dtiacc4
               b4i=i
            endif
            bdt=min(bdt,(1.d0/dtivel+1.d0/dtiacc+1.d0/dtiu
     $            +1.d0/dtiacc4)**(-1.d0))
            if(nblock.eq.1) dtpart(i)=min(dtpart(i),(1.d0/dtivel
     $           +1.d0/dtiacc+1.d0/dtiu+1.d0/dtiacc4)**(-1.d0))
         endif
      enddo
      rpart(1,ib)=bvl
      rpart(2,ib)=bvi
      rpart(3,ib)=bal
      rpart(4,ib)=bai
      rpart(5,ib)=bul
      rpart(6,ib)=bui
      rpart(7,ib)=b4l
      rpart(8,ib)=b4i
      rpart(9,ib)=bdt
      enddo
      do ib=1,nb
         if(rpart(1,ib).lt.mydtvel(1)) then
            mydtvel(1)=rpart(1,ib)
            mydtvel(2)=rpart(2,ib)
         endif
         if(rpart(3,ib).lt.mydtacc(1)) then
            mydtacc(1)=rpart(3,ib)
            mydtacc(2)=rpart(4,ib)
         endif
         if(rpart(5,ib).lt.mydtu(1)) then
            mydtu(1)=rpart(5,ib)
            mydtu(2)=rpart(6,ib)
         endif
         if(rpart(7,ib).lt.mydtacc4(1)) then
            mydtacc4(1)=rpart(7,ib)
            mydtacc4(2)=rpart(8,ib)
         endif
         mydt=min(mydt,rpart(9,ib))
      enddo


c     point particles (u=0): the pairwise limits cn5, cn6 (and the softening
c     term cn7) against every other particle.  Each rank takes the pairs whose
c     partner j it owns, so the work is shared instead of falling on the rank
c     that owns the point particle, and the minima are combined below.  hp of
c     a point particle is broadcast from its owner first, since otherwise only
c     the gravity processes are sure to hold its current value.
      npoint=0
!$omp parallel do schedule(static) private(i) reduction(+:npoint)
      do i=1,ntot
         if(u(i).eq.0.d0) npoint=npoint+1
      enddo
      do i=1,ntot
         if(npoint.eq.0) exit
         if(u(i).ne.0.d0) cycle
         do irank=0,nprocs-1
            if(i.gt.displs(irank+1) .and.
     $           i.le.displs(irank+1)+recvcounts(irank+1)) then
               call mpi_bcast(hp(i),1,mpi_double_precision,irank,
     $              mpi_comm_world,ierr)
               exit
            endif
         enddo
c     the partners j are shared among the threads in fixed blocks, with the
c     partial minima combined in block order
         call reserve_rpart(5,nb)
!$omp parallel do schedule(static) private(ib,j,rij,vij,vdotij,dtivel,dtiacc,bcvl,bcvi,bcal,bcai,bdtp)
         do ib=1,nb
            bcvl=1.d30
            bcvi=0.d0
            bcal=1.d30
            bcai=0.d0
            bdtp=1.d30
            do j=n_lower+(ib-1)*redblock_fine,
     $           min(n_lower+ib*redblock_fine-1,n_upper)
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
               if(dtivel.lt.bcvl) then
                  bcvl=dtivel
                  bcvi=j
               endif
               if(dtiacc.lt.bcal) then
                  bcal=dtiacc
                  bcai=j
               endif
               bdtp=min(bdtp,(1.d0/dtivel+1.d0/dtiacc)**(-1.d0))
c     block steps: the pairwise limit protects both partners, since gas
c     particles have no acceleration criterion of their own (cn2=cn4=1e30)
               if(nblock.eq.1) then
                  dtpart(j)=min(dtpart(j),
     $                    (1.d0/dtivel+1.d0/dtiacc)**(-1.d0))
               endif
            endif
            enddo
            rpart(1,ib)=bcvl
            rpart(2,ib)=bcvi
            rpart(3,ib)=bcal
            rpart(4,ib)=bcai
            rpart(5,ib)=bdtp
         enddo
         do ib=1,nb
            if(rpart(1,ib).lt.mydtvelco(1)) then
               mydtvelco(1)=rpart(1,ib)
               mydtvelco(2)=rpart(2,ib)
            endif
            if(rpart(3,ib).lt.mydtaccco(1)) then
               mydtaccco(1)=rpart(3,ib)
               mydtaccco(2)=rpart(4,ib)
            endif
            mydt=min(mydt,rpart(5,ib))
            if(nblock.eq.1) dtpart(i)=min(dtpart(i),rpart(5,ib))
         enddo
      enddo


c     mpi sync here dt should be min for all processes
      call mpi_allreduce(mydt,dt,1,mpi_double_precision,mpi_min, 
     $      mpi_comm_world,ierr)
      if(dtforce.gt.0.d0) dt=dtforce
      if(nblock.eq.1) call mpi_allreduce(mpi_in_place,dtpart,ntot,
     $     mpi_double_precision,mpi_min,mpi_comm_world,ierr)

c     start-up ramp (nstartramp=1).  A run that starts at t=0 from rest has
c     udot=0, so the u/udot limit (cn3) cannot act on the first step,
c     although du/dt grows as soon as the particles move.  Across an
c     initial discontinuity the first step would be ~10 times what cn3
c     allows one step later.  So the first step is 1/64 of what the limits
c     give, and each later step may at most double, until the limits
c     themselves are the smaller.  On by default only for new shock-tube
c     runs (see init.f), so that other runs are unchanged.  A restart has
c     t>0 at its first call and never ramps.
      if(nstartramp.eq.1 .and. dtforce.le.0.d0) then
         if(firstcall) then
            firstcall=.false.
            ramp=t.eq.0.d0
            if(ramp) dt=dt/64.d0
         elseif(ramp) then
            if(dt.gt.2.d0*dtlast) then
               dt=2.d0*dtlast
            else
               ramp=.false.
            endif
         endif
         dtlast=dt
         if(ramp .and. nblock.eq.1) then
            do i=1,ntot
               dtpart(i)=min(dtpart(i),dt)
            enddo
         endif
      endif

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
c     the particle that set it) and write them once per dtmaxblk.
         if(nblock.eq.1) then
            do k=1,6
               if(dtsnow(k).lt.dtsblk(k)) then
                  dtsblk(k)=dtsnow(k)
                  indblk(k)=indnow(k)
               endif
            enddo
            dtsblk(7)=min(dtsblk(7),dtsnow(7))
            if(.not.blkfull) return
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
