c *****************************************
c 
c starsmasher
c 
c primary developer: jamie lombardi (jamie.lombardi@allegheny.edu)
c 
c *****************************************

      include 'starsmasher.h'
      include 'mpif.h'

      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      integer i,ierr
      common /jumpcomm/ tjumpahead
      logical autotf
      common/autotfblock/autotf

      omeg=0.d0
      gonedynamic=.false.

c mpi initialization:
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)
      omega2=0
      call init

      if(myrank.eq.0) write(69,*) 'nprocs=',nprocs

c      if(myrank.eq.0) then
cc         open(72,file='energy.sph', status='unknown')
c         open(78,file='com.sph', status='unknown')
c      endif

c      write(69,*) "back in main", t, tf
c     main program loop:
      do while (t.le.tf)

         if(nrelax.ge.2 .and. .not.gonedynamic) then
            call getcoms
            do i=1,ntot
               if(dabs(x(i)).gt.dabs(xcm1)*2.5d0.and.
     $              dabs(x(i)).gt.dabs(xcm2)*2.5d0) then
                  write(69,*) 'particle',i,'at',
     $                 x(i),y(i),z(i),
     $                 'has overflowed an outer lagrangian point'
                  print *,'xcm1,xcm2,omeg=',xcm1,xcm2,omeg
                  write(69,*) 'will now stop run'
                  if(myrank.eq.0) then
                     close(22)
                     close(69)
                     if(reat.gt.0) close(43)
                  endif
                  call mpi_finalize(ierr)
                  stop
               endif
            enddo
         endif

         if(t.ge.treloff.and.nrelax.lt.3 .and.
     $        .not.gonedynamic .and. omeg.ne.0.d0)then
            if(myrank.eq.0) write(69,*)'main: going dynamic',dth
c            alpha=1.d0
c            beta=2.d0
            if(myrank.eq.0) write(69,*)'alpha,beta=',alpha,beta
            gonedynamic=.true.
            do i=1,ntot
               vx(i)=-omeg*(y(i)+omeg*x(i)*dth) !at half-timestep
               vy(i)= omeg*(x(i)-omeg*y(i)*dth) !at half-timestep
               vz(i)=0.d0
            enddo
            trelax=1.d30
             if(myrank.eq.0) write(69,*)'changing omeg from',omeg,'to 0'
            omega2=0.d0
            omeg=0.d0
         elseif(t.ge.treloff.and.nrelax.eq.3 .and.
     $           .not.gonedynamic)then
            if(myrank.eq.0) write(69,*)'main: going dynamic',dth
            gonedynamic=.true.
            trelax=1.d30
         elseif(t.ge.treloff.and.nrelax.eq.1 .and.
     $           .not.gonedynamic)then
            if(myrank.eq.0) write(69,*)'main: going dynamic',dth
            gonedynamic=.true.
            trelax=1.d30
         endif
         
         call mainit

         if(reat.ge.0.d0) call eatem

         if(t.gt.abs(tjumpahead)) then
            if(autotf) then
               if(myrank.eq.0) write(69,*) 'jumpping ahead at time t=',t
               call jumpahead
            else
               if(myrank.eq.0) write(69,*)
     $              'we might have jumpped ahead but decided not to'
            endif
            tjumpahead=1.d30
         endif

         if (t.ge.tf) then
            if(myrank.eq.0) then
               write (69,*) 'main: end of integration, t=',t,tf
               close(22)
               close(69)
               if(reat.gt.0) close(43)
            endif
            call mpi_finalize(ierr)
         endif
      enddo
c      call shiftboundmasstoorigin
      end
***********************************************************************
      subroutine mainit
c     main iteration
      include 'starsmasher.h'                                         
      include 'mpif.h'

      if(myrank.eq.nprocs-1)call cpu_time(time1)

c      if (mod(nit, 100) .eq. 0) then
c        write (6, '(a)')  'sorting particle in peano-hilber order'
c        call peano_sort(ntot, x, y, z, vx, vy, vz, hp, u)
c      end if
      nit=nit+1


      call advance


      if(myrank.eq.nprocs-1) then
         call cpu_time(time2)
         seconds = time2-time1
         write (6,'(a,f6.3,a)')
     $     'mainit: ',seconds
      endif

c     write checkpoint file (if necessary):
      call checkpt(1)
c     write results to out*.sph (if necessary):
      call output

      call flush(6)

      return
      end


*********************************************************
      subroutine gravquant
********************************************************
c
c   coppied from parallel code
c
c     release 1.0
c     calculate quantities for the parallel code, and check numbers
c     to make sure they work properly
c     called by init,lfstart,setup1em,setup1es,setup2cm,setup2cs
*****************************************************************

      include 'starsmasher.h'
      include 'mpif.h'

      integer stride,irank,nlow,nup
      integer gravstride,ngravlow,ngravup,i,mygravlength,logq,qtest
      integer comm_world, group_world, comm_worker, group_worker, ierr
      character*(MPI_MAX_PROCESSOR_NAME) name
      integer len
      real*8 range(ntot)
      real*8 mintime
      common/gravworkers/comm_worker
      integer ranks(nprocsmax)
      integer imax
      logical alreadyinitialized
      data alreadyinitialized /.false./
      save alreadyinitialized
      real*8 theta_angle
      real*8 numinteractions,totnuminteractions ! declare as real*8 so don't run into max integer problems

      if(myrank.eq.0)write(69,*)'called gravquant'
      if(nprocs.gt.nprocsmax)then
         write(69,*)'increase nprocsmax in starsmasher.h'
         stop
      endif
  
      if(ngravprocs.gt.nprocs)then
         if(myrank.eq.0)then
            write(69,*)'need to have ngravprocs <= nprocs'
            write(69,*)'setting ngravprocs = nprocs = ',nprocs
         endif
         ngravprocs=nprocs
      endif

      if(ngravprocs.lt.0)then
!     the following line assumes ppn cpu threads and abs(ngravprocs) gpu threads per node
         ngravprocs=min((nprocs+ppn-1)/ppn*abs(ngravprocs),nprocs)
         if(myrank.eq.0) then
            write(69,*)'because ngravprocs<0 in sph.input, we reset it:'
            write(69,*)'ngravprocs=',ngravprocs
         endif
      endif
      
      if(ngravprocs.gt.ngravprocsmax)then
         write(69,*)'increase ngravprocsmax in starsmasher.h'
         stop
      endif
     
      if(.true.)then
c      if(4*ngravprocs.ge.nprocs)then
c     let myrank<ngravprocs do both gravity and hydro:
         stride=(n-1)/nprocs+1
         n_lower=stride*myrank+1
         n_upper=min(n_lower+stride-1,n)
         do irank=0,nprocs-1
            displs(irank+1)=stride*irank
            nlow=displs(irank+1)+1
            nup=min(displs(irank+1)+stride,n)
            recvcounts(irank+1)=nup-nlow+1
            if(myrank.eq.0) write(69,*)'myrank,displ,recvcount=',
     $           irank,displs(irank+1),recvcounts(irank+1)
         enddo
      else 
c     let myrank<ngravprocs do only the gravity.
c     e.g. if ngravprocs=2 and nprocs=8 then myrank=0,1 do gravity (and no
c     hydro), while myrank=2,3,4,5,6,7 do hydro (and no gravity):
         stride=(n-1)/(nprocs-ngravprocs)+1
         if(myrank.lt.ngravprocs) then
            n_lower=1
            n_upper=0
         else
            n_lower=stride*(myrank-ngravprocs)+1
            n_upper=min(n_lower+stride-1,n)
         endif
         do irank=0,ngravprocs-1
            displs(irank+1)=0
            recvcounts(irank+1)=0
         enddo
         do irank=ngravprocs,nprocs-1
            nlow=stride*(irank-ngravprocs)+1
            nup=min(nlow+stride-1,n)
            displs(irank+1)=nlow-1
            recvcounts(irank+1)=nup-nlow+1
            if(myrank.eq.0) write(69,*)'myrank,displ,recvcount=',
     $           irank,displs(irank+1),recvcounts(irank+1)
         enddo         
      endif

!      if(myrank.lt.ngravprocs .and. computeexclusivemode.ne.1 .and. ngravprocs.gt.1 .and. .not. alreadyinitialized) then
!                                    ^^^^^^^^^^^^^                   ^^^^^^^^^^^^^^^
!                                                                     missing from below
!     the following line assumes 8 cpu threads per node
!         call gpu_init_dev(myrank/((nprocs+ppn-1)/ppn)) ! if the gpus are set up in device exclusive mode,
!                                   ! then you probably won't want to do this
!      endif

      call MPI_GET_PROCESSOR_NAME(name, len, ierr)
      if (ierr .ne. MPI_SUCCESS) then
         print *,'Error getting processor name. Terminating.'
         stop
      end if
      
      if(myrank.lt.ngravprocs .and. .not. alreadyinitialized .and. computeexclusivemode.ne.1 .and. nusegpus.eq.1) then
!          gpus must always be initialized, even if we use just 1 mpi process
!                                                                 
!     the following line assumes 8 cpu threads per node
         theta_angle = 0.3             ! 0.4 or 0.3 should be safe, experiment if 0.5 or larger works
         call gpu_init_dev(myrank/((nprocs+ppn-1)/ppn), theta_angle) ! if the gpus are set up in device exclusive mode,
                                   ! then you probably won't want to do this
         write(6,"('myrank=',I3,' is running on ',A,' with gpu',i3)")
     $        myrank,
     $        trim(name),
     $        myrank/((nprocs+ppn-1)/ppn)
      else
         write(6,"('myrank=',I3,' is running on ',A)")
     $        myrank,
     $        trim(name)
      endif
      
      if(ngravprocs.ne.0) then
         if(nusegpus.eq.1) then
            gravstride=(n-1)/ngravprocs+1 ! number of particles to be treated by a typical gravity process
            ngrav_lower=gravstride*myrank+1 ! smallest particle index treated by gravproc myrank
            ngrav_upper=min(ngrav_lower+gravstride-1,n) ! largest particle index treated by gravproc myrank
            do irank=0,ngravprocs-1
               gravdispls(irank+1)=gravstride*irank ! one less than smallest particle index treated by gravproc irank
               ngravlow=gravdispls(irank+1)+1 ! smallest particle index treated by gravproc irank
               ngravup=min(gravdispls(irank+1)+gravstride,n) ! largest particle index treated by gravproc irank
               gravrecvcounts(irank+1)=ngravup-ngravlow+1 ! number of particles treated by gravproc irank
               if(myrank.eq.0) write(69,*)
     $              'gpurank,gravdispl,gravrecvcount=',
     $              irank,gravdispls(irank+1),gravrecvcounts(irank+1)
            enddo
         else
            totnuminteractions=n*(n+1d0)/2 ! total number of gravity interactions that need followed.
                                           ! the '+' is because particle interacts with self

            if(myrank.eq.0) then
               write(69,*)'n,totnuminteractions,ngravporcs=',n,totnuminteractions,ngravprocs
               write(6,*)'n,totnuminteractions,ngravporcs=',n,totnuminteractions,ngravprocs
            endif
            irank=0
            gravdispls(irank+1)=0 ! one less than smallest particle index treated by gravproc irank=0
            numinteractions=0   ! initialize counter
            do i=1,n
               numinteractions=numinteractions+ n-i+1  ! particle i interacts iwth n-i+1 particles
               if(numinteractions+(n-i)/2.ge.(irank+1d0)*totnuminteractions/ngravprocs)then
                  ngravlow=gravdispls(irank+1)+1 ! smallest particle index treated by gravproc irank
                  ngravup=i     ! largest particle index treated by gravproc irank
                  if(myrank.eq.irank)then
                     ngrav_lower=ngravlow
                     ngrav_upper=ngravup
                  endif
                  gravrecvcounts(irank+1)=ngravup-ngravlow+1 ! number of particles treated by gravproc irank
                  if(myrank.eq.0) then
                     write(69,*)
     $                 'cpurank,gravdispl,gravrecvcount=',
     $                 irank,gravdispls(irank+1),gravrecvcounts(irank+1),numinteractions
                     write(6,*)
     $                 'cpurank,gravdispl,gravrecvcount=',
     $                 irank,gravdispls(irank+1),gravrecvcounts(irank+1),numinteractions
                  endif
                  irank=irank+1
                  gravdispls(irank+1)=ngravup
                  if(irank.eq.ngravprocs-1) goto 1235
               endif
            enddo
 1235       if(irank.ne.ngravprocs-1) then
               print *,'problem assigning particles to grav processes'
               stop
            endif
            ngravlow=gravdispls(irank+1)+1 ! smallest particle index treated by gravproc irank=ngravprocs-1
            ngravup=n           ! largest particle index treated by gravproc irank=ngravprocs-1
            if(myrank.eq.irank)then
               ngrav_lower=ngravlow
               ngrav_upper=ngravup
            endif
            gravrecvcounts(irank+1)=ngravup-ngravlow+1 ! number of particles treated by gravproc irank
            if(myrank.eq.0) then
               write(69,*)
     $              'cpurank,gravdispl,gravrecvcount=',
     $              irank,gravdispls(irank+1),gravrecvcounts(irank+1),totnuminteractions
               write(6,*)
     $              'cpurank,gravdispl,gravrecvcount=',
     $              irank,gravdispls(irank+1),gravrecvcounts(irank+1),totnuminteractions
            endif
         endif
         

         if(myrank.lt.ngravprocs .and. nusegpus.eq.1) then
            if(qthreads .lt. 0)then
               if(myrank.eq.0) write(69,*)'timing to find best q...'
               mygravlength=ngrav_upper-ngrav_lower+1
               do i=1,ntot
                  range(i) = 4*hp(i)**2;
                  if(u(i).eq.0.d0) range(i)=-range(i)
               enddo
               qtest=4
c     the first call to a gpu is usually slow, so let's not time this one:
               call firsthalf_grav_forces(ntot, ngrav_lower, mygravlength, x,y,z,
     $              am, range, qtest, nkernel)
               call lasthalf_grav_forces(ntot, gx, gy,gz, grpot)
               
               mintime=1d30
               imax=2
               do logq=0,2
                  qtest=2**logq
                  call cpu_time(time0)
                  do i=1,imax
                     call firsthalf_grav_forces(ntot, ngrav_lower, mygravlength, x,y,z,
     $                    am, range, qtest, nkernel)
                     call lasthalf_grav_forces(ntot, gx, gy,gz, grpot)
                  enddo
                  call cpu_time(time1)
                  write(6,*)myrank,'average time=',(time1-time0)/imax,
     $                 'for q=',qtest
                  if(time1-time0 .lt. mintime)then
                     mintime=time1-time0
                     q=qtest
                  endif
               enddo
            else if(qthreads .eq. 0)then
               if(myrank.eq.0) write(69,*)
     $              'automatically setting q (without timing)'
               mygravlength=ngrav_upper-ngrav_lower+1
               if(mygravlength .lt. 15000)then
                  q=8
               else if(mygravlength .lt. 60000)then
                  q=4
               else if(mygravlength .lt. 240000)then
                  q=2
               else
                  q=1
               endif
            else
               if(myrank.eq.0) write(69,*)
     $              'will set q=qthreads=',qthreads
               q=qthreads
            endif
            write(6,'(2(a,i4))')'myrank=',myrank,'will use q=',q  
         endif
         
         if(.not. alreadyinitialized) then
            comm_world = mpi_comm_world
            call mpi_comm_group(comm_world, group_world, ierr)
            
            if(ngravprocs.lt.nprocs) then
               do irank=1,nprocs-ngravprocs
                  ranks(irank)=irank+ngravprocs-1
               enddo
               call mpi_group_excl(group_world, nprocs-ngravprocs, ranks, group_worker, ierr)
               
               call mpi_comm_create(comm_world, group_worker, comm_worker, ierr)
            else
               comm_worker=mpi_comm_world
            endif
            alreadyinitialized=.true.
         endif
      endif

      end




