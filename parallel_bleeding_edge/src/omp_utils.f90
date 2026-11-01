!***********************************************************************
!     Helpers for the hybrid MPI x OpenMP build (make OPENMP=1).  Without
!     -fopenmp everything here compiles and behaves as a single thread, so
!     the same code serves both builds.
!
!     Deterministic reductions.  A floating-point sum depends on the order
!     of the additions, so a parallel sum whose partial sums follow the
!     threads would change with the number of threads.  Sums that are
!     parallel are therefore made of partial sums over fixed blocks of
!     redblock consecutive terms (a partition that does not depend on the
!     number of threads, nor on the schedule), added together serially in
!     block order.  Each partial sum is a serial loop in index order.  The
!     plain MPI build uses the same blocks, so both builds give the same
!     bits.  min and max do not depend on the order and use reduction().
!***********************************************************************
      module omp_utils
!$    use omp_lib
      implicit none
      integer, parameter :: redblock=4096
!     smaller blocks for the loops whose terms are expensive (tstep, enout): with 4096
!     the 10^5 particles of a typical run would be only 25 blocks to share among the threads
      integer, parameter :: redblock_fine=256
!     scratch for the partial sums of the blocks, rpart(k,block); sized by
!     reserve_rpart, which is called outside parallel regions
      real*8, allocatable :: rpart(:,:)
      logical :: openmp_active=.false.
      integer :: nthreads_used=1

      contains

!     number of blocks of redblock terms that n terms make
      integer function nredblocks(n)
      integer n
      nredblocks=(n+redblock-1)/redblock
      end function nredblocks

!     the same for blocks of b terms
      integer function nblk(n,b)
      integer n,b
      nblk=(n+b-1)/b
      end function nblk

!     make rpart at least k x nb
      subroutine reserve_rpart(k,nb)
      integer k,nb
      if(allocated(rpart)) then
         if(size(rpart,1).ge.k .and. size(rpart,2).ge.nb) return
         deallocate(rpart)
      endif
      allocate(rpart(max(k,8),max(nb,64)))
      end subroutine reserve_rpart

!     index of the calling thread and the number of threads a parallel region will have
      integer function thread_id()
      thread_id=0
!$    thread_id=omp_get_thread_num()
      end function thread_id

      integer function max_threads()
      max_threads=1
!$    max_threads=omp_get_max_threads()
      end function max_threads

!     MPI_Init, or MPI_Init_thread with MPI_THREAD_FUNNELED in the OpenMP
!     build: only the master thread calls MPI (and CUDA)
      subroutine omp_mpi_init(ierr)
      include 'mpif.h'
      integer ierr,provided
      openmp_active=.false.
!$    openmp_active=.true.
      if(openmp_active) then
         provided=-1
         call mpi_init_thread(mpi_thread_funneled,provided,ierr)
         if(provided.lt.mpi_thread_funneled) then
            write(6,*) 'StarSmasher: this MPI library does not provide MPI_THREAD_FUNNELED (provided=',provided,')'
            call mpi_abort(mpi_comm_world,1,ierr)
         endif
      else
         call mpi_init(ierr)
      endif
      end subroutine omp_mpi_init

!     Number of threads per rank.  If OMP_NUM_THREADS is not set and there
!     is more than one rank, each rank gets one thread: the default
!     (one thread per logical CPU) would make every rank use the whole node.
!     With one rank the OpenMP default stands.  An explicit OMP_NUM_THREADS
!     is always honoured.
      subroutine omp_setup(nprocs)
      integer nprocs
      character(len=64) :: val
      integer st
      nthreads_used=1
      st=0
      val=' '
!$    call get_environment_variable('OMP_NUM_THREADS',val,status=st)
!$    if(st.eq.1 .and. nprocs.gt.1) call omp_set_num_threads(1)
!$    nthreads_used=omp_get_max_threads()
      end subroutine omp_setup

!     one line with P, T and where the threads of rank 0 run, for fort.69
      subroutine omp_report(myrank,nprocs)
      integer myrank,nprocs
      integer it,nthr,np,pl
      integer, allocatable :: first(:)
      integer, allocatable :: ids(:)
      character(len=16) :: bindname
      nthr=nthreads_used
      if(myrank.ne.0) return
      if(.not.openmp_active) then
         write(69,*) 'OpenMP: not compiled in (OPENMP=0): P=',nprocs,' T=1'
         return
      endif
      allocate(first(0:nthr-1))
      first=-1
!$omp parallel private(it,pl,np,ids) num_threads(nthr)
!$    it=omp_get_thread_num()
!$    pl=omp_get_place_num()
!$    if(pl.ge.0) then
!$       np=omp_get_place_num_procs(pl)
!$       allocate(ids(np))
!$       call omp_get_place_proc_ids(pl,ids)
!$       first(it)=ids(1)
!$       deallocate(ids)
!$    endif
!$omp end parallel
      bindname='?'
!$    select case(omp_get_proc_bind())
!$    case(0)
!$       bindname='false'
!$    case(1)
!$       bindname='true'
!$    case(2)
!$       bindname='master'
!$    case(3)
!$       bindname='close'
!$    case(4)
!$       bindname='spread'
!$    end select
      write(69,'(a,i5,a,i5,a,a,a,i5,a)') ' OpenMP: P=',nprocs,' T=',nthr,' proc_bind=',trim(bindname),&
           ' places=',nplaces(),' (cpu of the first thread of each place of rank 0 below)'
      write(69,'(a,32(i4))') '   rank 0 threads -> cpu:',first(0:min(nthr-1,31))
      deallocate(first)
      contains
      integer function nplaces()
      nplaces=0
!$    nplaces=omp_get_num_places()
      end function nplaces
      end subroutine omp_report

      end module omp_utils

!***********************************************************************
      subroutine wall_time(t)
!     wall-clock seconds from an arbitrary origin, in double precision.
!     The timing lines of the log used cpu_time, which with threads adds
!     the CPU time of all of them.
      real*8 t
      integer*8 c,r
      call system_clock(c,r)
      t=dble(c)/dble(r)
      end
