!***********************************************************************
!     SPH on the GPUs (input options ngpusph and ngpusphcheck): the Fortran side of SPHgrav_lib2/sph_gpu.cu.
!
!     ngpusph is a bit mask; bit 1 (value 1) solves the smoothing lengths on the GPU and makes the neighbour lists, the
!     density and the other sums of the SPH (rho, omega, 0mega, psi, wn) there.  The results are the bits the CPU gives: the
!     same sequence of trial values and steps as solve_particle, the same sums in the same order (the neighbours in the order
!     of their index), in double precision without fused multiply-add.  A particle that the GPU does not do (more than 128
!     candidates) is solved by the CPU with the code that always did it.  Where the GPU cannot do what is asked the CPU does the
!     whole step, and a line of log0.sph says why (once, and again when the reason changes).  The reasons are fixed by the build
!     and the input, never by how long something takes:
!        1  this build has no GPU code (make cpu, the Metal build)
!        2  the CPU code of this build fuses multiplications and additions (MARCH, aarch64): its bits are not the GPU's
!        3  periodic box
!        4  dynhco>0 (the smoothing lengths of the point particles)
!        5  nblock=1 (block time steps)
!        6  this process has no GPU of its own (the node has too few)
!        7  an error of CUDA (the line before it says which); the CPU does the work from then on
!        8  the lists do not fit in list()
!        9  the first step of the run (dt=0: the guess of the smoothing length can be far off)
!     ngpusphcheck=1 makes the CPU do its work too and compares: the run stops at the first particle that differs, with the
!     quantity and the two values.
!***********************************************************************
      module gpusph_mod
      implicit none
!     0: the library has not been asked yet, 1: ready, 2: not usable (gs_why says why)
      integer, save :: gs_state=0, gs_why=0
!     the reason written last
      integer, save :: gs_said=-1
!     1 for the particles of the rank (i-n_lower+1) that the GPU did not do
      integer, allocatable, save :: gs_flag(:)
!     ngpusphcheck: what the GPU made, kept to compare with what the CPU makes
      real*8, allocatable, save :: ck_hp(:),ck_rho(:),ck_om(:),ck_om0(:),ck_psi(:),ck_wn(:),hp0(:)
      integer, allocatable, save :: ck_nn(:),ck_first(:),ck_list(:)
      integer, save :: ck_ltot=0

      contains

!     one line in log0.sph when the reason why the CPU does the work changes (0: the GPU does it)
      subroutine gpusph_say(what,why)
      include 'starsmasher.h'
      character(*) what
      integer why
      character(len=120) text
      if(why.eq.gs_said) return
      gs_said=why
      if(myrank.ne.0) return
      select case(why)
      case(0)
         text='on the GPU'
      case(1)
         text='on the CPUs: this build has no GPU code'
      case(2)
         text='on the CPUs: the CPU code of this build fuses multiplications and additions (MARCH)'
      case(3)
         text='on the CPUs: periodic box'
      case(4)
         text='on the CPUs: dynhco>0'
      case(5)
         text='on the CPUs: nblock=1'
      case(6)
         text='on the CPUs: this process has no GPU of its own'
      case(7)
         text='on the CPUs: error of CUDA (see the error log)'
      case(8)
         text='on the CPUs: the lists do not fit in list()'
      case(9)
         text='on the CPUs: the first step of the run'
      case default
         text='on the CPUs'
      end select
      write(69,'(a,a,a,a)') 'ngpusph: ',trim(what),' ',trim(text)
      call flush(69)
      end subroutine gpusph_say

!     why=0 if the GPU can be used for this step, else the number of the reason
      subroutine gpusph_usable(why)
      include 'starsmasher.h'
      integer why,ok,ierr
      if(gs_state.eq.0) then
         call gpu_sph_available(ok)
         if(ok.eq.0) then
            gs_state=2
            gs_why=1
         else if(ok.eq.2) then
            gs_state=2
            gs_why=2
         else
            call gpu_sph_setup(myrank/((nprocs+ppn-1)/ppn),ngpuperrank,ngravprocs,ierr)
            if(ierr.eq.0) call gpu_sph_tables(ntab,gtab,dgdhtab,wtab,dwdhtab,dphidhtab,ctab,ierr)
            if(ierr.eq.0) then
               gs_state=1
            else
               gs_state=2
               gs_why=7
               if(ierr.eq.1) gs_why=6
            endif
         endif
      endif
      why=0
      if(gs_state.eq.2) then
         why=gs_why
      else if(periodic) then
         why=3
      else if(dynhco.ne.0 .or. hcolim) then
         why=4
      else if(nblock.eq.1) then
         why=5
      else if(dt.le.0.d0) then
         why=9
      endif
      end subroutine gpusph_usable

!     The GPU found a problem it cannot get over: the CPU does the SPH from here on
      subroutine gpusph_fail(why)
      integer why
      gs_state=2
      gs_why=why
      end subroutine gpusph_fail

!***********************************************************************
!     The smoothing lengths, densities and neighbour lists of the particles of the rank, on the GPU (the first part of
!     rho_and_h).  done: they are in the arrays (hp, rho, bonet_omega, bonet_0mega, bonet_psi, bonet_wn, zeta, nn, first and
!     list(1:ltot)), except for the nflag particles with gs_flag=1, which the CPU solves.  With ngpusphcheck=1 the result is
!     kept and the arrays are put back as they were (done=.false.): the CPU then does the whole step, and gpusph_compare
!     compares (checking=.true.).
!***********************************************************************
      subroutine gpusph_rhoh(divv,done,nflag,ltot,checking)
      include 'starsmasher.h'
      real*8 divv(*)
      logical done,checking
      integer nflag,ltot
      integer why,ierr,nown,listcap,i,io
      character(len=32) what
      what='h, rho and neighbour lists'
      done=.false.
      checking=.false.
      nflag=0
      ltot=0
      call gpusph_usable(why)
      if(why.ne.0) then
         call gpusph_say(what,why)
         return
      endif
      nown=n_upper-n_lower+1
      if(nown.le.0) return
      if(allocated(gs_flag)) then
         if(size(gs_flag).lt.nown) deallocate(gs_flag)
      endif
      if(.not.allocated(gs_flag)) allocate(gs_flag(nown))
!     the guesses of h: they are what the CPU starts from if the GPU does not finish, and for the check
      if(allocated(hp0)) then
         if(size(hp0).lt.nown) deallocate(hp0)
      endif
      if(.not.allocated(hp0)) allocate(hp0(nown))
!$omp parallel do simd schedule(static)
      do io=1,nown
         hp0(io)=hp(n_lower+io-1)
      enddo
      listcap=nnmax*nmax
      call gpu_sph_rhoh(n,n_lower,n_upper,nmax,x,y,z,am,u,hp,divv,rho,bonet_omega,bonet_0mega,bonet_psi,bonet_wn,zeta,nn,first,&
           gs_flag,list,listcap,ltot,dt,hfloor,nnopt,nselfgravity,ierr)
      if(ierr.ne.0) then
         if(ierr.eq.2) then
            call gpusph_say(what,8)
         else
            call gpusph_fail(7)
            call gpusph_say(what,7)
         endif
         ltot=0
!        (the results may have been put in the arrays before the lists were found not to fit)
         hp(n_lower:n_upper)=hp0(1:nown)
         return
      endif
      call gpusph_say(what,0)
      nflag=0
      do io=1,nown
         if(gs_flag(io).ne.0) nflag=nflag+1
      enddo
      done=.true.
      if(ngpusphcheck.ne.0) then
!        keep what the GPU made and put back what the CPU starts from
         if(allocated(ck_hp)) deallocate(ck_hp,ck_rho,ck_om,ck_om0,ck_psi,ck_wn,ck_nn,ck_first,ck_list)
         allocate(ck_hp(nown),ck_rho(nown),ck_om(nown),ck_om0(nown),ck_psi(nown),ck_wn(nown),ck_nn(nown),ck_first(nown))
         allocate(ck_list(max(ltot,1)))
         do io=1,nown
            i=n_lower+io-1
            ck_hp(io)=hp(i)
            ck_rho(io)=rho(i)
            ck_om(io)=bonet_omega(i)
            ck_om0(io)=bonet_0mega(i)
            ck_psi(io)=bonet_psi(i)
            ck_wn(io)=bonet_wn(i)
            ck_nn(io)=nn(i)
            ck_first(io)=first(i)
         enddo
         ck_list(1:ltot)=list(1:ltot)
         ck_ltot=ltot
         hp(n_lower:n_upper)=hp0(1:nown)
         done=.false.
         checking=.true.
         ltot=0
      endif
      end subroutine gpusph_rhoh

!     true if the two numbers have the same bits
      logical function same8(a,b)
      real*8 a,b
      same8=transfer(a,0_8).eq.transfer(b,0_8)
      end function same8

!     ngpusphcheck: the result of the CPU (the arrays) against what the GPU made (ck_*, kept by gpusph_rhoh): the run stops at
!     the first particle, in the order of the index, that has a quantity with other bits, and says which
      subroutine gpusph_compare
      include 'starsmasher.h'
      include 'mpif.h'
      integer nown,io,i,k,ierr
      character(len=24) qty
      real*8 gval,cval
      logical bad
      nown=n_upper-n_lower+1
      do io=1,nown
         if(gs_flag(io).ne.0) cycle
         i=n_lower+io-1
         bad=.true.
         if(.not.same8(ck_hp(io),hp(i))) then
            qty='h'; gval=ck_hp(io); cval=hp(i)
         else if(.not.same8(ck_rho(io),rho(i))) then
            qty='rho'; gval=ck_rho(io); cval=rho(i)
         else if(.not.same8(ck_om(io),bonet_omega(i))) then
            qty='omega'; gval=ck_om(io); cval=bonet_omega(i)
         else if(.not.same8(ck_om0(io),bonet_0mega(i))) then
            qty='0mega'; gval=ck_om0(io); cval=bonet_0mega(i)
         else if(.not.same8(ck_psi(io),bonet_psi(i))) then
            qty='psi'; gval=ck_psi(io); cval=bonet_psi(i)
         else if(u(i).ne.0.d0 .and. .not.same8(ck_wn(io),bonet_wn(i))) then
            qty='wn'; gval=ck_wn(io); cval=bonet_wn(i)
         else if(ck_nn(io).ne.nn(i)) then
            qty='nn'; gval=dble(ck_nn(io)); cval=dble(nn(i))
         else
            bad=.false.
            do k=1,nn(i)
               if(ck_list(ck_first(io)+k).ne.list(first(i)+k)) then
                  bad=.true.
                  write(qty,'(a,i0)') 'list entry ',k
                  gval=dble(ck_list(ck_first(io)+k))
                  cval=dble(list(first(i)+k))
                  exit
               endif
            enddo
         endif
         if(bad) then
            write(6,*) 'ngpusphcheck: particle',i,' (rank',myrank,'): ',trim(qty),' differs between the GPU and the CPU'
            write(6,*) '   GPU',gval,'   CPU',cval
            write(6,*) '   x=',x(i),y(i),z(i),' u=',u(i),' h on entry=',hp0(io),'  nn GPU',ck_nn(io),' CPU',nn(i)
            if(myrank.eq.0) write(69,*) 'ngpusphcheck: particle',i,': ',trim(qty),' differs between the GPU and the CPU'
            call flush(6)
            call flush(69)
            call mpi_abort(mpi_comm_world,1,ierr)
            stop 'ngpusphcheck: the GPU and the CPU differ'
         endif
      enddo
      end subroutine gpusph_compare

      end module gpusph_mod
