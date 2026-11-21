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
!
!     Bit 2 (value 2) puts the pair sums of uvdots on the GPU: div v and curl v (getderivs), the accelerations, du/dt and
!     uijmax (uvchunk, uvpairs, Balsara's viscosity of balAV3.f).  It needs bit 1 (ngpusph is 0, 1 or 3): the lists and the
!     sums it reads are the ones of that step.  Every term of a pair is formed with the operations, tables and order of the
!     CPU, and the sums are made in the order of the CPU (the gather part in the order of the list, the scatter part chunk by
!     chunk), with no atomic operation on floating-point numbers, so that the results are the same bits in every run.  The
!     one difference is the root of the viscosity term, which the CPU takes with pow (r2**0.5d0) and the GPU with sqrt: the
!     bits of the CPU are given where no pair moves towards each other (the viscosity off, or nothing approaching), and
!     differ in the last bit of some terms otherwise.  Block steps (nblock=1) are done: the same particles as the CPU.  The
!     reasons for the CPU to do the whole step, for the forces (log line "ngpusph: the pair sums (forces)"):
!        1, 2, 3, 4, 6, 7  as above
!       10  nav=4 (the Riemann solver is on the CPU)
!       11  nintvar=3
!     Bit 4 (value 4, ngpusph=7) is the work of the third part, in steps.  The first is the Verlet lists: the candidates of the
!     first search of the smoothing length of each particle (every particle within a radius a fraction ngpusphskin larger than
!     the one asked for) are kept on the GPU from the call that made them and read from there while they are complete, that is
!     while the radius they cover, less twice the largest displacement of a particle since then, reaches the radius the search
!     needs; otherwise the search of that particle goes to the tree as before, and when more than 2% of the lists do not reach
!     the next call makes them again (a geometric criterion, never a time).  The candidates are the same, in the same order,
!     and the final filter r2<4h2 is the same: the same bits as ngpusph=3.  One GPU for the rank (with several the search is
!     as before).  The second is that the neighbour lists stay on the GPU (the forces read them there) and are not copied to list():
!     they are copied only if the CPU needs them (particles it solves itself, an error of CUDA in the forces).
!     ngpusphcheck=2 runs the CPU pair sums after the GPU ones on the same data, compares them, and writes one line to log0.sph
!     for each call with the differences (the run goes on with the CPU results); 3 is 1 and 2.
!***********************************************************************
      module gpusph_mod
      implicit none
!     0: the library has not been asked yet, 1: ready, 2: not usable (gs_why says why)
      integer, save :: gs_state=0, gs_why=0
!     the reason written last, for the smoothing lengths (1) and the forces (2)
      integer, save :: gs_said(2)=(/-1,-1/)
!     ngpusphcheck=2: calls of the forces compared
      integer, save :: gs_ncheck=0
!     the smoothing lengths, sums and lists of this step were just made on the GPU, all of them (the forces use what the
!     devices hold)
      logical, save :: gs_fresh=.false.
!     1 for the particles of the rank (i-n_lower+1) that the GPU did not do
      integer, allocatable, save :: gs_flag(:)
!     ngpusphcheck: what the GPU made, kept to compare with what the CPU makes
      real*8, allocatable, save :: ck_hp(:),ck_rho(:),ck_om(:),ck_om0(:),ck_psi(:),ck_wn(:),hp0(:)
      integer, allocatable, save :: ck_nn(:),ck_first(:),ck_list(:)
      integer, save :: ck_ltot=0

      contains

!     one line in log0.sph when the reason why the CPU does the work changes (0: the GPU does it)
      subroutine gpusph_say(what,why,part)
      include 'starsmasher.h'
      character(*) what
      integer why,part
      character(len=120) text
      if(why.eq.gs_said(part)) return
      gs_said(part)=why
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
      case(10)
         text='on the CPUs: nav=4 (Riemann solver)'
      case(11)
         text='on the CPUs: nintvar=3'
      case default
         text='on the CPUs'
      end select
      write(69,'(a,a,a,a)') 'ngpusph: ',trim(what),' ',trim(text)
      call flush(69)
      end subroutine gpusph_say

!     why=0 if the GPU can be used for this step, else the number of the reason; part 1: the smoothing lengths, 2: the forces
      subroutine gpusph_usable(why,part)
      include 'starsmasher.h'
      integer why,part,ok,ierr
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
            if(ierr.eq.0) call gpu_sph_ftables(ntab,dwtab,dgtab,ierr)
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
      else if(part.eq.2) then
         if(nav.eq.4) then
            why=10
         else if(nintvar.eq.3) then
            why=11
         endif
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
      integer why,why2,ierr,nown,listcap,i,io,skiplist
      real*8 skinarg
      logical needvdotsm
      external needvdotsm
      character(len=32) what
      what='h, rho and neighbour lists'
      done=.false.
      checking=.false.
      gs_fresh=.false.
      nflag=0
      ltot=0
      call gpusph_usable(why,1)
      if(why.ne.0) then
         call gpusph_say(what,why,1)
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
!     Verlet lists (bit 4 of ngpusph): the margin of the candidate lists kept between steps, as a fraction of the radius; 0 for none
      skinarg=0.d0
      if(iand(ngpusph,4).ne.0) skinarg=ngpusphskin
!     The lists stay on the GPU (bit 4) when the forces will be made there from them: shared steps, no check against the CPU, and
!     the forces can be done on the GPU in this run (the same reasons as in gpusph_forces).  If the CPU needs them they are copied;
!     it does when tstep calls vdotsm (cn2 or cn4 finite, see needvdotsm), which sums over list().
      skiplist=0
      if(iand(ngpusph,4).ne.0 .and. iand(ngpusph,2).ne.0 .and. ngpusphcheck.eq.0 .and. nblock.eq.0 .and. .not.needvdotsm()) then
         call gpusph_usable(why2,2)
         if(why2.eq.0) skiplist=1
      endif
      call gpu_sph_rhoh(n,n_lower,n_upper,nmax,x,y,z,am,u,hp,divv,rho,bonet_omega,bonet_0mega,bonet_psi,bonet_wn,zeta,nn,first,&
           gs_flag,list,listcap,ltot,dt,hfloor,nnopt,nselfgravity,skinarg,skiplist,ierr)
      if(ierr.ne.0) then
         if(ierr.eq.2) then
            call gpusph_say(what,8,1)
         else
            call gpusph_fail(7)
            call gpusph_say(what,7,1)
         endif
         ltot=0
!        (the results may have been put in the arrays before the lists were found not to fit)
         hp(n_lower:n_upper)=hp0(1:nown)
         return
      endif
      call gpusph_say(what,0,1)
      nflag=0
      do io=1,nown
         if(gs_flag(io).ne.0) nflag=nflag+1
      enddo
      done=.true.
      gs_fresh=(nflag.eq.0)
      if(iand(ngpusphcheck,1).ne.0) then
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
         gs_fresh=.false.
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

!***********************************************************************
!     The pair sums of uvdots on the GPU (bit 2 of ngpusph): the terms of getderivs and uvchunk for the particles of the work of
!     the rank (nloop places: the own range, or wmine with block steps; only the particles that take part in the substep,
!     actblk or refblk, are done), in mvx..mgz (the sums of the rank, over all the particles it has terms for), udot, uijmax and
!     divv, as uvpairs leaves them.  done=.false. (nothing changed): the CPU does it.  With ngpusphcheck=2 the CPU does it too,
!     after the GPU, on the same data; the two are compared (one line in log0.sph) and the arrays hold the results of the CPU.
!***********************************************************************
      subroutine gpusph_forces(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz,done)
      include 'starsmasher.h'
      include 'mpif.h'
      integer nloop
      real*8 alpha1,beta1
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*)
      logical done
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      double precision divv(nmax)
      common/commdivv/divv
      integer, allocatable, save :: gtp(:),gtiw(:)
      real*8, allocatable, save :: sv(:,:),sudot(:),suij(:),sdiv(:)
      integer why,ierr,nt,iw,i,k,ncontig,nlarge,blocks,relaxnoheat,udotother,listcap,useb
      logical allact
      character(len=32) what
      what='the pair sums (forces)'
      done=.false.
      call gpusph_usable(why,2)
      if(why.ne.0) then
         call gpusph_say(what,why,2)
!        (ngpusphcheck=2: the ranks that are not using their GPU still take part in the comparison, see gpusph_cmpidle)
         if(gpusph_cmpidle(why)) then
            call gpusph_fidle(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
            done=.true.
         endif
         return
      endif
      if(nloop.le.0) then
         done=.true.
         if(gpusph_cmpidle(0)) call gpusph_fidle(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
         return
      endif
      blocks=0
      if(nblock.eq.1) blocks=1
      if(allocated(gtp)) then
         if(size(gtp).lt.nloop) deallocate(gtp,gtiw)
      endif
      if(.not.allocated(gtp)) allocate(gtp(max(nloop,1)),gtiw(max(nloop,1)))
!     every particle of the range takes part (shared steps): the targets are consecutive and nothing is listed
      ncontig=0
      nt=0
      if(.not.blkdist) then
         allact=.true.
         do i=n_lower,n_upper
            if(.not.(actblk(i) .or. refblk(i))) then
               allact=.false.
               exit
            endif
         enddo
         if(allact) then
            ncontig=1
            nt=nloop
            gtp(1)=n_lower
            gtiw(1)=1
         endif
      endif
      if(ncontig.eq.0) then
         do iw=1,nloop
            if(blkdist) then
               i=wmine(iw)
            else
               i=n_lower+iw-1
            endif
            if(actblk(i) .or. refblk(i)) then
               nt=nt+1
               gtp(nt)=i
               gtiw(nt)=iw
            endif
         enddo
      endif
      relaxnoheat=0
      if(nrelax.eq.1 .and. t.lt.treloff) relaxnoheat=1
      udotother=0
      if(nrelax.eq.0 .or. t.ge.treloff) udotother=1
      listcap=nnmax*nmax
!     the devices hold the data of the smoothing-length kernel of this step (shared steps, every particle done on the GPU)
      useb=0
      if(gs_fresh .and. ncontig.eq.1 .and. blocks.eq.0) useb=1
      gs_fresh=.false.
      if(nt.gt.0) then
         call gpu_sph_forces(n,nmax,nt,nloop,ncontig,gtp,gtiw,x,y,z,vx,vy,vz,am,u,hp,rho,por2,bonet_omega,bonet_0mega,bonet_psi,&
              nn,first,list,listcap,actblk,blocks,alpha1,beta1,hfloor,gam,nav,nintvar,relaxnoheat,udotother,ngr,useb,&
              mvx,mvy,mvz,mgx,mgy,mgz,udot,uijmax,divv,nlarge,ierr)
         if(ierr.ne.0) then
            call gpusph_fail(7)
            call gpusph_say(what,7,2)
!           (the lists may have stayed on the GPU: the CPU forces need them)
            call gpu_sph_getlists(list,ierr)
            if(ierr.ne.0) then
               write(6,*) 'gpusph: the neighbour lists are not on the GPU either after an error of CUDA in the forces'
               call flush(6)
               call mpi_abort(mpi_comm_world,1,ierr)
            endif
            if(gpusph_cmpidle(7)) then
               call gpusph_fidle(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
               done=.true.
            endif
            return
         endif
      else
         nlarge=0
      endif
      call gpusph_say(what,0,2)
      done=.true.
!     the terms the viscosity gave a huge value (uvchunk writes one line for each and stops a run with a tiny time step)
      if(nlarge.gt.0) then
         write(69,*)'c warning large udbij: ',nlarge,' terms'
         if(dt.lt.1d-17) then
            write(69,*)'unreasonably small dt. abort run.'
            call flush(69)
            call mpi_abort(mpi_comm_world,1,ierr)
         endif
      endif
!     the factor of the variable that is integrated (uvchunk): nintvar=1 here (nintvar=3 is done on the CPU)
      if(nintvar.eq.1) then
!$omp parallel do schedule(static) private(k,i)
         do k=1,nt
            if(ncontig.eq.1) then
               i=n_lower+k-1
            else
               i=gtp(k)
            endif
            if(u(i).ne.0.d0) udot(i)=udot(i)*0.5d0*(gam-1.d0)*rho(i)**(1.d0-gam)
         enddo
      endif
      if(iand(ngpusphcheck,2).eq.0) return
!     ngpusphcheck=2: keep the GPU's results, let the CPU make its own and compare
      if(allocated(sv)) then
         if(size(sv,2).lt.n) deallocate(sv,sudot,suij,sdiv)
      endif
      if(.not.allocated(sv)) allocate(sv(6,n),sudot(n),suij(n),sdiv(n))
      do i=1,n
         sv(1,i)=mvx(i)
         sv(2,i)=mvy(i)
         sv(3,i)=mvz(i)
         sv(4,i)=mgx(i)
         sv(5,i)=mgy(i)
         sv(6,i)=mgz(i)
         sudot(i)=udot(i)
         suij(i)=uijmax(i)
         sdiv(i)=divv(i)
      enddo
      do i=1,n
         mvx(i)=0.d0
         mvy(i)=0.d0
         mvz(i)=0.d0
         mgx(i)=0.d0
         mgy(i)=0.d0
         mgz(i)=0.d0
         udot(i)=0.d0
         if(actblk(i)) uijmax(i)=0.d0
      enddo
      call uvpairs(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
      call gpusph_fcompare(mvx,mvy,mvz,mgx,mgy,mgz,sv,sudot,suij,sdiv)
      end subroutine gpusph_forces

!     ngpusphcheck=2 compares the sums of all the ranks (gpusph_fcompare: collective calls), so every rank has to go through it
!     in a step in which one does, also a rank that is not using its GPU: it has none of its own (more ranks than GPUs, reason
!     6), CUDA failed on it (7), or it has no particles to do.  The reasons that every rank has (the build, a periodic box,
!     dynhco, nav=4, nintvar=3: the GPU makes no forces anywhere) leave nobody comparing.
      logical function gpusph_cmpidle(why)
      include 'starsmasher.h'
      integer why
      gpusph_cmpidle=.false.
      if(iand(ngpusphcheck,2).eq.0 .or. nprocs.le.1) return
      if(why.ne.0 .and. why.ne.6 .and. why.ne.7) return
      if(periodic .or. dynhco.ne.0 .or. hcolim .or. nav.eq.4 .or. nintvar.eq.3) return
      gpusph_cmpidle=.true.
      end function gpusph_cmpidle

!     What a rank that is not using its GPU gives to the comparison: the CPU does its pair sums (uvpairs, as the caller would
!     do) and they are taken for the sums of the GPU too, so that this rank adds the same to both sides.
      subroutine gpusph_fidle(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
      include 'starsmasher.h'
      integer nloop
      real*8 alpha1,beta1
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*)
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      double precision divv(nmax)
      common/commdivv/divv
      real*8, allocatable :: sv(:,:),sudot(:),suij(:),sdiv(:)
      integer i
      if(nloop.gt.0) call uvpairs(nloop,alpha1,beta1,mvx,mvy,mvz,mgx,mgy,mgz)
      allocate(sv(6,n),sudot(n),suij(n),sdiv(n))
      do i=1,n
         sv(1,i)=mvx(i)
         sv(2,i)=mvy(i)
         sv(3,i)=mvz(i)
         sv(4,i)=mgx(i)
         sv(5,i)=mgy(i)
         sv(6,i)=mgz(i)
         sudot(i)=udot(i)
         suij(i)=uijmax(i)
         sdiv(i)=divv(i)
      enddo
      call gpusph_fcompare(mvx,mvy,mvz,mgx,mgy,mgz,sv,sudot,suij,sdiv)
      deallocate(sv,sudot,suij,sdiv)
      end subroutine gpusph_fidle

!     ngpusphcheck=2: the sums of the CPU (the arguments, and udot, uijmax and divv) against those of the GPU (sv, sudot, suij,
!     sdiv), over the active particles of the rank; one line in log0.sph with the largest and rms differences.  The
!     accelerations are compared as the sums of the hydrodynamic terms of all the ranks (the gather and scatter parts and the
!     gravity-correction parts added), relative to the size of a of the particle and to the largest a; du/dt, uijmax and div v
!     relative to the largest.
      subroutine gpusph_fcompare(cvx,cvy,cvz,cgx,cgy,cgz,sv,sudot,suij,sdiv)
      include 'starsmasher.h'
      include 'mpif.h'
      real*8 cvx(*),cvy(*),cvz(*),cgx(*),cgy(*),cgz(*)
      real*8 sv(6,*),sudot(*),suij(*),sdiv(*)
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      double precision divv(nmax)
      common/commdivv/divv
      real*8, allocatable :: ag(:,:),ac(:,:)
!     amx: 1 largest |da|/|a| of a particle, 2 largest |a|, 3 largest |da|, 4 largest |d udot|, 5 largest |udot|, 6 largest
!     |d uijmax|, 7 largest uijmax, 8 largest |d div v|, 9 largest |div v|; s2: sums of the squares of |da|/|a|, |da|, d udot
      real*8 amx(9),gmx(9),s2(3),gs2(3),cn(2),gcn(2)
      real*8 ax,ay,az,bx,by,bz,an,dn,r,ddiff
      integer i,k,ierr,nbit
      allocate(ag(6,n),ac(6,n))
      do i=1,n
         ag(1:6,i)=sv(1:6,i)
         ac(1,i)=cvx(i)
         ac(2,i)=cvy(i)
         ac(3,i)=cvz(i)
         ac(4,i)=cgx(i)
         ac(5,i)=cgy(i)
         ac(6,i)=cgz(i)
      enddo
      if(nprocs.gt.1) then
         call mpi_allreduce(mpi_in_place,ag,6*n,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
         call mpi_allreduce(mpi_in_place,ac,6*n,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      endif
      amx=0.d0
      s2=0.d0
      cn=0.d0
      do i=n_lower,n_upper
         if(.not.actblk(i)) cycle
         ax=ag(1,i)+ag(4,i)
         ay=ag(2,i)+ag(5,i)
         az=ag(3,i)+ag(6,i)
         bx=ac(1,i)+ac(4,i)
         by=ac(2,i)+ac(5,i)
         bz=ac(3,i)+ac(6,i)
         an=sqrt(bx*bx+by*by+bz*bz)
         dn=sqrt((ax-bx)**2+(ay-by)**2+(az-bz)**2)
         if(an.gt.0.d0) then
            r=dn/an
            amx(1)=max(amx(1),r)
            s2(1)=s2(1)+r*r
         endif
         amx(2)=max(amx(2),an)
         amx(3)=max(amx(3),dn)
         s2(2)=s2(2)+dn*dn
         cn(1)=cn(1)+1.d0
         ddiff=abs(sudot(i)-udot(i))
         amx(4)=max(amx(4),ddiff)
         amx(5)=max(amx(5),abs(udot(i)))
         s2(3)=s2(3)+ddiff*ddiff
         amx(6)=max(amx(6),abs(suij(i)-uijmax(i)))
         amx(7)=max(amx(7),abs(uijmax(i)))
         amx(8)=max(amx(8),abs(sdiv(i)-divv(i)))
         amx(9)=max(amx(9),abs(divv(i)))
         nbit=0
         do k=1,6
            if(transfer(ag(k,i),0_8).eq.transfer(ac(k,i),0_8)) nbit=nbit+1
         enddo
         if(nbit.eq.6 .and. transfer(sudot(i),0_8).eq.transfer(udot(i),0_8)) cn(2)=cn(2)+1.d0
      enddo
      call mpi_allreduce(amx,gmx,9,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(s2,gs2,3,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(cn,gcn,2,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      deallocate(ag,ac)
      gs_ncheck=gs_ncheck+1
      if(myrank.ne.0) return
      r=max(gcn(1),1.d0)
!     (key value pairs; the differences of du/dt, uijmax and div v are relative to the largest value)
      write(69,*)'ngpusphcheck forces: call',gs_ncheck,' t',t,' n',nint(gcn(1)),' bit_equal',nint(gcn(2)),&
           ' a_rel_max',gmx(1),' a_rel_rms',sqrt(gs2(1)/r),' a_max',gmx(2),' da_max',gmx(3),' da_rms',sqrt(gs2(2)/r),&
           ' du_max',gmx(4),' du_rms',sqrt(gs2(3)/r),' du_ref',gmx(5),' uij_max',gmx(6),' uij_ref',gmx(7),&
           ' div_max',gmx(8),' div_ref',gmx(9)
      call flush(69)
      end subroutine gpusph_fcompare

      end module gpusph_mod
