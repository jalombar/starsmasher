!***********************************************************************
!     The pair loop of uvdots (balAV3.f) adds each term to the particle i it
!     was computed for, which only one thread does, and to the neighbour j,
!     which any thread might.  The terms on j are therefore not added to the
!     accelerations as they come.  The own particles are taken in chunks of
!     pb_chunk consecutive indices (the same chunks for any number of
!     threads), and the terms on the j of a chunk are added up in a buffer
!     that belongs to the chunk and holds only the particles that the chunk
!     touched, in order of index.  When all the chunks are done, the buffers
!     are added to the accelerations: each particle receives the sum of its
!     chunk buffers in order of chunk, and each buffer is a sum in the order
!     of the particles and neighbour lists of its chunk.  So the bits of the
!     accelerations depend on neither the number of threads nor the schedule.
!     The particles of a chunk are close to each other, so the buffer of a
!     chunk of 128 has about a fifth as many entries as the chunk has pairs.
!***********************************************************************
      module pairbuf_mod
      use omp_utils
      implicit none
!     own particles per chunk, and targets per block of the final addition
      integer, parameter :: pb_chunk=128
      integer, parameter :: pb_tb=1024
      type pbrec
!        the particles that the chunk touched (their indices, ascending) and
!        the six sums on each (the three components of the acceleration and of its gc part)
         integer, allocatable :: tgt(:)
         real*8, allocatable :: val(:,:)
      end type pbrec
!     the working space of one thread: slot of each particle in the chunk being done (0 if it
!     is not touched yet), a bit for each particle touched (to list them in order of index
!     without sorting), and their sums in order of first touch
      type pbscr
         integer, allocatable :: smap(:)
         integer(8), allocatable :: bm(:)
         real*8, allocatable :: acc(:,:)
      end type pbscr
      type(pbrec), allocatable, save :: pb_rec(:)
      type(pbscr), allocatable, save :: pb_scr(:)
!     pb_bst(ic,tb): position in the buffer of chunk ic of the first particle of block tb of targets
      integer, allocatable, save :: pb_bst(:,:)
      integer, save :: pb_nchunk=0, pb_ntb=0
!     set by a thread that finds the pair loop in a state from which the run cannot go on
      integer, save :: pb_fail=0

      contains

!     prepare for nloop particles, targets 1..n and nthr threads (outside parallel regions)
      subroutine pb_begin(nloop,n,nthr)
      integer nloop,n,nthr,it
      pb_nchunk=(nloop+pb_chunk-1)/pb_chunk
      pb_ntb=(n+pb_tb-1)/pb_tb
      if(allocated(pb_rec)) then
         if(size(pb_rec).lt.pb_nchunk) deallocate(pb_rec)
      endif
      if(.not.allocated(pb_rec)) allocate(pb_rec(max(pb_nchunk,1)))
      if(allocated(pb_bst)) then
         if(size(pb_bst,1).lt.pb_nchunk .or. size(pb_bst,2).lt.pb_ntb+1) deallocate(pb_bst)
      endif
      if(.not.allocated(pb_bst)) allocate(pb_bst(max(pb_nchunk,1),pb_ntb+1))
      if(allocated(pb_scr)) then
         if(size(pb_scr).lt.nthr) deallocate(pb_scr)
      endif
      if(.not.allocated(pb_scr)) allocate(pb_scr(nthr))
      do it=1,nthr
         if(allocated(pb_scr(it)%smap)) then
            if(size(pb_scr(it)%smap).lt.n) deallocate(pb_scr(it)%smap,pb_scr(it)%bm,pb_scr(it)%acc)
         endif
         if(.not.allocated(pb_scr(it)%smap)) then
!           a chunk can touch every particle, so the space is for n of them; only what is used is touched
            allocate(pb_scr(it)%smap(n),pb_scr(it)%bm((n+63)/64+1),pb_scr(it)%acc(6,n))
            pb_scr(it)%smap=0
            pb_scr(it)%bm=0
         endif
      enddo
      end subroutine pb_begin

!     Chunk ic is done: cnt particles were touched (slots in smap, bits in bm, sums in acc),
!     the smallest index being jmin and the largest jmax.  Keep them in order of index and clear
!     smap and bm: the words of bm from jmin to jmax are read, and each bit set is a particle.
      subroutine pb_store(ic,cnt,jmin,jmax,smap,bm,acc)
      integer ic,cnt,jmin,jmax,smap(*)
      integer(8) bm(*)
      real*8 acc(6,*)
      integer k,kk,t,tb,slot,w,w0,w1,b
      integer(8) word
      if(cnt.eq.0) then
         pb_bst(ic,1:pb_ntb+1)=1
         return
      endif
      if(allocated(pb_rec(ic)%tgt)) then
         if(size(pb_rec(ic)%tgt).lt.cnt) deallocate(pb_rec(ic)%tgt,pb_rec(ic)%val)
      endif
      if(.not.allocated(pb_rec(ic)%tgt)) allocate(pb_rec(ic)%tgt(cnt+cnt/4+16),pb_rec(ic)%val(6,cnt+cnt/4+16))
      kk=0
      w0=ishft(jmin-1,-6)+1
      w1=ishft(jmax-1,-6)+1
!     (the arrays of the record are reached through local names: through pb_rec(ic)% every
!     access would reload the array descriptor)
      associate(rt=>pb_rec(ic)%tgt,rv=>pb_rec(ic)%val)
      do w=w0,w1
         word=bm(w)
         if(word.eq.0_8) cycle
         bm(w)=0_8
         do while(word.ne.0_8)
            b=trailz(word)
            word=iand(word,word-1_8)
            t=(w-1)*64+b+1
            slot=smap(t)
            kk=kk+1
            rt(kk)=t
            rv(1,kk)=acc(1,slot)
            rv(2,kk)=acc(2,slot)
            rv(3,kk)=acc(3,slot)
            rv(4,kk)=acc(4,slot)
            rv(5,kk)=acc(5,slot)
            rv(6,kk)=acc(6,slot)
            smap(t)=0
         enddo
      enddo
      k=1
      do tb=1,pb_ntb
         pb_bst(ic,tb)=k
         do while(k.le.kk)
            if(rt(k).gt.tb*pb_tb) exit
            k=k+1
         enddo
      enddo
      end associate
      pb_bst(ic,pb_ntb+1)=kk+1
      end subroutine pb_store

!     add the buffers of all the chunks to the sums mvx.. (the x,y,z of the acceleration) and
!     mgx.. (the gc part), each block of targets by one thread, the chunks in order
      subroutine pb_merge(mvx,mvy,mvz,mgx,mgy,mgz)
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*)
      integer tb,ic,k,t
!$omp parallel do schedule(dynamic,1) private(tb,ic,k,t)
      do tb=1,pb_ntb
         do ic=1,pb_nchunk
            do k=pb_bst(ic,tb),pb_bst(ic,tb+1)-1
               t=pb_rec(ic)%tgt(k)
               mvx(t)=mvx(t)+pb_rec(ic)%val(1,k)
               mvy(t)=mvy(t)+pb_rec(ic)%val(2,k)
               mvz(t)=mvz(t)+pb_rec(ic)%val(3,k)
               mgx(t)=mgx(t)+pb_rec(ic)%val(4,k)
               mgy(t)=mgy(t)+pb_rec(ic)%val(5,k)
               mgz(t)=mgz(t)+pb_rec(ic)%val(6,k)
            enddo
         enddo
      enddo
      end subroutine pb_merge

      end module pairbuf_mod
