!***********************************************************************
!     Parallel build of a binary k-d tree over points in three dimensions.
!
!     The kd tree of the neighbour search (kdtree2.f90) and the tree of the
!     gravity (gravtree.f90) are both made here.  The tree is returned as
!     arrays of nodes (struct of arrays), and the points as a permutation
!     and a copy of their positions in tree order, so that the points of
!     a node are the consecutive elements l..u of these.
!
!     A node is split along the dimension of its largest extent, at the mean
!     of the coordinate (rule 1, as kdtree2 always did) or at the middle of
!     the box (rule 2, as the gravity tree of Phantom); the points with a
!     coordinate not above the cut go to the left child.  A node with
!     leafmax points or fewer is a leaf.
!
!     How it is done, in the manner of build_top_parallel of Phantom (Price
!     et al. 2018, PASA 35, e031, src/main/kdtree.F90): the big nodes are
!     divided level by level, all the threads working on all the nodes of
!     a level.  The points of a node are cut into chunks of tbchunk points
!     (a fixed size: not a function of the number of threads), each chunk
!     counts its points on the left, and then writes them, in their old
!     order, into the other one of two buffers, at places that follow from
!     the counts of the chunks before it.  The bounding box and the sums of
!     the coordinates of the two halves are taken by the chunk that wrote
!     them, on data still in its cache, and added in the order of the
!     chunks.  The nodes that have fewer than tbtopmin points are built each
!     by one thread, in place, as a task of a dynamic loop.  Neither the
!     chunks nor the order of the additions depend on the number of threads,
!     so the tree is the same for any number of them.
!
!     The coordinates are carried along with the permutation, so that every
!     level reads and writes consecutive memory, and the positions in tree
!     order (what kdtree2 calls rearranged_data) come out of the build
!     itself instead of from a gather at the end.
!
!     Node numbers.  The nodes made while the big ones are divided are
!     numbered 1..ntop in order of creation (the root is 1, a node comes
!     before its children).  The nodes of a small subtree are numbered from
!     a base of its own, sub_base(j)+1..sub_base(j)+sub_used(j), parent
!     before children; the root of the subtree is one of the first ntop.
!     The bases are laid out in the order of the subtrees (2 numbers per
!     point of the subtree are kept for it) so the numbering is the same
!     whatever the threads do.  A pass from the bottom of the tree to the top
!     is a loop over the subtrees (each backwards over its numbers, then its
!     root) followed by a backwards loop over 1..ntop.
!***********************************************************************
module treebuild
!$ use omp_lib
  implicit none
  private
  public :: tbwork, tb_build
  integer, parameter :: dp = kind(1.d0)
  ! points per chunk of the parallel passes over the big nodes
  integer, parameter :: tbchunk = 2048
  ! a node of this many points or more is divided by all the threads; the smaller ones are subtrees
  integer, parameter :: tbtopmin = 2048

  type tbwork
     integer :: n = 0, ncap = 0
     ! nodes used: 1..ntop and the ranges of the subtrees
     integer :: ntop = 0, nsub = 0, nidmax = 0
     ! l..u: the points of the node; lc, rc: the children (0 for a leaf); cdim: the dimension of the cut (0 for a leaf)
     integer, allocatable :: nl(:), nu(:), lc(:), rc(:), cdim(:)
     ! bounding box of the points of the node
     real(dp), allocatable :: lo(:,:), hi(:,:)
     ! the subtrees: root, base of the numbers of the other nodes, how many of them are used
     integer, allocatable :: sub_id(:), sub_base(:), sub_used(:), sub_buf(:)
     real(dp), allocatable :: sub_lo(:,:), sub_hi(:,:), sub_sm(:,:)
     ! the second buffer
     integer, allocatable :: indb(:)
     real(dp), allocatable :: posb(:,:)
     ! lists of the big nodes of a level (current and next) and their chunks
     integer, allocatable :: wk_id(:), wk_l(:), wk_u(:), wk_c(:), wk_chf(:), wk_nch(:)
     real(dp), allocatable :: wk_cv(:), wk_lo(:,:), wk_hi(:,:), wk_sm(:,:)
     integer, allocatable :: nx_id(:), nx_l(:), nx_u(:)
     real(dp), allocatable :: nx_lo(:,:), nx_hi(:,:), nx_sm(:,:)
     integer, allocatable :: ch_k(:), ch_a(:), ch_b(:), ch_cnt(:)
     real(dp), allocatable :: ch_st(:,:)
  end type tbwork

contains

!     make the arrays of w big enough for n points
  subroutine tb_reserve(w, n)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: n
    integer :: ncap, nwk, nch, nsb
    ncap = 2*n + 4*(n/tbtopmin) + 64
    nwk = n/tbtopmin + 8
    nch = n/tbchunk + 2*nwk + 8
    nsb = 4*(n/tbtopmin) + 16
    if (w%ncap >= ncap .and. w%n >= n) return
    ! the part of the nodes that is not used is never touched, so the memory costs nothing until it is
    if (allocated(w%nl)) deallocate (w%nl, w%nu, w%lc, w%rc, w%cdim, w%lo, w%hi, w%sub_id, w%sub_base, w%sub_used, w%sub_buf, &
         w%sub_lo, w%sub_hi, w%sub_sm, w%indb, w%posb, w%wk_id, w%wk_l, w%wk_u, w%wk_c, w%wk_chf, w%wk_nch, w%wk_cv, &
         w%wk_lo, w%wk_hi, w%wk_sm, w%nx_id, w%nx_l, w%nx_u, w%nx_lo, w%nx_hi, w%nx_sm, w%ch_k, w%ch_a, w%ch_b, w%ch_cnt, w%ch_st)
    allocate (w%nl(ncap), w%nu(ncap), w%lc(ncap), w%rc(ncap), w%cdim(ncap), w%lo(3,ncap), w%hi(3,ncap))
    allocate (w%sub_id(nsb), w%sub_base(nsb), w%sub_used(nsb), w%sub_buf(nsb), w%sub_lo(3,nsb), w%sub_hi(3,nsb), w%sub_sm(3,nsb))
    allocate (w%indb(n), w%posb(3,n))
    allocate (w%wk_id(nwk), w%wk_l(nwk), w%wk_u(nwk), w%wk_c(nwk), w%wk_chf(nwk), w%wk_nch(nwk), w%wk_cv(nwk), &
         w%wk_lo(3,nwk), w%wk_hi(3,nwk), w%wk_sm(3,nwk))
    allocate (w%nx_id(nwk), w%nx_l(nwk), w%nx_u(nwk), w%nx_lo(3,nwk), w%nx_hi(3,nwk), w%nx_sm(3,nwk))
    allocate (w%ch_k(nch), w%ch_a(nch), w%ch_b(nch), w%ch_cnt(nch), w%ch_st(18,nch))
    w%ncap = ncap
    w%n = n
  end subroutine tb_reserve

!     Build the tree of the n points xyz(:,1:n).  ind_out gets the permutation (the point at place i of the tree
!     is xyz(:,ind_out(i))) and pos_out their positions in that order.
  subroutine tb_build(w, n, xyz, ind_out, pos_out, leafmax, rule)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: n, leafmax, rule
    real(dp), intent(in) :: xyz(3,n)
    integer, intent(out) :: ind_out(n)
    real(dp), intent(out) :: pos_out(3,n)
    integer :: nthr, nwk, nnx, nch, nid, nsub, src, nrootch
    integer :: ic, k, j, a, b, ntot, nl, offl, kl, kr, cntl, cntr, jj
    integer :: lcid, rcid, cnt, nbase, nidl
    real(dp) :: slo(3), shi(3), ssm(3), lolo(3), lohi(3), losm(3), rolo(3), rohi(3), rosm(3)

    call tb_reserve(w, n)
    nthr = 1
!$  if (.not. omp_in_parallel() .and. n >= 4*tbtopmin) nthr = omp_get_max_threads()
    nrootch = (n + tbchunk - 1)/tbchunk
    nwk = 0
    nnx = 0
    nch = 0
    nsub = 0
    nid = 1
    src = 1

!$omp parallel default(shared) num_threads(nthr) &
!$omp private(ic,k,j,a,b,ntot,nl,offl,kl,kr,cntl,cntr,jj,lcid,rcid,cnt,slo,shi,ssm,lolo,lohi,losm,rolo,rohi,rosm,nidl)

!     the root: copy the points into the second buffer and take their box and sums, chunk by chunk
!$omp do schedule(static)
    do ic = 1, nrootch
       a = (ic - 1)*tbchunk + 1
       b = min(n, ic*tbchunk)
       call copy_stats(xyz, w%posb, w%indb, a, b, w%ch_st(1:9,ic))
    end do
!$omp end do

!$omp single
    call combine(w%ch_st(1:9,1:nrootch), nrootch, slo, shi, ssm)
    w%nl(1) = 1
    w%nu(1) = n
    w%lo(:,1) = slo
    w%hi(:,1) = shi
    w%lc(1) = 0
    w%rc(1) = 0
    w%cdim(1) = 0
    if (n >= tbtopmin) then
       nwk = 1
       w%wk_id(1) = 1
       w%wk_l(1) = 1
       w%wk_u(1) = n
       w%wk_lo(:,1) = slo
       w%wk_hi(:,1) = shi
       w%wk_sm(:,1) = ssm
    else
       nsub = 1
       w%sub_id(1) = 1
       w%sub_buf(1) = 1
       w%sub_lo(:,1) = slo
       w%sub_hi(:,1) = shi
       w%sub_sm(:,1) = ssm
    end if
    call level_setup(w, nwk, nch, leafmax, rule)
!$omp end single

!     the big nodes, a level at a time
    do while (nwk > 0)
!$omp do schedule(static)
       do ic = 1, nch
          k = w%ch_k(ic)
          if (w%wk_c(k) == 0) then
             w%ch_cnt(ic) = 0
          else if (src == 1) then
             call chunk_count(w%posb, w%ch_a(ic), w%ch_b(ic), w%wk_c(k), w%wk_cv(k), w%ch_cnt(ic))
          else
             call chunk_count(pos_out, w%ch_a(ic), w%ch_b(ic), w%wk_c(k), w%wk_cv(k), w%ch_cnt(ic))
          end if
       end do
!$omp end do

!$omp do schedule(static)
       do ic = 1, nch
          k = w%ch_k(ic)
          if (w%wk_c(k) == 0) cycle
          a = w%ch_a(ic)
          b = w%ch_b(ic)
          ! the points on the left in the chunks before this one, and in the whole node
          offl = 0
          nl = 0
          do jj = w%wk_chf(k), w%wk_chf(k) + w%wk_nch(k) - 1
             if (jj < ic) offl = offl + w%ch_cnt(jj)
             nl = nl + w%ch_cnt(jj)
          end do
          ntot = w%wk_u(k) - w%wk_l(k) + 1
          if (nl == 0 .or. nl == ntot) cycle
          kl = w%wk_l(k) + offl
          kr = w%wk_l(k) + nl + (a - w%wk_l(k) - offl)
          cntl = w%ch_cnt(ic)
          cntr = b - a + 1 - cntl
          if (src == 1) then
             call chunk_move(w%posb, w%indb, pos_out, ind_out, a, b, w%wk_c(k), w%wk_cv(k), kl, kr)
             call stats(pos_out, kl, kl + cntl - 1, w%ch_st(1:9,ic), .true.)
             call stats(pos_out, kr, kr + cntr - 1, w%ch_st(10:18,ic), .true.)
          else
             call chunk_move(pos_out, ind_out, w%posb, w%indb, a, b, w%wk_c(k), w%wk_cv(k), kl, kr)
             call stats(w%posb, kl, kl + cntl - 1, w%ch_st(1:9,ic), .true.)
             call stats(w%posb, kr, kr + cntr - 1, w%ch_st(10:18,ic), .true.)
          end if
       end do
!$omp end do

!$omp single
!      the nodes of the level are done: make their children, and the lists of the next level
       nnx = 0
       do k = 1, nwk
          ntot = w%wk_u(k) - w%wk_l(k) + 1
          nl = 0
          do jj = w%wk_chf(k), w%wk_chf(k) + w%wk_nch(k) - 1
             nl = nl + w%ch_cnt(jj)
          end do
          if (w%wk_c(k) == 0 .or. nl == 0 .or. nl == ntot) then
             ! it cannot be divided (the points are all the same, or the cut is at an end): a leaf holding them all, where they are
             call grow_subs(w, nsub + 1, nsub)
             nsub = nsub + 1
             w%sub_id(nsub) = w%wk_id(k)
             w%sub_buf(nsub) = src
             w%sub_lo(:,nsub) = w%wk_lo(:,k)
             w%sub_hi(:,nsub) = w%wk_hi(:,k)
             w%sub_sm(:,nsub) = w%wk_sm(:,k)
             w%cdim(w%wk_id(k)) = -1
             cycle
          end if
          ! boxes and sums of the two halves, added in the order of the chunks
          call combine_lr(w%ch_st(:,w%wk_chf(k):w%wk_chf(k)+w%wk_nch(k)-1), w%wk_nch(k), lolo, lohi, losm, rolo, rohi, rosm)
          call grow_nodes(w, nid + 2, nid)
          lcid = nid + 1
          rcid = nid + 2
          nid = nid + 2
          w%lc(w%wk_id(k)) = lcid
          w%rc(w%wk_id(k)) = rcid
          w%cdim(w%wk_id(k)) = w%wk_c(k)
          w%nl(lcid) = w%wk_l(k)
          w%nu(lcid) = w%wk_l(k) + nl - 1
          w%nl(rcid) = w%wk_l(k) + nl
          w%nu(rcid) = w%wk_u(k)
          w%lo(:,lcid) = lolo
          w%hi(:,lcid) = lohi
          w%lo(:,rcid) = rolo
          w%hi(:,rcid) = rohi
          w%lc(lcid) = 0
          w%rc(lcid) = 0
          w%cdim(lcid) = 0
          w%lc(rcid) = 0
          w%rc(rcid) = 0
          w%cdim(rcid) = 0
          ! the children: more to divide, or a subtree (their points are now in the other buffer)
          if (nl >= tbtopmin) then
             nnx = nnx + 1
             w%nx_id(nnx) = lcid
             w%nx_l(nnx) = w%nl(lcid)
             w%nx_u(nnx) = w%nu(lcid)
             w%nx_lo(:,nnx) = lolo
             w%nx_hi(:,nnx) = lohi
             w%nx_sm(:,nnx) = losm
          else
             call grow_subs(w, nsub + 1, nsub)
             nsub = nsub + 1
             w%sub_id(nsub) = lcid
             w%sub_buf(nsub) = 1 - src
             w%sub_lo(:,nsub) = lolo
             w%sub_hi(:,nsub) = lohi
             w%sub_sm(:,nsub) = losm
          end if
          if (ntot - nl >= tbtopmin) then
             nnx = nnx + 1
             w%nx_id(nnx) = rcid
             w%nx_l(nnx) = w%nl(rcid)
             w%nx_u(nnx) = w%nu(rcid)
             w%nx_lo(:,nnx) = rolo
             w%nx_hi(:,nnx) = rohi
             w%nx_sm(:,nnx) = rosm
          else
             call grow_subs(w, nsub + 1, nsub)
             nsub = nsub + 1
             w%sub_id(nsub) = rcid
             w%sub_buf(nsub) = 1 - src
             w%sub_lo(:,nsub) = rolo
             w%sub_hi(:,nsub) = rohi
             w%sub_sm(:,nsub) = rosm
          end if
       end do
       nwk = nnx
       do k = 1, nwk
          w%wk_id(k) = w%nx_id(k)
          w%wk_l(k) = w%nx_l(k)
          w%wk_u(k) = w%nx_u(k)
          w%wk_lo(:,k) = w%nx_lo(:,k)
          w%wk_hi(:,k) = w%nx_hi(:,k)
          w%wk_sm(:,k) = w%nx_sm(:,k)
       end do
       src = 1 - src
       call level_setup(w, nwk, nch, leafmax, rule)
!$omp end single
    end do

!$omp single
!   the numbers of the nodes of the subtrees, laid out in the order of the subtrees
    w%ntop = nid
    nbase = nid
    do j = 1, nsub
       w%sub_base(j) = nbase
       cnt = w%nu(w%sub_id(j)) - w%nl(w%sub_id(j)) + 1
       if (cnt > leafmax .and. w%cdim(w%sub_id(j)) >= 0) nbase = nbase + 2*cnt
       w%sub_used(j) = 0
    end do
    w%nsub = nsub
    w%nidmax = nbase
!   the numbers up to nbase may be used: the arrays of the nodes must reach it
    call grow_nodes(w, nbase, w%ntop)
!$omp end single

!$omp do schedule(dynamic,1)
    do j = 1, nsub
       k = w%sub_id(j)
       nidl = w%sub_base(j)
       a = w%nl(k)
       b = w%nu(k)
       if (w%cdim(k) >= 0) then
          if (w%sub_buf(j) == 1) then
             call sub_build(w, k, nidl, leafmax, rule, w%indb, w%posb, w%sub_sm(:,j))
          else
             call sub_build(w, k, nidl, leafmax, rule, ind_out, pos_out, w%sub_sm(:,j))
          end if
       else
          w%cdim(k) = 0
       end if
       w%sub_used(j) = nidl - w%sub_base(j)
       ! the points of the subtree go to the final buffer if they are not in it
       if (w%sub_buf(j) == 1) then
          ind_out(a:b) = w%indb(a:b)
          pos_out(:,a:b) = w%posb(:,a:b)
       end if
    end do
!$omp end do
!$omp end parallel
  end subroutine tb_build

!     The sizes given by tb_reserve are estimates: how many nodes and how many subtrees there are depends on how the points are
!     distributed (a diffuse halo gives many thin big nodes).  These make the arrays larger when the build needs more.
!     Called by one thread, where the others wait.
  subroutine grow_nodes(w, need, nused)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: need, nused
    integer :: ncap
    if (need <= w%ncap) return
    ncap = max(need, w%ncap + w%ncap/2)
    call grow_i(w%nl, ncap, nused)
    call grow_i(w%nu, ncap, nused)
    call grow_i(w%lc, ncap, nused)
    call grow_i(w%rc, ncap, nused)
    call grow_i(w%cdim, ncap, nused)
    call grow_r(w%lo, ncap, nused)
    call grow_r(w%hi, ncap, nused)
    w%ncap = ncap
  end subroutine grow_nodes

!     room for need subtrees, of which nused are filled
  subroutine grow_subs(w, need, nused)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: need, nused
    integer :: ncap
    if (need <= size(w%sub_id)) return
    ncap = max(need, size(w%sub_id) + size(w%sub_id)/2)
    call grow_i(w%sub_id, ncap, nused)
    call grow_i(w%sub_base, ncap, nused)
    call grow_i(w%sub_used, ncap, nused)
    call grow_i(w%sub_buf, ncap, nused)
    call grow_r(w%sub_lo, ncap, nused)
    call grow_r(w%sub_hi, ncap, nused)
    call grow_r(w%sub_sm, ncap, nused)
  end subroutine grow_subs

  subroutine grow_i(a, ncap, nused)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: ncap, nused
    integer, allocatable :: t(:)
    allocate (t(ncap))
    t(1:nused) = a(1:nused)
    call move_alloc(t, a)
  end subroutine grow_i

  subroutine grow_r(a, ncap, nused)
    real(dp), allocatable, intent(inout) :: a(:,:)
    integer, intent(in) :: ncap, nused
    real(dp), allocatable :: t(:,:)
    allocate (t(3,ncap))
    t(:,1:nused) = a(:,1:nused)
    call move_alloc(t, a)
  end subroutine grow_r

!     The nodes of the next level: the cut of each (dimension of the largest extent, and the value) and the chunks
  subroutine level_setup(w, nwk, nch, leafmax, rule)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: nwk, leafmax, rule
    integer, intent(out) :: nch
    integer :: k, c, cnt, nc, j, a
    real(dp) :: ext(3)
    nch = 0
    do k = 1, nwk
       cnt = w%wk_u(k) - w%wk_l(k) + 1
       ext = w%wk_hi(:,k) - w%wk_lo(:,k)
       c = 1
       if (ext(2) > ext(c)) c = 2
       if (ext(3) > ext(c)) c = 3
       if (ext(c) <= 0.d0 .or. cnt <= leafmax) then
          w%wk_c(k) = 0
          w%wk_cv(k) = 0.d0
       else
          w%wk_c(k) = c
          if (rule == 1) then
             w%wk_cv(k) = w%wk_sm(c,k)/real(cnt, dp)
          else
             w%wk_cv(k) = 0.5d0*(w%wk_lo(c,k) + w%wk_hi(c,k))
          end if
       end if
       nc = (cnt + tbchunk - 1)/tbchunk
       w%wk_chf(k) = nch + 1
       w%wk_nch(k) = nc
       do j = 1, nc
          nch = nch + 1
          a = w%wk_l(k) + (j - 1)*tbchunk
          w%ch_k(nch) = k
          w%ch_a(nch) = a
          w%ch_b(nch) = min(w%wk_u(k), a + tbchunk - 1)
       end do
    end do
  end subroutine level_setup

!     copy the points a..b of x into the second buffer (with their numbers) and take their box and sums
  subroutine copy_stats(x, pos, ind, a, b, st)
    real(dp), intent(in) :: x(3,*)
    real(dp), intent(out) :: pos(3,*)
    integer, intent(out) :: ind(*)
    integer, intent(in) :: a, b
    real(dp), intent(out) :: st(9)
    integer :: i
    real(dp) :: x1, x2, x3, l1, l2, l3, h1, h2, h3, s1, s2, s3
    l1 = huge(1.d0); l2 = l1; l3 = l1
    h1 = -huge(1.d0); h2 = h1; h3 = h1
    s1 = 0.d0; s2 = 0.d0; s3 = 0.d0
    do i = a, b
       x1 = x(1,i); x2 = x(2,i); x3 = x(3,i)
       pos(1,i) = x1; pos(2,i) = x2; pos(3,i) = x3
       ind(i) = i
       l1 = min(l1, x1); l2 = min(l2, x2); l3 = min(l3, x3)
       h1 = max(h1, x1); h2 = max(h2, x2); h3 = max(h3, x3)
       s1 = s1 + x1; s2 = s2 + x2; s3 = s3 + x3
    end do
    st(1) = l1; st(2) = l2; st(3) = l3; st(4) = h1; st(5) = h2; st(6) = h3; st(7) = s1; st(8) = s2; st(9) = s3
  end subroutine copy_stats

!     box and sums of pos(:,a:b): st = lo(3), hi(3), sum(3)  (an empty range gives an empty box)
  subroutine stats(pos, a, b, st, wantsum)
    real(dp), intent(in) :: pos(3,*)
    integer, intent(in) :: a, b
    real(dp), intent(out) :: st(9)
    logical, intent(in) :: wantsum
    integer :: i
    real(dp) :: x1, x2, x3, l1, l2, l3, h1, h2, h3, s1, s2, s3
    l1 = huge(1.d0); l2 = l1; l3 = l1
    h1 = -huge(1.d0); h2 = h1; h3 = h1
    s1 = 0.d0; s2 = 0.d0; s3 = 0.d0
    if (wantsum) then
       do i = a, b
          x1 = pos(1,i); x2 = pos(2,i); x3 = pos(3,i)
          l1 = min(l1, x1); l2 = min(l2, x2); l3 = min(l3, x3)
          h1 = max(h1, x1); h2 = max(h2, x2); h3 = max(h3, x3)
          s1 = s1 + x1; s2 = s2 + x2; s3 = s3 + x3
       end do
    else
       do i = a, b
          x1 = pos(1,i); x2 = pos(2,i); x3 = pos(3,i)
          l1 = min(l1, x1); l2 = min(l2, x2); l3 = min(l3, x3)
          h1 = max(h1, x1); h2 = max(h2, x2); h3 = max(h3, x3)
       end do
    end if
    st(1) = l1; st(2) = l2; st(3) = l3; st(4) = h1; st(5) = h2; st(6) = h3; st(7) = s1; st(8) = s2; st(9) = s3
  end subroutine stats

!     the box and sums of the union of the chunks, added in the order of the chunks
  subroutine combine(st, nc, lo, hi, sm)
    integer, intent(in) :: nc
    real(dp), intent(in) :: st(9,nc)
    real(dp), intent(out) :: lo(3), hi(3), sm(3)
    integer :: j
    lo = st(1:3,1)
    hi = st(4:6,1)
    sm = st(7:9,1)
    do j = 2, nc
       lo = min(lo, st(1:3,j))
       hi = max(hi, st(4:6,j))
       sm = sm + st(7:9,j)
    end do
  end subroutine combine

  subroutine combine_lr(st, nc, lolo, lohi, losm, rolo, rohi, rosm)
    integer, intent(in) :: nc
    real(dp), intent(in) :: st(18,nc)
    real(dp), intent(out) :: lolo(3), lohi(3), losm(3), rolo(3), rohi(3), rosm(3)
    integer :: j
    lolo = st(1:3,1)
    lohi = st(4:6,1)
    losm = st(7:9,1)
    rolo = st(10:12,1)
    rohi = st(13:15,1)
    rosm = st(16:18,1)
    do j = 2, nc
       lolo = min(lolo, st(1:3,j))
       lohi = max(lohi, st(4:6,j))
       losm = losm + st(7:9,j)
       rolo = min(rolo, st(10:12,j))
       rohi = max(rohi, st(13:15,j))
       rosm = rosm + st(16:18,j)
    end do
  end subroutine combine_lr

!     how many of the points a..b have coordinate c not above cv
  subroutine chunk_count(pos, a, b, c, cv, cnt)
    real(dp), intent(in) :: pos(3,*)
    integer, intent(in) :: a, b, c
    real(dp), intent(in) :: cv
    integer, intent(out) :: cnt
    integer :: i, k
    k = 0
    do i = a, b
       if (pos(c,i) <= cv) k = k + 1
    end do
    cnt = k
  end subroutine chunk_count

!     write the points a..b of the source buffer into the destination, those with coordinate c not above cv from kl on and
!     the others from kr on, each group in its old order
  subroutine chunk_move(ps, is, pd, id, a, b, c, cv, kl, kr)
    real(dp), intent(in) :: ps(3,*)
    integer, intent(in) :: is(*)
    real(dp), intent(inout) :: pd(3,*)
    integer, intent(inout) :: id(*)
    integer, intent(in) :: a, b, c, kl, kr
    real(dp), intent(in) :: cv
    integer :: i, l, r
    l = kl
    r = kr
    do i = a, b
       if (ps(c,i) <= cv) then
          pd(1,l) = ps(1,i); pd(2,l) = ps(2,i); pd(3,l) = ps(3,i)
          id(l) = is(i)
          l = l + 1
       else
          pd(1,r) = ps(1,i); pd(2,r) = ps(2,i); pd(3,r) = ps(3,i)
          id(r) = is(i)
          r = r + 1
       end if
    end do
  end subroutine chunk_move

!     The subtree of node k, whose points are the elements l..u of ind and pos, built in place by one thread.
!     nid is the last number of a node used (the numbers of the nodes of the subtree follow it).
  recursive subroutine sub_build(w, k, nid, leafmax, rule, ind, pos, sm)
    type(tbwork), intent(inout) :: w
    integer, intent(in) :: k, leafmax, rule
    integer, intent(inout) :: nid
    integer, intent(inout) :: ind(*)
    real(dp), intent(inout) :: pos(3,*)
    real(dp), intent(in) :: sm(3)
    integer :: l, u, cnt, c, m, lcid, rcid, i, j, it
    real(dp) :: cv, ext(3), x1, x2, x3, st(9), lsm(3), rsm(3)
    l = w%nl(k)
    u = w%nu(k)
    cnt = u - l + 1
    w%lc(k) = 0
    w%rc(k) = 0
    w%cdim(k) = 0
    if (cnt <= leafmax) return
    ext = w%hi(:,k) - w%lo(:,k)
    c = 1
    if (ext(2) > ext(c)) c = 2
    if (ext(3) > ext(c)) c = 3
    if (ext(c) <= 0.d0) return
    if (rule == 1) then
       cv = sm(c)/real(cnt, dp)
    else
       cv = 0.5d0*(w%lo(c,k) + w%hi(c,k))
    end if
    ! the points not above the cut to the front, in place
    i = l
    j = u
    do while (i <= j)
       if (pos(c,i) <= cv) then
          i = i + 1
       else
          x1 = pos(1,i); x2 = pos(2,i); x3 = pos(3,i); it = ind(i)
          pos(1,i) = pos(1,j); pos(2,i) = pos(2,j); pos(3,i) = pos(3,j); ind(i) = ind(j)
          pos(1,j) = x1; pos(2,j) = x2; pos(3,j) = x3; ind(j) = it
          j = j - 1
       end if
    end do
    m = i - 1
    if (m < l .or. m >= u) return
    lcid = nid + 1
    rcid = nid + 2
    nid = nid + 2
    w%lc(k) = lcid
    w%rc(k) = rcid
    w%cdim(k) = c
    w%nl(lcid) = l
    w%nu(lcid) = m
    w%nl(rcid) = m + 1
    w%nu(rcid) = u
    call stats(pos, l, m, st, rule == 1 .and. m - l + 1 > leafmax)
    w%lo(:,lcid) = st(1:3)
    w%hi(:,lcid) = st(4:6)
    lsm = st(7:9)
    call stats(pos, m + 1, u, st, rule == 1 .and. u - m > leafmax)
    w%lo(:,rcid) = st(1:3)
    w%hi(:,rcid) = st(4:6)
    rsm = st(7:9)
    call sub_build(w, lcid, nid, leafmax, rule, ind, pos, lsm)
    call sub_build(w, rcid, nid, leafmax, rule, ind, pos, rsm)
  end subroutine sub_build

end module treebuild
