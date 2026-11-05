!***********************************************************************
!     Gravity with a tree on the cpus (ngravtree=1), after the tree code of
!     Phantom (Price et al. 2018, PASA 35, e031, src/main/kdtree.F90, GPLv3
!     like StarSmasher).
!
!     The particles are put in a binary k-d tree, divided at the middle of the
!     longest side of the box (treebuild.f90, the same builder as the one of
!     the neighbour search), with leaves of at most nleafgrav particles.
!     Each node holds its mass, centre of mass, the second moments of its mass
!     about the centre of mass (the quadrupole), the radius of the smallest
!     sphere about the centre of mass that holds all its particles (size) and
!     the largest softening radius 2h of its particles (hmax).
!
!     The force is found one leaf at a time, as in the "classic" walk of
!     Phantom.  The tree is walked from the root for the leaf L.  A node S is
!     well separated from L if
!        theta^2 r^2 > (size_L + size_S)^2   and   r > size_L + size_S + max(hmax_L, hmax_S)
!     with r the distance between their centres of mass: then no particle of
!     L is inside the softening of a particle of S (or the other way round), the
!     force between them is Newtonian, and the quadrupole expansion of S is
!     added to the Taylor series of the force about the centre of mass of L
!     (m2l_loop).  Otherwise, a leaf S is added
!     to the list of the particles that interact directly with those of L, and
!     for any other node the two children are visited, the left one first.  The
!     force on each particle of L is the Taylor series at its position plus
!     the direct sum over the list with the pair formula of get_gravity_using_cpus
!     (gravpair.f90), in double precision.
!
!     Nothing depends on the number of threads: the tree is the same (see
!     treebuild.f90), the moments of a node come from its children in a fixed
!     order, and a leaf is done by one thread, in an order of the nodes and of the
!     sums that is fixed by the tree.  With theta=0 no node is well separated
!     and the result is the direct sum.
!
!     ngravtree=2 (gt_gpu = .true.): the tree, its moments and the list of the
!     leaves are made here, on the cpus, as above, but the visit is done by the
!     gpus (SPHgrav_lib2/grav_tree_gpu.cu, which has the same physics: the same
!     criterion, the same Taylor series, the same pair formula).  gpu_visit puts
!     the particles and the nodes in the host buffers of the library (the nodes
!     are numbered again without the holes between the subtrees), lists the
!     leaves that have a particle to do, and calls the library.
!***********************************************************************
module gravtree
!$ use omp_lib
  use iso_c_binding
  use treebuild
  use gravpair
  implicit none
  private
  public :: gt_forces, gt_gpu, gt_start, gt_finish, gt_pend
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: lenf = 20
  ! node data: 1:3 centre of mass, 4 size, 5 hmax, 6 mass, 7:12 second moments (xx,xy,xz,yy,yz,zz)
  integer, parameter :: nfd = 12
  integer, parameter :: stackmax = 1024
  ! leaves handed to a thread at a time
  integer, parameter :: leafchunk = 4
  ! ngravtree=2: the visit of the tree on the gpus (set by the caller before gt_forces)
  logical, save :: gt_gpu = .false.
  ! ngravoverlap=1 with ngravtree=2: gt_start has made the tree and queued the visit on the gpus, and gt_finish has not
  ! collected it yet
  logical, save :: gt_pend = .false.
  integer, save :: ps_n = 0, ps_ilow = 0, ps_ihigh = 0, ps_nwork = 0, ps_nnode = 0
  integer(c_signed_char), pointer, contiguous, save :: ps_needf(:) => null()
  real(c_double), pointer, contiguous, save :: ps_out(:,:) => null()

  type(tbwork), save :: gtw
  integer :: nsave = 0
  real(dp), allocatable, save :: gxyz(:,:), gpos(:,:)
  integer, allocatable, save :: gind(:)
  ! the particles in tree order
  real(dp), allocatable, save :: xt(:), yt(:), zt(:), mt(:), ht(:), h2t(:), ft(:)
  real(dp), allocatable, save :: nd(:,:)
  integer, allocatable, save :: leaves(:), leafcnt(:), leafoff(:)

  ! the lists of a thread: mb nodes (dx dy dz q0 and the second moments), nb sources that are Newtonian for the whole leaf
  ! (x y z m), cb the other sources (x y z m 2h (2h)^2 flag r^2), sb the pairs inside a softening of one particle
  type walkbuf
     integer :: capm = 0, capn = 0, capc = 0
     real(dp), allocatable :: mb(:,:), nb(:,:), cb(:,:), sb(:,:)
     integer, allocatable :: ci(:)
     integer, allocatable :: stack(:)
  end type walkbuf
  type(walkbuf), allocatable, save :: wb(:)

  ! the library of the gpus (SPHgrav_lib2/grav_tree_gpu.cu; stubs in cpu_grav.f)
  interface
     subroutine gpu_tree_reserve(n, nnode, pa, pb, pneed, pnf, pni, pwork, pout) bind(C, name='gpu_tree_reserve_')
       import :: c_int, c_ptr
       integer(c_int), intent(in) :: n, nnode
       type(c_ptr), intent(out) :: pa, pb, pneed, pnf, pni, pwork, pout
     end subroutine gpu_tree_reserve
     subroutine gpu_tree_launch(n, nnode, nwork, nkernel, theta2) bind(C, name='gpu_tree_launch_')
       import :: c_int, c_double
       integer(c_int), intent(in) :: n, nnode, nwork, nkernel
       real(c_double), intent(in) :: theta2
     end subroutine gpu_tree_launch
     subroutine gpu_tree_wait(ierr) bind(C, name='gpu_tree_wait_')
       import :: c_int
       integer(c_int), intent(out) :: ierr
     end subroutine gpu_tree_wait
  end interface
  ! work arrays of gpu_visit: the first particle of the leaf that starts at a place (the number of the node plus 1), the
  ! subtrees' offsets
  ! in the numbering without holes, the leaves of each block of places
  integer, allocatable, save :: lmark(:), soff(:), bcnt(:), boff(:)

contains

!     The tree of the n particles, its moments, and the particles in tree order (xt.. ft).
  subroutine make_tree(n, x, y, z, am, hp, u, nleaf)
    integer, intent(in) :: n, nleaf
    real(dp), intent(in) :: x(n), y(n), z(n), am(n), hp(n), u(n)
    integer :: i, j
    call reserve(n, nleaf)
!$omp parallel do simd schedule(static) private(i)
    do i = 1, n
       gxyz(1,i) = x(i)
       gxyz(2,i) = y(i)
       gxyz(3,i) = z(i)
    end do
    call tb_build(gtw, n, gxyz, gind, gpos, nleaf, 2)
    ! the builder may have needed more node numbers than the estimate of reserve
    if (size(nd,2) < gtw%nidmax) then
       deallocate (nd)
       allocate (nd(nfd,gtw%nidmax))
    end if
!$omp parallel do simd schedule(static) private(i, j)
    do i = 1, n
       j = gind(i)
       xt(i) = gpos(1,i)
       yt(i) = gpos(2,i)
       zt(i) = gpos(3,i)
       mt(i) = am(j)
       ht(i) = 2*hp(j)
       h2t(i) = ht(i)*ht(i)
       if (u(j) /= 0.d0) then
          ft(i) = 1.d0
       else
          ft(i) = 0.d0
       end if
    end do
    call moments()
  end subroutine make_tree

!     Gravity on the particles ilow..ihigh (and, with block steps, only the active ones: the others get zero).
  subroutine gt_forces(n, x, y, z, am, hp, u, act, blk, ilow, ihigh, nk, theta, nleaf, gx, gy, gz, gpot)
    integer, intent(in) :: n, ilow, ihigh, nk, nleaf
    real(dp), intent(in) :: x(n), y(n), z(n), am(n), hp(n), u(n), theta
    logical, intent(in) :: act(n), blk
    real(dp), intent(inout) :: gx(n), gy(n), gz(n), gpot(n)
    integer :: i, id, nlf
    real(dp) :: theta2
    theta2 = theta*theta
!$omp parallel do simd schedule(static) private(i)
    do i = ilow, ihigh
       gx(i) = 0.d0
       gy(i) = 0.d0
       gz(i) = 0.d0
       gpot(i) = 0.d0
    end do
    call make_tree(n, x, y, z, am, hp, u, nleaf)
    if (gt_gpu) then
       ! the visit on the gpus
       call gpu_visit(n, ilow, ihigh, act, blk, nk, theta2, gx, gy, gz, gpot)
       return
    end if
    call leaflist(nlf)
!$omp parallel default(shared) private(i,id)
    call walk_setup()
!$omp do schedule(dynamic,leafchunk)
    do i = 1, nlf
       id = leaves(i)
       call walk_leaf(id, nk, theta2, ilow, ihigh, act, blk, gx, gy, gz, gpot, n)
    end do
!$omp end do
!$omp end parallel
  end subroutine gt_forces

!     ngravtree=2 with ngravoverlap=1: the tree is made and the visit is queued on the gpus (gt_start), the cpus do something
!     else, and
!     the results are collected later (gt_finish).  The same calculation as gt_forces with gt_gpu = .true., in two halves: the
!     bits are
!     the same.  Between the two the positions, the masses, the smoothing lengths and u must not change.
  subroutine gt_start(n, x, y, z, am, hp, u, act, blk, ilow, ihigh, nk, theta, nleaf)
    integer, intent(in) :: n, ilow, ihigh, nk, nleaf
    real(dp), intent(in) :: x(n), y(n), z(n), am(n), hp(n), u(n), theta
    logical, intent(in) :: act(n), blk
    call make_tree(n, x, y, z, am, hp, u, nleaf)
    call gpu_visit_start(n, ilow, ihigh, act, blk, nk, theta*theta)
  end subroutine gt_start

!     the second half: wait for the visit queued by gt_start and put the forces of the particles ilow..ihigh of gt_start in gx,
!     gy, gz, gpot
  subroutine gt_finish(n, gx, gy, gz, gpot)
    integer, intent(in) :: n
    real(dp), intent(inout) :: gx(n), gy(n), gz(n), gpot(n)
    call gpu_visit_finish(gx, gy, gz, gpot)
  end subroutine gt_finish

!     arrays for n particles, and the buffers of the threads
  subroutine reserve(n, nleaf)
    integer, intent(in) :: n, nleaf
    integer :: nthr, ncap
    ncap = 2*n + 4*(n/2048) + 64
    if (nsave < n .or. .not. allocated(gxyz)) then
       if (allocated(gxyz)) deallocate (gxyz, gpos, gind, xt, yt, zt, mt, ht, h2t, ft, nd, leaves, leafcnt, leafoff)
       allocate (gxyz(3,n), gpos(3,n), gind(n), xt(n), yt(n), zt(n), mt(n), ht(n), h2t(n), ft(n))
       allocate (nd(nfd,ncap))
       allocate (leaves(n), leafcnt(4*(n/2048)+64), leafoff(4*(n/2048)+64))
       nsave = n
    end if
    nthr = 1
!$  nthr = omp_get_max_threads()
    if (.not. allocated(wb)) allocate (wb(0:nthr-1))
    if (size(wb) < nthr) then
       deallocate (wb)
       allocate (wb(0:nthr-1))
    end if
  end subroutine reserve

!     the buffer of the calling thread
  subroutine walk_setup()
    integer :: it
    it = 0
!$  it = omp_get_thread_num()
    if (wb(it)%capm == 0) then
       wb(it)%capm = 2048
       wb(it)%capn = 4096
       wb(it)%capc = 4096
       allocate (wb(it)%mb(2056,10), wb(it)%nb(4104,4), wb(it)%cb(4104,8), wb(it)%ci(4104), wb(it)%sb(4104,6))
       allocate (wb(it)%stack(stackmax))
    end if
  end subroutine walk_setup

!     The moments of all the nodes, from the bottom of the tree up: the subtrees, then the nodes above them.
  subroutine moments()
    integer :: j, id, k, ntop
    ntop = gtw%ntop
!$omp parallel default(shared) private(j,id,k)
!$omp do schedule(dynamic,1)
    do j = 1, gtw%nsub
       do k = gtw%sub_base(j) + gtw%sub_used(j), gtw%sub_base(j) + 1, -1
          call node_moments(k)
       end do
       call node_moments(gtw%sub_id(j))
    end do
!$omp end do
!$omp end parallel
    ! the nodes made while the big ones were divided: the children have the larger numbers
    do id = ntop, 1, -1
       if (gtw%lc(id) /= 0) call node_moments(id)
    end do
  end subroutine moments

!     the moments of node id: from its particles if it is a leaf, from its two children otherwise
  subroutine node_moments(id)
    integer, intent(in) :: id
    integer :: i, l, u, a, b
    real(dp) :: m, sx, sy, sz, cx, cy, cz, dx, dy, dz, r2max, hm
    real(dp) :: qxx, qxy, qxz, qyy, qyz, qzz, ma, mb, dax, day, daz, dbx, dby, dbz, sa, sb
    if (gtw%lc(id) == 0) then
       l = gtw%nl(id)
       u = gtw%nu(id)
       m = 0.d0
       sx = 0.d0
       sy = 0.d0
       sz = 0.d0
       hm = 0.d0
       do i = l, u
          m = m + mt(i)
          sx = sx + mt(i)*xt(i)
          sy = sy + mt(i)*yt(i)
          sz = sz + mt(i)*zt(i)
          hm = max(hm, ht(i))
       end do
       cx = sx/m
       cy = sy/m
       cz = sz/m
       qxx = 0.d0
       qxy = 0.d0
       qxz = 0.d0
       qyy = 0.d0
       qyz = 0.d0
       qzz = 0.d0
       r2max = 0.d0
       do i = l, u
          dx = xt(i) - cx
          dy = yt(i) - cy
          dz = zt(i) - cz
          qxx = qxx + mt(i)*dx*dx
          qxy = qxy + mt(i)*dx*dy
          qxz = qxz + mt(i)*dx*dz
          qyy = qyy + mt(i)*dy*dy
          qyz = qyz + mt(i)*dy*dz
          qzz = qzz + mt(i)*dz*dz
          r2max = max(r2max, dx*dx + dy*dy + dz*dz)
       end do
       nd(1,id) = cx
       nd(2,id) = cy
       nd(3,id) = cz
       nd(4,id) = sqrt(r2max)
       nd(5,id) = hm
       nd(6,id) = m
       nd(7,id) = qxx
       nd(8,id) = qxy
       nd(9,id) = qxz
       nd(10,id) = qyy
       nd(11,id) = qyz
       nd(12,id) = qzz
    else
       a = gtw%lc(id)
       b = gtw%rc(id)
       ma = nd(6,a)
       mb = nd(6,b)
       m = ma + mb
       cx = (ma*nd(1,a) + mb*nd(1,b))/m
       cy = (ma*nd(2,a) + mb*nd(2,b))/m
       cz = (ma*nd(3,a) + mb*nd(3,b))/m
       dax = nd(1,a) - cx
       day = nd(2,a) - cy
       daz = nd(3,a) - cz
       dbx = nd(1,b) - cx
       dby = nd(2,b) - cy
       dbz = nd(3,b) - cz
       nd(1,id) = cx
       nd(2,id) = cy
       nd(3,id) = cz
       sa = sqrt(dax*dax + day*day + daz*daz) + nd(4,a)
       sb = sqrt(dbx*dbx + dby*dby + dbz*dbz) + nd(4,b)
       nd(4,id) = max(sa, sb)
       nd(5,id) = max(nd(5,a), nd(5,b))
       nd(6,id) = m
       nd(7,id) = nd(7,a) + nd(7,b) + ma*dax*dax + mb*dbx*dbx
       nd(8,id) = nd(8,a) + nd(8,b) + ma*dax*day + mb*dbx*dby
       nd(9,id) = nd(9,a) + nd(9,b) + ma*dax*daz + mb*dbx*dbz
       nd(10,id) = nd(10,a) + nd(10,b) + ma*day*day + mb*dby*dby
       nd(11,id) = nd(11,a) + nd(11,b) + ma*day*daz + mb*dby*dbz
       nd(12,id) = nd(12,a) + nd(12,b) + ma*daz*daz + mb*dbz*dbz
    end if
  end subroutine node_moments

!     the leaves, in the order of the subtrees (then the ones above them): leaves(1:nlf)
  subroutine leaflist(nlf)
    integer, intent(out) :: nlf
    integer :: j, k, c, nsub
    nsub = gtw%nsub
    if (size(leafcnt) < nsub + 1) then
       deallocate (leafcnt, leafoff)
       allocate (leafcnt(nsub+64), leafoff(nsub+64))
    end if
!$omp parallel default(shared) private(j,k,c)
!$omp do schedule(static)
    do j = 1, nsub
       c = 0
       if (gtw%lc(gtw%sub_id(j)) == 0) c = 1
       do k = gtw%sub_base(j) + 1, gtw%sub_base(j) + gtw%sub_used(j)
          if (gtw%lc(k) == 0) c = c + 1
       end do
       leafcnt(j) = c
    end do
!$omp end do
!$omp end parallel
    nlf = 0
    do j = 1, nsub
       leafoff(j) = nlf
       nlf = nlf + leafcnt(j)
    end do
!$omp parallel default(shared) private(j,k,c)
!$omp do schedule(static)
    do j = 1, nsub
       c = leafoff(j)
       if (gtw%lc(gtw%sub_id(j)) == 0) then
          c = c + 1
          leaves(c) = gtw%sub_id(j)
       end if
       do k = gtw%sub_base(j) + 1, gtw%sub_base(j) + gtw%sub_used(j)
          if (gtw%lc(k) == 0) then
             c = c + 1
             leaves(c) = k
          end if
       end do
    end do
!$omp end do
!$omp end parallel
  end subroutine leaflist

!     ngravtree=2.  The visit of the tree on the gpus, for the particles ilow..ihigh (and, with block steps, the active ones).
!
!     What the library needs (see grav_tree_gpu.cu), in host buffers that it gives: the particles in the order of the tree as
!     pairs of
!     floats (hi, lo) of their position relative to the centre of mass of the system, their mass and 2h (negative for a point
!     particle); a flag for the ones whose force is wanted; the nodes, without the holes between the subtrees (the numbers of the
!     nodes of the subtrees have bases of their own, two for each particle): the top nodes first, then the nodes of the subtrees in
!     the order of the subtrees, the root being node 0; and the leaves to do, in the order of their first particle.  The result
!     comes
!     back for the particles in the order of the tree.
  subroutine gpu_visit(n, ilow, ihigh, act, blk, nk, theta2, gx, gy, gz, gpot)
    integer, intent(in) :: n, ilow, ihigh, nk
    real(dp), intent(in) :: theta2
    logical, intent(in) :: act(n), blk
    real(dp), intent(inout) :: gx(n), gy(n), gz(n), gpot(n)
    call gpu_visit_start(n, ilow, ihigh, act, blk, nk, theta2)
    call gpu_visit_finish(gx, gy, gz, gpot)
  end subroutine gpu_visit

!     The first half: the buffers of the library are filled and the visit is queued on the gpus (the call returns before it is
!     done).
  subroutine gpu_visit_start(n, ilow, ihigh, act, blk, nk, theta2)
    integer, intent(in) :: n, ilow, ihigh, nk
    real(dp), intent(in) :: theta2
    logical, intent(in) :: act(n), blk
    integer :: ntop, nsub, nnode, nwork, nblk, i, j, k, b, c, p, id, ig, hi, base, newb
    real(dp) :: xr, yr, zr, d
    type(c_ptr) :: cpa, cpb, cpneed, cpnf, cpni, cpwork, cpout
    real(c_float), pointer, contiguous :: pa(:,:), pb(:,:), nf(:,:)
    integer(c_signed_char), pointer, contiguous :: needf(:)
    integer(c_int), pointer, contiguous :: ni(:,:), work(:)
    real(c_double), pointer, contiguous :: out(:,:)
    real(c_float) :: fh

    if (gt_pend) call gpu_visit_discard()
    ntop = gtw%ntop
    nsub = gtw%nsub
    if (.not. allocated(soff)) then
       allocate (soff(nsub+64))
    else if (size(soff) < nsub + 1) then
       deallocate (soff)
       allocate (soff(nsub+64))
    end if
    ! the nodes of subtree j (they follow its root, which is one of the top nodes) are the numbers soff(j) .. soff(j+1)-1 (from 0)
    soff(1) = ntop
    do j = 1, nsub
       soff(j+1) = soff(j) + gtw%sub_used(j)
    end do
    nnode = soff(nsub+1)

    call gpu_tree_reserve(int(n, c_int), int(nnode, c_int), cpa, cpb, cpneed, cpnf, cpni, cpwork, cpout)
    call c_f_pointer(cpa, pa, [4, n])
    call c_f_pointer(cpb, pb, [4, n])
    call c_f_pointer(cpneed, needf, [n])
    call c_f_pointer(cpnf, nf, [16, nnode])
    call c_f_pointer(cpni, ni, [4, nnode])
    call c_f_pointer(cpwork, work, [n])
    call c_f_pointer(cpout, out, [4, n])
    if (.not. allocated(lmark)) then
       allocate (lmark(n))
    else if (size(lmark) < n) then
       deallocate (lmark)
       allocate (lmark(n))
    end if
    xr = nd(1,1)
    yr = nd(2,1)
    zr = nd(3,1)

!$omp parallel default(shared) private(i, j, k, ig, fh, d, id, base, newb)
    ! the particles: position (hi, lo), mass, 2h, and whether the force is wanted
!$omp do schedule(static)
    do i = 1, n
       lmark(i) = 0
       d = xt(i) - xr
       fh = real(d, c_float)
       pa(1,i) = fh
       pb(1,i) = real(d - real(fh, dp), c_float)
       d = yt(i) - yr
       fh = real(d, c_float)
       pa(2,i) = fh
       pb(2,i) = real(d - real(fh, dp), c_float)
       d = zt(i) - zr
       fh = real(d, c_float)
       pa(3,i) = fh
       pb(3,i) = real(d - real(fh, dp), c_float)
       pa(4,i) = real(mt(i), c_float)
       if (ft(i) /= 0.d0) then
          pb(4,i) = real(ht(i), c_float)
       else
          pb(4,i) = -real(ht(i), c_float)
       end if
       ig = gind(i)
       needf(i) = 0
       if (ig >= ilow .and. ig <= ihigh) then
          if (.not. blk) then
             needf(i) = 1
          else if (act(ig)) then
             needf(i) = 1
          end if
       end if
    end do
!$omp end do

    ! the top nodes whose children are top nodes too (or that are leaves)
!$omp do schedule(static)
    do id = 1, ntop
       if (gtw%lc(id) == 0) then
          call put_node(id, id - 1, -1, -1, xr, yr, zr, nf, ni)
       else if (gtw%lc(id) <= ntop) then
          call put_node(id, id - 1, gtw%lc(id) - 1, gtw%rc(id) - 1, xr, yr, zr, nf, ni)
       end if
    end do
!$omp end do

    ! the subtrees: the root if it was divided, and the nodes that follow it
!$omp do schedule(dynamic,1)
    do j = 1, nsub
       base = gtw%sub_base(j)
       newb = soff(j) - base - 1
       id = gtw%sub_id(j)
       if (gtw%lc(id) /= 0) call put_node(id, id - 1, gtw%lc(id) + newb, gtw%rc(id) + newb, xr, yr, zr, nf, ni)
       do k = base + 1, base + gtw%sub_used(j)
          if (gtw%lc(k) == 0) then
             call put_node(k, k + newb, -1, -1, xr, yr, zr, nf, ni)
          else
             call put_node(k, k + newb, gtw%lc(k) + newb, gtw%rc(k) + newb, xr, yr, zr, nf, ni)
          end if
       end do
    end do
!$omp end do
!$omp end parallel

    ! the leaves to do, in the order of their first particle: blocks of 4096 places, counted, then written
    nblk = (n + 4095)/4096
    if (.not. allocated(bcnt)) then
       allocate (bcnt(nblk+64), boff(nblk+64))
    else if (size(bcnt) < nblk + 1) then
       deallocate (bcnt, boff)
       allocate (bcnt(nblk+64), boff(nblk+64))
    end if
!$omp parallel do schedule(static) private(b, c, p, hi)
    do b = 1, nblk
       c = 0
       do p = (b - 1)*4096 + 1, min(n, b*4096)
          if (lmark(p) > 0) then
             hi = ni(4, lmark(p)) + 1
             if (any(needf(p:hi) /= 0)) then
                c = c + 1
             else
                lmark(p) = -lmark(p)
             end if
          end if
       end do
       bcnt(b) = c
    end do
    nwork = 0
    do b = 1, nblk
       boff(b) = nwork
       nwork = nwork + bcnt(b)
    end do
!$omp parallel do schedule(static) private(b, c, p)
    do b = 1, nblk
       c = boff(b)
       do p = (b - 1)*4096 + 1, min(n, b*4096)
          if (lmark(p) > 0) then
             c = c + 1
             work(c) = lmark(p) - 1
          end if
       end do
    end do

    ! the visit, queued
    if (nwork > 0) &
         call gpu_tree_launch(int(n, c_int), int(nnode, c_int), int(nwork, c_int), int(nk, c_int), real(theta2, c_double))
    ps_n = n
    ps_ilow = ilow
    ps_ihigh = ihigh
    ps_nwork = nwork
    ps_nnode = nnode
    ps_needf => needf
    ps_out => out
    gt_pend = .true.
  end subroutine gpu_visit_start

!     The second half: wait for the gpus and put the results in gx, gy, gz, gpot (those of ilow..ihigh are zero but for the
!     particles
!     that were done).
  subroutine gpu_visit_finish(gx, gy, gz, gpot)
    real(dp), intent(inout) :: gx(*), gy(*), gz(*), gpot(*)
    integer :: i, ig, ierr
    if (.not. gt_pend) stop 'gravtree: gpu_visit_finish without a visit'
    if (ps_nwork > 0) then
       call gpu_tree_wait(ierr)
       if (ierr == 1) then
          write (6, '(a)') 'gravtree: the stack of the tree visit on the gpu had no room (see TREE_STK in grav_tree_gpu.cu)'
          stop 'gravtree: no room in the lists of the visit on the gpu'
       else if (ierr /= 0) then
          write (6, '(a)') 'gravtree: the pool of the chunks of the leaves with a long list had no room on the gpu'
          stop 'gravtree: no room in the pool of the visit on the gpu'
       end if
    end if
    gt_pend = .false.

    ! the results back to the order of the particles
!$omp parallel do schedule(static) private(i)
    do i = ps_ilow, ps_ihigh
       gx(i) = 0.d0
       gy(i) = 0.d0
       gz(i) = 0.d0
       gpot(i) = 0.d0
    end do
!$omp parallel do schedule(static) private(i, ig)
    do i = 1, ps_n
       if (ps_needf(i) /= 0) then
          ig = gind(i)
          gx(ig) = ps_out(1,i)
          gy(ig) = ps_out(2,i)
          gz(ig) = ps_out(3,i)
          gpot(ig) = ps_out(4,i)
       end if
    end do
  end subroutine gpu_visit_finish

!     A visit that was queued is not wanted (a new one starts): wait for it.
  subroutine gpu_visit_discard()
    integer :: ierr
    if (ps_nwork > 0) call gpu_tree_wait(ierr)
    gt_pend = .false.
  end subroutine gpu_visit_discard

!     the record of node old (its number from 0 in the numbering without holes: nw; its children, from 0, or -1 for a leaf), in
!     the buffers of the library: 16 reals (centre of mass relative to (xr,yr,zr) as hi and lo, size, hmax, mass, second moments), 4
!     integers (children, first and last particle from 0)
  subroutine put_node(old, nw, lcn, rcn, xr, yr, zr, nf, ni)
    integer, intent(in) :: old, nw, lcn, rcn
    real(dp), intent(in) :: xr, yr, zr
    real(c_float), intent(inout) :: nf(16,*)
    integer(c_int), intent(inout) :: ni(4,*)
    real(dp) :: ex, ey, ez
    real(c_float) :: hx, hy, hz
    integer :: q
    ex = nd(1,old) - xr
    ey = nd(2,old) - yr
    ez = nd(3,old) - zr
    hx = real(ex, c_float)
    hy = real(ey, c_float)
    hz = real(ez, c_float)
    nf(1,nw+1) = hx
    nf(2,nw+1) = hy
    nf(3,nw+1) = hz
    nf(4,nw+1) = real(nd(4,old), c_float)
    nf(5,nw+1) = real(ex - real(hx, dp), c_float)
    nf(6,nw+1) = real(ey - real(hy, dp), c_float)
    nf(7,nw+1) = real(ez - real(hz, dp), c_float)
    nf(8,nw+1) = real(nd(5,old), c_float)
    do q = 6, 12
       nf(q+3,nw+1) = real(nd(q,old), c_float)
    end do
    nf(16,nw+1) = 0.0
    ni(1,nw+1) = lcn
    ni(2,nw+1) = rcn
    ni(3,nw+1) = gtw%nl(old) - 1
    ni(4,nw+1) = gtw%nu(old) - 1
    if (lcn < 0) lmark(gtw%nl(old)) = nw + 1
  end subroutine put_node

!     The force on the particles of leaf lf that are in ilow..ihigh (and active).
!
!     The tree is walked once for the leaf.  The nodes that are well separated go to the list m (their quadrupoles are summed
!     into the Taylor coefficients in a vectorised loop); the particles of the leaves that are too near go to the list n if
!     they cannot be inside a softening of any particle of this leaf (their distance from the centre of the leaf is more than
!     size_L + max(hmax_L, 2h_j)), to the list c otherwise.  The force on a particle of the leaf is then the Taylor series, a
!     branch-free sum over n, and a sum over c in which the pairs inside a softening are done together in a vectorised loop.
  subroutine walk_leaf(lf, nk, theta2, ilow, ihigh, act, blk, gx, gy, gz, gpot, n)
    integer, intent(in) :: lf, nk, ilow, ihigh, n
    real(dp), intent(in) :: theta2
    logical, intent(in) :: act(n), blk
    real(dp), intent(inout) :: gx(n), gy(n), gz(n), gpot(n)
    integer :: it, l, u, i, ig, s, ns, a, b, jj, nm, nn, nc, nmp, nnp, ncp
    real(dp) :: xl, yl, zl, szl, hml, dx, dy, dz, r2, sz, rc, fnode(lenf), thr, rj2
    real(dp) :: xi, yi, zi, h2i, hi, ax, ay, az, pt
    logical :: need
    it = 0
!$  it = omp_get_thread_num()
    l = gtw%nl(lf)
    u = gtw%nu(lf)
    ! is there a particle to do in this leaf?
    need = .false.
    do i = l, u
       ig = gind(i)
       if (ig >= ilow .and. ig <= ihigh) then
          if (.not. blk) then
             need = .true.
          else if (act(ig)) then
             need = .true.
          end if
       end if
    end do
    if (.not. need) return
    xl = nd(1,lf)
    yl = nd(2,lf)
    zl = nd(3,lf)
    szl = nd(4,lf)
    hml = nd(5,lf)
    nm = 0
    nn = 0
    nc = 0
    ns = 1
    wb(it)%stack(1) = 1
    do while (ns > 0)
       s = wb(it)%stack(ns)
       ns = ns - 1
       dx = xl - nd(1,s)
       dy = yl - nd(2,s)
       dz = zl - nd(3,s)
       r2 = dx*dx + dy*dy + dz*dz
       sz = szl + nd(4,s)
       rc = sz + max(hml, nd(5,s))
       if (theta2*r2 > sz*sz .and. r2 > rc*rc) then
          ! well separated: to the list of the nodes
          if (nm + 1 > wb(it)%capm) call grow_m(it, nm + 1)
          nm = nm + 1
          wb(it)%mb(nm,1) = dx
          wb(it)%mb(nm,2) = dy
          wb(it)%mb(nm,3) = dz
          wb(it)%mb(nm,4) = nd(6,s)
          wb(it)%mb(nm,5) = nd(7,s)
          wb(it)%mb(nm,6) = nd(8,s)
          wb(it)%mb(nm,7) = nd(9,s)
          wb(it)%mb(nm,8) = nd(10,s)
          wb(it)%mb(nm,9) = nd(11,s)
          wb(it)%mb(nm,10) = nd(12,s)
       else if (gtw%lc(s) == 0) then
          ! a leaf too near: its particles are summed directly
          a = gtw%nl(s)
          b = gtw%nu(s)
          if (nn + (b - a + 1) > wb(it)%capn) call grow_n(it, nn + (b - a + 1))
          if (nc + (b - a + 1) > wb(it)%capc) call grow_c(it, nc + (b - a + 1))
          do jj = a, b
             dx = xt(jj) - xl
             dy = yt(jj) - yl
             dz = zt(jj) - zl
             rj2 = dx*dx + dy*dy + dz*dz
             thr = (szl + max(hml, ht(jj)))*(1.d0 + 1.d-9)
             if (rj2 >= thr*thr) then
                nn = nn + 1
                wb(it)%nb(nn,1) = xt(jj)
                wb(it)%nb(nn,2) = yt(jj)
                wb(it)%nb(nn,3) = zt(jj)
                wb(it)%nb(nn,4) = mt(jj)
             else
                nc = nc + 1
                wb(it)%cb(nc,1) = xt(jj)
                wb(it)%cb(nc,2) = yt(jj)
                wb(it)%cb(nc,3) = zt(jj)
                wb(it)%cb(nc,4) = mt(jj)
                wb(it)%cb(nc,5) = ht(jj)
                wb(it)%cb(nc,6) = h2t(jj)
                wb(it)%cb(nc,7) = ft(jj)
                wb(it)%ci(nc) = jj
             end if
          end do
       else
          if (ns + 2 > stackmax) stop 'gravtree: the stack of the walk of the tree is too small'
          ns = ns + 1
          wb(it)%stack(ns) = gtw%rc(s)
          ns = ns + 1
          wb(it)%stack(ns) = gtw%lc(s)
       end if
    end do
    ! the Taylor coefficients of the far field
    nmp = 8*((nm + 7)/8)
    if (nm > 0) then
       if (nmp > wb(it)%capm) call grow_m(it, nmp)
       do jj = nm + 1, nmp
          wb(it)%mb(jj,1) = 1.d0
          wb(it)%mb(jj,2) = 0.d0
          wb(it)%mb(jj,3) = 0.d0
          wb(it)%mb(jj,4:10) = 0.d0
       end do
    end if
    call m2l_loop(nmp, size(wb(it)%mb,1), wb(it)%mb, fnode)
    ! sources without mass, far away, to fill the last group of 8 of each list
    nnp = 8*((nn + 7)/8)
    ncp = 8*((nc + 7)/8)
    if (nnp > wb(it)%capn) call grow_n(it, nnp)
    if (ncp > wb(it)%capc) call grow_c(it, ncp)
    do jj = nn + 1, nnp
       wb(it)%nb(jj,1) = 1.d100
       wb(it)%nb(jj,2) = 0.d0
       wb(it)%nb(jj,3) = 0.d0
       wb(it)%nb(jj,4) = 0.d0
    end do
    do jj = nc + 1, ncp
       wb(it)%cb(jj,1) = 1.d100
       wb(it)%cb(jj,2) = 0.d0
       wb(it)%cb(jj,3) = 0.d0
       wb(it)%cb(jj,4) = 0.d0
       wb(it)%cb(jj,5) = 1.d0
       wb(it)%cb(jj,6) = 0.d0
       wb(it)%cb(jj,7) = 0.d0
       wb(it)%ci(jj) = 0
    end do
    do i = l, u
       ig = gind(i)
       if (ig < ilow .or. ig > ihigh) cycle
       if (blk) then
          if (.not. act(ig)) cycle
       end if
       xi = xt(i)
       yi = yt(i)
       zi = zt(i)
       hi = ht(i)
       h2i = h2t(i)
       call expand_taylor(fnode, xi - xl, yi - yl, zi - zl, ax, ay, az, pt)
       call near_sum(it, nn, nc, nnp, ncp, i, xi, yi, zi, h2i, hi, nk, ax, ay, az, pt)
       gx(ig) = ax
       gy(ig) = ay
       gz(ig) = az
       gpot(ig) = pt
    end do
  end subroutine walk_leaf

!     The direct sum for the particle i over the lists n (all the pairs are Newtonian) and c of the thread it.
  subroutine near_sum(it, nn, nc, nnp, ncp, i, xi, yi, zi, h2i, hi, nk, ax, ay, az, pt)
    integer, intent(in) :: it, nn, nc, nnp, ncp, i, nk
    real(dp), intent(in) :: xi, yi, zi, h2i, hi
    real(dp), intent(inout) :: ax, ay, az, pt
    integer :: jj, ns, nsp
    real(dp) :: sx, sy, sz, sp, r2, hm, gacc, gpt, selff
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    sp = 0.d0
    call newton_all(nnp, size(wb(it)%nb,1), xi, yi, zi, wb(it)%nb, sx, sy, sz, sp)
    call newton_masked(ncp, size(wb(it)%cb,1), xi, yi, zi, h2i, wb(it)%cb, sx, sy, sz, sp)
    ! the pairs inside a softening: collected, then done together
    ns = 0
    do jj = 1, nc
       r2 = wb(it)%cb(jj,8)
       hm = max(h2i, wb(it)%cb(jj,6))
       if (r2 < hm) then
          if (r2 <= 0.d0 .or. nk /= 2) then
             ! at zero distance (the particle itself), and the kernels other than the Wendland C4, one by one
             selff = 0.5d0
             if (wb(it)%ci(jj) == i) selff = 1.d0
             call grav_soft(nk, r2, hi, wb(it)%cb(jj,5), wb(it)%cb(jj,4), wb(it)%cb(jj,7), selff, gacc, gpt)
             sx = sx + gacc*(wb(it)%cb(jj,1) - xi)
             sy = sy + gacc*(wb(it)%cb(jj,2) - yi)
             sz = sz + gacc*(wb(it)%cb(jj,3) - zi)
             sp = sp + gpt
          else
             ns = ns + 1
             wb(it)%sb(ns,1) = r2
             wb(it)%sb(ns,2) = wb(it)%cb(jj,1) - xi
             wb(it)%sb(ns,3) = wb(it)%cb(jj,2) - yi
             wb(it)%sb(ns,4) = wb(it)%cb(jj,3) - zi
             wb(it)%sb(ns,5) = wb(it)%cb(jj,5)
             wb(it)%sb(ns,6) = wb(it)%cb(jj,4)
          end if
       end if
    end do
    if (ns > 0) then
       nsp = 8*((ns + 7)/8)
       do jj = ns + 1, nsp
          wb(it)%sb(jj,1) = 1.d0
          wb(it)%sb(jj,2) = 1.d0
          wb(it)%sb(jj,3) = 0.d0
          wb(it)%sb(jj,4) = 0.d0
          wb(it)%sb(jj,5) = 1.d0
          wb(it)%sb(jj,6) = 0.d0
       end do
       call soft_loop(nsp, size(wb(it)%sb,1), h2i, hi, wb(it)%sb, sx, sy, sz, sp)
    end if
    ax = ax + sx
    ay = ay + sy
    az = az + sz
    pt = pt + sp
  end subroutine near_sum

!     the Newtonian sums over sources b(:,1:4) = x y z m (the number a multiple of 8), in 8 lanes added in a fixed order,
!     so that they do not depend on how the compiler peels or splits the loop, or on where the arrays are in memory
  subroutine newton_all(nlp, ld, xi, yi, zi, b, sx, sy, sz, sp)
    integer, intent(in) :: nlp, ld
    real(dp), intent(in) :: xi, yi, zi
    real(dp), intent(in) :: b(ld,4)
    real(dp), intent(inout) :: sx, sy, sz, sp
    integer :: jb, l, jj
    real(dp) :: dx, dy, dz, r2, rinv, mr1, mr3
    real(dp) :: accx(8), accy(8), accz(8), accp(8)
    accx = 0.d0
    accy = 0.d0
    accz = 0.d0
    accp = 0.d0
    do jb = 1, nlp, 8
!$omp simd private(jj,dx,dy,dz,r2,rinv,mr1,mr3)
       do l = 1, 8
          jj = jb + l - 1
          dx = b(jj,1) - xi
          dy = b(jj,2) - yi
          dz = b(jj,3) - zi
          r2 = dx*dx + dy*dy + dz*dz
          rinv = 1.d0/sqrt(r2)
          mr1 = rinv*b(jj,4)
          mr3 = rinv*rinv*mr1
          accx(l) = accx(l) + mr3*dx
          accy(l) = accy(l) + mr3*dy
          accz(l) = accz(l) + mr3*dz
          accp(l) = accp(l) - mr1
       end do
    end do
    sx = sx + (((accx(1) + accx(2)) + (accx(3) + accx(4))) + ((accx(5) + accx(6)) + (accx(7) + accx(8))))
    sy = sy + (((accy(1) + accy(2)) + (accy(3) + accy(4))) + ((accy(5) + accy(6)) + (accy(7) + accy(8))))
    sz = sz + (((accz(1) + accz(2)) + (accz(3) + accz(4))) + ((accz(5) + accz(6)) + (accz(7) + accz(8))))
    sp = sp + (((accp(1) + accp(2)) + (accp(3) + accp(4))) + ((accp(5) + accp(6)) + (accp(7) + accp(8))))
  end subroutine newton_all

!     r^2 of the pairs with the sources c(:,1:7) = x y z m 2h (2h)^2 flag, kept in c(:,8), and the Newtonian sums over the pairs
!     that are not inside a softening
  subroutine newton_masked(nlp, ld, xi, yi, zi, h2i, c, sx, sy, sz, sp)
    integer, intent(in) :: nlp, ld
    real(dp), intent(in) :: xi, yi, zi, h2i
    real(dp), intent(inout) :: c(ld,8)
    real(dp), intent(inout) :: sx, sy, sz, sp
    integer :: jb, l, jj
    real(dp) :: dx, dy, dz, r2, rinv, mr1, mr3, hm
    real(dp) :: accx(8), accy(8), accz(8), accp(8)
    logical :: nwt
    accx = 0.d0
    accy = 0.d0
    accz = 0.d0
    accp = 0.d0
    do jb = 1, nlp, 8
!$omp simd private(jj,dx,dy,dz)
       do l = 1, 8
          jj = jb + l - 1
          dx = c(jj,1) - xi
          dy = c(jj,2) - yi
          dz = c(jj,3) - zi
          c(jj,8) = dx*dx + dy*dy + dz*dz
       end do
!$omp simd private(jj,dx,dy,dz,r2,rinv,mr1,mr3,nwt,hm)
       do l = 1, 8
          jj = jb + l - 1
          dx = c(jj,1) - xi
          dy = c(jj,2) - yi
          dz = c(jj,3) - zi
          r2 = c(jj,8)
          rinv = 1.d0/sqrt(r2)
          mr1 = rinv*c(jj,4)
          mr3 = rinv*rinv*mr1
          hm = merge(h2i, c(jj,6), h2i > c(jj,6))
          nwt = (r2 >= hm)
          if (nwt) then
             accx(l) = accx(l) + mr3*dx
             accy(l) = accy(l) + mr3*dy
             accz(l) = accz(l) + mr3*dz
             accp(l) = accp(l) - mr1
          end if
       end do
    end do
    sx = sx + (((accx(1) + accx(2)) + (accx(3) + accx(4))) + ((accx(5) + accx(6)) + (accx(7) + accx(8))))
    sy = sy + (((accy(1) + accy(2)) + (accy(3) + accy(4))) + ((accy(5) + accy(6)) + (accy(7) + accy(8))))
    sz = sz + (((accz(1) + accz(2)) + (accz(3) + accz(4))) + ((accz(5) + accz(6)) + (accz(7) + accz(8))))
    sp = sp + (((accp(1) + accp(2)) + (accp(3) + accp(4))) + ((accp(5) + accp(6)) + (accp(7) + accp(8))))
  end subroutine newton_masked

!     The pairs inside a softening, Wendland C4 kernel (nkernel=2), for the particle with 2h = hi, 8 at a time:
!     s(:,1:6) = r^2, dx, dy, dz, 2h_j, m_j.  The same expressions as grav_soft (gravpair.f90), without the branches.
  subroutine soft_loop(nsp, ld, h2i, hi, s, sx, sy, sz, sp)
    integer, intent(in) :: nsp, ld
    real(dp), intent(in) :: h2i, hi
    real(dp), intent(in) :: s(ld,6)
    real(dp), intent(inout) :: sx, sy, sz, sp
    integer :: jb, l, jj
    real(dp) :: r2, rinv, mrinv1, mrinv3, qq1, qq2, q21, q22, acc1, acc2, pot1, pot2, hj, h2j, mj
    real(dp) :: mj11, mj12, mj21, mj22, gacc, gpot, dx, dy, dz
    logical :: in1, in2
    real(dp) :: accx(8), accy(8), accz(8), accp(8)
    accx = 0.d0
    accy = 0.d0
    accz = 0.d0
    accp = 0.d0
    do jb = 1, nsp, 8
!$omp simd private(jj,r2,rinv,mrinv1,mrinv3,qq1,qq2,q21,q22,acc1,acc2,pot1,pot2,hj,h2j,mj,mj11,mj12,mj21,mj22,gacc,gpot, &
!$omp dx,dy,dz,in1,in2)
       do l = 1, 8
          jj = jb + l - 1
          r2 = s(jj,1)
          dx = s(jj,2)
          dy = s(jj,3)
          dz = s(jj,4)
          hj = s(jj,5)
          mj = s(jj,6)
          h2j = hj*hj
          rinv = 1.d0/sqrt(r2)
          mrinv1 = rinv*mj
          mrinv3 = rinv*rinv*mrinv1
          qq1 = 1.d0/(rinv*hi)
          qq2 = 1.d0/(rinv*hj)
          q21 = qq1*qq1
          q22 = qq2*qq2
          mj11 = mj/hi
          mj12 = mj/hj
          mj21 = mj11/h2i
          mj22 = mj12/h2j
          in1 = (r2 <= h2i)
          in2 = (r2 <= h2j)
          acc1 = 20.625d0 + q21*((-115.5d0) + q21*(618.75d0 + qq1*((-1155.d0) + qq1*(962.5d0 + qq1*((-396.d0) + 65.625d0*qq1)))))
          acc2 = 20.625d0 + q22*((-115.5d0) + q22*(618.75d0 + qq2*((-1155.d0) + qq2*(962.5d0 + qq2*((-396.d0) + 65.625d0*qq2)))))
          pot1 = (-3.4375d0) + q21*(10.3125d0 + q21*((-28.875d0) + q21*(103.125d0 + qq1*((-165.d0) &
               + qq1*(120.3125d0 + qq1*((-44.d0) + 6.5625d0*qq1))))))
          pot2 = (-3.4375d0) + q22*(10.3125d0 + q22*((-28.875d0) + q22*(103.125d0 + qq2*((-165.d0) &
               + qq2*(120.3125d0 + qq2*((-44.d0) + 6.5625d0*qq2))))))
          gacc = 0.5d0*(merge(mj21*acc1, mrinv3, in1) + merge(mj22*acc2, mrinv3, in2))
          gpot = 0.5d0*(merge(mj11*pot1, -mrinv1, in1) + merge(mj12*pot2, -mrinv1, in2))
          accx(l) = accx(l) + gacc*dx
          accy(l) = accy(l) + gacc*dy
          accz(l) = accz(l) + gacc*dz
          accp(l) = accp(l) + gpot
       end do
    end do
    sx = sx + (((accx(1) + accx(2)) + (accx(3) + accx(4))) + ((accx(5) + accx(6)) + (accx(7) + accx(8))))
    sy = sy + (((accy(1) + accy(2)) + (accy(3) + accy(4))) + ((accy(5) + accy(6)) + (accy(7) + accy(8))))
    sz = sz + (((accz(1) + accz(2)) + (accz(3) + accz(4))) + ((accz(5) + accz(6)) + (accz(7) + accz(8))))
    sp = sp + (((accp(1) + accp(2)) + (accp(3) + accp(4))) + ((accp(5) + accp(6)) + (accp(7) + accp(8))))
  end subroutine soft_loop

!     The Taylor coefficients of the force at the centre of the leaf from the nodes of the list m(:,1:10) = dx dy dz q0 qxx qxy qxz
!     qyy qyz qzz (the number a multiple of 8), 8 nodes at a time with the sums kept in 8 lanes.
!     (Adapted from Phantom (Price et al. 2018, PASA 35, e031), src/main/kdtree.F90, routine compute_M2L, after Dehnen 2002; GPLv3.)
!     fnode(1:3) the acceleration, (4:9) its derivatives xx xy xz yy yz zz, (10:19) its second derivatives
!     xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz, (20) the potential.
  subroutine m2l_loop(nmp, ld, m, fnode)
    integer, intent(in) :: nmp, ld
    real(dp), intent(in) :: m(ld,10)
    real(dp), intent(inout) :: fnode(lenf)
    integer :: jb, l, jj, k
    real(dp) :: dx, dy, dz, q0, qxx, qxy, qxz, qyy, qyz, qzz, dx2, dx3, dy2, dy3, dz2, dz3, r2, dr1
    real(dp) :: dr12, g0, g1, g2, g3, g2dx, g2dy, g2dz
    real(dp) :: d31, d32, d33, d34, d35, d36, d37, d38, d39, d310, d21, d22, d23, d24, d25, d26, d11, d12, d13
    real(dp) :: f(8,lenf)
    f = 0.d0
    do jb = 1, nmp, 8
!$omp simd private(jj,dx,dy,dz,q0,qxx,qxy,qxz,qyy,qyz,qzz,dx2,dx3,dy2,dy3,dz2,dz3,r2,dr1,dr12,g0,g1,g2,g3,g2dx,g2dy,g2dz, &
!$omp d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,d21,d22,d23,d24,d25,d26,d11,d12,d13)
       do l = 1, 8
          jj = jb + l - 1
          dx = m(jj,1)
          dy = m(jj,2)
          dz = m(jj,3)
          q0 = m(jj,4)
          qxx = m(jj,5)
          qxy = m(jj,6)
          qxz = m(jj,7)
          qyy = m(jj,8)
          qyz = m(jj,9)
          qzz = m(jj,10)
          r2 = dx*dx + dy*dy + dz*dz
          dr1 = 1.d0/sqrt(r2)
          dr12 = dr1*dr1
          dx2 = dx*dx
          dx3 = dx*dx2
          dy2 = dy*dy
          dy3 = dy*dy2
          dz2 = dz*dz
          dz3 = dz*dz2
          g0 = -dr1
          g1 = dr12*g0
          g2 = -3.d0*dr12*g1
          g3 = -5.d0*dr12*g2
          g2dx = g2*dx
          g2dy = g2*dy
          g2dz = g2*dz
          d31 = 3.d0*g2dx + g3*dx3
          d32 = g2dy + g3*dx2*dy
          d33 = g2dz + g3*dx2*dz
          d34 = g2dx + g3*dy2*dx
          d35 = g3*dx*dy*dz
          d36 = g2dx + g3*dz2*dx
          d37 = 3.d0*g2dy + g3*dy3
          d38 = g2dz + g3*dy2*dz
          d39 = g2dy + g3*dz2*dy
          d310 = 3.d0*g2dz + g3*dz3
          d21 = g1 + g2*dx2
          d22 = g2dx*dy
          d23 = g2dx*dz
          d24 = g1 + g2*dy2
          d25 = g2dy*dz
          d26 = g1 + g2*dz2
          d11 = g1*dx
          d12 = g1*dy
          d13 = g1*dz
          f(l,1) = f(l,1) + (d11*q0 + 0.5d0*(d31*qxx + 2.d0*(d32*qxy + d33*qxz + d35*qyz) + d34*qyy + d36*qzz))
          f(l,2) = f(l,2) + (d12*q0 + 0.5d0*(d32*qxx + 2.d0*(d34*qxy + d35*qxz + d38*qyz) + d37*qyy + d39*qzz))
          f(l,3) = f(l,3) + (d13*q0 + 0.5d0*(d33*qxx + 2.d0*(d35*qxy + d36*qxz + d39*qyz) + d38*qyy + d310*qzz))
          f(l,4) = f(l,4) + d21*q0
          f(l,5) = f(l,5) + d22*q0
          f(l,6) = f(l,6) + d23*q0
          f(l,7) = f(l,7) + d24*q0
          f(l,8) = f(l,8) + d25*q0
          f(l,9) = f(l,9) + d26*q0
          f(l,10) = f(l,10) + d31*q0
          f(l,11) = f(l,11) + d32*q0
          f(l,12) = f(l,12) + d33*q0
          f(l,13) = f(l,13) + d34*q0
          f(l,14) = f(l,14) + d35*q0
          f(l,15) = f(l,15) + d36*q0
          f(l,16) = f(l,16) + d37*q0
          f(l,17) = f(l,17) + d38*q0
          f(l,18) = f(l,18) + d39*q0
          f(l,19) = f(l,19) + d310*q0
          f(l,20) = f(l,20) + g0*q0 - 0.5d0*(d21*qxx + d24*qyy + d26*qzz + 2*(d22*qxy + d23*qxz + d25*qyz))
       end do
    end do
    do k = 1, lenf
       fnode(k) = (((f(1,k) + f(2,k)) + (f(3,k) + f(4,k))) + ((f(5,k) + f(6,k)) + (f(7,k) + f(8,k))))
    end do
  end subroutine m2l_loop

!     more room in the lists of thread it
  subroutine grow_m(it, need)
    integer, intent(in) :: it, need
    integer :: nc
    real(dp), allocatable :: t(:,:)
    nc = max(need, 2*wb(it)%capm)
    allocate (t(nc+8,10))
    t(1:wb(it)%capm,:) = wb(it)%mb(1:wb(it)%capm,:)
    call move_alloc(t, wb(it)%mb)
    wb(it)%capm = nc
  end subroutine grow_m

  subroutine grow_n(it, need)
    integer, intent(in) :: it, need
    integer :: nc
    real(dp), allocatable :: t(:,:)
    nc = max(need, 2*wb(it)%capn)
    allocate (t(nc+8,4))
    t(1:wb(it)%capn,:) = wb(it)%nb(1:wb(it)%capn,:)
    call move_alloc(t, wb(it)%nb)
    wb(it)%capn = nc
  end subroutine grow_n

!     the list c, and the list of the pairs inside a softening, which is never longer than c
  subroutine grow_c(it, need)
    integer, intent(in) :: it, need
    integer :: nc
    real(dp), allocatable :: t(:,:)
    integer, allocatable :: ti(:)
    nc = max(need, 2*wb(it)%capc)
    allocate (t(nc+8,8))
    t(1:wb(it)%capc,:) = wb(it)%cb(1:wb(it)%capc,:)
    call move_alloc(t, wb(it)%cb)
    allocate (ti(nc+8))
    ti(1:wb(it)%capc) = wb(it)%ci(1:wb(it)%capc)
    call move_alloc(ti, wb(it)%ci)
    deallocate (wb(it)%sb)
    allocate (wb(it)%sb(nc+8,6))
    wb(it)%capc = nc
  end subroutine grow_c

!     The force and the potential at (dx,dy,dz) from the centre of the leaf from the Taylor series of the force
!     (adapted from Phantom, src/main/kdtree.F90, routine expand_fgrav_in_taylor_series; GPLv3).  The potential is the integral
!     of the series: Phantom's expression for it has a coefficient 1/2 where the cubic term has 1/6.
  pure subroutine expand_taylor(f, dx, dy, dz, fxi, fyi, fzi, poti)
    real(dp), intent(in) :: f(lenf), dx, dy, dz
    real(dp), intent(out) :: fxi, fyi, fzi, poti
    real(dp) :: fx0, fy0, fz0, dfxx, dfxy, dfxz, dfyy, dfyz, dfzz
    real(dp) :: d2fxxx, d2fxxy, d2fxxz, d2fxyy, d2fxyz, d2fxzz, d2fyyy, d2fyyz, d2fyzz, d2fzzz
    real(dp) :: ux, uy, uz, t3
    fx0 = f(1)
    fy0 = f(2)
    fz0 = f(3)
    dfxx = f(4)
    dfxy = f(5)
    dfxz = f(6)
    dfyy = f(7)
    dfyz = f(8)
    dfzz = f(9)
    d2fxxx = f(10)
    d2fxxy = f(11)
    d2fxxz = f(12)
    d2fxyy = f(13)
    d2fxyz = f(14)
    d2fxzz = f(15)
    d2fyyy = f(16)
    d2fyyz = f(17)
    d2fyzz = f(18)
    d2fzzz = f(19)
    ! u_i = (dx . H_i)(dx) / 2, the second order term of component i
    ux = 0.5d0*(dx*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz) + dy*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz) &
         + dz*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
    uy = 0.5d0*(dx*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz) + dy*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz) &
         + dz*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
    uz = 0.5d0*(dx*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz) + dy*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz) &
         + dz*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))
    fxi = fx0 + dx*dfxx + dy*dfxy + dz*dfxz + ux
    fyi = fy0 + dx*dfxy + dy*dfyy + dz*dfyz + uy
    fzi = fz0 + dx*dfxz + dy*dfyz + dz*dfzz + uz
    ! phi = phi0 - f0.d - (1/2) d.J.d - (1/6) H d d d
    t3 = dx*ux + dy*uy + dz*uz
    poti = f(20) - (fx0*dx + fy0*dy + fz0*dz) &
         - 0.5d0*(dx*(dx*dfxx + dy*dfxy + dz*dfxz) + dy*(dx*dfxy + dy*dfyy + dz*dfyz) + dz*(dx*dfxz + dy*dfyz + dz*dfzz)) &
         - t3/3.d0
  end subroutine expand_taylor

end module gravtree

!***********************************************************************
!     The call from the fixed-form code: the gravity of the particles ngrav_lower..ngrav_upper of this rank in gx, gy, gz, grpot,
!     as lasthalf_grav_forces gives them for the gpu.  Used when ngravtree=1 (the visit on the cpus) and 2 (the visit on the gpus).
!***********************************************************************
subroutine tree_grav_forces(ntot_in, gxo, gyo, gzo, grpoto)
  use gravtree
  include 'starsmasher.h'
  integer ntot_in
  real*8 gxo(ntot_in), gyo(ntot_in), gzo(ntot_in), grpoto(ntot_in)
  logical blk
  if (gt_pend) then
     ! the visit was queued by tree_gpu_start
     call gt_finish(ntot_in, gxo, gyo, gzo, grpoto)
     return
  end if
  blk = (nblock.eq.1)
  if (blk) blk = (count(actblk(1:ntot_in)).lt.ntot_in)
  gt_gpu = (ngravtree.eq.2)
  call gt_forces(ntot_in, x, y, z, am, hp, u, actblk, blk, ngrav_lower, ngrav_upper, nkernel, treeacc, nleafgrav, &
       gxo, gyo, gzo, grpoto)
  return
end subroutine tree_grav_forces

!***********************************************************************
!     ngravtree=2 with ngravoverlap=1, to be called where the direct sum on the gpus is started (gravforce): the tree is made
!     and the
!     visit queued on the gpus; grav_collect finishes it.
!***********************************************************************
subroutine tree_gpu_start
  use gravtree
  include 'starsmasher.h'
  logical blk
  if (ngravtree.ne.2 .or. ngravoverlap.eq.0) return
  if (myrank.ge.ngravprocs) return
  blk = (nblock.eq.1)
  if (blk) blk = (count(actblk(1:ntot)).lt.ntot)
  call gt_start(ntot, x, y, z, am, hp, u, actblk, blk, ngrav_lower, ngrav_upper, nkernel, treeacc, nleafgrav)
  return
end subroutine tree_gpu_start

!***********************************************************************
!     The gravity of the particles ngrav_lower..ngrav_upper of this rank, to be called where the gpu's is collected
!     (lasthalf_grav_forces): from the tree when ngravtree=1 or 2, from the library of the gpu otherwise.
!***********************************************************************
subroutine grav_collect(ntot_in, gxo, gyo, gzo, grpoto)
  include 'starsmasher.h'
  integer ntot_in
  real*8 gxo(ntot_in), gyo(ntot_in), gzo(ntot_in), grpoto(ntot_in)
  if (ngravtree.ne.0) then
     call tree_grav_forces(ntot_in, gxo, gyo, gzo, grpoto)
  else
     call lasthalf_grav_forces(ntot_in, gxo, gyo, gzo, grpoto)
  end if
  return
end subroutine grav_collect
