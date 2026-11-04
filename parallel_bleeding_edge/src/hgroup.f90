!***********************************************************************
!     Neighbour search of the smoothing-length solve (rho_and_h in
!     advance.f90) for groups of particles.
!
!     A group is a few consecutive particles in the order of the leaves of
!     the neighbour tree, so that they lie close together.  One search of
!     the tree finds every particle that is within reach of any of them:
!     those closer to the centre of the group than the radius of its box
!     plus the largest 2h to be tried.  They are put in the order of their
!     index (a bitmap of the indices, read from the lowest word), and their
!     coordinates are copied into arrays of their own.  Each particle of the
!     group then takes from that list, in one pass that the compiler
!     vectorizes, the candidates it needs: those within the radius of its
!     own search, in the order of the index, with their squared distances.
!     That is the list that kdtree2_r_nearest_around_point and
!     sort_candidates would give for the particle, element by element and bit
!     by bit: the squared distance is the sum (dx*dx+dy*dy)+dz*dz of the
!     tree (this file is compiled without fused multiply-add, see HGFLAGS
!     in the Makefile), and the set is every particle at a squared distance
!     not above the radius.  Where the radius of a search goes beyond what
!     the list of the group holds the caller searches the tree for the
!     particle itself.
!
!     The bits of the results therefore do not depend on the groups, nor on
!     the number of threads, provided that the tree search forms the squared
!     distance without fused multiply-add too and the sums of rho_and_h are
!     made the same way on both paths, as the default build for x86-64 does.
!     With a build for the machine it runs on (MARCH) or on aarch64 the
!     compiler may fuse multiplications and additions differently in the two
!     paths, and then a result can differ in its last bit.
!
!     This file is compiled for AVX2 on x86-64 (HGFLAGS in the Makefile).  A
!     processor without it (older than 2013, or a virtual machine that hides
!     it) would stop on the first instruction of the file, so rho_and_h asks
!     hg_cpu_ok (hgcpu.f90, compiled for any x86-64) before it calls anything
!     here, and without avx2 goes on without groups (hgroup=0), which gives
!     the same bits.
!
!     The groups are the input parameter hgroup (particles per group; 0 keeps them off).  A group is halved, again and
!     again, when its members lie far apart (the box is wider than hg_split times the largest radius) or when their
!     smoothing lengths differ by more than a factor hg_hsplit, because the list reaches as far as the member with the
!     largest radius and every member reads all of it.  The threads take the groups hg_chunk at a time, or fewer when
!     the rank has few of them, and the neighbour lists are laid out in the order of the particles (hg_listord), which is
!     the order the loop over the pairs reads them in.
!***********************************************************************
module hgroup_mod
  use kdtree2_module
  implicit none
  private
  public :: hgscratch, hg_tune, hg_order, hg_cands, hg_dist2, hg_compact, hg_listsub
  public :: hg_ngs, hg_gpad, hg_on, hg_ngsmax, hg_chunk, hg_split, hg_hsplit, hg_listord
  integer, parameter :: dp = kind(1.d0)

!     particles per group and whether groups are used (hg_tune)
  integer, save :: hg_ngs = 32
  logical, save :: hg_on = .false.
!     the most groups a thread takes at a time (the number used is the groups of the rank over 16 times the threads, but not more)
  integer, parameter :: hg_chunk = 64
!     a group whose box is wider than hg_split times its largest radius is halved
  real(dp), parameter :: hg_split = 1.5d0
!     a group whose members' largest first guess of h is more than hg_hsplit times the smallest is halved as well (0: no such rule)
  real(dp), parameter :: hg_hsplit = 1.5d0
!     the neighbour lists are laid out in the order of the particles (instead of in the order in which the threads solve them)
  logical, parameter :: hg_listord = .true.
!     the fraction by which the radius of the list of a group exceeds the first search of each member
  real(dp), parameter :: hg_gpad = 0.03d0
!     the most particles in a group
  integer, parameter :: hg_ngsmax = 256

!     the working space of one thread (the records of the threads lie one after the other in an array; the padding keeps two
!     of them from sharing a cache line)
  type hgscratch
     integer :: pad1(32) = 0
!    the candidates of the group, in the order of the index: index, coordinates, squared distance to the particle being solved
     integer :: ncg = 0
     integer, allocatable :: cidx(:)
     real(dp), allocatable :: cx(:), cy(:), cz(:), cd2(:)
!    their masses and gas flags (the u of the particle), copied once for the group
     real(dp), allocatable :: cam(:), cu(:)
!    the candidates of the particle being solved (hg_compact): the place of each in the arrays of the group
     integer, allocatable :: kloc(:)
!    the neighbours of the particle (hg_listsub): squared distance, mass, gas flag, in the order of the list
     integer :: nl = 0
     real(dp), allocatable :: ldis(:), lam(:), lu(:)
!    bit map of the indices (one bit per particle), found by the tree search, and a second map with a bit for each word of the
!    first that is not zero (so that the words to read are found without a pass over all of them)
     integer(8), allocatable :: bm(:), bm1(:)
     type(kdtree2_result), allocatable :: tres(:)
!    the squared distances of the candidates have been found for the particle being solved
     logical :: cd2ok = .false.
!    the candidates of the particle being solved (s%res) come from the list of the group (they are those of the last search, which
!    may have been one of the tree for the particle itself)
     logical :: fromgrp = .false.
!    the list of the group is valid for any search of the particle out to this squared radius
     real(dp) :: gcov2 = 0.d0
     logical :: gactive = .false.
     integer :: pad2(32) = 0
  end type hgscratch

contains

!     the input parameter hgroup: particles per group (at most hg_ngsmax, which init checks), 0 for no groups
  subroutine hg_tune(hgin)
    integer, intent(in) :: hgin
    if (hgin > 0) then
       hg_ngs = min(hgin, hg_ngsmax)
       hg_on = .true.
    else
       hg_on = .false.
    endif
  end subroutine hg_tune

!     The indices ind(1:nt) of the points of a tree in the order of its leaves; gord gets those between lo and hi (inclusive) whose
!     mask(ind+off) is not zero, plus off, in that order, ngl of them.  (Outside the parallel regions: it makes its own.)
  subroutine hg_order(ind, nt, lo, hi, off, mask, gord, ngl)
    integer, intent(in) :: ind(*), nt, lo, hi, off
    integer(1), intent(in) :: mask(*)
    integer, allocatable, intent(inout) :: gord(:)
    integer, intent(out) :: ngl
    integer, parameter :: blk = 8192
    integer :: nb, ib, k, c
    integer, allocatable :: cnt(:)
    nb = (nt + blk - 1)/blk
    allocate(cnt(nb + 1))
    if (allocated(gord)) then
       if (size(gord) < nt) deallocate(gord)
    endif
    if (.not. allocated(gord)) allocate(gord(max(nt, 1)))
!$omp parallel do schedule(static) private(ib, k, c)
    do ib = 1, nb
       c = 0
       do k = (ib - 1)*blk + 1, min(ib*blk, nt)
          if (ind(k) >= lo .and. ind(k) <= hi) then
             if (mask(ind(k) + off) /= 0_1) c = c + 1
          endif
       enddo
       cnt(ib) = c
    enddo
    c = 0
    do ib = 1, nb
       k = cnt(ib)
       cnt(ib) = c
       c = c + k
    enddo
    ngl = c
!$omp parallel do schedule(static) private(ib, k, c)
    do ib = 1, nb
       c = cnt(ib)
       do k = (ib - 1)*blk + 1, min(ib*blk, nt)
          if (ind(k) >= lo .and. ind(k) <= hi) then
             if (mask(ind(k) + off) /= 0_1) then
                c = c + 1
                gord(c) = ind(k) + off
             endif
          endif
       enddo
    enddo
    deallocate(cnt)
  end subroutine hg_order

!     The candidates of the group: the particles of the tree within sqrt(r2g) of the point qc, in the order of the index (global:
!     through map(1:nsub) if it is present, which gives the global index of each point of the tree), with their coordinates.
!     n is the number of particles (indices are 1..n), x, y, z their coordinates.
  subroutine hg_cands(tree, qc, r2g, n, x, y, z, am, uu, g, map)
    type(kdtree2), pointer :: tree
    real(dp), intent(in) :: qc(3), r2g
    integer, intent(in) :: n
    real(dp), intent(in) :: x(*), y(*), z(*), am(*), uu(*)
    type(hgscratch), intent(inout) :: g
    integer, intent(in), optional :: map(*)
    integer :: nf, k, j, w, b, nw, nw1, v, v0, v1, wb
    integer :: nc
    logical :: ovf
    integer(8) :: word, word1
    real(dp) :: qv(3)
    qv = qc
    if (.not. allocated(g%tres)) allocate(g%tres(4096))
    do
       call kdtree2_r_nearest(tp=tree, qv=qv, r2=r2g, nfound=nf, nalloc=size(g%tres), results=g%tres, overflow=ovf)
       if (.not. ovf) exit
       deallocate(g%tres)
       allocate(g%tres(2*nf))
    enddo
    nw = (n + 63)/64 + 1
    nw1 = (nw + 63)/64 + 1
    if (allocated(g%bm)) then
       if (size(g%bm) < nw) deallocate(g%bm, g%bm1)
    endif
    if (.not. allocated(g%bm)) then
       allocate(g%bm(nw), g%bm1(nw1))
       g%bm = 0_8
       g%bm1 = 0_8
    endif
    v0 = nw1
    v1 = 0
    if (present(map)) then
       do k = 1, nf
          j = map(g%tres(k)%idx) - 1
          w = ishft(j, -6) + 1
          g%bm(w) = ior(g%bm(w), ishft(1_8, iand(j, 63)))
          v = ishft(w - 1, -6) + 1
          g%bm1(v) = ior(g%bm1(v), ishft(1_8, iand(w - 1, 63)))
          v0 = min(v0, v)
          v1 = max(v1, v)
       enddo
    else
       do k = 1, nf
          j = g%tres(k)%idx - 1
          w = ishft(j, -6) + 1
          g%bm(w) = ior(g%bm(w), ishft(1_8, iand(j, 63)))
          v = ishft(w - 1, -6) + 1
          g%bm1(v) = ior(g%bm1(v), ishft(1_8, iand(w - 1, 63)))
          v0 = min(v0, v)
          v1 = max(v1, v)
       enddo
    endif
    if (.not. allocated(g%cidx)) then
       allocate(g%cidx(max(nf + 64, 2048)), g%cx(max(nf + 64, 2048)), g%cy(max(nf + 64, 2048)), g%cz(max(nf + 64, 2048)), &
            g%cd2(max(nf + 64, 2048)), g%cam(max(nf + 64, 2048)), g%cu(max(nf + 64, 2048)), g%kloc(max(nf + 64, 2048)))
    else if (size(g%cidx) < nf + 8) then
       deallocate(g%cidx, g%cx, g%cy, g%cz, g%cd2, g%cam, g%cu, g%kloc)
       allocate(g%cidx(2*nf + 64), g%cx(2*nf + 64), g%cy(2*nf + 64), g%cz(2*nf + 64), g%cd2(2*nf + 64), &
            g%cam(2*nf + 64), g%cu(2*nf + 64), g%kloc(2*nf + 64))
    endif
    nc = 0
    do v = v0, v1
       word1 = g%bm1(v)
       if (word1 == 0_8) cycle
       g%bm1(v) = 0_8
       do while (word1 /= 0_8)
          wb = trailz(word1)
          word1 = iand(word1, word1 - 1_8)
          w = (v - 1)*64 + wb + 1
          word = g%bm(w)
          g%bm(w) = 0_8
          do while (word /= 0_8)
             b = trailz(word)
             word = iand(word, word - 1_8)
             nc = nc + 1
             g%cidx(nc) = (w - 1)*64 + b + 1
          enddo
       enddo
    enddo
    g%ncg = nc
    do k = 1, nc
       j = g%cidx(k)
       g%cx(k) = x(j)
       g%cy(k) = y(j)
       g%cz(k) = z(j)
       g%cam(k) = am(j)
       g%cu(k) = uu(j)
    enddo
  end subroutine hg_cands

!     the squared distances from (xi,yi,zi) to the candidates, as the tree finds them: (dx*dx+dy*dy)+dz*dz
  subroutine hg_dist2(g, xi, yi, zi)
    type(hgscratch), intent(inout) :: g
    real(dp), intent(in) :: xi, yi, zi
    call hg_dist2_loop(g%ncg, g%cx, g%cy, g%cz, xi, yi, zi, g%cd2)
  end subroutine hg_dist2

  subroutine hg_dist2_loop(nc, cx, cy, cz, xi, yi, zi, d2)
    integer, intent(in) :: nc
    real(dp), intent(in) :: cx(*), cy(*), cz(*), xi, yi, zi
    real(dp), intent(out) :: d2(*)
    integer :: k
    real(dp) :: dx, dy, dz
!$omp simd
    do k = 1, nc
       dx = cx(k) - xi
       dy = cy(k) - yi
       dz = cz(k) - zi
       d2(k) = dx*dx + dy*dy + dz*dz
    enddo
  end subroutine hg_dist2_loop

!     The candidates at a squared distance not above r2 from the particle (hg_dist2 has been called for it), in the order of
!     the index, put in res(1:nc) as the tree search of kdtree2 gives them: res(k)%idx, res(k)%dis.  res must have room for
!     ncg+1 elements.
  subroutine hg_compact(g, r2, res, nc)
    type(hgscratch), intent(inout) :: g
    real(dp), intent(in) :: r2
    type(kdtree2_result), intent(inout) :: res(*)
    integer, intent(out) :: nc
    integer :: k, m
    m = 0
    do k = 1, g%ncg
       res(m + 1)%idx = g%cidx(k)
       res(m + 1)%dis = g%cd2(k)
       g%kloc(m + 1) = k
       if (g%cd2(k) <= r2) m = m + 1
    enddo
    nc = m
  end subroutine hg_compact

!     The neighbours of the particle: the candidates res(1:ncand) (hg_compact) with a squared distance below r2, in the same order,
!     with that distance, the mass and the gas flag of each, in g%ldis, g%lam and g%lu (g%nl of them).  The final sums of the
!     smoothing-length solve read these three arrays, which lie together, instead of the coordinates, masses and flags of the
!     particles themselves, which lie all over the arrays of the code.
  subroutine hg_listsub(g, res, ncand, r2)
    type(hgscratch), intent(inout) :: g
    type(kdtree2_result), intent(in) :: res(*)
    integer, intent(in) :: ncand
    real(dp), intent(in) :: r2
    integer :: k, m
    if (.not. allocated(g%ldis)) then
       allocate(g%ldis(max(ncand + 16, 512)), g%lam(max(ncand + 16, 512)), g%lu(max(ncand + 16, 512)))
    else if (size(g%ldis) < ncand + 1) then
       deallocate(g%ldis, g%lam, g%lu)
       allocate(g%ldis(2*ncand + 16), g%lam(2*ncand + 16), g%lu(2*ncand + 16))
    endif
    m = 0
    do k = 1, ncand
       if (res(k)%dis < r2) then
          m = m + 1
          g%ldis(m) = res(k)%dis
          g%lam(m) = g%cam(g%kloc(k))
          g%lu(m) = g%cu(g%kloc(k))
       endif
    enddo
    g%nl = m
  end subroutine hg_listsub

end module hgroup_mod
