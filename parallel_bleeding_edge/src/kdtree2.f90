!
!(c) matthew kennel, institute for nonlinear science (2004)
!
! licensed under the academic free license version 1.1 found in file license
! with additional provisions found in that same file.
!
module kdtree2_precision_module
  
  integer, parameter :: sp = kind(0.0)
  integer, parameter :: dp = kind(0.0d0)

  private :: sp, dp

  !
  ! you must comment out exactly one
  ! of the two lines.  if you comment
  ! out kdkind = sp then you get single precision
  ! and if you comment out kdkind = dp 
  ! you get double precision.
  !

  !integer, parameter :: kdkind = sp  
  integer, parameter :: kdkind = dp  
  public :: kdkind

end module kdtree2_precision_module

module kdtree2_priority_queue_module
  use kdtree2_precision_module
  !
  ! maintain a priority queue (pq) of data, pairs of 'priority/payload', 
  ! implemented with a binary heap.  this is the type, and the 'dis' field
  ! is the priority.
  !
  type kdtree2_result
      ! a pair of distances, indexes
      real(kdkind)    :: dis!=0.0
      integer :: idx!=-1   initializers cause some bugs in compilers.
  end type kdtree2_result
  !
  ! a heap-based priority queue lets one efficiently implement the following
  ! operations, each in log(n) time, as opposed to linear time.
  !
  ! 1)  add a datum (push a datum onto the queue, increasing its length) 
  ! 2)  return the priority value of the maximum priority element 
  ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
  !     the size of the queue. 
  ! 4)  replace the datum with the maximum priority with a supplied datum
  !     (of either higher or lower priority), maintaining the size of the
  !     queue. 
  !
  !
  ! in the k-d tree case, the 'priority' is the square distance of a point in
  ! the data set to a reference point.   the goal is to keep the smallest m
  ! distances to a reference point.  the tree algorithm searches terminal
  ! nodes to decide whether to add points under consideration.
  !
  ! a priority queue is useful here because it lets one quickly return the
  ! largest distance currently existing in the list.  if a new candidate
  ! distance is smaller than this, then the new candidate ought to replace
  ! the old candidate.  in priority queue terms, this means removing the
  ! highest priority element, and inserting the new one.
  !
  ! algorithms based on cormen, leiserson, rivest, _introduction
  ! to algorithms_, 1990, with further optimization by the author.
  !
  ! originally informed by a c implementation by sriranga veeraraghavan.
  !
  ! this module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type pq
      !
      ! the priority queue consists of elements
      ! priority(1:heap_size), with associated payload(:).
      !
      ! there are heap_size active elements. 
      ! assumes the allocation is always sufficient.  will not increase it
      ! to match.
      integer :: heap_size = 0
      type(kdtree2_result), pointer :: elems(:) 
  end type pq

  public :: kdtree2_result

  public :: pq
  public :: pq_create
  public :: pq_delete, pq_insert
  public :: pq_extract_max, pq_max, pq_replace_max, pq_maxpri
  private

contains


  function pq_create(results_in) result(res)
    !
    ! create a priority queue from already allocated
    ! array pointers for storage.  note! it will not
    ! add any alements to the heap, i.e. any existing
    ! data in the input arrays will not be used and may
    ! be overwritten.
    ! 
    ! usage:
    !    real(kdkind), pointer :: x(:)
    !    integer, pointer :: k(:)
    !    allocate(x(1000),k(1000))
    !    pq => pq_create(x,k)
    !
    type(kdtree2_result), target:: results_in(:) 
    type(pq) :: res
    !
    !
    integer :: nalloc

    nalloc = size(results_in,1)
    if (nalloc .lt. 1) then
       write (*,*) 'pq_create: error, input arrays must be allocated.'
    end if
    res%elems => results_in
    res%heap_size = 0
    return
  end function pq_create

  !
  ! operations for getting parents and left + right children
  ! of elements in a binary heap.
  !

!
! these are written inline for speed.
!    
!  integer function parent(i)
!    integer, intent(in) :: i
!    parent = (i/2)
!    return
!  end function parent

!  integer function left(i)
!    integer, intent(in) ::i
!    left = (2*i)
!    return
!  end function left

!  integer function right(i)
!    integer, intent(in) :: i
!    right = (2*i)+1
!    return
!  end function right

!  logical function compare_priority(p1,p2)
!    real(kdkind), intent(in) :: p1, p2
!
!    compare_priority = (p1 .gt. p2)
!    return
!  end function compare_priority

  subroutine heapify(a,i_in)
    !
    ! take a heap rooted at 'i' and force it to be in the
    ! heap canonical form.   this is performance critical 
    ! and has been tweaked a little to reflect this.
    !
    type(pq),pointer   :: a
    integer, intent(in) :: i_in
    !
    integer :: i, l, r, largest

    real(kdkind)    :: pri_i, pri_l, pri_r, pri_largest


    type(kdtree2_result) :: temp

    i = i_in

bigloop:  do
       l = 2*i ! left(i)
       r = l+1 ! right(i)
       ! 
       ! set 'largest' to the index of either i, l, r
       ! depending on whose priority is largest.
       !
       ! note that l or r can be larger than the heap size
       ! in which case they do not count.


       ! does left child have higher priority? 
       if (l .gt. a%heap_size) then
          ! we know that i is the largest as both l and r are invalid.
          exit 
       else
          pri_i = a%elems(i)%dis
          pri_l = a%elems(l)%dis 
          if (pri_l .gt. pri_i) then
             largest = l
             pri_largest = pri_l
          else
             largest = i
             pri_largest = pri_i
          endif

          !
          ! between i and l we have a winner
          ! now choose between that and r.
          !
          if (r .le. a%heap_size) then
             pri_r = a%elems(r)%dis
             if (pri_r .gt. pri_largest) then
                largest = r
             endif
          endif
       endif

       if (largest .ne. i) then
          ! swap data in nodes largest and i, then heapify

          temp = a%elems(i)
          a%elems(i) = a%elems(largest)
          a%elems(largest) = temp 
          ! 
          ! canonical heapify() algorithm has tail-ecursive call: 
          !
          !        call heapify(a,largest)   
          ! we will simulate with cycle
          !
          i = largest
          cycle bigloop ! continue the loop 
       else
          return   ! break from the loop
       end if
    enddo bigloop
    return
  end subroutine heapify

  subroutine pq_max(a,e) 
    !
    ! return the priority and its payload of the maximum priority element
    ! on the queue, which should be the first one, if it is 
    ! in heapified form.
    !
    type(pq),pointer :: a
    type(kdtree2_result),intent(out)  :: e

    if (a%heap_size .gt. 0) then
       e = a%elems(1) 
    else
       write (*,*) 'pq_max: error, heap_size < 1'
       stop
    endif
    return
  end subroutine pq_max
  
  real(kdkind) function pq_maxpri(a)
    type(pq), pointer :: a

    if (a%heap_size .gt. 0) then
       pq_maxpri = a%elems(1)%dis
    else
       write (*,*) 'pq_max_pri: error, heapsize < 1'
       stop
    endif
    return
  end function pq_maxpri

  subroutine pq_extract_max(a,e)
    !
    ! return the priority and payload of maximum priority
    ! element, and remove it from the queue.
    ! (equivalent to 'pop()' on a stack)
    !
    type(pq),pointer :: a
    type(kdtree2_result), intent(out) :: e
    
    if (a%heap_size .ge. 1) then
       !
       ! return max as first element
       !
       e = a%elems(1) 
       
       !
       ! move last element to first
       !
       a%elems(1) = a%elems(a%heap_size) 
       a%heap_size = a%heap_size-1
       call heapify(a,1)
       return
    else
       write (*,*) 'pq_extract_max: error, attempted to pop non-positive pq'
       stop
    end if
    
  end subroutine pq_extract_max


  real(kdkind) function pq_insert(a,dis,idx) 
    !
    ! insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    type(pq),pointer  :: a
    real(kdkind), intent(in) :: dis
    integer, intent(in) :: idx
    !    type(kdtree2_result), intent(in) :: e
    !
    integer :: i, isparent
    real(kdkind)    :: parentdis
    !

    !    if (a%heap_size .ge. a%max_elems) then
    !       write (*,*) 'pq_insert: error, attempt made to insert element on full pq'
    !       stop
    !    else
    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while (i .gt. 1)
       isparent = int(i/2)
       parentdis = a%elems(isparent)%dis
       if (dis .gt. parentdis) then
          ! move what was in i's parent into i.
          a%elems(i)%dis = parentdis
          a%elems(i)%idx = a%elems(isparent)%idx
          i = isparent
       else
          exit
       endif
    end do

    ! insert the element at the determined position
    a%elems(i)%dis = dis
    a%elems(i)%idx = idx

    pq_insert = a%elems(1)%dis 
    return
    !    end if

  end function pq_insert

  subroutine pq_adjust_heap(a,i)
    type(pq),pointer  :: a
    integer, intent(in) :: i
    !
    ! nominally arguments (a,i), but specialize for a=1
    !
    ! this routine assumes that the trees with roots 2 and 3 are already heaps, i.e.
    ! the children of '1' are heaps.  when the procedure is completed, the
    ! tree rooted at 1 is a heap.
    real(kdkind) :: prichild
    integer :: parent, child, n

    type(kdtree2_result) :: e

    e = a%elems(i) 

    parent = i
    child = 2*i
    n = a%heap_size
    
    do while (child .le. n)
       if (child .lt. n) then
          if (a%elems(child)%dis .lt. a%elems(child+1)%dis) then
             child = child+1
          endif
       endif
       prichild = a%elems(child)%dis
       if (e%dis .ge. prichild) then
          exit 
       else
          ! move child into parent.
          a%elems(parent) = a%elems(child) 
          parent = child
          child = 2*parent
       end if
    end do
    a%elems(parent) = e
    return
  end subroutine pq_adjust_heap
    

  real(kdkind) function pq_replace_max(a,dis,idx) 
    !
    ! replace the extant maximum priority element
    ! in the pq with (dis,idx).  return
    ! the new maximum priority, which may be larger
    ! or smaller than the old one.
    !
    type(pq),pointer         :: a
    real(kdkind), intent(in) :: dis
    integer, intent(in) :: idx
!    type(kdtree2_result), intent(in) :: e
    ! not tested as well!

    integer :: parent, child, n
    real(kdkind)    :: prichild, prichildp1

    type(kdtree2_result) :: etmp
    
    if (.true.) then
       n=a%heap_size
       if (n .ge. 1) then
          parent =1
          child=2

          loop: do while (child .le. n)
             prichild = a%elems(child)%dis

             !
             ! posibly child+1 has higher priority, and if
             ! so, get it, and increment child.
             !

             if (child .lt. n) then
                prichildp1 = a%elems(child+1)%dis
                if (prichild .lt. prichildp1) then
                   child = child+1
                   prichild = prichildp1
                endif
             endif

             if (dis .ge. prichild) then
                exit loop  
                ! we have a proper place for our new element, 
                ! bigger than either children's priority.
             else
                ! move child into parent.
                a%elems(parent) = a%elems(child) 
                parent = child
                child = 2*parent
             end if
          end do loop
          a%elems(parent)%dis = dis
          a%elems(parent)%idx = idx
          pq_replace_max = a%elems(1)%dis
       else
          a%elems(1)%dis = dis
          a%elems(1)%idx = idx
          pq_replace_max = dis
       endif
    else
       !
       ! slower version using elementary pop and push operations.
       !
       call pq_extract_max(a,etmp) 
       etmp%dis = dis
       etmp%idx = idx
       pq_replace_max = pq_insert(a,dis,idx)
    endif
    return
  end function pq_replace_max

  subroutine pq_delete(a,i)
    ! 
    ! delete item with index 'i'
    !
    type(pq),pointer :: a
    integer           :: i

    if ((i .lt. 1) .or. (i .gt. a%heap_size)) then
       write (*,*) 'pq_delete: error, attempt to remove out of bounds element.'
       stop
    endif

    ! swap the item to be deleted with the last element
    ! and shorten heap by one.
    a%elems(i) = a%elems(a%heap_size) 
    a%heap_size = a%heap_size - 1

    call heapify(a,i)

  end subroutine pq_delete

end module kdtree2_priority_queue_module


module kdtree2_module
  use kdtree2_precision_module
  use kdtree2_priority_queue_module
  use treebuild
!$ use omp_lib
  ! k-d tree routines in fortran 90 by matt kennel.
  ! original program was written in sather by steve omohundro and
  ! matt kennel.  only the euclidean metric is supported. 
  !
  !
  ! this module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! in otherwords for an embedding of n d-dimensional vectors, the
  ! data file is here, in natural fortran order  data(1:d, 1:n)
  ! because fortran lays out columns first,
  !
  ! whereas conventionally (c-style) it is data(1:n,1:d)
  ! as in the original kd_tree module. 
  !
  !-------------data type, creation, deletion---------------------
  public :: kdkind
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy
  !---------------------------------------------------------------
  !-------------------search routines-----------------------------
  public :: kdtree2_n_nearest,kdtree2_n_nearest_around_point
  ! return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window. 
  !
  public :: kdtree2_r_nearest, kdtree2_r_nearest_around_point
  ! return points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_sort_results
  ! sort, in order of increasing distance, rseults from above.
  !
  public :: kdtree2_r_count, kdtree2_r_count_around_point 
  ! count points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force
  ! brute force of kdtree2_[n|r]_nearest
  !----------------------------------------------------------------


  integer, parameter :: bucket_size = 12
  ! the maximum number of points to keep in a terminal node.

  ! The tree is made by the parallel builder of treebuild.f90 (the mean of the coordinate as the cut, as it always was);
  ! the nodes it numbers are then copied into the pointer nodes below, in a pool kept for the next tree.
  type(tbwork), save, target :: tbw

  type interval
      real(kdkind) :: lower,upper
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(kdkind) :: cut_val
      ! where to cut the dimension
      real(kdkind) :: cut_val_left, cut_val_right  
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (tree_node), pointer :: left, right
      type(interval) :: box(3)
      ! child pointers
      ! points included in this node are indexes[k] with k \in [l,u] 


  end type tree_node

  type :: kdtree2
      ! global information about the tree, one per tree
      integer :: dimen=0, n=0
      ! dimensionality and total # of points
      real(kdkind), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array 
      ! 
      !  important note:  it is dimensioned   the_data(1:d,1:n)
      !  which may be opposite of what may be conventional.
      !  this is, because in fortran, the memory layout is such that
      !  the first dimension is in sequential order.  hence, with
      !  (1:d,1:n), all components of the vector will be in consecutive
      !  memory locations.  the search time is dominated by the
      !  evaluation of distances in the terminal nodes.  putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable 
      !  vectorization on some processors and compilers. 

      integer, pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical       :: sort = .false.
      ! do we always sort output results?
      logical       :: rearrange = .false. 
      real(kdkind), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer :: root => null()
      ! root pointer of the tree
      type (tree_node), pointer :: pool(:) => null()
      ! the nodes, root = pool(1)
  end type kdtree2

  ! The arrays of the last tree that was destroyed are kept for the next one of the same size (a tree is made at every step,
  ! and the memory of a big one would be new, and slow to touch, every time).
  integer, pointer, save :: spare_ind(:) => null()
  real(kdkind), pointer, save :: spare_pos(:,:) => null()
  type (tree_node), pointer, save :: spare_pool(:) => null()

  type :: tree_search_record
      !
      ! one of these is created for each search.
      !
      private
      ! 
      ! many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer           :: dimen   
      integer           :: nn, nfound
      real(kdkind)      :: ballsize
      integer           :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer           :: nalloc  ! how much allocated for results(:)?
      logical           :: rearrange  ! are the data rearranged or original? 
      ! did the # of points found overflow the storage provided?
      logical           :: overflow
      real(kdkind), pointer :: qv(:)  ! query vector
      type(kdtree2_result), pointer :: results(:) ! results
      type(pq) :: pq
      real(kdkind), pointer :: data(:,:)  ! temp pointer to data
      integer, pointer      :: ind(:)     ! temp pointer to indexes
  end type tree_search_record

  private
  ! everything else is private.

  type(tree_search_record), save, target :: sr   ! a global variable for search
  ! every thread searches with its own record: the searches only read the tree
!$omp threadprivate(sr)

contains

  function kdtree2_create(input_data,dim,sort,rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! note, input data is input_data(1:d,1:n), not the other way around.
    ! this is the reverse of the previous version of this module.
    ! the reason for it is cache friendliness, improving performance.
    !
    ! optional arguments:  if 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from size(input_data,1).  (Only three dimensions
    !                      are built now.)
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !                      
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness. 
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. function return cut_value ..
    type (kdtree2), pointer :: mr
    integer, intent(in), optional      :: dim
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    ! ..
    ! .. array arguments ..
    real(kdkind), target :: input_data(:,:)
    !
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       write (*,*) 'kd_tree_trans: likely user error.'
       write (*,*) 'kd_tree_trans: you passed in matrix with d=',mr%dimen
       write (*,*) 'kd_tree_trans: and n=',mr%n
       write (*,*) 'kd_tree_trans: note, that new format is data(1:d,1:n)'
       write (*,*) 'kd_tree_trans: with usually n >> d.   if n =approx= d, then a k-d tree'
       write (*,*) 'kd_tree_trans: is not an appropriate data structure.'
       stop
    end if
    if (mr%dimen /= 3 .or. size(input_data,1) < 3) then
       write (*,*) 'kdtree2_create: only trees in three dimensions are built (dimen=', mr%dimen, ')'
       stop
    end if

    call build_tree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (.not. mr%rearrange) then
       ! the positions in tree order were only the way to the permutation
       if (associated(mr%rearranged_data)) deallocate (mr%rearranged_data)
       nullify (mr%rearranged_data)
    endif

  end function kdtree2_create

    subroutine build_tree(tp)
      ! the permutation and the positions in tree order come from the builder; the nodes are copied from its arrays
      type (kdtree2), pointer :: tp
      ! ..
      integer :: n, j, id, ncap, k
      n = tp%n
      ! memory of the tree destroyed last, if it fits
!$omp critical (kdtree2_spare)
      if (associated(spare_ind)) then
         if (size(spare_ind) == n) then
            tp%ind => spare_ind
         else
            deallocate (spare_ind)
         end if
         nullify (spare_ind)
      end if
      if (associated(spare_pos)) then
         if (size(spare_pos,2) == n) then
            tp%rearranged_data => spare_pos
         else
            deallocate (spare_pos)
         end if
         nullify (spare_pos)
      end if
!$omp end critical (kdtree2_spare)
      if (.not. associated(tp%ind)) allocate (tp%ind(n))
      if (.not. associated(tp%rearranged_data)) allocate (tp%rearranged_data(3,n))
      call tb_build(tbw, n, tp%the_data, tp%ind, tp%rearranged_data, bucket_size + 1, 1)
      ncap = tbw%nidmax
!$omp critical (kdtree2_spare)
      if (associated(spare_pool)) then
         if (size(spare_pool) >= ncap) then
            tp%pool => spare_pool
         else
            deallocate (spare_pool)
         end if
         nullify (spare_pool)
      end if
!$omp end critical (kdtree2_spare)
      if (.not. associated(tp%pool)) allocate (tp%pool(max(ncap,tbw%ncap)))
      ! the nodes of the top of the tree, then those of the subtrees
!$omp parallel default(shared) private(id,j,k) if(n.gt.8*2048 .and. .not.omp_in_parallel())
!$omp do schedule(static)
      do id = 1, tbw%ntop
         call copy_node(tp, id)
      end do
!$omp end do nowait
!$omp do schedule(dynamic,4)
      do j = 1, tbw%nsub
         do k = tbw%sub_base(j) + 1, tbw%sub_base(j) + tbw%sub_used(j)
            call copy_node(tp, k)
         end do
      end do
!$omp end do
!$omp end parallel
      tp%root => tp%pool(1)
    end subroutine build_tree

    subroutine copy_node(tp, id)
      type (kdtree2), pointer :: tp
      integer, intent(in) :: id
      type (tree_node), pointer :: nd
      integer :: c, i, lc, rc
      nd => tp%pool(id)
      nd%l = tbw%nl(id)
      nd%u = tbw%nu(id)
      do i = 1, 3
         nd%box(i)%lower = tbw%lo(i,id)
         nd%box(i)%upper = tbw%hi(i,id)
      end do
      lc = tbw%lc(id)
      if (lc == 0) then
         ! a terminal node
         nd%cut_dim = 0
         nd%cut_val = 0.0
         nd%cut_val_left = 0.0
         nd%cut_val_right = 0.0
         nullify (nd%left)
         nullify (nd%right)
      else
         rc = tbw%rc(id)
         c = tbw%cdim(id)
         nd%cut_dim = c
         nd%left => tp%pool(lc)
         nd%right => tp%pool(rc)
         nd%cut_val_left = tbw%hi(c,lc)
         nd%cut_val_right = tbw%lo(c,rc)
         nd%cut_val = (nd%cut_val_left + nd%cut_val_right)/2
      end if
    end subroutine copy_node

  subroutine kdtree2_destroy(tp)
    ! deallocates all memory for the tree, except input data matrix
    ! (the arrays of the tree go to the spare ones if there are none, to serve the next tree)
    ! .. structure arguments ..
    type (kdtree2), pointer :: tp
    ! ..
!$omp critical (kdtree2_spare)
    if (.not. associated(spare_ind)) then
       spare_ind => tp%ind
    else
       deallocate (tp%ind)
    end if
    nullify (tp%ind)
    if (associated(tp%rearranged_data)) then
       if (.not. associated(spare_pos)) then
          spare_pos => tp%rearranged_data
       else
          deallocate (tp%rearranged_data)
       end if
       nullify (tp%rearranged_data)
    end if
    if (associated(tp%pool)) then
       if (.not. associated(spare_pool)) then
          spare_pool => tp%pool
       else
          deallocate (tp%pool)
       end if
       nullify (tp%pool)
    end if
!$omp end critical (kdtree2_spare)
    nullify (tp%root)
    deallocate(tp)
    return
  end subroutine kdtree2_destroy

  subroutine kdtree2_n_nearest(tp,qv,nn,results)
    ! find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (kdtree2), pointer      :: tp
    real(kdkind), target, intent (in)    :: qv(:)
    integer, intent (in)         :: nn
    type(kdtree2_result), target :: results(:)


    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false. 

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif
    sr%dimen = tp%dimen

    call validate_query_storage(nn) 
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
!    deallocate(sr%pqp)
    return
  end subroutine kdtree2_n_nearest

  subroutine kdtree2_n_nearest_around_point(tp,idxin,correltime,nn,results)
    ! find the 'nn' vectors in the tree nearest to point 'idxin',
    ! with correlation window 'correltime', returing results in
    ! results(:), which must be pre-allocated upon entry.
    type (kdtree2), pointer        :: tp
    integer, intent (in)           :: idxin, correltime, nn
    type(kdtree2_result), target   :: results(:)

    ! the query is a column of the tree's data: point at it instead of
    ! allocating and copying a vector for every search
    sr%qv => tp%the_data(:,idxin)
    sr%ballsize = huge(1.0)       ! the largest real(kdkind) number
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%nn = nn
    sr%nfound = 0

    sr%dimen = tp%dimen
    sr%nalloc = nn

    sr%results => results

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif

    call validate_query_storage(nn)
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
    return
  end subroutine kdtree2_n_nearest_around_point

  subroutine kdtree2_r_nearest(tp,qv,r2,nfound,nalloc,results) 
    ! find the nearest neighbors to point 'idxin', within squared
    ! euclidean distance 'r2'.   upon entry, nalloc must be the
    ! size of memory allocated for results(1:nalloc).  upon
    ! exit, nfound is the number actually found within the ball. 
    !
    !  note that if nfound .gt. nalloc then more neighbors were found
    !  than there were storage to store.  the resulting list is not
    !  the smallest ball inside norm r^2 
    !
    ! results are not sorted unless tree was created with sort option.
    type (kdtree2), pointer      :: tp
    real(kdkind), target, intent (in)    :: qv(:)
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    integer, intent (in)         :: nalloc
    type(kdtree2_result), target :: results(:)

    !
    sr%qv => qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    sr%results => results

    call validate_query_storage(nalloc)
    sr%nalloc = nalloc
    sr%overflow = .false. 
    sr%ind => tp%ind
    sr%rearrange= tp%rearrange

    if (tp%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound, results)
    endif

    if (sr%overflow) then
       write (*,*) 'kd_tree_trans: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'kd_tree_trans: than storage was provided for.  answer is not smallest ball'
       write (*,*) 'kd_tree_trans: with that number of neighbors!  i.e. it is wrong.'
    endif

    return
  end subroutine kdtree2_r_nearest

  subroutine kdtree2_r_nearest_around_point(tp,idxin,correltime,r2,&
   nfound,nalloc,results,overflow)
    !
    ! like kdtree2_r_nearest, but around a point 'idxin' already existing
    ! in the data set. 
    ! 
    ! results are not sorted unless tree was created with sort option.
    !
    type (kdtree2), pointer      :: tp
    integer, intent (in)         :: idxin, correltime, nalloc
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    type(kdtree2_result), target :: results(:)
    ! if present, the search reports in overflow that results(:) was too small (nfound is then the
    ! number of points found) instead of printing a warning, so that the caller can search again
    logical, intent(out), optional :: overflow
    ! ..
    ! .. intrinsic functions ..
    intrinsic huge
    ! ..
    ! the query is a column of the tree's data: point at it instead of
    ! allocating and copying a vector for every search
    sr%qv => tp%the_data(:,idxin)
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%results => results

    sr%nalloc = nalloc
    sr%overflow = .false.

    call validate_query_storage(nalloc)

    !    sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !    sr%il = -1               ! set to invalid indexes

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (tp%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif
    sr%rearrange = tp%rearrange
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound,results)
    endif

    if (present(overflow)) then
       overflow = sr%overflow
    else if (sr%overflow) then
       write (*,*) 'kd_tree_trans: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'kd_tree_trans: than storage was provided for.  answer is not smallest ball'
       write (*,*) 'kd_tree_trans: with that number of neighbors!  i.e. it is wrong.'
    endif

    return
  end subroutine kdtree2_r_nearest_around_point

  function kdtree2_r_count(tp,qv,r2) result(nfound)
    ! count the number of neighbors within square distance 'r2'. 
    type (kdtree2), pointer   :: tp
    real(kdkind), target, intent (in) :: qv(:)
    real(kdkind), intent(in)          :: r2
    integer                   :: nfound
    ! ..
    ! .. intrinsic functions ..
    intrinsic huge
    ! ..
    sr%qv => qv
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    
    nullify(sr%results) ! for some reason, ftn 95 chokes on '=> null()'

    sr%nalloc = 0            ! we do not allocate any storage but that's ok
                             ! for counting.
    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count

  function kdtree2_r_count_around_point(tp,idxin,correltime,r2) &
   result(nfound)
    ! count the number of neighbors within square distance 'r2' around
    ! point 'idxin' with decorrelation time 'correltime'.
    !
    type (kdtree2), pointer :: tp
    integer, intent (in)    :: correltime, idxin
    real(kdkind), intent(in)        :: r2
    integer                 :: nfound
    ! ..
    ! ..
    ! .. intrinsic functions ..
    intrinsic huge
    ! ..
    ! the query is a column of the tree's data: point at it instead of
    ! allocating and copying a vector for every search
    sr%qv => tp%the_data(:,idxin)
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    nullify(sr%results)

    sr%nalloc = 0            ! we do not allocate any storage but that's ok
                             ! for counting.

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%data => tp%rearranged_data
    else
       sr%data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count_around_point


  subroutine validate_query_storage(n)
    !
    ! make sure we have enough storage for n
    !
    integer, intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'kd_tree_trans:  you did not provide enough storage for results(1:n)'
       stop
       return
    endif

    return
  end subroutine validate_query_storage

  function square_distance(d, iv,qv) result (res)
    ! distance between iv[1:n] and qv[1:n] 
    ! .. function return value ..
    ! re-implemented to improve vectorization.
    real(kdkind) :: res
    ! ..
    ! ..
    ! .. scalar arguments ..
    integer :: d
    ! ..
    ! .. array arguments ..
    real(kdkind) :: iv(:),qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:d)-qv(1:d))**2 )
  end function square_distance
  
  recursive subroutine search(node)
    !
    ! this is the innermost core routine of the kd-tree search.  along
    ! with "process_terminal_node", it is the performance bottleneck. 
    !
    ! this version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (tree_node), pointer          :: node
    ! ..
    type(tree_node),pointer            :: ncloser, nfarther
    !
    integer                            :: cut_dim, i
    ! ..
    real(kdkind)                               :: qval, dis
    real(kdkind)                               :: ballsize
    real(kdkind), pointer           :: qv(:)
    type(interval), pointer :: box(:) 

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call process_terminal_node_fixedball(node)
       else
          call process_terminal_node(node)
       endif
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
       endif

       if (associated(ncloser)) call search(ncloser)

       ! we may need to search the second node. 
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false. 
             !
             box => node%box(1:)
             do i=1,sr%dimen
                if (i .ne. cut_dim) then
                   dis = dis + dis2_from_bnd(qv(i),box(i)%lower,box(i)%upper)
                   if (dis > ballsize) then
                      return
                   endif
                endif
             end do
             
             !
             ! if we are still here then we need to search mroe.
             !
             call search(nfarther)
          endif
       endif
    end if
  end subroutine search


  real(kdkind) function dis2_from_bnd(x,amin,amax) result (res)
    real(kdkind), intent(in) :: x, amin,amax

    if (x > amax) then
       res = (x-amax)**2;
       return
    else
       if (x < amin) then
          res = (amin-x)**2;
          return
       else
          res = 0.0
          return
       endif
    endif
    return
  end function dis2_from_bnd

  logical function box_in_search_range(node, sr) result(res)
    !
    ! return the distance from 'qv' to the closest corner of node's
    ! bounding box
    ! for all coordinates outside the box.   coordinates inside the box
    ! contribute nothing to the distance.
    !
    type (tree_node), pointer :: node
    type (tree_search_record), pointer :: sr

    integer :: dimen, i
    real(kdkind)    :: dis, ballsize
    real(kdkind)    :: l, u

    dimen = sr%dimen
    ballsize = sr%ballsize
    dis = 0.0
    res = .true.
    do i=1,dimen
       l = node%box(i)%lower
       u = node%box(i)%upper
       dis = dis + (dis2_from_bnd(sr%qv(i),l,u))
       if (dis > ballsize) then
          res = .false.
          return
       endif
    end do
    res = .true.
    return
  end function box_in_search_range


  subroutine process_terminal_node(node)
    !
    ! look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (tree_node), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: dimen, i, indexofi, k, centeridx, correltime
    real(kdkind)                   :: ballsize, sd, newpri
    logical                :: rearrange
    type(pq), pointer      :: pqp 
    !
    ! copy values from sr to local variables
    !
    !
    ! notice, making local pointers with an explicit lower bound
    ! seems to generate faster code.
    ! why?  i don't know.
    qv => sr%qv(1:) 
    pqp => sr%pq
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%data(1:,1:)     
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! do we have a decorrelation window? 
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       endif


       ! 
       ! two choices for any point.  the list so far is either undersized,
       ! or it is not.
       !
       ! if it is undersized, then add the point and its distance
       ! unconditionally.  if the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! if the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1 
          newpri = pq_insert(pqp,sd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the pq. 

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list. 
          ! hence we replace that with the current one.
          !
          ballsize = pq_replace_max(pqp,sd,indexofi)
       endif
    end do mainloop
    !
    ! reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize 

  end subroutine process_terminal_node

  subroutine process_terminal_node_fixedball(node)
    !
    ! look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (tree_node), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: nfound
    integer                :: dimen, i, indexofi, k
    integer                :: centeridx, correltime, nn
    real(kdkind)                   :: ballsize, sd
    logical                :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       ! 
       ! two choices for any point.  the list so far is either undersized,
       ! or it is not.
       !
       ! if it is undersized, then add the point and its distance
       ! unconditionally.  if the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! if the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       ! which index to the point do we use? 

       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       endif

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%idx = indexofi
       endif
    end do mainloop
    !
    ! reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine process_terminal_node_fixedball

  subroutine kdtree2_n_nearest_brute_force(tp,qv,nn,results) 
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is slow!  the
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(kdkind), intent (in)       :: qv(:)
    integer, intent (in)    :: nn
    type(kdtree2_result)    :: results(:) 

    integer :: i, j, k
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
       results(i)%dis =  huge(1.0)
       results(i)%idx = -1
    end do
    do i = 1, tp%n
       if (all_distances(i)<results(nn)%dis) then
          ! insert it somewhere on the list
          do j = 1, nn
             if (all_distances(i)<results(j)%dis) exit
          end do
          ! now we know 'j'
          do k = nn - 1, j, -1
             results(k+1) = results(k)
          end do
          results(j)%dis = all_distances(i)
          results(j)%idx = i
       end if
    end do
    deallocate (all_distances)
  end subroutine kdtree2_n_nearest_brute_force
  

  subroutine kdtree2_r_nearest_brute_force(tp,qv,r2,nfound,results) 
    ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
    ! only use this subroutine for testing, as it is slow!  the
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(kdkind), intent (in)       :: qv(:)
    real(kdkind), intent (in)       :: r2
    integer, intent(out)    :: nfound
    type(kdtree2_result)    :: results(:) 

    integer :: i, nalloc
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do
    
    nfound = 0
    nalloc = size(results,1)

    do i = 1, tp%n
       if (all_distances(i)< r2) then
          ! insert it somewhere on the list
          if (nfound .lt. nalloc) then
             nfound = nfound+1
             results(nfound)%dis = all_distances(i)
             results(nfound)%idx = i
          endif
       end if
    enddo
    deallocate (all_distances)

    call kdtree2_sort_results(nfound,results)


  end subroutine kdtree2_r_nearest_brute_force

  subroutine kdtree2_sort_results(nfound,results)
    !  use after search to sort results(1:nfound) in order of increasing 
    !  distance.
    integer, intent(in)          :: nfound
    type(kdtree2_result), target :: results(:) 
    !
    !

    !this is buggy with intel fortran
    !    if (nfound .gt. 1) call heapsort(results(1:nfound)%dis,results(1:nfound)%ind,nfound)
    !
    if (nfound .gt. 1) call heapsort_struct(results,nfound)

    return
  end subroutine kdtree2_sort_results

  subroutine heapsort(a,ind,n)
    !
    ! sort a(1:n) in ascending order, permuting ind(1:n) similarly.
    ! 
    ! if ind(k) = k upon input, then it will give a sort index upon output.
    !
    integer,intent(in)          :: n
    real(kdkind), intent(inout)         :: a(:) 
    integer, intent(inout)      :: ind(:)

    !
    !
    real(kdkind)        :: value   ! temporary for a value from a()
    integer     :: ivalue  ! temporary for a value from ind()

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft); ivalue=ind(ileft)
       else
          value=a(iright); ivalue=ind(iright)
          a(iright)=a(1); ind(iright)=ind(1)
          iright=iright-1
          if (iright == 1) then
             a(1)=value;ind(1)=ivalue
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j) < a(j+1)) j=j+1
          endif
          if(value < a(j)) then
             a(i)=a(j); ind(i)=ind(j)
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value; ind(i)=ivalue
    end do
  end subroutine heapsort

  subroutine heapsort_struct(a,n)
    !
    ! sort a(1:n) in ascending order
    ! 
    !
    integer,intent(in)                 :: n
    type(kdtree2_result),intent(inout) :: a(:)

    !
    !
    type(kdtree2_result) :: value ! temporary value

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          endif
          if(value%dis < a(j)%dis) then
             a(i)=a(j); 
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value
    end do
  end subroutine heapsort_struct

end module kdtree2_module

