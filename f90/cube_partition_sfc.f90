!-------------------------------------------------------------------------------
   module cube_partition_sfc
!-------------------------------------------------------------------------------
!
! abstract : partitioning the cubed-sphere using the space-filling curves
!
! history log :
!   2013-09-08  ki-hwan kim  initial work with python
!   2017-04-06  ki-hwan kim  convert to f90
!   2017-05-17  ki-hwan kim  bug fix the intervals in make_cube_rank()
!   2017-06-16  ki-hwan kim  add make_elem_coord()
!
!-------------------------------------------------------------------------------
!
   implicit none
!
   private
   integer, dimension(4,4)            :: direction_vector_table
   integer, dimension(2,2,4), target  :: hilbert
   integer, dimension(3,3,4), target  :: peano
   integer, dimension(5,5,4), target  :: cinco
!
   type sfc_t
     integer                            :: p
     integer, dimension(:,:,:), pointer :: sfc
   end type sfc_t
!
!
! subroutines
!
   public :: rot, inv_x, inv_y
   public :: make_sfcs
   public :: find_size_factors
   public :: find_factors
   public :: make_panel_sfc
   public :: make_global_sfc
   public :: make_cube_rank
   public :: make_elem_coord
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rot(n, num_rot, arr, ret)
!-------------------------------------------------------------------------------
! rotate a matrix by multiple of 90 degrees in anti-clockwise direction
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                 intent(in   ) :: n, num_rot
   integer, dimension(n,n), intent(in   ) :: arr
   integer, dimension(n,n), intent(  out) :: ret
!
   integer :: i, j
!-------------------------------------------------------------------------------
!
   do j=1,n
     do i=1,n
       if (num_rot .eq. 0) then
         ret(i,j) = arr(i,j)
!
       else if (num_rot .eq. 1) then  ! 90 degrees
         ret(i,j) = arr(j,n-i+1)
!
       else if (num_rot .eq. 2) then  ! 180 degrees
         ret(i,j) = arr(n-i+1,n-j+1)
!
       else if (num_rot .eq. 3) then  ! 270 degrees
         ret(i,j) = arr(n-j+1,i)
!
       else
         stop "the 'num_rot' should be one of 0, 1, 2, 3"
!
       end if
     end do
   end do
!
   end subroutine rot
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inv_x(n, arr, ret)
!-------------------------------------------------------------------------------
! inversion of a matrix in x-direction
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                 intent(in   ) :: n
   integer, dimension(n,n), intent(in   ) :: arr
   integer, dimension(n,n), intent(  out) :: ret
!
   integer :: i, j
!-------------------------------------------------------------------------------
!
   do j=1,n
     do i=1,n
       ret(n-i+1,j) = arr(i,j)
     end do
   end do
!
   end subroutine inv_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inv_y(n, arr, ret)
!-------------------------------------------------------------------------------
! inversion of a matrix in y-direction
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                 intent(in   ) :: n
   integer, dimension(n,n), intent(in   ) :: arr
   integer, dimension(n,n), intent(  out) :: ret
!
   integer :: i, j
!-------------------------------------------------------------------------------
!
   do j=1,n
     do i=1,n
       ret(i,n-j+1) = arr(i,j)
     end do
   end do
!
   end subroutine inv_y
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_sfcs(hilbert, peano, cinco)
!-------------------------------------------------------------------------------
! make the space-filling curves along the direction vectors
!
! edge coordinates
!   4 <-- 3
!   ^     |
!   |     V
!   1 --> 2
!
! direction vector
!   1: 1->2
!   2: 1->4
!   3: 3->2
!   4: 3->4
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, dimension(2,2,4), intent(inout) :: hilbert
   integer, dimension(3,3,4), intent(inout) :: peano
   integer, dimension(5,5,4), intent(inout) :: cinco
!
   integer, dimension(2,2) :: tmp2
   integer, dimension(3,3) :: tmp3
   integer, dimension(5,5) :: tmp5
!-------------------------------------------------------------------------------
!
   hilbert(:,:,1) = reshape((/ 1, 4, &
                               2, 3/), (/2,2/))
!
   peano(:,:,1)   = reshape((/ 1, 8, 9, &
                               2, 7, 6, &
                               3, 4, 5/), (/3,3/))
!
   cinco(:,:,1)   = reshape((/ 1, 2, 3,24,25, &
                               8, 7, 4,23,22, &
                               9, 6, 5,20,21, &
                              10,13,14,19,18, &
                              11,12,15,16,17/), (/5,5/))
!
!   print *, hilbert(1,1,1), hilbert(2,1,1), hilbert(1,2,1), hilbert(2,2,1)
!
   call rot(2, 1, hilbert(:,:,1), tmp2)  ! rotate by 270 degrees
   call rot(3, 1,   peano(:,:,1), tmp3)
   call rot(5, 1,   cinco(:,:,1), tmp5)
   call inv_x(2, tmp2, hilbert(:,:,2))
   call inv_x(3, tmp3,   peano(:,:,2))
   call inv_x(5, tmp5,   cinco(:,:,2))
!
   call rot(2, 1, hilbert(:,:,1), tmp2)  ! rotate by 90 degrees
   call rot(3, 1,   peano(:,:,1), tmp3)
   call rot(5, 1,   cinco(:,:,1), tmp5)
   call inv_y(2, tmp2, hilbert(:,:,3))
   call inv_y(3, tmp3,   peano(:,:,3))
   call inv_y(5, tmp5,   cinco(:,:,3))
!
   call rot(2, 2, hilbert(:,:,1), hilbert(:,:,4))  ! rotate by 180 degrees
   call rot(3, 2,   peano(:,:,1),   peano(:,:,4))
   call rot(5, 2,   cinco(:,:,1),   cinco(:,:,4))
!
   end subroutine make_sfcs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function find_size_factors(ne) result(size_factors)
!-------------------------------------------------------------------------------
! initialize the factors of ne
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, intent(in   ) :: ne
!
   integer :: n, p, i
   integer :: size_factors
   integer, dimension(3) :: prime_numbers = (/2,3,5/)
!-------------------------------------------------------------------------------
!
   n = ne
   size_factors = 0
   do i=1,3
     p = prime_numbers(i)
     do while (mod(n,p) .eq. 0)
       n = n/p
       size_factors = size_factors + 1
     end do
   end do
   if (n .ne. 1) stop "The 'ne' should be a composite number of 2, 3, 5."
!
   end function find_size_factors
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine find_factors(ne, nproc, size_factors, factors)
!-------------------------------------------------------------------------------
! initialize the factors of ne
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                          intent(in   ) :: ne, nproc, size_factors
   integer, dimension(size_factors), intent(  out) :: factors
!
   integer :: i, k
   integer :: n, p
   integer, dimension(3) :: prime_numbers
!-------------------------------------------------------------------------------
!
! take the case where square division is possible
!
!   if (mod(nproc,6*2*2) .eq. 0) then
!     prime_numbers(:) = (/3,5,2/)
!!
!   else if (mod(nproc,6*3*3) .eq. 0) then
!     prime_numbers(:) = (/2,5,3/)
!!
!   else
!     prime_numbers(:) = (/2,3,5/)
!!
!   end if
   prime_numbers(:) = (/2,3,5/)
!
!
! make the factors
!
   factors(:) = 0
   n = ne
   k = 1
   do i=1,3
     p = prime_numbers(i)
     do while (mod(n,p) .eq. 0)
       factors(k) = p
       n = n/p
       k = k + 1
     end do
   end do
!
   end subroutine find_factors
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_direction_vector_table(dv_table)
!-------------------------------------------------------------------------------
! make the direction vector table
!
! neighbor indices
!       4
!     ----
!   1|    |
!    |    |2
!     ----
!      3
!
! edge coordinates
!   4 <-- 3
!   ^     |
!   |     V
!   1 --> 2
!
! direction vector
!   1: 1->2
!   2: 1->4
!   3: 3->2
!   4: 3->4
!
! relation diagram with the previous element
!
!            4 ---- 3
!             |    |
!             |    |
!            1 ---- 2
!
!  4 ---- 3  4 ---- 3  4 ---- 3
!   |    |    |prev|    |    |
!   |    |    |    |    |    |
!  1 ---- 2  1 ---- 2  1 ---- 2
!
!            4 ---- 3
!             |    |
!             |    |
!            1 ---- 2
!
!
! dv_table(prev nbr index, next nbr index) = direction vector index
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, dimension(4,4), intent(  out) :: dv_table
!-------------------------------------------------------------------------------
!
   dv_table(1,1) = -1
   dv_table(1,2) = 1
   dv_table(1,3) = 1
   dv_table(1,4) = 2
!
   dv_table(2,1) = 4
   dv_table(2,2) = -1
   dv_table(2,3) = 3
   dv_table(2,4) = 4
!
   dv_table(3,1) = 2
   dv_table(3,2) = 1
   dv_table(3,3) = -1
   dv_table(3,4) = 2
!
   dv_table(4,1) = 4
   dv_table(4,2) = 3
   dv_table(4,3) = 3
   dv_table(4,4) = -1
!
   end subroutine make_direction_vector_table
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine find_direction_vector(nlev, lev, uplev_seq, prev_nbrs, next_nbrs, curves, dvs)
!-------------------------------------------------------------------------------
! inversion of a matrix in y-direction
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                      intent(in   ) :: nlev, lev, uplev_seq
   type(sfc_t), dimension(nlev), intent(in   ) :: curves
   integer,     dimension(nlev), intent(inout) :: prev_nbrs, next_nbrs
   integer,     dimension(nlev), intent(inout) :: dvs
!
   integer :: next_nbr=-1, prev_nbr
   integer :: uplev_size
   integer, dimension(4) :: next2prev
   integer, dimension(2) :: ij1, ij2, dij
!-------------------------------------------------------------------------------
!
   uplev_size = curves(lev+1)%p * curves(lev+1)%p
!
   if (uplev_seq .eq. uplev_size) then
     next_nbr = next_nbrs(lev+1)
!
   else
     ij1(:) = minloc(abs(curves(lev+1)%sfc(:,:,dvs(lev+1))-uplev_seq))
     ij2(:) = minloc(abs(curves(lev+1)%sfc(:,:,dvs(lev+1))-uplev_seq-1))
     dij(:) = ij2 - ij1
!
     if ((dij(1) .eq. -1) .and. (dij(2) .eq. 0)) then
       next_nbr = 1
     else if ((dij(1) .eq. 1) .and. (dij(2) .eq. 0)) then
       next_nbr = 2
     else if ((dij(1) .eq. 0) .and. (dij(2) .eq. -1)) then
       next_nbr = 3
     else if ((dij(1) .eq. 0) .and. (dij(2) .eq. 1)) then
       next_nbr = 4
     end if
   end if
!
   next2prev(:) = (/2,1,4,3/)
   prev_nbr = next2prev(next_nbrs(lev))
!
   next_nbrs(lev) = next_nbr
   prev_nbrs(lev) = prev_nbr
   dvs(lev) = direction_vector_table(prev_nbr, next_nbr)
!
   end subroutine find_direction_vector
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_panel_sfc(ne, nproc, panel_sfc)
!-------------------------------------------------------------------------------
! make a space-filling curve in a panel of cubed-sphere
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                   intent(in   ) :: ne, nproc
   integer, dimension(ne,ne), intent(  out) :: panel_sfc
!
   integer :: nlev, lev
   integer :: seq, uplev_seq
   integer :: dx, gi, gj
   integer :: i, p
   integer, allocatable :: factors(:)                  ! prime numbers of ne
   integer, allocatable :: dvs(:)                      ! direction vector index
   integer, allocatable :: prev_nbrs(:), next_nbrs(:)  ! neighbor index
   integer, allocatable :: sizes(:)                    ! size of curve matrix
   integer, allocatable :: intervals(:)
   type(sfc_t), allocatable :: curves(:)
!-------------------------------------------------------------------------------
!
! initialize
!
   panel_sfc(:,:) = 0
!
   nlev = find_size_factors(ne)
   allocate(factors(nlev))
   allocate(dvs(nlev))
   allocate(prev_nbrs(nlev))
   allocate(next_nbrs(nlev))
   allocate(sizes(nlev))
   allocate(intervals(nlev))
   allocate(curves(nlev))

   call make_sfcs(hilbert, peano, cinco)
   call make_direction_vector_table(direction_vector_table)
   call find_factors(ne, nproc, nlev, factors)
!
   dvs(:) = 1        ! direction vector index
   prev_nbrs(:) = 1  ! neighbor index
   next_nbrs(:) = 2  ! neighbor index
!
   do i=1,nlev
     p = factors(i)
     sizes(i) = p*p
     curves(i)%p = p
     if (p .eq. 2) then
       curves(i)%sfc => hilbert
     else if (p .eq. 3) then
       curves(i)%sfc => peano
     else if (p .eq. 5) then
       curves(i)%sfc => cinco
     end if
   end do
!
   intervals(1) = 1
   do i=2,nlev
     intervals(i) = intervals(i-1) * sizes(i)
   end do
!
   dx = factors(1)
   gi = 1-dx
   gj = 1
!
!
! make overlapped curves
!
   do seq=0,ne*ne/sizes(1)-1
     do lev=nlev-1,1,-1
       if (mod(seq,intervals(lev)) .eq. 0) then
         uplev_seq = mod(seq/intervals(lev), sizes(lev+1)) + 1
         call find_direction_vector(nlev, lev, uplev_seq, prev_nbrs, next_nbrs, curves, dvs)
       end if
     end do
!
     if (prev_nbrs(1) .eq. 1) then
       gi = gi + dx
     else if (prev_nbrs(1) .eq. 2) then
       gi = gi - dx
     else if (prev_nbrs(1) .eq. 3) then
       gj = gj + dx
     else if (prev_nbrs(1) .eq. 4) then
       gj = gj - dx
     end if
!
     panel_sfc(gi:gi+dx-1,gj:gj+dx-1) = curves(1)%sfc(:,:,dvs(1)) + seq*sizes(1)
   end do
!
!
! free allocated
!
   deallocate(factors)
   deallocate(dvs)
   deallocate(prev_nbrs)
   deallocate(next_nbrs)
   deallocate(sizes)
   deallocate(intervals)
   deallocate(curves)
!
   end subroutine make_panel_sfc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_global_sfc(ne, nproc, cube_gid)
!-------------------------------------------------------------------------------
! numbering the elements on the cubed-sphere along the space-filling curves
!
! direction vectors on the cube (HOMME style)
!           <---
!          | 6  |
!          |(3) |
!           ---- 
! 
!   ----    --->    --->    ---- 
!  | 4  |  | 1  |  | 2  |  | 3  |
!  V(4) |  |(1) |  |(2) |  |(6) |
!   ----    ----    ----    ---> 
! 
!           ----
!          | 5  |
!          |(5) |
!           --->
!-------------------------------------------------------------------------------
!
   implicit none
   integer,                     intent(in   ) :: ne, nproc
   integer, dimension(ne,ne,6), intent(  out) :: cube_gid
!
   integer, dimension(ne,ne) :: panel_sfc
   integer, dimension(ne,ne) :: tmp
!-------------------------------------------------------------------------------
!
   call make_panel_sfc(ne, nproc, panel_sfc)
!
   call inv_y(ne, panel_sfc, tmp) 
   cube_gid(:,:,1) = tmp
   cube_gid(:,:,2) = tmp + ne*ne
!
   cube_gid(:,:,3) = panel_sfc + 5*ne*ne
!   
   call rot(ne, 3, panel_sfc, tmp)  ! rotate by 260 degrees
   cube_gid(:,:,4) = tmp + 3*ne*ne 
!
   cube_gid(:,:,5) = panel_sfc + 4*ne*ne
!
   call rot(ne, 2, panel_sfc, tmp)  ! rotate by 180 degrees
   cube_gid(:,:,6) = tmp + 2*ne*ne
!
   end subroutine make_global_sfc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_cube_rank(ne, nproc, nelems, cube_rank, cube_lid)
!-------------------------------------------------------------------------------
! assign the process numbers to elements on the cubed-sphere
! along the space-filling curves
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                     intent(in   ) :: ne, nproc
   integer, dimension(nproc),   intent(  out) :: nelems
   integer, dimension(ne,ne,6), intent(  out) :: cube_rank
   integer, dimension(ne,ne,6), intent(  out) :: cube_lid
!
   integer :: i, ei, ej, p
   integer :: remain_nelem, gid, proc
   integer, dimension(nproc)   :: lids
   integer, dimension(nproc+1) :: accum_nelems
   integer, dimension(ne,ne,6) :: global_elem_id
!-------------------------------------------------------------------------------
!
   call make_global_sfc(ne, nproc, global_elem_id)
!
   remain_nelem = mod(ne*ne*6, nproc)
   do i=1,nproc
     nelems(i) = ne*ne*6/nproc
     if (i .le. remain_nelem) nelems(i) = nelems(i) + 1
   end do
!
   accum_nelems(1) = 0
   do i=2,nproc+1
     accum_nelems(i) = accum_nelems(i-1) + nelems(i-1)
   end do
!
   lids(:) = 1
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         gid = global_elem_id(ei,ej,p)
!
         do proc=1,nproc
           if ((gid .gt. accum_nelems(proc)) .and. &
               (gid .le. accum_nelems(proc+1))) then
             cube_rank(ei,ej,p) = proc-1
             cube_lid(ei,ej,p) = lids(proc)
             lids(proc) = lids(proc) + 1
             exit
           end if
         end do
!
       end do
     end do
   end do
!
   end subroutine make_cube_rank
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_elem_coord(ne, iproc, nelem, cube_rank, cube_lid, elem_coord)
!-------------------------------------------------------------------------------
! coordinates of the elements in my rank
! elem_coord(3,nelem)->(ei,ej,panel)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,                     intent(in   ) :: ne, iproc, nelem
   integer, dimension(ne,ne,6), intent(in   ) :: cube_rank
   integer, dimension(ne,ne,6), intent(in   ) :: cube_lid
   integer, dimension(3,nelem), intent(  out) :: elem_coord
!
   integer :: ei, ej, p
   integer :: rank, lid
!-------------------------------------------------------------------------------
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         rank = cube_rank(ei,ej,p)
         lid  = cube_lid(ei,ej,p)
!
         if (iproc .eq. rank) then
           elem_coord(:,lid) = (/ei,ej,p/)
         end if
!
       end do
     end do
   end do
!
   end subroutine make_elem_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
end module cube_partition_sfc
!-------------------------------------------------------------------------------
