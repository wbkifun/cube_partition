!-------------------------------------------------------------------------------
   module cube_partition_band
!-------------------------------------------------------------------------------
!
! abstract : partitioning the cubed-sphere using the band search method
!
! history log :
!   2018-03-06  ki-hwan kim  start
!
!-------------------------------------------------------------------------------
!
   use cube_neighbor, only : convert_nbr_eij
!
   implicit none
   logical, parameter :: debug=.false.
!
   private
!
   public :: calc_perimeter_ratio
   public :: find_optimal_band
   public :: band_partition
   public :: make_cube_rank
   public :: make_elem_coord
   public :: global_perimeter_ratio
   public :: global_communication_ratio
   public :: make_cube_color
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function calc_perimeter_ratio(nx, ny, nproc,                                &
       start_rank, end_rank, nelems, i12, box) result(mean_perimeter_ratio) 
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: nx, ny
   integer, intent(in   ) :: nproc
   integer, intent(in   ) :: start_rank  ! start from zero
   integer, intent(in   ) :: end_rank
   integer, intent(in   ) :: nelems(0:nproc-1)
   integer, intent(inout) :: i12(4)  ! (i1, i2, band_elem, remain_elem)
   integer, intent(inout) :: box(nx,ny)
   real(8) :: mean_perimeter_ratio
!
   integer :: i, j, k
   integer :: i1, i2, j1, j2
   integer :: seq
   integer :: rank, myrank
   integer :: band_elem, required_elem, remain_elem
   integer :: same_sides
   integer :: num_nbrs(2,end_rank-start_rank+1)
   real(8) :: sum_perimeter_ratio
!-------------------------------------------------------------------------------
!
! determine the i2: band interval
!
   i1 = i12(1)
   i2 = i12(2)
   band_elem = i12(3)
!
   if (debug) print *, 'i2, start_rank, end_rank',                             &
       i2, start_rank, end_rank
!
   required_elem = sum(nelems(start_rank:end_rank))
   do while (i2.lt.nx .and. required_elem.gt.band_elem)
     i2 = i2 + 1
     band_elem = band_elem + count(box(i2,:) .eq. -1)
     if (i2 .eq. nx) exit
   end do
   if (debug) print *, 'i2, required_elem, band_elem',                         &
       i2, required_elem, band_elem
!
!
! exit condition
!
   if (i2.eq.nx .and. required_elem.gt.band_elem) then
     mean_perimeter_ratio = -1.D0
     if (debug) print *, 'i2 exceeded the nx (exit)'
!
   else if (i2.ne.nx .and. i1.eq.i2) then
     mean_perimeter_ratio = 4.D0
     if (debug) print *, 'i1=i2 single-line band (exit)'
!
   else
     remain_elem = band_elem - required_elem
     if (remain_elem .ne. 0) then
       do j=ny,1,-1
         if (box(i2,j) .eq. -1) exit
       end do
       j2 = j
       j1 = j - remain_elem + 1
       box(i2,j1:j2) = -3  ! temporary mask
     end if
     if (debug) print *, 'remain_elem', remain_elem
!
     !
     ! partitioning in the band
     !
     rank = start_rank
     seq = 1
     do j=1,ny
       do i=i2,i1,-1
         if (box(i,j) .eq. -1) then
           box(i,j) = rank
           if (seq .eq. nelems(rank)) then
             rank = rank + 1
             seq = 1
           else
             seq = seq + 1
           end if
         else if (box(i,j) .eq. -3) then
           box(i,j) = -1
         end if
       end do
     end do 
     if (rank-1 .ne. end_rank) stop 'The rank is greater than the end_rank in calc_perimeter_ratio()'
!
     ! 
     ! compare the perimeter ratio
     ! 
     num_nbrs(:,:) = 0
     do j=1,ny
       do i=i1,i2
         myrank = box(i,j)
         same_sides = 0  ! contact with an element having the same rank number
!
         if (myrank.ge.start_rank .and. myrank.le.end_rank) then
           if (i .gt. 1) then
             if (box(i-1,j) .eq. myrank) same_sides = same_sides + 1
           end if
           if (i .lt. nx) then
             if (box(i+1,j) .eq. myrank) same_sides = same_sides + 1
           end if
           if (j .gt. 1) then
             if (box(i,j-1) .eq. myrank) same_sides = same_sides + 1
           end if
           if (j .lt. ny) then
             if (box(i,j+1) .eq. myrank) same_sides = same_sides + 1
           end if
!
         k = myrank - start_rank + 1
         num_nbrs(1,k) = num_nbrs(1,k) + 1
         num_nbrs(2,k) = num_nbrs(2,k) + 4 - same_sides 
         end if
       end do
     end do
     sum_perimeter_ratio = sum(num_nbrs(2,:)*1.D0/num_nbrs(1,:))
     mean_perimeter_ratio = sum_perimeter_ratio/(end_rank-start_rank+1)
     i12(:) = (/i1, i2, band_elem, remain_elem/)
!
   end if
!
   end function calc_perimeter_ratio
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine find_optimal_band(nx, ny, nproc, start_rank, start_i, nelems,    &
       box, ret)
!-------------------------------------------------------------------------------
! find the optimal band partitioning with minimum perimeter ratio
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: nx, ny
   integer, intent(in   ) :: nproc
   integer, intent(in   ) :: start_rank  ! start from zero
   integer, intent(in   ) :: start_i
   integer, intent(in   ) :: nelems(0:nproc-1)
   integer, intent(inout) :: box(nx,ny)
   integer, intent(  out) :: ret(2)  ! (rank, i2)
!
   integer :: end_rank
   integer :: next_i2
   real(8) :: perimeter_ratio, prev_perimeter_ratio
   integer :: tmp_box(nx,ny), prev_box(nx,ny)
   integer :: i12(4), prev_i12(4)  ! (i1, i2, band_elem, remain_elem)
!-------------------------------------------------------------------------------
!
   i12(:) = (/start_i, start_i, count(box(start_i,:) .eq. -1), -1/)
   end_rank = start_rank !+ int(i12(3)/nelems(start_rank))
   prev_perimeter_ratio = 4.D0  ! max perimeter ratio
   prev_i12(:) = i12(:)
   prev_box(:,:) = box(:,:)
!
   search_loop: do
     if (debug) print *, ''
     tmp_box(:,:) = box(:,:)
     perimeter_ratio = calc_perimeter_ratio(                                   &
         nx, ny, nproc, start_rank, end_rank, nelems, i12, tmp_box) 
     if (debug) print *, 'preimeter_ratio', prev_perimeter_ratio, perimeter_ratio
!
     if (perimeter_ratio .lt. 0.D0 .or.                                        &
         perimeter_ratio .gt. prev_perimeter_ratio) then
       box(:,:) = prev_box(:,:)
       exit search_loop
!
     else
       prev_perimeter_ratio = perimeter_ratio
       end_rank = end_rank + 1
       prev_i12(:) = i12(:)
!
       if (end_rank .eq. nproc) then
         box(:,:) = tmp_box(:,:)
         exit search_loop
       else
         prev_box(:,:) = tmp_box(:,:)
       end if
     end if
!
   end do search_loop
   if (debug) print *, '---------- end search_loop ----------'
!
   next_i2 = prev_i12(2)
   if (prev_i12(4) .eq. 0) next_i2 = next_i2 + 1
   ret(:) = (/end_rank, next_i2/)
!
   if (debug) print *, '(rank,i2)', ret
!
   end subroutine find_optimal_band
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine band_partition(ne, nproc, nelems, cube_rank)
!-------------------------------------------------------------------------------
! assign the rank number to elements on the cubed-sphere
! band partitioning the cubed-sphere
! Caution: The nproc must be greater than 3.
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: ne, nproc
   integer, intent(in   ) :: nelems(nproc)
   integer, intent(  out) :: cube_rank(ne,ne,6)
!
   integer :: i, j, k
   integer :: start_rank, start_i
   integer :: prev_start_rank, prev_start_i
   integer :: box2(2*ne,ne), tmp_box2(2*ne,ne)
   integer :: box4(2*ne,2*ne)
   integer :: ret(2)  ! (rank, i2)
!
   integer :: num_block
   integer :: band1_start_rank, band2_start_rank
   integer :: box(ne,ne), tmp_box(ne,ne)
   integer :: nelem1, nelem2
   integer :: i12(4)  ! (i1, i2, band_elem, remain_elem)
   real(8) :: ratio1, ratio2
   real(8), allocatable :: ratios(:)
   integer, allocatable :: boxes(:,:,:)
!-------------------------------------------------------------------------------
!
   box4(ne+1:2*ne,ne+1:2*ne) = -2  ! permanant mask
!
!
! panel 6
!
   box2(:,:) = -1
   tmp_box2(:,:) = -1
   start_rank = 0
   start_i = 1
   do 
     prev_start_rank = start_rank
     prev_start_i = start_i
     call find_optimal_band(2*ne, ne, nproc, start_rank, start_i, nelems,      &
         box2, ret)
     start_rank = ret(1)
     start_i = ret(2)
!
     if (start_i .gt. ne) then
       if (count(box2(ne+1:2*ne,:).ne.-1) .gt.                                 &
           count(tmp_box2(1:ne,:).eq.-1)) then
         do j=1,2*ne
           do i=1,ne
             box4(i,j) = tmp_box2(2*ne-j+1,i)
           end do
         end do
         start_rank = prev_start_rank
         start_i = prev_start_i
       else
         do j=1,2*ne
           do i=1,ne
             box4(i,j) = box2(2*ne-j+1,i)
           end do
         end do
       end if
       exit
     else
       tmp_box2(:,:) = box2(:,:)
     end if
   end do
   if (debug) print *, '========== end panel 6 =========='
!
!
! panel 1
!
   box4(ne+1:2*ne,1:ne) = -1 ! candidate mask
   start_i = 1
   do
     call find_optimal_band(2*ne, 2*ne, nproc, start_rank, start_i, nelems,    &
         box4, ret)
     start_rank = ret(1)
     start_i = ret(2)
!
     if (start_i .gt. ne) then
       box2(1:ne,:) = box4(ne+1:2*ne,1:ne)
       cube_rank(:,:,6) = box4(1:ne,ne+1:2*ne)
       cube_rank(:,:,1) = box4(1:ne,1:ne)
       exit
     end if
   end do
   if (debug) print *, '========== end panel 1 =========='
!
!
! panel 2
!
   box2(ne+1:2*ne,:) = -1  ! candidate mask
   start_i = start_i - ne
   do
     call find_optimal_band(2*ne, ne, nproc, start_rank, start_i, nelems,      &
         box2, ret)
     start_rank = ret(1)
     start_i = ret(2)
!
     if (start_i .gt. ne) then
       do j=1,ne
         do i=1,ne
           tmp_box2(i,j) = box2(ne+i,ne-j+1)
           cube_rank(i,j,2) = box2(i,j)
         end do
       end do
       exit
     end if
   end do
   if (debug) print *, '========== end panel 2 =========='
!
!
! panel 3
!
   box2(1:ne,:) = tmp_box2(1:ne,:)
   box2(ne+1:2*ne,:) = -1  ! candidate mask
   start_i = start_i - ne
   do 
     prev_start_rank = start_rank
     prev_start_i = start_i
     call find_optimal_band(2*ne, ne, nproc, start_rank, start_i, nelems,      &
         box2, ret)
     start_rank = ret(1)
     start_i = ret(2)
!
     if (start_i .gt. ne) then
       if (count(box2(ne:2*ne,:) .ne. -1) .gt.                                 &
           count(tmp_box2(1:ne,:) .eq. -1)) then
         do j=1,2*ne
           do i=1,ne
             box4(i,j) = tmp_box2(2*ne-j+1,i)
           end do
         end do
         start_rank = prev_start_rank
         start_i = prev_start_i
       else
         do j=1,2*ne
           do i=1,ne
             box4(i,j) = box2(2*ne-j+1,i)
           end do
         end do
       end if
       exit
     else
       tmp_box2(:,:) = box2(:,:)
     end if
   end do
   if (debug) print *, '========== end panel 3 =========='
!
!
! panel 4,5
!
   box4(ne+1:2*ne,1:ne) = -1  ! candidate mask
   start_i = 1
   do
     prev_start_rank = start_rank
     call find_optimal_band(2*ne, 2*ne, nproc, start_rank, start_i, nelems,    &
         box4, ret)
     start_rank = ret(1)
     start_i = ret(2)
!
     if (count(box4.eq.-1) .eq. 0) then
       do j=1,ne
       do i=start_i,2*ne
         if (box4(i,j).eq.-1) box4(i,j) = 0
       end do
       end do
       do j=1,ne
         do i=1,ne
           cube_rank(i,j,3) = box4(ne-j+1,2*ne-i+1)
           cube_rank(i,j,4) = box4(ne-j+1,ne-i+1)
           cube_rank(i,j,5) = box4(ne+i,ne-j+1)
         end do
       end do
       exit
     end if
   end do
   if (debug) print *, '========== end panel 4,5 =========='
   if (any(cube_rank.eq.-1)) stop 'cube_rank has -1 rank number'
!
!
! rearrange last two bands
!
   band1_start_rank = box4(2*ne,1)
   do i=2*ne-1,ne+1,-1
     band2_start_rank = box4(i,1)
     if (band2_start_rank .ne. band1_start_rank) exit
   end do
   num_block = nproc - band2_start_rank
   if (debug) print *, 'num_block', num_block
!
   box(:,:) = box4(ne+1:2*ne,1:ne)
   start_i = 1
   do i=ne,1,-1
     do j=ne,1,-1
       if (box(i,j) .ge. band2_start_rank) then
         box(i,j) = -1
         start_i = i
       end if
     end do
   end do
   if (debug) print *, 'start_i', start_i
!
   if (count(box(:,:).eq.-1) .lt. ne*ne) then
     allocate(ratios(num_block))
     allocate(boxes(ne,ne,num_block))
!
     tmp_box(:,:) = box(:,:)
     i12(:) = (/start_i, start_i, count(tmp_box(start_i,:).eq.-1), -1/)
     ratio2 = calc_perimeter_ratio(ne, ne, nproc,                            &
         band2_start_rank, nproc-1, nelems, i12, tmp_box) 
     ratios(num_block) = ratio2
     boxes(:,:,num_block) = tmp_box(:,:)
!
     do k=1,num_block-1
       tmp_box(:,:) = box(:,:)
       i12(:) = (/start_i, start_i, count(tmp_box(start_i,:).eq.-1), -1/)
       ratio2 = calc_perimeter_ratio(ne, ne, nproc,                            &
           band2_start_rank, nproc-k-1, nelems, i12, tmp_box) 
       nelem2 = sum(nelems(band2_start_rank:nproc-k-1))
!
       i12(:) = (/i12(2), i12(2), count(tmp_box(i12(2),:).eq.-1), -1/)
       ratio1 = calc_perimeter_ratio(ne, ne, nproc,                            &
           nproc-k, nproc-1, nelems, i12, tmp_box) 
       nelem1 = sum(nelems(nproc-k:nproc-1))
!
       ratios(k) = (ratio2*nelem2 + ratio1*nelem1)/(nelem2 + nelem1)
       boxes(:,:,k) = tmp_box(:,:)
     end do
!
     k = minloc(ratios, dim=1)
     do j=1,ne
       do i=1,ne
         cube_rank(i,j,5) = boxes(i,ne-j+1,k)
       end do
     end do
!
     deallocate(boxes)
     deallocate(ratios)
   end if
   if (debug) print *, '========== end rearrange last two bands =========='
!
   if (any(cube_rank.eq.-1)) stop 'cube_rank has -1 rank number'
!
   end subroutine band_partition
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_cube_rank(ne, nproc, nelems, cube_rank, cube_lid)
!-------------------------------------------------------------------------------
! assign the process numbers to elements on the cubed-sphere
! using the band-partitioning
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
   integer :: remain_elem
   integer :: proc
   integer, dimension(0:nproc-1) :: lids
!-------------------------------------------------------------------------------
!
! set the nelems
!
   remain_elem = mod(ne*ne*6, nproc)
   do i=1,nproc
     nelems(i) = ne*ne*6/nproc
     if (i .gt. nproc-remain_elem) nelems(i) = nelems(i) + 1
   end do
!
!
! partitioning
!
   if (nproc .eq. 1) then
     cube_rank(:,:,:) = 0
   else if (nproc .eq. 2) then
     cube_rank(:,:,6) = 0
     cube_rank(:,:,1:2) = 0
     cube_rank(:,:,3:5) = 1
   else if (nproc .eq. 3) then
     cube_rank(:,:,6) = 0
     cube_rank(:,:,1) = 0
     cube_rank(:,:,2:3) = 1
     cube_rank(:,:,4:5) = 2
   else
     call band_partition(ne, nproc, nelems, cube_rank)
   end if
!
!
! local numbering
!
   lids(:) = 1
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         proc = cube_rank(ei,ej,p)
         cube_lid(ei,ej,p) = lids(proc)
         lids(proc) = lids(proc) + 1
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
   integer, intent(in   ) :: ne, iproc, nelem
   integer, intent(in   ) :: cube_rank(ne,ne,6)
   integer, intent(in   ) :: cube_lid(ne,ne,6)
   integer, intent(  out) :: elem_coord(3,nelem)
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
   function global_perimeter_ratio(ne, nproc, cube_rank, num_nbrs)             &
       result(perimeter_ratio)
!-------------------------------------------------------------------------------
! count neighbor elements which have same rank
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, intent(in   ) :: ne
   integer, intent(in   ) :: nproc
   integer, intent(in   ) :: cube_rank(ne,ne,6)
   integer, intent(  out) :: num_nbrs(2,0:nproc-1)  !(num elem, diff nbrs)
   real(8) :: perimeter_ratio
!
   integer :: ei, ej, p
   integer :: myrank, w_rank, e_rank, s_rank, n_rank
   integer :: eij(4)  ! (ei,ej,panel,rot)
   integer :: diff_sides
!-------------------------------------------------------------------------------
!
   num_nbrs(:,:) = 0
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         myrank = cube_rank(ei,ej,p)
!
         diff_sides = 0  ! contact with an element having the different rank
!
         call convert_nbr_eij(ne, ei-1, ej, p, eij)
         w_rank = cube_rank(eij(1),eij(2),eij(3))
         if (w_rank .ne. myrank) diff_sides = diff_sides + 1
!
         call convert_nbr_eij(ne, ei+1, ej, p, eij)
         e_rank = cube_rank(eij(1),eij(2),eij(3))
         if (e_rank .ne. myrank) diff_sides = diff_sides + 1
!
         call convert_nbr_eij(ne, ei, ej-1, p, eij)
         s_rank = cube_rank(eij(1),eij(2),eij(3))
         if (s_rank .ne. myrank) diff_sides = diff_sides + 1
!
         call convert_nbr_eij(ne, ei, ej+1, p, eij)
         n_rank = cube_rank(eij(1),eij(2),eij(3))
         if (n_rank .ne. myrank) diff_sides = diff_sides + 1
!
         num_nbrs(1,myrank) = num_nbrs(1,myrank) + 1
         num_nbrs(2,myrank) = num_nbrs(2,myrank) + diff_sides
       end do
     end do
   end do
!
   perimeter_ratio = sum(num_nbrs(2,:)*1.D0/num_nbrs(1,:))/nproc
!
   end function global_perimeter_ratio
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function global_communication_ratio(ne, np, nproc, cube_rank, num_pts)      &
       result(comm_ratio)
!-------------------------------------------------------------------------------
! communication points/computation points ratio
! of the spectral-element method
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, intent(in   ) :: ne, np
   integer, intent(in   ) :: nproc
   integer, intent(in   ) :: cube_rank(ne,ne,6)
   integer, intent(  out) :: num_pts(2,0:nproc-1)  !(comp points, comm points)
   real(8) :: comm_ratio
!
   integer :: ei, ej, p
   integer :: myrank, nbr_rank
   integer :: eij(4)  ! (ei,ej,panel,rot)
   integer :: comm_pts
!-------------------------------------------------------------------------------
!
   num_pts(:,:) = 0
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         myrank = cube_rank(ei,ej,p)
!
         comm_pts = 0  ! communication points
!
!
! side
!
         call convert_nbr_eij(ne, ei-1, ej, p, eij)
         nbr_rank = cube_rank(eij(1),eij(2),eij(3))
         if (nbr_rank .ne. myrank) comm_pts = comm_pts + np
!
         call convert_nbr_eij(ne, ei+1, ej, p, eij)
         nbr_rank = cube_rank(eij(1),eij(2),eij(3))
         if (nbr_rank .ne. myrank) comm_pts = comm_pts + np
!
         call convert_nbr_eij(ne, ei, ej-1, p, eij)
         nbr_rank = cube_rank(eij(1),eij(2),eij(3))
         if (nbr_rank .ne. myrank) comm_pts = comm_pts + np
!
         call convert_nbr_eij(ne, ei, ej+1, p, eij)
         nbr_rank = cube_rank(eij(1),eij(2),eij(3))
         if (nbr_rank .ne. myrank) comm_pts = comm_pts + np
!
!
! corner
!
         call convert_nbr_eij(ne, ei-1, ej-1, p, eij)
         if (eij(3) .ne. -1) then
           nbr_rank = cube_rank(eij(1),eij(2),eij(3))
           if (nbr_rank .ne. myrank) comm_pts = comm_pts + 1
         end if
!
         call convert_nbr_eij(ne, ei+1, ej-1, p, eij)
         if (eij(3) .ne. -1) then
           nbr_rank = cube_rank(eij(1),eij(2),eij(3))
           if (nbr_rank .ne. myrank) comm_pts = comm_pts + 1
         end if
!
         call convert_nbr_eij(ne, ei-1, ej+1, p, eij)
         if (eij(3) .ne. -1) then
           nbr_rank = cube_rank(eij(1),eij(2),eij(3))
           if (nbr_rank .ne. myrank) comm_pts = comm_pts + 1
         end if
!
         call convert_nbr_eij(ne, ei+1, ej+1, p, eij)
         if (eij(3) .ne. -1) then
           nbr_rank = cube_rank(eij(1),eij(2),eij(3))
           if (nbr_rank .ne. myrank) comm_pts = comm_pts + 1
         end if
!
         num_pts(1,myrank) = num_pts(1,myrank) + np*np
         num_pts(2,myrank) = num_pts(2,myrank) + comm_pts
       end do
     end do
   end do
!
   comm_ratio = sum(num_pts(2,:)*1.D0/num_pts(1,:))/nproc
!
   end function global_communication_ratio
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_cube_color(ne, nproc, cube_rank, cube_color)
!-------------------------------------------------------------------------------
! assign a different number to other adjacent domains for coloring
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, intent(in   ) :: ne
   integer, intent(in   ) :: nproc
   integer, intent(in   ) :: cube_rank(ne,ne,6)
   integer, intent(  out) :: cube_color(ne,ne,6)
!
   integer, parameter :: max_nbr=20  ! empirical
   integer :: k
   integer :: ei, ej, p
   integer :: myrank, nbr_rank
   integer :: a, b, eij(4)  ! (ei,ej,panel,rot)
   integer :: proc_links(max_nbr,0:nproc-1)
   integer :: proc_colors(0:nproc-1)  ! color index 1~6
   integer :: nbr_colors(max_nbr)
!-------------------------------------------------------------------------------
!
   proc_links(:,:) = -1
!
!
! find neighbor ranks
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         myrank = cube_rank(ei,ej,p)
!
         do b=-1,1
           do a=-1,1
             call convert_nbr_eij(ne, ei+a, ej+b, p, eij)
             if (eij(3) .eq. -1) cycle
             nbr_rank = cube_rank(eij(1),eij(2),eij(3))
             if (nbr_rank .ne. myrank) then
               do k=1,max_nbr
                 if (proc_links(k,myrank) .eq. nbr_rank) exit
                 if (proc_links(k,myrank) .eq. -1) then
                   proc_links(k,myrank) = nbr_rank
                   exit
                 end if
               end do
               !if (k .eq. max_nbr) stop 'max_nbr must be increased'
             end if
           end do
         end do
!
       end do
     end do
   end do
!
!
! color indexing
!
   proc_colors(:) = -1
   proc_colors(0) = 1
   do myrank=1,nproc-1
     nbr_colors(:) = -1
     do k=1,max_nbr
       if (proc_links(k,myrank) .ne. -1) then
         nbr_colors(k) = proc_colors(proc_links(k,myrank))
       end if
     end do
!
     k = 1
     do
       if (all(nbr_colors.ne.k)) then
         proc_colors(myrank) = k
         exit
       end if
       k = k + 1
     end do
   end do
!
   do p=1,6
     do ej=1,ne
       do ei=1,ne
         myrank = cube_rank(ei,ej,p)
         cube_color(ei,ej,p) = proc_colors(myrank)
       end do
     end do
   end do
!
   end subroutine make_cube_color
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
end module cube_partition_band
!-------------------------------------------------------------------------------
