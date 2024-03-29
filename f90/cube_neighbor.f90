!-------------------------------------------------------------------------------
   module cube_neighbor
!-------------------------------------------------------------------------------
!
! abstract : find a neighbor element/point on the cubed-sphere grid
!
! history log :
!   2017-04-15  ki-hwan kim  start
!   2017-05-09  ki-hwan kim  remove convert_neighbor_ij(), find_nbr_elem_ij()
!                            add convert_nbr_elem_ij()
!   2017-05-15  ki-hwan kim  convert_nbr_elem_ij() -> convert_nbr_eij()
!                            add convert_nbr_gij()
!   2017-06-16  ki-hwan kim  init_nbr_panels() -> to init_nbr_panels()
!
!-------------------------------------------------------------------------------
!
   implicit none
!
   private
   integer, dimension(2,4,6) :: nbr_panels
!
   public :: quotient
   public :: init_nbr_panels
   public :: convert_rotated_ij
   public :: convert_nbr_eij
   public :: convert_nbr_gij
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function quotient(n, i) result(q)
!-------------------------------------------------------------------------------
! example: n=3
! i   : -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7
! i/n : -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2
!-------------------------------------------------------------------------------
!
   integer, intent(in   ) :: n, i
!
   integer :: q
!-------------------------------------------------------------------------------
!
   if (i .lt. 0) then
     q = (i-n+1)/n
   else
     q = i/n
   end if
!
   end function quotient
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine init_nbr_panels()
!-------------------------------------------------------------------------------
! initialize the nbr_panels
! shape: (panel, rotation) x (4 directions) x (6 panels)
!
! sequence of neighbor panels:
!    ---4---
!   |       |
!   1       2
!   |       |
!    ---3---
!-------------------------------------------------------------------------------
!
   ! panel 1
   nbr_panels(:,1,1) = (/4,0/)
   nbr_panels(:,2,1) = (/2,0/)
   nbr_panels(:,3,1) = (/5,0/)
   nbr_panels(:,4,1) = (/6,0/)
   ! panel 2
   nbr_panels(:,1,2) = (/1,0/)
   nbr_panels(:,2,2) = (/3,0/)
   nbr_panels(:,3,2) = (/5,1/)
   nbr_panels(:,4,2) = (/6,3/)
   ! panel 3
   nbr_panels(:,1,3) = (/2,0/)
   nbr_panels(:,2,3) = (/4,0/)
   nbr_panels(:,3,3) = (/5,2/)
   nbr_panels(:,4,3) = (/6,2/)
   ! panel 4
   nbr_panels(:,1,4) = (/3,0/)
   nbr_panels(:,2,4) = (/1,0/)
   nbr_panels(:,3,4) = (/5,3/)
   nbr_panels(:,4,4) = (/6,1/)
   ! panel 5
   nbr_panels(:,1,5) = (/4,1/)
   nbr_panels(:,2,5) = (/2,3/)
   nbr_panels(:,3,5) = (/3,2/)
   nbr_panels(:,4,5) = (/1,0/)
   ! panel 6
   nbr_panels(:,1,6) = (/4,3/)
   nbr_panels(:,2,6) = (/2,1/)
   nbr_panels(:,3,6) = (/1,0/)
   nbr_panels(:,4,6) = (/3,2/)
!
   end subroutine init_nbr_panels
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_rotated_ij(n, i, j, rot, ret)
!-------------------------------------------------------------------------------
! return rotated (i,j) in anti-clockwise rotation
! rotation indices: (0, 1, 2, 3)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,               intent(in   ) :: n, i, j, rot
   integer, dimension(2), intent(  out) :: ret
!
   !integer :: ri, rj
!-------------------------------------------------------------------------------
!
   !ri = (abs(rot-2)-1)*i - (abs(rot-1)-1)*j + (n+1)*(rot/2)
   !rj = (abs(rot-1)-1)*i + (abs(rot-2)-1)*j - (n+1)*int(abs(rot-1.5)-1.5)
   !ret(:) = (/ri,rj/)
!
   select case (rot)
     case(0)
       ret(:) = (/i,j/)
     case(1)
       ret(:) = (/j,n-i+1/)
     case(2)
       ret(:) = (/n-i+1,n-j+1/)
     case(3)
       ret(:) = (/n-j+1,i/)
   end select
!
   end subroutine convert_rotated_ij
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_nbr_eij(ne, ei, ej, panel, ret)
!-------------------------------------------------------------------------------
! find a neighbor element in any relative location
! input : (ei+a, ej+b, panel)
! return: (ei2, ej2, panel2, rotation)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,               intent(in   ) :: ne, ei, ej, panel
   integer, dimension(4), intent(  out) :: ret
!
   integer :: qei, qej  ! quotient
   integer :: rei, rej  ! remainder
   integer :: p2, rot, dir
   integer :: eij(2)
!-------------------------------------------------------------------------------
!
   p2 = -1; rot = -1
!
   qei = quotient(ne, ei-1)
   qej = quotient(ne, ej-1)
!
   if ((abs(qei) .ge. 2) .or. (abs(qej) .ge. 2)) then
     stop "out of bounds of neighbors"
!
   else if ((qei .eq. 0) .and. (qej .eq. 0)) then
     ret(:) = (/ei, ej, panel, 0/)
!
   else if (qej .eq. 0) then
     rei = modulo(ei-1,ne) + 1
     dir = int(0.5*(qei+3))
     p2  = nbr_panels(1,dir,panel)
     rot = nbr_panels(2,dir,panel)
     call convert_rotated_ij(ne, rei, ej, rot, eij)
     ret(:) = (/eij(1), eij(2), p2, rot/)
!
   else if (qei .eq. 0) then
     rej = modulo(ej-1,ne) + 1
     dir = int(0.5*(qej+7))
     p2  = nbr_panels(1,dir,panel)
     rot = nbr_panels(2,dir,panel)
     call convert_rotated_ij(ne, ei, rej, rot, eij)
     ret(:) = (/eij(1), eij(2), p2, rot/)
!
   else
     ret(:) = (/-1, -1, -1, -1/)
   end if
!
   end subroutine convert_nbr_eij
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_nbr_gij(ne, np, gi, gj, ei, ej, panel, ret)
!-------------------------------------------------------------------------------
! find a neighbor Gauss-quadrature point in any relative location
! input : (gi+a, gj+b, ei, ej, panel)
! return: (gi2, gj2, ei2, ej2, panel2)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer,               intent(in   ) :: ne, np, gi, gj, ei, ej, panel
   integer, dimension(5), intent(  out) :: ret
!
   integer :: qgi, qgj  ! quotient
   integer :: rgi, rgj  ! remainder
   integer :: eij(4)    ! (ei2, ej2, panel2, rot)
   integer :: rij(2)    ! rotated (gi, gj)
!-------------------------------------------------------------------------------
!
   qgi = quotient(np, gi-1)
   qgj = quotient(np, gj-1)
!
   call convert_nbr_eij(ne, ei+qgi, ej+qgj, panel, eij)
!
   if (eij(3) .ne. -1) then
     rgi = modulo(gi-1,np) + 1
     rgj = modulo(gj-1,np) + 1
     call convert_rotated_ij(np, rgi, rgj, eij(4), rij)
     ret(:) = (/rij(1), rij(2), eij(1), eij(2), eij(3)/)
   else
     ret(:) = (/-1, -1, -1, -1, -1/)
   end if
!
   end subroutine convert_nbr_gij
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
end module cube_neighbor
!-------------------------------------------------------------------------------
