!> cleaned
module mesh
  use kinds, only : DP
  implicit none

  integer, allocatable :: vertex(:)
  logical, allocatable :: ifdfrm(:) !>!< is the element deformed?
  logical, allocatable :: iffast(:) !>!< can we use a fast method on the element?
  logical :: ifsolv !>!< are ifdfrm and iffast up to date?
  integer :: niterhm
  logical :: if_ortho = .true.

  real(DP) :: start_x(3)
  real(DP) :: end_x(3)
  integer  :: shape_x(3)



  contains

  subroutine init_mesh()
    use size_m, only : lelt, ldim
    implicit none

    ifsolv = .false.
    allocate(vertex((2**ldim)*lelt))
    allocate(ifdfrm(LELT), iffast(lelt))
    
  end subroutine init_mesh

  function ieg_to_xyz(ieg) result(xyz)
    integer, intent(in) :: ieg
    integer :: xyz(3)
 
    xyz(1) = mod(ieg - 1, shape_x(1))
    xyz(2) = mod((ieg-1)/shape_x(1), shape_x(2))
    xyz(3) = mod((ieg-1)/(shape_x(1)*shape_x(2)), shape_x(3))

  end function ieg_to_xyz

  integer function xyz_to_ieg(xyz)
    implicit none
    integer, intent(in) :: xyz(3)
    xyz_to_ieg = 1 + xyz(1) + xyz(2)*shape_x(1) + xyz(3)*shape_x(1)*shape_x(2)
    return
  end function xyz_to_ieg

end module mesh
