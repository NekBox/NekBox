!> cleaned
module mesh
  use kinds, only : DP
  implicit none

  integer, allocatable :: vertex(:)
  logical, allocatable :: ifdfrm(:) !>!< is the element deformed?
  logical, allocatable :: iffast(:) !>!< can we use a fast method on the element?
  logical :: ifsolv !>!< are ifdfrm and iffast up to date?
  integer :: niterhm

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

end module mesh
