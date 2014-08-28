module mesh
  implicit none

  integer, allocatable :: vertex(:)
  logical, allocatable :: ifdfrm(:) !>!< is the element deformed?
  logical, allocatable :: iffast(:) !>!< can we use a fast method on the element?
  logical :: ifsolv !>!< are ifdfrm and iffast up to date?

  contains

  subroutine init_mesh()
    use size_m, only : lelt, ldim
    implicit none

    allocate(vertex((2**ldim)*lelt))
    allocate(ifdfrm(LELT), iffast(lelt))
    
  end subroutine init_mesh

end module mesh
