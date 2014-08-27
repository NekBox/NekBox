module mesh
  implicit none

  integer, allocatable :: vertex(:)

  contains

  subroutine init_mesh()
    use size_m, only : lelt, ldim
    implicit none

    allocate(vertex((2**ldim)*lelt))
    
  end subroutine init_mesh

end module mesh
