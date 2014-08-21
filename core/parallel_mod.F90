module parallel
!     Communication information
!     NOTE: NID is stored in 'SIZE' for greater accessibility
  implicit none

  INTEGER ::        NODE,PID,NP,NULLPID,NODE0

!     Maximum number of elements (limited to 2**31/12, at least for now)
  integer, parameter :: NELGT_MAX = 178956970

  integer, allocatable :: nelg(:), lglel(:), gllel(:), gllnid(:)
  integer :: nvtot, nelgv, nelgt

  LOGICAL :: IFGPRNT
  INTEGER :: WDSIZE,ISIZE,LSIZE,CSIZE,WDSIZI
  LOGICAL :: IFDBLAS

!     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solve)

  integer :: cr_h, gsh
  integer, allocatable :: gsh_fld(:), xxth(:)

  contains

  subroutine init_parallel()
    use size_m
    implicit none

    allocate(NELG(0:LDIMT1), LGLEL(LELT), GLLEL(LELG), GLLNID(LELG))
    allocate(gsh_fld(0:ldimt3), xxth(ldimt3))

  end subroutine init_parallel

end module parallel
