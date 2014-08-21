module opctr
!     OPCTR is a set of arrays for tracking the number of operations,
!     and number of calls for a particular subroutine
  use kinds, only : DP
  implicit none

  integer, parameter :: maxrts = 1000
  character(6) ::     rname(MAXRTS)

  real(DP) :: dcount,dct(MAXRTS),rct(MAXRTS)
  integer :: ncall(MAXRTS), nrout

  integer, save :: myrout,isclld
  data    myrout /0/
  data    isclld /0/

end module opctr
