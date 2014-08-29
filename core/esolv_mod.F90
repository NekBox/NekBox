module esolv
  use kinds, only : DP
  implicit none

  integer :: iesolv
  logical, allocatable :: ifalgn(:), ifrsxy(:)
  real(DP), allocatable :: volel(:)

  contains

  subroutine init_esolv()
    use size_m
    implicit none

!    allocate(ifalgn(lelv), ifrsxy(lelv))
    allocate(volel(lelv)) ! verified
  end subroutine init_esolv

end module esolv
