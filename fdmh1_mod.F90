module fdmh1
!>     'FDMH1'
  use kinds, only : DP
  implicit none

  real(DP), allocatable :: dd(:,:),fds(:,:), fdst(:,:)
  real(DP), allocatable :: elsize(:,:), bhalf(:,:,:,:)

  integer :: kfldfdm
  integer, allocatable :: ktype(:,:,:)

  logical ::         ifbhalf

  contains

  subroutine init_fdmh1()
    use size_m
    implicit none

    allocate(ktype(lelt, 3, 0:4))
    allocate(dd(lx1,9),fds(lx1*lx1,9),fdst(lx1*lx1,9) &
    , elsize(3,lelt) &
    , bhalf(lx1,ly1,lz1,lelt) &
    )

  end subroutine init_fdmh1

end module fdmh1
