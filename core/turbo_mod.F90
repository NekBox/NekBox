!> cleaned
module turbo
!       Common block for turbulence model
  use kinds, only : DP
  implicit none

#if 0
  real(DP), allocatable, dimension(:,:,:,:) :: VTURB, TURBL, UWALL, ZWALL
  real(DP), allocatable, dimension(:,:,:,:) :: TWX, TWY, TWZ

  real(DP) :: CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT &
    , BETA1,BETA2 &
    , CMI,SKI,SEI,VKI,BTI,STI &
    , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
#endif

  integer :: IFLDK,IFLDTK,IFLDE,IFLDTE
  logical :: IFSWALL,IFCWUZ
!  logical, allocatable :: IFTWSH(:,:)

  contains

  subroutine init_turbo
    use size_m
    implicit none

#if 0
    allocate(VTURB (LX1M,LY1M,LZ1M,LELV) &
    , TURBL (LX1M,LY1M,LZ1M,LELV) &
    , UWALL (LX1M,LZ1M,6,LELV) &
    , ZWALL (LX1M,LZ1M,6,LELV) &
    , TWX   (LX1M,LZ1M,6,LELV) &
    , TWY   (LX1M,LZ1M,6,LELV) &
    , TWZ   (LX1M,LZ1M,6,LELV) )

    allocate(IFTWSH(6,LELV))
#endif

  end subroutine init_turbo

end module turbo
