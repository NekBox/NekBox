!> cleaned
module mass
  use kinds, only : DP
  implicit none

  real(DP), allocatable, dimension(:,:,:,:) :: &
    BM1, BM2, BINVM1, BINTM1, BM2INV, BAXM1, YINVM1
  real(DP), allocatable :: BM1LAG(:,:,:,:,:)

  real(DP) :: VOLVM1,VOLVM2,VOLTM1,VOLTM2

  contains

  subroutine init_mass
    use size_m
    implicit none

    allocate( BM1(LX1,LY1,LZ1,LELT),  BM2(LX2,LY2,LZ2,LELV) &
    ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT) &
    ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT) &
    ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1) & ! verified
    ,YINVM1(LX1,LY1,LZ1,LELT))
    bm1lag = 0._dp

  end subroutine init_mass

end module mass
