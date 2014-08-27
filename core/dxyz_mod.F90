
!     Elemental derivative operators
module dxyz
  use kinds, only : DP
  implicit none

  real(DP), allocatable, dimension(:,:) :: &
     DXM1,  DXM12 &
    ,DYM1,  DYM12 &
    ,DZM1,  DZM12 &
    ,DXTM1, DXTM12 &
    ,DYTM1, DYTM12 &
    ,DZTM1, DZTM12 &
    ,DXM3,  DXTM3 &
    ,DYM3,  DYTM3 &
    ,DZM3,  DZTM3 &
    ,DCM1,  DCTM1 &
    ,DCM3,  DCTM3 &
    ,DCM12, DCTM12 &
    ,DAM1,  DATM1 &
    ,DAM12, DATM12 &
    ,DAM3,  DATM3

  real(DP), allocatable :: WDDX(:,:), WDDYT(:,:), WDDZT(:,:) 

  contains

  subroutine init_dxyz()
    use size_m
    implicit none

    allocate( &
     DXM1(LX1,LX1),  DXM12(LX2,LX1) &
    ,DYM1(LY1,LY1),  DYM12(LY2,LY1) &
    ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1) &
    ,DXTM1(LX1,LX1), DXTM12(LX1,LX2) &
    ,DYTM1(LY1,LY1), DYTM12(LY1,LY2) &
    ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2) &
    ,DXM3(LX3,LX3),  DXTM3(LX3,LX3) &
    ,DYM3(LY3,LY3),  DYTM3(LY3,LY3) &
    ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3) &
    ,DCM1(LY1,LY1),  DCTM1(LY1,LY1) &
    ,DCM3(LY3,LY3),  DCTM3(LY3,LY3) &
    ,DCM12(LY2,LY1), DCTM12(LY1,LY2) &
    ,DAM1(LY1,LY1),  DATM1(LY1,LY1) &
    ,DAM12(LY2,LY1), DATM12(LY1,LY2) &
    ,DAM3(LY3,LY3),  DATM3(LY3,LY3) &
            )

    allocate(WDDX(LX1,LY1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1))

  end subroutine init_dxyz

end module dxyz
