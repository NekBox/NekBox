
!> Geometry arrays
module geom
  use kinds, only : DP
  implicit none

  real(DP), allocatable, dimension(:,:,:,:) :: &
     XM1  &
    ,YM1 &
    ,ZM1 &
    ,XM2 &
    ,YM2 &
    ,ZM2

  real(DP), allocatable, dimension(:,:,:,:) :: &
     RXM1 &
    ,SXM1 &
    ,TXM1 &
    ,RYM1 &
    ,SYM1 &
    ,TYM1 &
    ,RZM1 &
    ,SZM1 &
    ,TZM1 &
    ,JACM1 
  real(DP), allocatable :: jacmi (:,:)

  real(DP), allocatable, dimension(:,:,:,:) :: &
     RXM2 &
    ,SXM2 &
    ,TXM2 &
    ,RYM2 &
    ,SYM2 &
    ,TYM2 &
    ,RZM2 &
    ,SZM2 &
    ,TZM2 &
    ,JACM2 

  real(DP), allocatable :: rx(:,:,:)

  real(DP), allocatable, dimension(:,:,:,:) :: &
     G1M1 &
    ,G2M1 &
    ,G3M1 &
    ,G4M1 &
    ,G5M1 &
    ,G6M1

  real(DP), allocatable, dimension(:,:,:,:) :: &
     UNX &
    ,UNY &
    ,UNZ &
    ,T1X &
    ,T1Y &
    ,T1Z &
    ,T2X &
    ,T2Y &
    ,T2Z &
    ,AREA 
  real(DP) :: DLAM

  real(DP), allocatable, dimension(:,:,:,:) :: &
     VNX &
    ,VNY &
    ,VNZ &
    ,V1X &
    ,V1Y &
    ,V1Z &
    ,V2X &
    ,V2Y &
    ,V2Z 

  logical :: IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO, IFBCOR 
  logical, allocatable :: IFRZER(:), IFQINP(:,:), IFEPPM(:,:)
  logical, allocatable :: IFLMSF(:), IFLMSE(:), IFLMSC(:)
  logical, allocatable :: IFMSFC(:,:,:), IFMSEG(:,:,:),IFMSCR(:,:,:)
  logical, allocatable :: IFNSKP(:,:) 

  contains

  subroutine init_geom()
    use size_m
    implicit none

    allocate( &
     XM1 (LX1,LY1,LZ1,LELT) &
    ,YM1 (LX1,LY1,LZ1,LELT) &
    ,ZM1 (LX1,LY1,LZ1,LELT) &
    ,XM2 (LX2,LY2,LZ2,LELV) &
    ,YM2 (LX2,LY2,LZ2,LELV) &
    ,ZM2 (LX2,LY2,LZ2,LELV) &
    )

    allocate( &
     RXM1 (LX1,LY1,LZ1,LELT) &
    ,SXM1 (LX1,LY1,LZ1,LELT) &
    ,TXM1 (LX1,LY1,LZ1,LELT) &
    ,RYM1 (LX1,LY1,LZ1,LELT) &
    ,SYM1 (LX1,LY1,LZ1,LELT) &
    ,TYM1 (LX1,LY1,LZ1,LELT) &
    ,RZM1 (LX1,LY1,LZ1,LELT) &
    ,SZM1 (LX1,LY1,LZ1,LELT) &
    ,TZM1 (LX1,LY1,LZ1,LELT) &
    ,JACM1 (LX1,LY1,LZ1,LELT) &
    ,jacmi (lx1*ly1*lz1,lelt) &
    )

    allocate( &
     RXM2 (LX2,LY2,LZ2,LELV) &
    ,SXM2 (LX2,LY2,LZ2,LELV) &
    ,TXM2 (LX2,LY2,LZ2,LELV) &
    ,RYM2 (LX2,LY2,LZ2,LELV) &
    ,SYM2 (LX2,LY2,LZ2,LELV) &
    ,TYM2 (LX2,LY2,LZ2,LELV) &
    ,RZM2 (LX2,LY2,LZ2,LELV) &
    ,SZM2 (LX2,LY2,LZ2,LELV) &
    ,TZM2 (LX2,LY2,LZ2,LELV) &
    ,JACM2 (LX2,LY2,LZ2,LELV) &
    )

    allocate(rx(lxd*lyd*lzd,ldim*ldim,lelv))

    allocate( &
     G1M1 (LX1,LY1,LZ1,LELT) &
    ,G2M1 (LX1,LY1,LZ1,LELT) &
    ,G3M1 (LX1,LY1,LZ1,LELT) &
    ,G4M1 (LX1,LY1,LZ1,LELT) &
    ,G5M1 (LX1,LY1,LZ1,LELT) &
    ,G6M1 (LX1,LY1,LZ1,LELT) &
    )

    allocate( &
     UNX (LX1,LZ1,6,LELT) &
    ,UNY (LX1,LZ1,6,LELT) &
    ,UNZ (LX1,LZ1,6,LELT) &
    ,T1X (LX1,LZ1,6,LELT) &
    ,T1Y (LX1,LZ1,6,LELT) &
    ,T1Z (LX1,LZ1,6,LELT) &
    ,T2X (LX1,LZ1,6,LELT) &
    ,T2Y (LX1,LZ1,6,LELT) &
    ,T2Z (LX1,LZ1,6,LELT) &
    ,AREA (LX1,LZ1,6,LELT) &
    )

    allocate( &
     VNX (LX1M,LY1M,LZ1M,LELT) &
    ,VNY (LX1M,LY1M,LZ1M,LELT) &
    ,VNZ (LX1M,LY1M,LZ1M,LELT) &
    ,V1X (LX1M,LY1M,LZ1M,LELT) &
    ,V1Y (LX1M,LY1M,LZ1M,LELT) &
    ,V1Z (LX1M,LY1M,LZ1M,LELT) &
    ,V2X (LX1M,LY1M,LZ1M,LELT) &
    ,V2Y (LX1M,LY1M,LZ1M,LELT) &
    ,V2Z (LX1M,LY1M,LZ1M,LELT) &
    )

    allocate( &
     IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV) &
    ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1) &
    ,IFMSFC(6,LELT,0:1) &
    ,IFMSEG(12,LELT,0:1) &
    ,IFMSCR(8,LELT,0:1) &
    ,IFNSKP(8,LELT) &
    )

  end subroutine init_geom

end module geom
