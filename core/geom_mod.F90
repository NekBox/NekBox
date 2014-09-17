!> Geometry arrays
module geom
  use kinds, only : DP
  implicit none

  real(DP), allocatable, target, dimension(:,:,:,:) :: XM1, YM1, ZM1
  real(DP), pointer,             dimension(:,:,:,:) :: XM2, YM2, ZM2

  real(DP), allocatable :: jacmi (:,:,:,:)

  real(DP), allocatable, target, dimension(:,:,:,:) :: &
     RXM1, SXM1, TXM1, RYM1, SYM1, TYM1, RZM1, SZM1, TZM1, JACM1 
  real(DP), pointer,             dimension(:,:,:,:) :: &
     RXM2, SXM2, TXM2, RYM2, SYM2, TYM2, RZM2, SZM2, TZM2, JACM2

  real(DP), allocatable :: rx(:,:,:) !>!< saved data for lx1=>lxd interpolation

  real(DP), allocatable, dimension(:,:,:,:) :: &
     G1M1, G2M1, G3M1, G4M1, G5M1, G6M1 

  real(DP), allocatable, dimension(:,:,:,:) :: &
     UNX, UNY, UNZ, AREA 

  real(DP) :: DLAM

  real(DP), allocatable, dimension(:,:,:,:) :: &
     VNX, VNY, VNZ, V1X, V1Y, V1Z, V2X, V2Y, V2Z 

  real(DP), allocatable, target :: BM1(:,:,:,:) 
  real(DP), pointer             :: BM2(:,:,:,:) 

  real(DP), allocatable, dimension(:,:,:,:) :: &
    BINVM1, BINTM1, BM2INV, BAXM1, YINVM1

  real(DP), allocatable :: BM1LAG(:,:,:,:,:)

  real(DP) :: VOLVM1,VOLVM2,VOLTM1,VOLTM2

  logical :: IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO, IFBCOR 
  logical, allocatable :: IFRZER(:), IFQINP(:,:), IFEPPM(:,:)
  logical, allocatable :: IFLMSF(:), IFLMSE(:), IFLMSC(:)
  logical, allocatable :: IFMSFC(:,:,:), IFMSEG(:,:,:),IFMSCR(:,:,:)
  logical, allocatable :: IFNSKP(:,:) 

  logical :: bm1_compress !>!< are bm1's elements identical?

  contains

  subroutine init_geom()
    use size_m
    implicit none

    allocate( XM1(LX1,LY1,LZ1,LELT), YM1(LX1,LY1,LZ1,LELT), ZM1(LX1,LY1,LZ1,LELT) )

    allocate( RXM1(LX1,LY1,LZ1,LELT), SXM1(LX1,LY1,LZ1,LELT), TXM1(LX1,LY1,LZ1,LELT) &
            , RYM1(LX1,LY1,LZ1,LELT), SYM1(LX1,LY1,LZ1,LELT), TYM1(LX1,LY1,LZ1,LELT) &
            , RZM1(LX1,LY1,LZ1,LELT), SZM1(LX1,LY1,LZ1,LELT), TZM1(LX1,LY1,LZ1,LELT) )
    allocate( JACM1(LX1,LY1,LZ1,LELT), jacmi(lx1,ly1,lz1,lelt) )


    allocate(rx(lxd*lyd*lzd,ldim*ldim,lelv)) !verified

    allocate( G1M1(LX1,LY1,LZ1,LELT), G2M1(LX1,LY1,LZ1,LELT), G3M1(LX1,LY1,LZ1,LELT) &
            , G4M1(LX1,LY1,LZ1,LELT), G5M1(LX1,LY1,LZ1,LELT), G6M1(LX1,LY1,LZ1,LELT) )

!> \todo remove un{x,y,z} for regular geometry case
    allocate( UNX(LX1,LZ1,6,LELT), UNY(LX1,LZ1,6,LELT), UNZ(LX1,LZ1,6,LELT) )
    allocate( AREA(LX1,LZ1,6,LELT) )

    allocate( &
     IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV) &
    ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1) &
    ,IFMSFC(6,LELT,0:1) &
    ,IFMSEG(12,LELT,0:1) &
    ,IFMSCR(8,LELT,0:1) &
    ,IFNSKP(8,LELT) &
    )

    allocate( BM1(LX1,LY1,LZ1,LELT), BINVM1(LX1,LY1,LZ1,LELV) )
    IFLMSF = .False.

  end subroutine init_geom

  subroutine compress_geom()
    use size_m, only : lx1,ly1,lz1,nelt, nid
    implicit none

    integer, parameter :: lxyz = lx1*ly1*lz1
    real(DP) :: thresh
    real(DP), external :: dnrm2
    integer :: ie

    bm1_compress = .true.
    thresh = 1.d-4 * dnrm2(lxyz*nelt, bm1, 1) / nelt
    
    do ie = 2, nelt
      if (dnrm2(lxyz,bm1(:,:,:,ie) - bm1(:,:,:,ie-1),1) > thresh) then
        bm1_compress = .false.
        if (nid == 0) then
        write(*,*) bm1(:,1,1,ie)
        write(*,*) bm1(:,1,1,ie-1)
        endif
      endif
    enddo
    write(*,*) "MAX: bm1_compress =", bm1_compress

  end subroutine compress_geom

end module geom
