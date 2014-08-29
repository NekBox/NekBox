!> cleaned
module mvgeom
!     Moving mesh data
  use kinds, only : DP

  real(DP), allocatable :: WX(:,:,:,:) 
  real(DP), allocatable :: WY(:,:,:,:) 
  real(DP), allocatable :: WZ(:,:,:,:) 

  real(DP), allocatable :: WXLAG (:,:,:,:,:)
  real(DP), allocatable :: WYLAG (:,:,:,:,:)
  real(DP), allocatable :: WZLAG (:,:,:,:,:)

  real(DP), allocatable :: W1MASK(:,:,:,:) 
  real(DP), allocatable :: W2MASK(:,:,:,:) 
  real(DP), allocatable :: W3MASK(:,:,:,:) 
  real(DP), allocatable :: WMULT(:,:,:,:) 

  real(DP), allocatable :: EV1(:,:,:,:)
  real(DP), allocatable :: EV2(:,:,:,:)
  real(DP), allocatable :: EV3(:,:,:,:)

  contains

  subroutine init_mvgeom()
    use size_m
    implicit none

!    allocate( &
!        WX    (LX1M,LY1M,LZ1M,LELT) &
!    ,   WY    (LX1M,LY1M,LZ1M,LELT) &
!    ,   WZ    (LX1M,LY1M,LZ1M,LELT) )
!    allocate( WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1) &
!    ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1) &
!    ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1) )
!    allocate( W1MASK(LX1M,LY1M,LZ1M,LELT) &
!    ,   W2MASK(LX1M,LY1M,LZ1M,LELT) &
!    ,   W3MASK(LX1M,LY1M,LZ1M,LELT) &
!    ,   WMULT (LX1M,LY1M,LZ1M,LELT) )
!    allocate( EV1   (LX1M,LY1M,LZ1M,LELV) &
!    , EV2   (LX1M,LY1M,LZ1M,LELV) &
!    , EV3   (LX1M,LY1M,LZ1M,LELV) )

  end subroutine init_mvgeom

end module mvgeom
