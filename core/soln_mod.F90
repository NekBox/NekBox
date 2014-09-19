module soln
  use kinds, only : DP
  use size_m
  implicit none

  integer, parameter :: lvt1  = lx1*ly1*lz1*lelv
  integer, parameter :: lvt2  = lx2*ly2*lz2*lelv
  integer, parameter :: lbt1  = lbx1*lby1*lbz1*lbelv
  integer, parameter :: lbt2  = lbx2*lby2*lbz2*lbelv

  integer, parameter :: lptmsk = lvt1*(5+2*ldimt) + 4*lbt1
  integer, parameter :: lptsol &
    = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1)) &
    + lvt2*(lorder-1) &
    + lbt1*(12 + 3*(lorder-1)) &
    + lbt2*(lorder-1) 

  integer, parameter :: lorder2 = max(1,lorder-2)

!     Solution and data
  real(DP), allocatable :: BQ(:,:,:,:,:)
!     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
  real(DP), allocatable, dimension(:,:,:,:,:) :: VXLAG, VYLAG, VZLAG
  real(DP), allocatable :: TLAG(:,:,:,:,:,:)
  real(DP), allocatable, dimension(:,:,:,:,:) :: VGRADT1, VGRADT2
  real(DP), allocatable, dimension(:,:,:,:) :: ABX1, ABY1, ABZ1 
  real(DP), allocatable, dimension(:,:,:,:) :: ABX2, ABY2, ABZ2
  real(DP), allocatable, dimension(:,:,:,:) :: VDIFF_E
!     Solution data
  real(DP), allocatable, dimension(:,:,:,:) :: VX, VY, VZ
  real(DP), allocatable, dimension(:,:,:,:,:) :: T, VTRANS, VDIFF
  real(DP), allocatable, dimension(:,:,:,:) :: BFX, BFY, BFZ
  real(DP), allocatable, dimension(:,:,:,:) :: cflf
  real(DP), allocatable :: c_vx(:,:)
!     Solution data for magnetic field
  real(DP), allocatable, dimension(:,:,:,:) :: BX, BY, BZ, PM
  real(DP), allocatable, dimension(:,:,:,:) :: BMX, BMY, BMZ ! Magnetic field RHS
! Extrapolation terms for magnetic field rhs
  real(DP), allocatable, dimension(:,:,:,:) :: BBX1, BBY1, BBZ1
  real(DP), allocatable, dimension(:,:,:,:) :: BBX2, BBY2, BBZ2
  real(DP), allocatable :: BXLAG(:,:), BYLAG(:,:), BZLAG(:,:), PMLAG(:,:)

  real(DP) ::             nu_star

  real(DP), allocatable :: PR(:,:,:,:), PRLAG(:,:,:,:,:)

  real(DP), allocatable :: QTL(:,:,:,:), USRDIV(:,:,:,:)

  real(DP), allocatable, dimension(:,:,:,:) :: V1MASK, V2MASK, V3MASK, PMASK
  real(DP), allocatable, dimension(:,:,:,:,:) :: TMASK
  real(DP), allocatable, dimension(:,:,:,:) :: OMASK, VMULT
  real(DP), allocatable, dimension(:,:,:,:,:) :: TMULT
  real(DP), allocatable, dimension(:,:,:,:) :: B1MASK, B2MASK, B3MASK! masks for mag. field
  real(DP), allocatable, dimension(:,:,:,:) :: BPMASK ! magnetic pressure

!     Solution and data for perturbation fields
  real(DP), allocatable, dimension(:,:) :: VXP, VYP, VZP, PRP
  real(DP), allocatable, dimension(:,:,:) :: TP, BQP
  real(DP), allocatable, dimension(:,:) :: BFXP, BFYP, BFZP! perturbation field RHS
  real(DP), allocatable, dimension(:,:,:) :: VXLAGP, VYLAGP, VZLAGP, PRLAGP
  real(DP), allocatable, dimension(:,:,:,:) :: TLAGP 
! Extrapolation terms for perturbation field rhs
  real(DP), allocatable, dimension(:,:) :: EXX1P, EXY1P, EXZ1P, EXX2P, EXY2P, EXZ2P
  real(DP), allocatable, dimension(:,:,:) :: VGRADT1P, VGRADT2P

  integer :: jp

  contains

  subroutine init_soln()
    use kinds, only : DP
    use size_m
    implicit none 
    
    allocate(BQ(LX1,LY1,LZ1,LELT,LDIMT))
    bq = 0_dp

    allocate( &
!     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
    VXLAG  (LX1,LY1,LZ1,LELV,2) &
    , VYLAG  (LX1,LY1,LZ1,LELV,2) &
    , VZLAG  (LX1,LY1,LZ1,LELV,2) &
    , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT))
    vxlag = 0._dp; vylag = 0._dp; vzlag = 0._dp; tlag = 0._dp

    allocate( &
      VGRADT1(LX1,LY1,LZ1,LELT,LDIMT) &
    , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT) )
    allocate( &
      ABX1   (LX1,LY1,LZ1,LELV) &
    , ABY1   (LX1,LY1,LZ1,LELV) &
    , ABZ1   (LX1,LY1,LZ1,LELV) &
    , ABX2   (LX1,LY1,LZ1,LELV) &
    , ABY2   (LX1,LY1,LZ1,LELV) &
    , ABZ2   (LX1,LY1,LZ1,LELV) &
!    , VDIFF_E(LX1,LY1,LZ1,LELT) &
!     Solution data
    , VX     (LX1,LY1,LZ1,LELV) &
    , VY     (LX1,LY1,LZ1,LELV) &
    , VZ     (LX1,LY1,LZ1,LELV) &
    , T      (LX1,LY1,LZ1,LELT,LDIMT) &
    , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1) &
    , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1) &
    , BFX    (LX1,LY1,LZ1,LELV) &
    , BFY    (LX1,LY1,LZ1,LELV) &
    , BFZ    (LX1,LY1,LZ1,LELV) &
!    , cflf   (lx1,ly1,lz1,lelv) &
!    , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1)  & ! characteristics
!     Solution data for magnetic field
    , BX     (LBX1,LBY1,LBZ1,LBELV) &
    , BY     (LBX1,LBY1,LBZ1,LBELV) &
    , BZ     (LBX1,LBY1,LBZ1,LBELV) &
    , PM     (LBX2,LBY2,LBZ2,LBELV) &
    , BMX    (LBX1,LBY1,LBZ1,LBELV)   & ! Magnetic field RHS
    , BMY    (LBX1,LBY1,LBZ1,LBELV) &
    , BMZ    (LBX1,LBY1,LBZ1,LBELV) &
    , BBX1   (LBX1,LBY1,LBZ1,LBELV)  & ! Extrapolation terms for
    , BBY1   (LBX1,LBY1,LBZ1,LBELV)  & ! magnetic field rhs
    , BBZ1   (LBX1,LBY1,LBZ1,LBELV) &
    , BBX2   (LBX1,LBY1,LBZ1,LBELV) &
    , BBY2   (LBX1,LBY1,LBZ1,LBELV) &
    , BBZ2   (LBX1,LBY1,LBZ1,LBELV) &
    , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1) &
    , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1) &
    , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1) &
    , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2) )

    allocate(PR(LX2,LY2,LZ2,LELV))
!    allocate(PRLAG(LX2,LY2,LZ2,LELV,LORDER2)) 

    !> \todo Is qtl ever non-zero for .not. iflomach?
    allocate(QTL(LX2,LY2,LZ2,LELT))
!    allocate(USRDIV(LX2,LY2,LZ2,LELT))

    allocate(V1MASK (LX1,LY1,LZ1,LELV) &
    , V2MASK (LX1,LY1,LZ1,LELV) &
    , V3MASK (LX1,LY1,LZ1,LELV) &
    , PMASK  (LX1,LY1,LZ1,LELV) &
    , TMASK  (LX1,LY1,LZ1,LELT,LDIMT) &
!    , OMASK  (LX1,LY1,LZ1,LELT) &
    , VMULT  (LX1,LY1,LZ1,LELV) &
    , TMULT  (LX1,LY1,LZ1,LELT,LDIMT) &
    , B1MASK (LBX1,LBY1,LBZ1,LBELV)   & ! masks for mag. field
    , B2MASK (LBX1,LBY1,LBZ1,LBELV) &
    , B3MASK (LBX1,LBY1,LBZ1,LBELV) &
    , BPMASK (LBX1,LBY1,LBZ1,LBELV) )  ! magnetic pressure

!     Solution and data for perturbation fields
    allocate(VXP(LPX1*LPY1*LPZ1*LPELV,lpert) &
    , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert) &
    , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert) &
    , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert) &
    , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)   & ! perturbation field RHS
    , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert) &
    , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert) &
    , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert) &
    , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert) &
    , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert) &
    , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert)  & ! Extrapolation terms for
    , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert)  & ! perturbation field rhs
    , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert) &
    , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert) &
    ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert) &
    ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert) )

  end subroutine init_soln

end module soln
