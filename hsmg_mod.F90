!> Module containing data for HSMG
!

!> Module containing data for HSMG
!
module hsmg
  use kinds, only : DP
  use size_m
  implicit none

!     Allocate MHD memory only if lbx1==lx1
  integer, parameter :: lmg_mhd=1-(lx1-lbx1)/(lx1-1) !1 if MHD is true, 0 otherwise

  integer, parameter :: lmgs=1 + lmg_mhd         ! max number of multigrid solvers
  integer, parameter :: lmgn=3                   ! max number of multigrid levels
  integer, parameter :: lmgx=lmgn+1              ! max number of mg index levels
  integer, parameter :: lxm=lx2+2,lym=lxm,lzm=lz2+2*(ldim-2) ! mgrid sizes
  integer, parameter :: lmg_rwt=2*lxm*lzm        ! restriction weight max size
  integer, parameter :: lmg_fasts=2*lxm*lxm      ! FDM S max size
  integer, parameter :: lmg_fastd=2*lxm*lym*lzm  ! FDM D max size
  integer, parameter :: lmg_swt=2*lxm*lzm        ! schwarz weight max size
  integer, parameter :: lmg_g=2*lx2*ly2*lz2      ! metrics max size
  integer, parameter :: lmg_solve=2*lxm*lym*lzm  ! solver r,e max size

  integer :: mg_lmax !>!< number of multigrid levels
  integer :: mg_fld !> @var active mg field
  integer, allocatable :: mg_nx(:), mg_ny(:), mg_nz(:) !level poly order
  integer, allocatable :: mg_nh(:), mg_nhz(:) !number of 1d nodes
  integer, allocatable :: mg_gsh_schwarz_handle(:,:) !dssum schwarz handles
  integer, allocatable :: mg_gsh_handle(:,:) !dssum handle
  integer, allocatable :: mg_rstr_wt_index(:,:), mg_mask_index(:,:)
  integer, allocatable :: mg_fast_s_index(:,:), mg_fast_d_index(:,:)
  integer, allocatable :: mg_solve_index(:,:)
  integer, allocatable :: mg_g_index(:,:), mg_schwarz_wt_index(:,:)

  real, allocatable :: mg_jh(:,:) !c-to-f interpolation matrices
  real, allocatable :: mg_jht(:,:) !transpose of mg_jh
  real, allocatable :: mg_jhfc(:,:) !c-to-f interpolation matrices
  real, allocatable :: mg_jhfct(:,:) !transpose of mg_jh

  real, allocatable :: mg_ah(:,:) !A hat matrices
  real, allocatable :: mg_bh(:,:) !B hat matrices
  real, allocatable :: mg_ch(:,:) !C hat matrices
  real, allocatable :: mg_dh(:,:) !D hat matrices
  real, allocatable :: mg_dht(:,:) !D hat transpose matrices
  real, allocatable :: mg_zh(:,:) !Nodal coordinates
  real, allocatable :: mg_rstr_wt(:) !restriction wt
  real, allocatable :: mg_mask(:) !b.c. mask (Max: might not be used)
  real, allocatable :: mg_fast_s(:), mg_fast_d(:)
  real, allocatable :: mg_schwarz_wt(:)
  real, allocatable :: mg_solve_e(:), mg_solve_r(:)
  real, allocatable, dimension(:) :: mg_h1,mg_h2,mg_b
  real, allocatable :: mg_g(:) ! metrics matrices

! must be able to hold two lower level extended schwarz arrays
  real, allocatable :: mg_work(:),mg_work2(:),mg_worke(:,:) 

  integer, allocatable :: mg_imask(:) ! For h1mg, mask is a ptr

!     Specific to h1 multigrid:
  integer :: mg_h1_lmax
  integer, allocatable :: mg_h1_n(:,:)
  integer, allocatable, dimension(:,:) :: p_mg_h1,p_mg_h2,p_mg_b,p_mg_g,p_mg_msk

  real(DP), allocatable, dimension(:) :: lr,ls,lt &
                          , llr,lls,llt &
                          , lmr,lms,lmt &
                          , lrr,lrs,lrt

  contains

  subroutine init_hsmg()
    use size_m
    implicit none
    allocate(mg_nx(lmgn), mg_ny(lmgn), mg_nz(lmgn))
    allocate(mg_nh(lmgn), mg_nhz(lmgn))
    allocate(mg_gsh_schwarz_handle(lmgn,lmgs), mg_gsh_handle(lmgn,lmgs))
    allocate(mg_rstr_wt_index(lmgx,0:lmgs), mg_mask_index(lmgx,0:lmgs))
    allocate(mg_solve_index(lmgx,0:lmgs))
    allocate(mg_fast_s_index(lmgx,0:lmgs), mg_fast_d_index(lmgx,0:lmgs))
    allocate(mg_schwarz_wt_index(lmgx,0:lmgs), mg_g_index(lmgx,0:lmgs))


    allocate(mg_jh(lxm*lxm,lmgn) &
    , mg_jht(lxm*lxm,lmgn)      & 
    , mg_jhfc (lxm*lxm,lmgn)    & 
    , mg_jhfct(lxm*lxm,lmgn)    &  ! verified
    , mg_ah(lxm*lxm,lmgn)       & 
    , mg_bh(lxm,lmgn)           & 
    , mg_dh(lxm*lxm,lmgn)       & 
    , mg_dht(lxm*lxm,lmgn)      & 
    , mg_zh(lxm,lmgn)           & 
    , mg_rstr_wt   (0:lmgs*lmg_rwt*2*ldim*lelt-1)  & !restriction wt
    , mg_mask      (0:lmgs*lmg_rwt*4*ldim*lelt-1)    & !b.c. mask
    , mg_fast_s    (0:lmgs*lmg_fasts*2*ldim*lelt-1) &
    , mg_fast_d    (0:lmgs*lmg_fastd*lelt-1) & ! verified
    , mg_schwarz_wt(0:lmgs*lmg_swt*4*ldim*lelt-1) & ! verified
!    , mg_solve_e   (0:lmg_solve*lelt-1) &
!    , mg_solve_r   (0:lmg_solve*lelt-1) &
    , mg_h1        (0:lmg_g*lelt-1) &
    , mg_h2        (0:lmg_g*lelt-1) &
    , mg_b         (0:lmg_g*lelt-1) & ! verified
    , mg_g         (0:lmg_g*((ldim-1)*3)*lelt-1)  & !metrics matrices (verified)

    , mg_work      (2*lxm*lym*lzm*lelt)  & ! verified
!    , mg_work2     (lxm*lym*lzm*lelt)    & ! two lower level extended
    , mg_worke     (lxm*lym*lzm,6) &      ! schwarz arrays
    )

    allocate(mg_imask(0:lmgs*lmg_rwt*4*ldim*lelt-1)) ! verified

!     Specific to h1 multigrid:
    allocate(mg_h1_n  (lmgx,ldimt1) &
    , p_mg_h1  (lmgx,ldimt1),p_mg_h2(lmgx,ldimt1) &
    , p_mg_b   (lmgx,ldimt1),p_mg_g (lmgx,ldimt1) &
    , p_mg_msk (lmgx,ldimt1) )

    allocate(lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4) &
  , llr(lelt),lls(lelt),llt(lelt) &
  , lmr(lelt),lms(lelt),lmt(lelt) &
  , lrr(lelt),lrs(lelt),lrt(lelt) )

  end subroutine init_hsmg

end module hsmg
