module zper
!     Eigenvalue arrays and pointers for Global Tensor Product
!     parameter (lelg_sm=2)
!     parameter (ltfdm2 =2)
!     parameter (lelg_sm=lelg)
!     parameter (ltfdm2=2*lx2*ly2*lz2*lelt)
!     parameter (leig2=2*lx2*lx2*(lelx*lelx+lely*lely+lelz*lelz))
!     parameter (leig =2*lx2*(lelx+lely+lelz))
  use kinds, only : DP
  use size_m
  implicit none

  integer, parameter :: lfdm0 = 1-lfdm
  integer, parameter :: lelg_sm=lfdm0+lfdm*lelg
  integer, parameter :: ltfdm2 =lfdm0+lfdm*2*lx2*ly2*lz2*lelt
  integer, parameter :: leig2  =lfdm0+lfdm*2*lx2*lx2 &
    *(lelx*lelx+lely*lely+lelz*lelz)
  integer, parameter :: leig   =lfdm0+lfdm*2*lx2*(lelx+lely+lelz)

  integer :: neigx, neigy, neigz, pvalx ,pvaly ,pvalz &
    , pvecx ,pvecy ,pvecz

  real(DP), allocatable :: sp(:), spt(:), eigp(:), wavep(:)
  real(DP) :: msp(3,2),mlp(3,2)


!     Logical, array and geometry data for tensor-product box
  logical ::          ifycrv,ifzper,ifgfdm,ifgtp,ifemat

  integer :: nelx,nely,nelz,nelxy &
    ,lex2pst(3),pst2lex(3) &
    ,ngfdm_p(3),ngfdm_v(3,2)

!     Complete exchange arrays for pressure
!     common /gfdmcx/  part_in(0:lp),part_out(0:lp)
  integer, parameter :: lp_small = 256
  integer :: part_in(0:lp_small),part_out(0:lp_small), msg_id(0:lp_small,2), mcex

!     Permutation arrays for gfdm pressure solve
  integer, allocatable :: tpn1(:), tpn2(:), tpn3(:), ind23(:)

  integer, parameter :: lfdx =lfdm0+lfdm*lx2*lelx
  integer, parameter :: lfdy =lfdm0+lfdm*ly2*lely
  integer, parameter :: lfdz =lfdm0+lfdm*lz2*lelz
  real(DP), allocatable :: xgtp(:), ygtp(:), zgtp(:)
  real(DP), allocatable :: xmlt(:), ymlt(:), zmlt(:)

!     Metrics for 2D x tensor-product solver
  real(DP), allocatable, dimension(:,:,:) :: rx2, ry2, sx2, sy2, w2d, bxyi

  character(3) ::      gtp_cbc(6,0:ldimt1+1)

  contains

  subroutine init_zper
    use size_m
    implicit none

    allocate(sp(leig2),spt(leig2),eigp(leig), wavep(ltfdm2))

    allocate(tpn1(ltfdm2),tpn2(ltfdm2), tpn3(ltfdm2),ind23(ltfdm2))

    allocate(xgtp(0:lelx),ygtp(0:lely),zgtp(0:lelz))
    allocate(xmlt(lfdx),ymlt(lfdy),zmlt(lfdz))


!     Metrics for 2D x tensor-product solver
    allocate(rx2(lx2,ly2,lelv) &
    , ry2  (lx2,ly2,lelv) &
    , sx2  (lx2,ly2,lelv) &
    , sy2  (lx2,ly2,lelv) &
    , w2d  (lx2,ly2,lelv) &
    , bxyi (lx1,ly1,lelv) )

  end subroutine init_zper

end module zper
