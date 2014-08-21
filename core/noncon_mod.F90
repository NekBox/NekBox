module noncon
  use kinds, only : DP
  implicit none
! 34567

  real(DP), allocatable :: umult(:), Jmat(:,:,:,:)
  real(DP), allocatable, dimension(:,:) :: xsp, ysp, zsp, xch, ych, zch
  real(DP), allocatable :: rtwid(:), stwid(:), dtrk(:), rs(:,:,:)

  integer, allocatable :: noncon_f(:), noncon_e(:), noncon_ip(:)
  integer, allocatable :: mortar(:,:), imin(:,:)
  integer :: mort_m

  logical :: ifnc, ifhalf
  logical, allocatable :: ifJt(:)

  contains

  subroutine init_noncon()
    use size_m
    implicit none

    allocate(umult(lx1*ly1*lz1*lelt) &
    ,  Jmat(lx1,lx1,2,maxmor) &
    ,  xsp(lx1,lz1),ysp(lx1,lz1),zsp(lx1,lz1) &
    ,  xch(lx1,lz1),ych(lx1,lz1),zch(lx1,lz1) &
    ,  rtwid(lx1),stwid(lx1) &
    ,  dtrk(ldim) &
    ,  rs( 2 , 2 , 2 ) )

    allocate(noncon_f(maxmor) &
    ,  noncon_e(maxmor) &
    ,  noncon_ip(maxmor) &
    ,  mortar(6,lelt) &
    ,  imin(3,2) )

    allocate(ifJt(maxmor))

  end subroutine init_noncon

end module noncon
