module dealias
  use kinds, only : DP
  implicit none
  real(DP), allocatable :: vxd(:,:,:,:)
  real(DP), allocatable :: vyd(:,:,:,:)
  real(DP), allocatable :: vzd(:,:,:,:)

  real(DP), allocatable :: imd1(:,:)
  real(DP), allocatable :: imd1t(:,:)
  real(DP), allocatable :: im1d(:,:)
  real(DP), allocatable :: im1dt(:,:)

  real(DP), allocatable :: pmd1(:,:)
  real(DP), allocatable :: pmd1t(:,:)

  contains

  subroutine init_dealias()
    use size_m
    implicit none
    allocate(vxd(lxd, lyd, lzd, lelv)) !verified
    allocate(vyd(lxd, lyd, lzd, lelv)) !verified
    allocate(vzd(lxd, lyd, lzd, lelv)) !verified

    allocate(imd1(lx1, lxd))
    allocate(imd1t(lxd, lx1))
    allocate(im1d(lxd, lx1))
    allocate(im1dt(lx1, lxd))

    allocate(pmd1(lx1,lxd))
    allocate(pmd1t(lxd,lx1))

  end subroutine init_dealias

end module dealias
