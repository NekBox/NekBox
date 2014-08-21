module scratch
  use kinds, only : DP
  USE, INTRINSIC :: ISO_C_BINDING
!     Scratch arrays used in CONNECT and associated subroutines.

!     ring pass arrays

  real(DP), allocatable :: rmxs(:), rmax(:)
  real(DP), allocatable :: xcg(:), ycg(:), zcg(:)
  real(DP), allocatable :: xgs(:), ygs(:), zgs(:)
  real(DP), allocatable :: xml(:,:,:,:), xms(:,:,:,:)
  real(DP), allocatable :: yml(:,:,:,:), yms(:,:,:,:)
  real(DP), allocatable :: zml(:,:,:,:), zms(:,:,:,:)

  real(DP), allocatable :: side(:,:,:), sides(:,:,:)
  real(DP), allocatable :: flag(:,:,:,:)
  real(DP), allocatable :: tmp(:,:,:,:), tmp2(:,:,:,:)
  real(DP), allocatable :: lmult(:,:,:,:)
  real(DP), allocatable :: bcs(:,:,:), xyz(:,:,:)

!     nested dissection arrays

  character(3), allocatable :: cbcs(:,:)
  integer, allocatable :: ibrnch(:), nbrnch(:), list(:), list1(:), list2(:)
  logical, allocatable :: ifcnst(:,:)

  real(DP), allocatable :: xyzl(:,:,:), cg(:,:)

!max    dimension xyzl(3,8,lelt),cg(3,lelt)
!max    equivalence (xyzl,xms)
!max    equivalence (cg,xgs)

  contains

  subroutine init_scratch()
    use size_m
    implicit none

    allocate( rmxs(lelt),rmax(lelt) &
    ,xcg(lelt),ycg(lelt),zcg(lelt) &
    ,xgs(lelt),ygs(lelt),zgs(lelt) &
    ,xml(3,3,lzl,lelt),xms(3,3,lzl,lelt) &
    ,yml(3,3,lzl,lelt),yms(3,3,lzl,lelt) &
    ,zml(3,3,lzl,lelt),zms(3,3,lzl,lelt) )
    allocate(side(4,6,lelt),sides(4,6,lelt))
    allocate(flag(3,3,lzl,lelt),tmp2(3,3,lzl,lelt) &
    ,lmult(3,3,lzl,lelt),bcs(5,6,lelt) &
    ,xyz(3,8,lelt)) 

!     nested dissection arrays

    allocate(cbcs(6,lelt))
    allocate(ibrnch(lelt),nbrnch(lelt) &
    ,list(lelt),list1(lelt) ,list2(lelt) &
    ,ifcnst(6,lelt) )

!max    call C_F_POINTER (C_LOC(xms), xyzl, [3,8,lelt])
!max    call C_F_POINTER (C_LOC(xgs), cg, [3,lelt])
    allocate(xyzl(3,8,lelt), cg(3,lelt))

  end subroutine init_scratch

end module scratch
