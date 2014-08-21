module semhat
  use kinds, only : DP
  use size_m

  implicit none
  integer, parameter :: lr2 = 2*lx1*lx1 
  integer, parameter :: l3 = lx1*(lx1+1)*(lx1+2)/3
  integer, parameter :: l2 =lx1*(lx1+1)/2

  real(DP), allocatable, dimension(:) :: ah,  bh,  ch,  dh
  real(DP), allocatable, dimension(:) :: dph, jph, zh,  wh
  real(DP), allocatable, dimension(:) :: bgl, zgl, dgl, jgl
  real(DP), allocatable, dimension(:) :: dd, zp, ww

  contains

  subroutine init_semhat()
    implicit none
    
    allocate(ah(lr2), bh(lr2), ch(lr2), dh(lr2))
    allocate(dph(lr2), jph(lr2))
    allocate(zh(lr2),wh (lr2), bgl(lr2), zgl(lr2), dgl(lr2), jgl(lr2))
    allocate(dd(l3), zp(l2), ww(l2))

  end subroutine init_semhat
  
end module semhat
