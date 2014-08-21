module wzf
! Points (z) and weights (w) on velocity, pressure

!     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
!     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp = nx-2)
  use kinds, only : DP
  use size_m
  implicit none

!     parameter (lxm = lx1)
  integer, parameter :: lxq = lx2

  real(DP) ::  zgl(lx1),wgl(lx1) &
    ,   zgp(lx1),wgp(lxq)

!     Tensor- (outer-) product of 1D weights   (for volumetric integration)
  real(DP) :: wgl1(lx1*lx1),wgl2(lxq*lxq) &
    ,  wgli(lx1*lx1)


!    Frequently used derivative matrices:

!    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
!    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)

!    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
!                   (currently the same as D1 and D1t...)


  real(DP) ::  d1    (lx1*lx1) , d1t    (lx1*lx1) &
    ,  d2    (lx1*lx1) , b2p    (lx1*lx1) &
    ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1) &
    ,  da    (lx1*lx1) , dat    (lx1*lx1) &
    ,  iggl  (lx1*lxq) , igglt  (lx1*lxq) &
    ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq) &
    ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)

end module wzf
