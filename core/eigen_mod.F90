
!     Eigenvalues
module eigen
  use kinds, only : DP
  implicit none   

  real(DP) :: EIGAA, EIGAS, EIGAST, EIGAE &
  ,EIGGA, EIGGS, EIGGST, EIGGE

  LOGICAL :: IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST

end module eigen
