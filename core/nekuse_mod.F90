!> cleaned
module nekuse
  use kinds, only : DP
  use size_m
  implicit none

  real(DP) :: x, y, z, r, theta
  real(DP) :: ux, uy, uz, un, u1, u2
  real(DP) :: trx, try, trz, trn, tr1, tr2
  real(DP) :: PA
  real(DP) :: ffx, ffy, ffz
  real(DP) :: temp, flux, hc, hrad, tinf, qvol
  real(DP) :: udiff, utrans
  real(DP) :: si2, si3, sigma
  real(DP) :: turbk, turbe
  real(DP) :: ps(ldimt)

  character(len=3) :: cbu

end module nekuse
