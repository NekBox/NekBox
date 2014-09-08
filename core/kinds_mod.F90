module kinds
  implicit none
  integer, parameter :: DP = selected_real_kind(14,200)
  real(DP), parameter :: DP_eps = 1.d-15
end module kinds
