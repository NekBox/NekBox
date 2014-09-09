module kinds
  implicit none
  integer, parameter :: DP = selected_real_kind(14,200)
  integer, parameter :: i8 = 8
  integer, parameter :: i4 = 4
  integer, parameter :: r4 = 4
  real(DP), parameter :: DP_eps = 1.d-15
end module kinds
