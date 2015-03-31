module kinds
  implicit none
  integer, parameter :: SP = 4 !selected_real_kind(6,37)
  integer, parameter :: DP = 8 !selected_real_kind(15,307)
  integer, parameter :: QP = 16 !selected_real_kind(30,600)
  integer, parameter :: i8 = selected_int_kind(15)
  integer, parameter :: i4 = 4
  integer, parameter :: r4 = 4
  real(DP), parameter :: DP_eps = 1.d-15
end module kinds
