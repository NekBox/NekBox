module fftw3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f'
  integer, parameter :: plan_kind = 8
end module fftw3
