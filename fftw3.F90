module fftw3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  integer, parameter :: plan_kind = 8
  !integer, external :: fftw_mpi_local_size_3d_transposed
end module fftw3
