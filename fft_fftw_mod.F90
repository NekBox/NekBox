!> \file fft_fftw_mod.F90
!! \brief Parallel FFT
!! \date December 2015
!! \author Max Hutchinson
!!
!! This module provides a 3D parallel FFT through 1D FFTs and parallel
!! transposes.  It also provides access to the wave-number of each mode,
!! given its index and the type of transform used for that axis.
!!
!! This implementation uses FFTW for the 1D FFT and DCTs and 
!! MPI_AllToAll for the transposes.
module fft
  use, intrinsic :: iso_c_binding
  use fftw3, only : FFTW_R2HC, FFTW_HC2R, FFTW_REDFT00, FFTW_REDFT10, FFTW_REDFT01
  implicit none

  public :: fft_r2r, transpose_grid, wavenumber
  public :: W_FORWARD, W_BACKWARD, P_FORWARD, P_BACKWARD 

  private

  integer, parameter :: nplans = 10
  integer :: n_r2r_plans = 0
  type(C_PTR) :: r2r_plans(nplans)
  integer :: r2r_plan_lengths(nplans)
  integer :: r2r_plan_nums(nplans)
  integer(C_INT) :: r2r_plan_kinds(nplans)

  integer, parameter :: W_FORWARD = FFTW_REDFT10
  integer, parameter :: W_BACKWARD = FFTW_REDFT01
  integer, parameter :: P_FORWARD = FFTW_R2HC
  integer, parameter :: P_BACKWARD = FFTW_HC2R

contains

!> \brief Compute a number of 1D FFTs or DCTs 
!!
!! The transforms are performed inplace.  Rescale
!! is rescaled by the normalization factor of the
!! transform.
subroutine fft_r2r(u, length, num, kind, rescale)
  use kinds, only : DP
  use fftw3, only : FFTW_EXHAUSTIVE, FFTW_ESTIMATE
  use fftw3, only : fftw_plan_many_r2r
  use fftw3, only : fftw_execute_r2r
  use parallel, only : nid

  real(DP), intent(inout) :: u(:,:,:) !>!< data, leading dim is transformed
  integer,  intent(in)    :: length  !>!< length of leading dimension
  integer,  intent(in)    :: num     !>!< number of transforms
  integer(C_INT)          :: kind    !>!< Type of transform.  See pulic parameters
  real(DP), intent(inout) :: rescale !>!< Normalization factor

  integer :: i, plan_idx
  real(DP), allocatable :: proxy(:)

  ! Look for an existing plan
  plan_idx = -1
  do i = 1, n_r2r_plans
    if ( &
      length == r2r_plan_lengths(i) .and. &
      num    == r2r_plan_nums(i)    .and. &
      kind   == r2r_plan_kinds(i) &
       ) then
      plan_idx = i
      exit
    endif
  enddo

  ! If there's no plan, make a new one
  if (plan_idx < 0) then
    if (n_r2r_plans == nplans) then
      write(*,*) "Ran out of plans!"
    endif

    n_r2r_plans = n_r2r_plans + 1
    plan_idx = n_r2r_plans
    allocate(proxy(num*length))
    r2r_plans(plan_idx) = fftw_plan_many_r2r(1, &
                            (/length/), num, &
                            proxy, (/length/), 1, length, &
                            proxy, (/length/), 1, length, &
                            (/kind/), FFTW_EXHAUSTIVE) 
    deallocate(proxy)
    r2r_plan_lengths(plan_idx) = length
    r2r_plan_nums(plan_idx)    = num
    r2r_plan_kinds(plan_idx)   = kind
  endif

  ! execute the plan
  call fftw_execute_r2r(r2r_plans(plan_idx), u, u)

  ! rescale the normalization factor
  if (kind == FFTW_REDFT00) then
    rescale = rescale * sqrt(2.*real(length-1, kind=DP))
  else if (kind == W_FORWARD .or. kind == W_BACKWARD) then
    rescale = rescale * sqrt(2.*real(length, kind=DP))
  else if (kind == FFTW_R2HC .or. kind == FFTW_HC2R) then
    rescale = rescale * sqrt(real(length, kind=DP))
  else
    write(*,*) "Don't know how to deal with FFT kind", kind
  endif

end subroutine fft_r2r

!> \brief Get the wavenumber of the ith mode
real(DP) function wavenumber(i, N, L, kind)
  use kinds, only : DP
  implicit none
  integer,  intent(in) :: i    !>!< index of mode
  integer,  intent(in) :: N    !>!< number of modes in transform
  real(DP), intent(in) :: L    !>!< Length of domain in transform
  integer,  intent(in) :: kind !>!< Type of transform

  real(DP), parameter :: pi = 4.*atan(1.)

  if (kind == P_FORWARD) then
    if (i <= N / 2) then
      wavenumber = 2*pi*i/(L)
    else
      wavenumber = 2*pi*(N - i)/(L)
    endif
  else if (kind == W_FORWARD) then
    wavenumber = pi*i/(L)
  else
    write(*,*) "Don't know how to deal with FFT kind", kind
  endif

end function wavenumber

!> \brief Global (parallel) transform of grid data
subroutine transpose_grid(grid, grid_t, shape_x, idx, idx_t, comm)
  use kinds, only : DP
  use fftw3, only : FFTW_EXHAUSTIVE, FFTW_ESTIMATE
  use fftw3, only : fftw_mpi_plan_many_transpose
  use fftw3, only : fftw_mpi_execute_r2r
  use parallel, only : nid, nekreal

  real(DP), intent(inout) :: grid(0:,0:,0:)
  real(DP), intent(out)   :: grid_t(0:,0:,0:)
  integer,  intent(in)    :: shape_x(3)
  integer, intent(in)     :: idx
  integer, intent(in)     :: idx_t
  integer, intent(in)     :: comm

  real(DP), allocatable :: tmp(:,:)
  real(DP), allocatable :: tmp_t(:,:)
  integer :: block0, block1, num
  integer :: i, j, k, ierr, comm_size


  if (size(grid) < 1) return

  call MPI_Comm_size(comm, comm_size, ierr)
 
  if (idx == 1 .or. idx_t == 1) then
    block0 = size(grid,2)
    block1 = size(grid,1) / comm_size
    num    = size(grid,3) 
  else if (idx == 3 .or. idx_t == 3) then
    block0 = size(grid,3)
    block1 = size(grid,1) / comm_size
    num    = size(grid,2) 
  endif
  allocate(tmp(0:block0-1,   0:size(grid,1)-1))
  allocate(tmp_t(0:block0-1, 0:size(grid,1)-1))

  if (idx == 1 .or. idx_t == 1) then
    do i = 0, num - 1
      tmp = transpose(grid(:,:,i))
      call mpi_alltoall(tmp,   block0*block1, nekreal, &
                        tmp_t, block0*block1, nekreal, comm, ierr)
      if (ierr /= 0) write(*,*) "alltoall errored", ierr, nid
      do j = 0, size(grid_t,1) - 1
        do k = 0, size(grid_t,2) - 1
          grid_t(j,k,i) = tmp_t( mod(j,int(block0)), (j / block0) * block1 + k )
        enddo
      enddo
    enddo
  else if (idx == 3 .or. idx_t == 3) then
    do i = 0, num - 1
      tmp = transpose(grid(:,i,:))
      call mpi_alltoall(tmp,   block0*block1, nekreal, &
                        tmp_t, block0*block1, nekreal, comm, ierr)
      if (ierr /= 0) write(*,*) "alltoall errored", ierr, nid
      do j = 0, size(grid_t,1) - 1
        do k = 0, size(grid_t,3) - 1
          grid_t(j,i,k) = tmp_t( mod(j,int(block0)), (j / block0) * block1 + k )
        enddo
      enddo
    enddo
  else
    write(*,*) "Something went wrong in transpose", nid
  endif

end subroutine transpose_grid

end module
