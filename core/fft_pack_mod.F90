!> \file fft_pack_mod.F90
!! \brief fftpack-based implementation of the fft module
!!
!! \note This is using a DCT-I, when it should be using
!! a DCT-II.  It doesn't seem to make a huge difference 
!! (an extra gmres iteration here and there), but it is
!! wrong, so that's annoying.

module fft
  use, intrinsic :: iso_c_binding
  use kinds, only : DP 
  implicit none

  public :: fft_r2r, transpose_grid, wavenumber
  public :: W_FORWARD, W_BACKWARD, P_FORWARD, P_BACKWARD 

  private

  integer, parameter :: nplans = 10
  integer, parameter :: lensav = 256
  integer :: n_r2r_plans = 0
  real(DP) :: r2r_plans(lensav, nplans)
  integer :: r2r_plan_lengths(nplans)
  integer :: r2r_plan_kinds(nplans)

  integer, parameter :: W_FORWARD = 0
  integer, parameter :: W_BACKWARD = 1 
  integer, parameter :: P_FORWARD = 2 
  integer, parameter :: P_BACKWARD = 3 

contains

subroutine fft_r2r(u, length, num, kind, rescale)
  use kinds, only : DP
  use parallel, only : nid

  real(DP), intent(inout) :: u(length,*)
  integer, intent(in) :: length, num
  real(DP), intent(inout) :: rescale
  integer(C_INT) :: kind
  integer :: i, plan_idx, ierr
  real(DP), allocatable :: work(:)

  plan_idx = -1
  do i = 1, n_r2r_plans
    if ( &
      length == r2r_plan_lengths(i) .and. &
      kind/2   == r2r_plan_kinds(i) &
       ) then
      plan_idx = i
      exit
    endif
  enddo

  if (plan_idx < 0) then
    if (n_r2r_plans == nplans) then
      write(*,*) "Ran out of plans!"
    endif

    n_r2r_plans = n_r2r_plans + 1
    plan_idx = n_r2r_plans

    if (kind/2 == 0) then
      call dcost1i( length, r2r_plans(:,plan_idx), lensav, ierr )
    else if (kind/2 == 1) then
      call dfft1i( length, r2r_plans(:,plan_idx), lensav, ierr )
    endif

    r2r_plan_lengths(plan_idx) = length
    r2r_plan_kinds(plan_idx)   = kind / 2
  endif

  if (kind == W_FORWARD) then
    allocate(work(length))
    do i = 1, num
      call dcost1f(length, 1, u(:,i), length, r2r_plans(:,plan_idx), lensav, work, length, ierr)
    enddo
  else if (kind == W_BACKWARD) then
    allocate(work(length))
    do i = 1, num
      call dcost1b(length, 1, u(:,i), length, r2r_plans(:,plan_idx), lensav, work, length, ierr)
    enddo
  else if (kind == P_FORWARD) then
    allocate(work(length))
    do i = 1, num
      call dfft1f(length, 1, u(:,i), length, r2r_plans(:,plan_idx), lensav, work, length, ierr)
    enddo
  else if (kind == P_BACKWARD) then
    allocate(work(length))
    do i = 1, num
      call dfft1b(length, 1, u(:,i), length, r2r_plans(:,plan_idx), lensav, work, length, ierr)
    enddo
  endif

end subroutine fft_r2r

real(DP) function wavenumber(i, N, L, kind)
  use kinds, only : DP
  implicit none
  integer,  intent(in) :: i, N
  real(DP), intent(in) :: L
  integer,  intent(in) :: kind

  real(DP), parameter :: pi = 4.*atan(1.)

  if (kind == P_FORWARD) then
#if 0
    if (i <= N / 2) then
      wavenumber = 2*pi*i/(L)
    else
      wavenumber = 2*pi*(N - i)/(L)
    endif
#else
    wavenumber = 2*pi*i/(L)
#endif
  else if (kind == W_FORWARD) then
    wavenumber = pi*i*N/L
  else
    write(*,*) "Don't know how to deal with FFT kind", kind
  endif

end function wavenumber

subroutine transpose_grid(grid, grid_t, shape_x, idx, idx_t, comm)
  use kinds, only : DP
  use parallel, only : nid, nekreal

  real(DP), intent(inout) :: grid(0:,0:,0:)
  real(DP), intent(out) :: grid_t(0:,0:,0:)
  integer,  intent(in) :: shape_x(3)
  integer, intent(in) :: idx
  integer, intent(in) :: idx_t
  integer, intent(in) :: comm

  real(DP), allocatable :: tmp(:,:)
  real(DP), allocatable :: tmp_t(:,:)
  integer :: block0, block1, num
  integer :: i, j, k, ierr

  if (size(grid) < 1) return
 
  if (idx == 1 .or. idx_t == 1) then
    block0 = size(grid,2)
    block1 = size(grid_t,2)
    num    = size(grid,3) 
  else if (idx == 3 .or. idx_t == 3) then
    block0 = size(grid,3)
    block1 = size(grid_t,3)
    num    = size(grid,2) 
  endif
  allocate(tmp(0:block0-1,   0:size(grid,1)-1))
  allocate(tmp_t(0:block1-1, 0:size(grid_t,1)-1))

  if (idx == 1 .or. idx_t == 1) then
    do i = 0, num - 1
      tmp = transpose(grid(:,:,i))
      call mpi_alltoall(tmp,   block0*block1, nekreal, &
                        tmp_t, block0*block1, nekreal, comm, ierr)
      if (ierr /= 0) write(*,*) "alltoall errored", ierr, nid
      do j = 0, size(grid_t,1) - 1
        do k = 0, size(grid_t,2) - 1
          grid_t(j,k,i) = tmp_t( mod(j,int(block1)), (j / block1) * block1 + k )
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
        do k = 0, size(grid_t,2) - 1
          grid_t(j,i,k) = tmp_t( mod(j,int(block1)), (j / block1) * block1 + k )
        enddo
      enddo
    enddo
  else
    write(*,*) "Something went wrong in transpose", nid
  endif

end subroutine transpose_grid

end module
