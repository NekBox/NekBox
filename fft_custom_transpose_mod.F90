module fft
  use, intrinsic :: iso_c_binding
  use fftw3, only : FFTW_R2HC, FFTW_HC2R, FFTW_REDFT00, FFTW_REDFT10, FFTW_REDFT01
  implicit none

  public :: fft_r2r, transpose_grid
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

  integer, parameter :: max_transpose_plans = 100
  integer :: n_transpose_plans = 0
  type(C_PTR) :: transpose_plans(max_transpose_plans)
  integer :: transpose_comm(max_transpose_plans)
  integer :: transpose_num(max_transpose_plans)
  integer :: transpose_n0(max_transpose_plans)
  integer :: transpose_n1(max_transpose_plans)
  integer :: transpose_b0(max_transpose_plans)
  integer :: transpose_b1(max_transpose_plans)

contains

subroutine fft_r2r(u, length, num, kind, rescale)
  use kinds, only : DP
  use fftw3, only : FFTW_EXHAUSTIVE, FFTW_ESTIMATE
  use fftw3, only : fftw_plan_many_r2r
  use fftw3, only : fftw_execute_r2r
  use parallel, only : nid

  real(DP), intent(inout) :: u(:,:,:)
  integer, intent(in) :: length, num
  real(DP), intent(inout) :: rescale
  integer(C_INT) :: kind
  integer :: i, plan_idx
  real(DP), allocatable :: proxy(:)

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

  call fftw_execute_r2r(r2r_plans(plan_idx), u, u)

  if (kind == FFTW_REDFT00) then
    rescale = rescale * sqrt(2.*real(length-1, kind=DP))
  else if (kind == W_FORWARD .or. kind == W_BACKWARD) then
    rescale = rescale * sqrt(2.*real(length, kind=DP))
  else if (kind == FFTW_R2HC .or. kind == FFTW_HC2R) then
    rescale = rescale * sqrt(real(length, kind=DP))
  else
    write(*,*) "Don't know how to deal with FFTW kind", kind
  endif

end subroutine fft_r2r

subroutine transpose_grid(grid, grid_t, shape_x, idx, idx_t, comm)
  use kinds, only : DP
  use fftw3, only : FFTW_EXHAUSTIVE, FFTW_ESTIMATE
  use fftw3, only : fftw_mpi_plan_many_transpose
  use fftw3, only : fftw_mpi_execute_r2r
  use parallel, only : nid, nekreal

  real(DP), intent(inout) :: grid(0:,0:,0:)
  real(DP), intent(out) :: grid_t(0:,0:,0:)
  integer,  intent(in) :: shape_x(3)
  integer, intent(in) :: idx
  integer, intent(in) :: idx_t
  integer, intent(in) :: comm

  real(DP), allocatable :: tmp(:,:)
  real(DP), allocatable :: tmp_t(:,:)
  integer :: shape_c(3), block0, block1, n0, n1, num
  integer, parameter :: one = 1
  type(C_PTR) :: transpose_plan
  integer :: i, j, k, plan_idx, ierr

  shape_c = shape_x

  n0 = shape_c(idx_t)
  n1 = shape_c(idx)

  if (size(grid) < 1) return
!  if (nid == 0) write(*,*) shape_x, idx, idx_t
 
  if (idx == 1 .or. idx_t == 1) then
    block0 = size(grid,2)
    block1 = size(grid_t,2)
    num    = size(grid,3) 
    allocate(tmp(0:block0-1,   0:size(grid,1)-1))
    allocate(tmp_t(0:block1-1, 0:size(grid_t,1)-1))
  else if (idx == 3 .or. idx_t == 3) then
    block0 = size(grid,3)
    block1 = size(grid_t,3)
    num    = size(grid,2) 
    allocate(tmp(0:block0-1,   0:size(grid,1)-1))
    allocate(tmp_t(0:block1-1, 0:size(grid_t,1)-1))
  endif

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

  !call nekgsync()
  !if (nid == 0) write(*,*) "executed"

end subroutine transpose_grid

end module
