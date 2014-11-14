module fft
  use, intrinsic :: iso_c_binding
  implicit none

  public fft_r2r
  private
  integer, parameter :: nplans = 10
  integer :: n_r2r_plans = 0
  type(C_PTR) :: r2r_plans(nplans)
  integer :: r2r_plan_lengths(nplans)
  integer :: r2r_plan_nums(nplans)
  integer(C_INT) :: r2r_plan_kinds(nplans)

contains

subroutine fft_r2r(u, length, num, kind, rescale)
  use kinds, only : DP
  use fftw3, only : FFTW_EXHAUSTIVE, FFTW_ESTIMATE
  use fftw3, only : FFTW_R2HC, FFTW_HC2R, FFTW_REDFT00
  use fftw3, only : fftw_plan_many_r2r
  use fftw3, only : fftw_execute_r2r

  real(DP), intent(inout) :: u(*)
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
    write(*,*) "Made a plan", length, num, kind
  endif
  call fftw_execute_r2r(r2r_plans(plan_idx), u, u)

  if (kind == FFTW_REDFT00) then
    rescale = rescale * sqrt(2.*real(length-1, kind=DP))
  else if (kind == FFTW_R2HC .or. kind == FFTW_HC2R) then
    rescale = rescale * sqrt(real(length, kind=DP))
  else
    write(*,*) "Don't know how to deal with FFTW kind", kind
  endif

end subroutine fft_r2r

end module
