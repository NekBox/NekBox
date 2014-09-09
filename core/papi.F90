subroutine nek_flops(flops,mflops)
  use kinds, only : r4, i8
  implicit none
  real(r4) :: mflops
  integer(i8) :: flops

  call getflops_papi(flops,mflops)

  return
end subroutine nek_flops

subroutine getflops_papi(flops,mflops)
  use kinds, only : r4, i8
  implicit none
  real(r4) :: mflops
  integer(i8) :: flops
#ifdef PAPI
  real(r4) :: rtime,ptime
  include 'f77papi.h'

  call papif_flops(rtime,ptime,flops,mflops,ierr)
  if(ierr > 0) then
      flops = -1
      mflops = -1
  endif
#endif
     
  return
end subroutine getflops_papi
