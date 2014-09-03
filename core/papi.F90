subroutine nek_flops(flops,mflops)
  implicit none
  real*4 :: rtime,ptime,mflops
  integer*8 :: flops

  call getflops_papi(flops,mflops)

  return
end subroutine nek_flops

subroutine getflops_papi(flops,mflops)
  implicit none
  real*4 :: rtime,ptime,mflops
  integer*8 :: flops
#ifdef PAPI
  include 'f77papi.h'

  call papif_flops(rtime,ptime,flops,mflops,ierr)
  if(ierr > 0) then
      flops = -1
      mflops = -1
  endif
#endif
     
  return
end subroutine getflops_papi
