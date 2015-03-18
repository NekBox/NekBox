!> \file mxm_wrapper.F90
!! \brief Matrix-multiply wrappers

!> \brief Compute matrix-matrix product C = A*B
!!
!! A, B, and C must be contiguously packed
subroutine mxm(a,n1,b,n2,c,n3)
  use kinds, only : DP, i8
  implicit none
    
  integer, intent(in) :: n1, n2, n3
  real(DP), intent(in) :: a(n1,n2),b(n2,n3)
  real(DP), intent(out) :: c(n1,n3)
  integer :: aligned
  integer :: K10_mxm

#ifdef BGQ
    integer(i8) :: tt = 32

    if (n2 == 8 .and. mod(n1,4) == 0 &
!        .and. MOD(LOC(a),tt)==0 & 
!        .and. MOD(LOC(b),tt)==0 & 
!        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_8(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 12 .and. mod(n1,4) == 0 &
!        .and. MOD(LOC(a),tt)==0 & 
!        .and. MOD(LOC(b),tt)==0 & 
!        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_12(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 16 .and. mod(n1,4) == 0 &
!        .and. MOD(LOC(a),tt)==0 & 
!        .and. MOD(LOC(b),tt)==0 & 
!        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_16(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 10 .and. mod(n1,4) == 0 .and. mod(n3,2) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_10(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 6 .and. mod(n1,4) == 0 .and. mod(n3,2) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_6(a, n1, b, n2, c, n3)  
     return
   endif
#endif

!#define BLAS_MXM
#ifdef BLAS_MXM
    call dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
    return
#endif
     
#ifdef BG
  call bg_aligned3(a,b,c,aligned)
  if (n2 == 2) then
      call mxm44_2(a,n1,b,n2,c,n3)
  else if ((aligned == 1) .AND. &
      (n1 >= 8) .AND. (n2 >= 8) .AND. (n3 >= 8) .AND. &
      (modulo(n1,2) == 0) .AND. (modulo(n2,2) == 0) ) then
      if (modulo(n3,4) == 0) then
          call bg_mxm44(a,n1,b,n2,c,n3)
      else
          call bg_mxm44_uneven(a,n1,b,n2,c,n3)
      endif
  else if((aligned == 1) .AND. &
      (modulo(n1,6) == 0) .AND. (modulo(n3,6) == 0) .AND. &
      (n2 >= 4) .AND. (modulo(n2,2) == 0) ) then
      call bg_mxm3(a,n1,b,n2,c,n3)
  else
      call mxm44_0(a,n1,b,n2,c,n3)
  endif
  return
#endif

#ifdef K10_MXM
  ! fow now only supported for lx1=8
  ! tuned for AMD K10
  ierr = K10_mxm(a,n1,b,n2,c,n3)
  if (ierr > 0) call mxmf2(a,n1,b,n2,c,n3)
  return
#endif

  call mxmf2(a,n1,b,n2,c,n3)

  return
end subroutine mxm
