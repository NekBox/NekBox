!> \file mxm_wrapper.F90
!! \brief Matrix-multiply wrappers

!> \brief Compute matrix-matrix product C = A*B
!!
!! A, B, and C must be contiguously packed
!#define USE_LIBXSMM


#ifdef USE_LIBXSMM
!#define XSMM_DIRECT
#define XSMM_DISPATCH
#include "libxsmm.f90"
#endif

subroutine mxm(a,n1,b,n2,c,n3)
  use kinds, only : DP, i8
#ifdef USE_LIBXSMM
#ifdef XSMM_DIRECT
  use iso_c_binding
  use libxsmm, only : libxsmm_mm
#else
  use libxsmm
#endif
#endif
  implicit none
    
  integer, intent(in) :: n1, n2, n3
  real(DP), intent(in) :: a(n1,n2),b(n2,n3)
  real(DP), intent(out) :: c(n1,n3)
  integer :: aligned
  integer :: K10_mxm

#ifdef XSMM_DISPATCH
  INTEGER, PARAMETER :: T = LIBXSMM_DOUBLE_PRECISION
  PROCEDURE(LIBXSMM_DMM_FUNCTION), POINTER :: dmm
  TYPE(C_FUNPTR) :: f
#endif

#ifdef XSMM_DIRECT
  interface libxsmm_dmm_8_8_8
    subroutine libxsmm_dmm_8_8_8(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_64_8_8
    subroutine libxsmm_dmm_64_8_8(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_8_64_8
    subroutine libxsmm_dmm_8_64_8(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_4_4_4
    subroutine libxsmm_dmm_4_4_4(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_4_16_4
    subroutine libxsmm_dmm_4_16_4(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_16_4_4
    subroutine libxsmm_dmm_16_4_4(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_10_10_10
    subroutine libxsmm_dmm_10_10_10(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_10_100_10
    subroutine libxsmm_dmm_10_100_10(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_100_10_10
    subroutine libxsmm_dmm_100_10_10(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface

  interface libxsmm_dmm_144_12_8
    subroutine libxsmm_dmm_144_12_8(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface



  interface libxsmm_dmm_12_12_12
    subroutine libxsmm_dmm_12_12_12(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_144_12_12
    subroutine libxsmm_dmm_144_12_12(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
  interface libxsmm_dmm_12_144_12
    subroutine libxsmm_dmm_12_144_12(a,b,c) BIND(C)
      use iso_c_binding, only : c_ptr
      type(c_ptr), value :: a, b, c
    end subroutine
  end interface
#endif

#ifdef BGQ
    integer(i8) :: tt = 32

    if (n2 == 8 .and. mod(n1,4) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_8(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 12 .and. mod(n1,4) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
       ) then
     call mxm_bgq_12(a, n1, b, n2, c, n3)  
     return
   endif
   if (n2 == 16 .and. mod(n1,4) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
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

#ifdef USE_LIBXSMM


#ifdef XSMM_DISPATCH
  f = libxsmm_dispatch(1.0_dp, 0.0_dp, n1, n3, n2)
  if (C_ASSOCIATED(f)) then
      CALL C_F_PROCPOINTER(f, dmm)
      CALL dmm(1.0_dp, 0.0_dp, a, b, c)
      return
  endif
#endif

#ifdef XSMM_DIRECT
    if (n2 == 8) then
      if (n1 == 8 .and. n3 == 8) then
        call libxsmm_dmm_8_8_8(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 64 .and. n3 == 8) then
        call libxsmm_dmm_64_8_8(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 8 .and. n3 == 64) then
        call libxsmm_dmm_8_64_8(c_loc(a),c_loc(b),c_loc(c))
        return
      endif
      else if (n1 == 144 .and. n3 == 12) then
        call libxsmm_dmm_144_12_8(c_loc(a),c_loc(b),c_loc(c))
        return
      endif
    else if (n2 == 12) then
      if (n1 == 12 .and. n3 == 12) then
        call libxsmm_dmm_12_12_12(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 144 .and. n3 == 12) then
        call libxsmm_dmm_144_12_12(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 12 .and. n3 == 144) then
        call libxsmm_dmm_12_144_12(c_loc(a),c_loc(b),c_loc(c))
        return
      endif
    else if (n2 == 10) then
      if (n1 == 10 .and. n3 == 10) then
        call libxsmm_dmm_10_10_10(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 100 .and. n3 == 10) then
        call libxsmm_dmm_100_10_10(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 10 .and. n3 == 100) then
        call libxsmm_dmm_10_100_10(c_loc(a),c_loc(b),c_loc(c))
        return
      endif
    else if (n2 == 4) then
      if (n1 == 4 .and. n3 == 4) then
        call libxsmm_dmm_4_4_4(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 16 .and. n3 == 4) then
        call libxsmm_dmm_16_4_4(c_loc(a),c_loc(b),c_loc(c))
        return
      else if (n1 == 4 .and. n3 == 16) then
        call libxsmm_dmm_4_16_4(c_loc(a),c_loc(b),c_loc(c))
        return
      endif
    endif
#endif
#ifndef XSMM_DISPATCH
    call libxsmm_mm(n1, n3, n2, a, b, c)
    return
#endif
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
