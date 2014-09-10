    subroutine mxm(a,n1,b,n2,c,n3)

!     Compute matrix-matrix product C = A*B
!     for contiguously packed matrices A,B, and C.


    use size_m
    use opctr
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use geom
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    real :: a(n1,n2),b(n2,n3),c(n1,n3)
    integer :: aligned
    integer :: K10_mxm

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'mxm   '
    endif
    isbcnt = n1*n3*(2*n2-1)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! define BLAS_MXM
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
