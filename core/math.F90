!> \brief local inner product, with weight
real(DP) FUNCTION VLSC3(X,Y,B,N)
  use kinds, only : DP
  use opctr, only : isclld, nrout, myrout, rname, dct, ncall, dcount
  implicit none

  real(DP) :: X(1),Y(1),B(1)
  integer :: n

  REAL(DP) :: DT, T
  integer :: isbcnt, i

  if (isclld == 0) then
      isclld=1
      nrout=nrout+1
      myrout=nrout
      rname(myrout) = 'VLSC3 '
  endif
  isbcnt = 3*n
  dct(myrout) = dct(myrout) + float(isbcnt)
  ncall(myrout) = ncall(myrout) + 1
  dcount      =      dcount + float(isbcnt)

  DT = 0.0
  DO 10 I=1,N
      T = X(I)*Y(I)*B(I)
      DT = DT+T
  10 END DO
  T=DT
  VLSC3 = T
  RETURN
END FUNCTION VLSC3

!-----------------------------------------------------------------------
!> \brief blank a string
SUBROUTINE BLANK(A,N)
  implicit none
  CHARACTER(1) :: A(*)
  integer :: n
  CHARACTER(1) :: BLNK = ' '
  integer :: i

  DO 10 I=1,N
      A(I)=BLNK
  10 END DO
  RETURN
END SUBROUTINE BLANK

!-----------------------------------------------------------------------
    subroutine copy(a,b,n)
    real :: a(1),b(1)

    do i=1,n
        a(i)=b(i)
    enddo

    return
    end subroutine copy
!-----------------------------------------------------------------------
    subroutine chcopy(a,b,n)
    CHARACTER(1) :: A(1), B(1)

    DO 100 I = 1, N
        A(I) = B(I)
    100 END DO
    return
    end subroutine chcopy

    subroutine icopy(a,b,n)
    INTEGER :: A(1), B(1)

    DO 100 I = 1, N
        A(I) = B(I)
    100 END DO
    return
    end subroutine icopy
!-----------------------------------------------------------------------
    subroutine i8copy(a,b,n)
    INTEGER*8 :: A(1), B(1)

    DO 100 I = 1, N
        A(I) = B(I)
    100 END DO
    return
    end subroutine i8copy
!-----------------------------------------------------------------------
    subroutine chsign(a,n)
    REAL :: A(1)

    DO 100 I=1,N
        A(I) = -A(I)
    100 END DO
    return
    end subroutine chsign

!-----------------------------------------------------------------------
    real function vlmin(vec,n)
    REAL :: VEC(1)
    TMIN = 99.0E20

    DO 100 I=1,N
        TMIN = MIN(TMIN,VEC(I))
    100 END DO
    VLMIN = TMIN
    return
    end function vlmin
!-----------------------------------------------------------------------
    integer function ivlmin(vec,n)
    integer :: vec(1),tmin
    if (n == 0) then
        ivlmin=0
        return
    endif
    tmin = 8888888
    do i=1,n
        tmin = min(tmin,vec(i))
    enddo
    ivlmin = tmin
    return
    end function ivlmin
!-----------------------------------------------------------------------
    integer function ivlmax(vec,n)
    integer :: vec(1),tmax
    if (n == 0) then
        ivlmax=0
        return
    endif
    TMAX =-8888888
    do i=1,n
        TMAX = MAX(TMAX,VEC(I))
    enddo
    Ivlmax = tmax
    return
    end function ivlmax
!-----------------------------------------------------------------------
    real function vlmax(vec,n)
    REAL :: VEC(1)
    TMAX =-99.0E20
    do i=1,n
        TMAX = MAX(TMAX,VEC(I))
    enddo
    VLMAX = TMAX
    return
    end function vlmax
!-----------------------------------------------------------------------
    real function vlamax(vec,n)
    REAL :: VEC(1)
    TAMAX = 0.0

    DO 100 I=1,N
        TAMAX = MAX(TAMAX,ABS(VEC(I)))
    100 END DO
    VLAMAX = TAMAX
    return
    end function vlamax
!-----------------------------------------------------------------------
    real function vlsum(vec,n)
    use opctr
    REAL :: VEC(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'vlsum '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    SUM = 0.

    DO 100 I=1,N
        SUM=SUM+VEC(I)
    100 END DO
    VLSUM = SUM
    return
    end function vlsum
!-----------------------------------------------------------------------
    subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)

!     Compute a Cartesian vector cross product.

    DIMENSION U1(1),U2(1),U3(1)
    DIMENSION V1(1),V2(1),V3(1)
    DIMENSION W1(1),W2(1),W3(1)


    DO 100 I=1,N
        U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
        U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
        U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
    100 END DO
    return
    end subroutine vcross
!-----------------------------------------------------------------------
    subroutine vdot2 (dot,u1,u2,v1,v2,n)

!     Compute a Cartesian vector dot product. 2-d version

    DIMENSION DOT(1)
    DIMENSION U1(1),U2(1)
    DIMENSION V1(1),V2(1)


    DO 100 I=1,N
        DOT(I) = U1(I)*V1(I) + U2(I)*V2(I)
    100 END DO
    return
    end subroutine vdot2
!-----------------------------------------------------------------------
    subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)

!     Compute a Cartesian vector dot product. 3-d version

    DIMENSION DOT(1)
    DIMENSION U1(1),U2(1),U3(1)
    DIMENSION V1(1),V2(1),V3(1)


    DO 100 I=1,N
        DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
    100 END DO
    return
    end subroutine vdot3
!-----------------------------------------------------------------------
    function mod1(i,n)

!     Yields MOD(I,N) with the exception that if I=K*N, result is N.

    MOD1=0
    IF (I == 0) THEN
        return
    ENDIF
    IF (N == 0) THEN
        WRITE(6,*) &
        'WARNING:  Attempt to take MOD(I,0) in function mod1.'
        return
    ENDIF
    II = I+N-1
    MOD1 = MOD(II,N)+1
    return
    end function mod1
!-----------------------------------------------------------------------
    integer function log2(k)
    RK=(K)
    RLOG=LOG10(RK)
    RLOG2=LOG10(2.0)
    RLOG=RLOG/RLOG2+0.5
    LOG2=INT(RLOG)
    return
    end function log2
!-----------------------------------------------------------------------
    subroutine iswap(b,ind,n,temp)
    INTEGER :: B(1),IND(1),TEMP(1)
!***
!***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
!***  INTO ITEM(I), WHERE JJ=IND(I).
!***
    DO 20 I=1,N
        JJ=IND(I)
        TEMP(I)=B(JJ)
    20 END DO
    DO 30 I=1,N
        B(I)=TEMP(I)
    30 END DO
    return
    end subroutine iswap
!-----------------------------------------------------------------------
    real function vlsc2(x,y,n)
    use size_m
    use parallel
    use opctr
    REAL :: X(1),Y(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'VLSC2 '
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    s = 0.
    do i=1,n
        s = s + x(i)*y(i)
    enddo
    vlsc2=s
    return
    end function vlsc2
!-----------------------------------------------------------------------

!----------------------------------------------------------------------------

!     Vector reduction routines which require communication
!     on a parallel machine. These routines must be substituted with
!     appropriate routines which take into account the specific architecture.

!----------------------------------------------------------------------------


    function glsc3(a,b,mult,n)

!     Perform inner-product in double precision

    use opctr
    REAL :: A(1),B(1),MULT(1)
    REAL :: TMP,WORK(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'glsc3 '
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    TMP = 0.0
    DO 10 I=1,N
        TMP = TMP + A(I)*B(I)*MULT(I)
    10 END DO
    CALL GOP(TMP,WORK,'+  ',1)
    GLSC3 = TMP
    return
    end function glsc3
!-----------------------------------------------------------------------
    function glsc2(x,y,n)

!     Perform inner-product in double precision

    use opctr

    real :: x(1), y(1)
    real :: tmp,work(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'glsc2 '
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    tmp=0.0
    do 10 i=1,n
        tmp = tmp+ x(i)*y(i)
    10 END DO
    CALL GOP(TMP,WORK,'+  ',1)
    GLSC2 = TMP
    return
    end function glsc2
!-----------------------------------------------------------------------
    function glsc23(x,y,z,n)

!     Perform inner-product  x*x*y*z

    real :: x(1), y(1),z(1)
    real :: tmp,work(1)

    ds = 0.0
    do 10 i=1,n
        ds=ds+x(i)*x(i)*y(i)*z(i)
    10 END DO
    tmp=ds
    call gop(tmp,work,'+  ',1)
    glsc23 = tmp
    return
    end function glsc23
!-----------------------------------------------------------------------
    function glsum (x,n)
    DIMENSION X(1)
    DIMENSION TMP(1),WORK(1)
    TSUM = 0.
    DO 100 I=1,N
        TSUM = TSUM+X(I)
    100 END DO
    TMP(1)=TSUM
    CALL GOP(TMP,WORK,'+  ',1)
    GLSUM = TMP(1)
    return
    end function glsum
!-----------------------------------------------------------------------
    real function glamax(a,n)
    REAL :: A(1)
    DIMENSION TMP(1),WORK(1)
    TMAX = 0.0
    DO 100 I=1,N
        TMAX = MAX(TMAX,ABS(A(I)))
    100 END DO
    TMP(1)=TMAX
    CALL GOP(TMP,WORK,'M  ',1)
    GLAMAX=ABS(TMP(1))
    return
    end function glamax
!-----------------------------------------------------------------------
    function iglmin(a,n)
    integer :: a(1),tmin
    integer :: tmp(1),work(1)
    tmin=  999999999
    do i=1,n
        tmin=min(tmin,a(i))
    enddo
    tmp(1)=tmin
    call igop(tmp,work,'m  ',1)
    iglmin=tmp(1)
    return
    end function iglmin
!-----------------------------------------------------------------------
    function iglmax(a,n)
    integer :: a(1),tmax
    integer :: tmp(1),work(1)
    tmax= -999999999
    do i=1,n
        tmax=max(tmax,a(i))
    enddo
    tmp(1)=tmax
    call igop(tmp,work,'M  ',1)
    iglmax=tmp(1)
    return
    end function iglmax
!-----------------------------------------------------------------------
    function iglsum(a,n)
    integer :: a(1),tsum
    integer :: tmp(1),work(1)
    tsum= 0
    do i=1,n
        tsum=tsum+a(i)
    enddo
    tmp(1)=tsum
    call igop(tmp,work,'+  ',1)
    iglsum=tmp(1)
    return
    end function iglsum
!-----------------------------------------------------------------------
    integer*8 function i8glsum(a,n)
    integer*8 :: a(1),tsum
    integer*8 :: tmp(1),work(1)
    tsum= 0
    do i=1,n
        tsum=tsum+a(i)
    enddo
    tmp(1)=tsum
    call i8gop(tmp,work,'+  ',1)
    i8glsum=tmp(1)
    return
    END function
!-----------------------------------------------------------------------
    function glmax(a,n)
    REAL :: A(1)
    DIMENSION TMP(1),WORK(1)
    TMAX=-99.0e20
    DO 100 I=1,N
        TMAX=MAX(TMAX,A(I))
    100 END DO
    TMP(1)=TMAX
    CALL GOP(TMP,WORK,'M  ',1)
    GLMAX=TMP(1)
    return
    end function glmax
!-----------------------------------------------------------------------
    function glmin(a,n)
    REAL :: A(1)
    DIMENSION TMP(1),WORK(1)
    TMIN=99.0e20
    DO 100 I=1,N
        TMIN=MIN(TMIN,A(I))
    100 END DO
    TMP(1)=TMIN
    CALL GOP(TMP,WORK,'m  ',1)
    GLMIN = TMP(1)
    return
    end function glmin
!-----------------------------------------------------------------------
    subroutine gllog(la,lb)

!     If ANY LA=LB, then ALL LA=LB.

    LOGICAL :: LA,LB
    DIMENSION TMP(1),WORK(1)

    TMP(1)=1.0
    IF (LB) THEN
        IF (LA) TMP(1)=0.0
    ELSE
        IF ( .NOT. LA) TMP(1)=0.0
    ENDIF
    CALL GOP(TMP,WORK,'*  ',1)
    IF (TMP(1) == 0.0) LA=LB
    return
    end subroutine gllog
!-----------------------------------------------------------------------

!========================================================================
!     Double precision matrix and vector routines
!========================================================================

!-----------------------------------------------------------------------
    subroutine chswapr(b,L,ind,n,temp)
    INTEGER :: IND(1)
    CHARACTER(6) :: B(1),TEMP(1)
!***
!***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
!***  INTO ITEM(I), WHERE JJ=IND(I).
!***
    DO 20 I=1,N
        JJ=IND(I)
        TEMP(I)=B(JJ)
    20 END DO
    DO 30 I=1,N
        B(I)=TEMP(I)
    30 END DO
    return
    end subroutine chswapr
!-----------------------------------------------------------------------
    subroutine drcopy(r,d,N)
    real*8 ::    d(1)
    dimension r(1)
    do 10 i=1,n
        r(i)=d(i)
    10 END DO
    return
    end subroutine drcopy
!-----------------------------------------------------------------------
    subroutine rrcopy(r,d,N)
    real*4 :: d(1)
    real*4 :: r(1)
    do 10 i=1,n
        r(i)=d(i)
    10 END DO
    return
    end subroutine rrcopy
!-----------------------------------------------------------------------
    subroutine isort(a,ind,n)

!     Use Heap Sort (p 231 Num. Rec., 1st Ed.)

    integer :: a(1),ind(1)
    integer :: aa

    dO 10 j=1,n
        ind(j)=j
    10 END DO

    if (n <= 1) return
    L=n/2+1
    ir=n
    100 continue
    if (l > 1) then
        l=l-1
        aa  = a  (l)
        ii  = ind(l)
    else
        aa =   a(ir)
        ii = ind(ir)
        a(ir) =   a( 1)
        ind(ir) = ind( 1)
        ir=ir-1
        if (ir == 1) then
            a(1) = aa
            ind(1) = ii
            return
        endif
    endif
    i=l
    j=l+l
    200 continue
    if (j <= ir) then
        if (j < ir) then
            if ( a(j) < a(j+1) ) j=j+1
        endif
        if (aa < a(j)) then
            a(i) = a(j)
            ind(i) = ind(j)
            i=j
            j=j+j
        else
            j=ir+1
        endif
        GOTO 200
    endif
    a(i) = aa
    ind(i) = ii
    GOTO 100
    end subroutine isort
    subroutine sort(a,ind,n)

!     Use Heap Sort (p 231 Num. Rec., 1st Ed.)

    real :: a(1),aa
    integer :: ind(1)

    dO 10 j=1,n
        ind(j)=j
    10 END DO

    if (n <= 1) return
    L=n/2+1
    ir=n
    100 continue
    if (l > 1) then
        l=l-1
        aa  = a  (l)
        ii  = ind(l)
    else
        aa =   a(ir)
        ii = ind(ir)
        a(ir) =   a( 1)
        ind(ir) = ind( 1)
        ir=ir-1
        if (ir == 1) then
            a(1) = aa
            ind(1) = ii
            return
        endif
    endif
    i=l
    j=l+l
    200 continue
    if (j <= ir) then
        if (j < ir) then
            if ( a(j) < a(j+1) ) j=j+1
        endif
        if (aa < a(j)) then
            a(i) = a(j)
            ind(i) = ind(j)
            i=j
            j=j+j
        else
            j=ir+1
        endif
        GOTO 200
    endif
    a(i) = aa
    ind(i) = ii
    GOTO 100
    end subroutine sort
!-----------------------------------------------------------------------
    integer*8 function i8glmax(a,n)
    integer*8 :: a(1),tmax
    integer*8 :: tmp(1),work(1)
    tmax= -999999
    do i=1,n
        tmax=max(tmax,a(i))
    enddo
    tmp(1)=tmax
    call i8gop(tmp,work,'M  ',1)
    i8glmax=tmp(1)
    if (i8glmax == -999999) i8glmax=0
    return
    END function
!-----------------------------------------------------------------------
subroutine ident(a,n)
  use kinds, only : DP
  implicit none
  integer :: n, i
  real(DP) ::  a(n,n)

  a = 0._dp
  do i=1,n
      a(i,i) = 1._dp
  enddo
  return
end subroutine ident
!-----------------------------------------------------------------------
subroutine iswapt_ip(x,p,n)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: x(n)
  integer, intent(inout) :: p(n)

  integer :: j, k, loop_start, next, nextp, t1, t2

!   In-place permutation: x'(p) = x

  do k=1,n
      if (p(k) > 0) then   ! not swapped
          loop_start = k
          next       = p(loop_start)
          t1         = x(loop_start)
          do j=1,n
              if (next < 0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt
              elseif (next == loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
              else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
              endif
          enddo
          10 continue
      endif
  enddo

  do k=1,n
      p(k) = -p(k)
  enddo
  return
end subroutine iswapt_ip
