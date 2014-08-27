!-----------------------------------------------------------------------
    SUBROUTINE BLANK(A,N)
    CHARACTER(1) :: A(1)
    CHARACTER(1) :: BLNK
    SAVE        BLNK
    DATA        BLNK /' '/

    DO 10 I=1,N
        A(I)=BLNK
    10 END DO
    RETURN
    END SUBROUTINE BLANK
!-----------------------------------------------------------------------
    SUBROUTINE VSQ (A,N)
    use opctr
    DIMENSION  A(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'vsq   '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I = 1, N
        A(I) = A(I)**2
    100 END DO
    RETURN
    END SUBROUTINE VSQ
!-----------------------------------------------------------------------
    SUBROUTINE VSQRT(A,N)
    use opctr
    DIMENSION  A(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'vsqrt '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I = 1, N
        A(I) = SQRT(A(I))
    100 END DO
    RETURN
    END SUBROUTINE VSQRT
!-----------------------------------------------------------------------
    subroutine invers2(a,b,n)
    use opctr
    REAL :: A(1),B(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'inver2'
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=1./B(I)
    100 END DO
    return
    end subroutine invers2
!-----------------------------------------------------------------------
    subroutine invcol1(a,n)
    use opctr
    REAL :: A(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'invcl1'
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=1./A(I)
    100 END DO
    return
    end subroutine invcol1
!-----------------------------------------------------------------------
    subroutine invcol2(a,b,n)
    use ctimer
    use opctr

    REAL :: A(1),B(1)

#ifndef NOTIMER
    if (icalld == 0) tinvc=0.0
    icalld=icalld+1
    ninvc=icalld
    etime1=dnekclock()



    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'invcl2'
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)/B(I)
    100 END DO
#ifndef NOTIMER
    tinvc=tinvc+(dnekclock()-etime1)
#endif
    return
    end subroutine invcol2
!-----------------------------------------------------------------------
    subroutine invcol3(a,b,c,n)
    use ctimer
    use opctr
    REAL :: A(1),B(1),C(1)


#ifndef NOTIMER
    if (icalld == 0) tinv3=0.0
    icalld=icalld+1
    ninv3=icalld
    etime1=dnekclock()


    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'invcl3'
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=B(I)/C(I)
    100 END DO
#ifndef NOTIMER
    tinv3=tinv3+(dnekclock()-etime1)
#endif
    return
    end subroutine invcol3
!-----------------------------------------------------------------------
    subroutine col4(a,b,c,d,n)
    use opctr
    REAL :: A(1),B(1),C(1),D(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'col4  '
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=B(I)*C(I)*D(I)
    100 END DO
    return
    end subroutine col4
!-----------------------------------------------------------------------
    subroutine Xaddcol3(a,b,c,n)
    use opctr
    REAL :: A(1),B(1),C(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'addcl3'
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)+B(I)*C(I)
    100 END DO
    return
    end subroutine Xaddcol3
!-----------------------------------------------------------------------
    subroutine addcol4(a,b,c,d,n)
    use opctr
    REAL :: A(1),B(1),C(1),D(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'addcl4'
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)+B(I)*C(I)*D(I)
    100 END DO
    return
    end subroutine addcol4
!-----------------------------------------------------------------------
    subroutine ascol5 (a,b,c,d,e,n)
    use opctr
    REAL :: A(1),B(1),C(1),D(1),E(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'ascol5'
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I) = B(I)*C(I)-D(I)*E(I)
    100 END DO
    return
    end subroutine ascol5
!-----------------------------------------------------------------------
    subroutine sub2(a,b,n)
    use opctr
    REAL :: A(1),B(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'sub2  '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)-B(I)
    100 END DO
    return
    end subroutine sub2
!-----------------------------------------------------------------------
    subroutine sub3(a,b,c,n)
    use opctr
    REAL :: A(1),B(1),C(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'sub3  '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=B(I)-C(I)
    100 END DO
    return
    end subroutine sub3
!-----------------------------------------------------------------------
    subroutine subcol3(a,b,c,n)
    use opctr
    REAL :: A(1),B(1),C(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'subcl3'
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)-B(I)*C(I)
    100 END DO
    return
    end subroutine subcol3
!-----------------------------------------------------------------------
    subroutine subcol4(a,b,c,d,n)
    use opctr
    REAL :: A(1),B(1),C(1),D(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'subcl4'
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)-B(I)*C(I)*D(I)
    100 END DO
    return
    end subroutine subcol4
!-----------------------------------------------------------------------
    subroutine rzero(a,n)
    DIMENSION  A(1)
    DO 100 I = 1, N
        A(I ) = 0.0
    100 END DO
    return
    end subroutine rzero
!-----------------------------------------------------------------------
    subroutine izero(a,n)
    INTEGER :: A(1)

    DO 100 I = 1, N
        A(I ) = 0
    100 END DO
    return
    end subroutine izero
!-----------------------------------------------------------------------
    subroutine ione(a,n)
    INTEGER ::   A(1)
    DO 100 I = 1, N
        A(I ) = 1
    100 END DO
    return
    end subroutine ione
!-----------------------------------------------------------------------
    subroutine rone(a,n)
    DIMENSION  A(1)
    DO 100 I = 1, N
        A(I ) = 1.0
    100 END DO
    return
    end subroutine rone
!-----------------------------------------------------------------------
    subroutine cfill(a,b,n)
    DIMENSION  A(1)

    DO 100 I = 1, N
        A(I) = B
    100 END DO
    return
    end subroutine cfill
!-----------------------------------------------------------------------
    subroutine ifill(ia,ib,n)
    DIMENSION IA(1)

    DO 100 I = 1, N
        IA(I) = IB
    100 END DO
    return
    end subroutine ifill
!-----------------------------------------------------------------------
    subroutine copy(a,b,n)
    real :: a(1),b(1)

    do i=1,n
        a(i)=b(i)
    enddo

    return
    end subroutine copy
!-----------------------------------------------------------------------
    subroutine copyi4(a,b,n)
    integer :: a(1)
    real ::    b(1)

    do i=1,n
        a(i)=b(i)
    enddo

    return
    end subroutine copyi4
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
    subroutine cmult(a,const,n)
    use opctr
    REAL :: A(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'cmult '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)*CONST
    100 END DO
    return
    end subroutine cmult
!-----------------------------------------------------------------------
    subroutine cadd(a,const,n)
    use opctr
    REAL :: A(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'cadd  '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)+CONST
    100 END DO
    return
    end subroutine cadd
!-----------------------------------------------------------------------
    subroutine iadd(i1,iscal,n)
    DIMENSION I1(1)

    DO 10 I=1,N
        I1(I)=I1(I)+ISCAL
    10 END DO
    return
    end subroutine iadd
!-----------------------------------------------------------------------
    subroutine cadd2(a,b,const,n)
    use opctr
    REAL :: A(1),B(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'cadd2 '
    endif
    isbcnt = n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=B(I)+CONST
    100 END DO
    return
    end subroutine cadd2
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
    subroutine addtnsr(s,h1,h2,h3,nx,ny,nz)

!     Map and add to S a tensor product form of the three functions H1,H2,H3.
!     This is a single element routine used for deforming geometry.

    DIMENSION H1(1),H2(1),H3(1)
    DIMENSION S(NX,NY,NZ)

    DO 200 IZ=1,NZ
        DO 200 IY=1,NY
            HH = H2(IY)*H3(IZ)
            DO 100 IX=1,NX
                S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
            100 END DO
    200 END DO
    return
    end subroutine addtnsr

    function ltrunc(string,l)
    CHARACTER(1) :: STRING(L)
    CHARACTER(1) ::   BLNK
    DATA BLNK/' '/

    DO 100 I=L,1,-1
        L1=I
        IF (STRING(I) /= BLNK) GOTO 200
    100 END DO
    L1=0
    200 CONTINUE
    LTRUNC=L1
    return
    end function ltrunc
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
    subroutine iflip(i1,n)
    DIMENSION I1(1)
    N1=N+1
    N2=N/2
    DO 10 I=1,N2
        ILAST=N1-I
        ITMP=I1(ILAST)
        I1(ILAST)=I1(I)
        I1(I)=ITMP
    10 END DO
    return
    end subroutine iflip
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
    subroutine col2(a,b,n)
    use opctr
    real :: a(1),b(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'col2  '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! bm* unroll (10)
    do i=1,n
        a(i)=a(i)*b(i)
    enddo

    return
    end subroutine col2
!-----------------------------------------------------------------------
    subroutine col2c(a,b,c,n)
    real :: a(1),b(1),c

    do i=1,n
        a(i)=a(i)*b(i)*c
    enddo

    return
    end subroutine col2c
!-----------------------------------------------------------------------
    subroutine col3(a,b,c,n)
    use opctr
    real :: a(1),b(1),c(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'col3  '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! bm* unroll (10)
    do i=1,n
        a(i)=b(i)*c(i)
    enddo
    return
    end subroutine col3
!-----------------------------------------------------------------------
    subroutine add2(a,b,n)
    use opctr
    real :: a(1),b(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'ADD2  '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! bm* unroll (10)
    do i=1,n
        a(i)=a(i)+b(i)
    enddo
    return
    end subroutine add2
!-----------------------------------------------------------------------
    subroutine add3(a,b,c,n)
    use opctr
    real :: a(1),b(1),c(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'ADD3  '
    endif
    isbcnt = N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! bm* unroll (10)
    do i=1,n
        a(i)=b(i)+c(i)
    enddo
    return
    end subroutine add3
!-----------------------------------------------------------------------
    subroutine addcol3(a,b,c,n)
    use opctr
    real :: a(1),b(1),c(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'addcl3'
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

! bm* unroll (10)
    do i=1,n
        a(i)=a(i)+b(i)*c(i)
    enddo
    return
    end subroutine addcol3
!-----------------------------------------------------------------------
    subroutine add2s1(a,b,c1,n)
    use opctr
    real :: a(1),b(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'add2s1'
    endif
    isbcnt = 2*N
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
    100 END DO
    return
    end subroutine add2s1

!-----------------------------------------------------------------------
    subroutine add2s2(a,b,c1,n)
    use opctr
    real :: a(1),b(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'add2s2'
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
    100 END DO
    return
    end subroutine add2s2

!-----------------------------------------------------------------------
    subroutine add3s2(a,b,c,c1,c2,n)
    use opctr
    real :: a(1),b(1),c(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'add3s2'
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
    100 END DO
    return
    end subroutine add3s2

!-----------------------------------------------------------------------
    subroutine add4(a,b,c,d,n)
    use opctr
    REAL :: A(1),B(1),C(1),D(1)


#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'add4  '
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    DO 100 I=1,N
        A(I)=B(I)+C(I)+D(I)
    100 END DO
    return
    end subroutine add4
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
    real function vlsc21(x,y,n)
    use size_m
    use parallel
    use opctr
    real :: x(1),y(1)

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'VLSC21'
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    s = 0.
    do i=1,n
        s = s + x(i)*x(i)*y(i)
    enddo
    vlsc21=s
    return
    end function vlsc21


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
!> \brief global 2-norm
real function gl2norm(a,n)
  use kinds, only : DP
  use mass, only : bm1, volvm1
  implicit none

  real(DP) :: a(1)
  integer :: n
  real(DP), external :: glsum
  real(DP) :: w1(n)

  call col3 (w1,a,a,n)
  call col2 (w1,bm1,n)
  gl2norm = sqrt(glsum (w1,n)/volvm1)

  return
end function gl2norm
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
    real function glamin(a,n)
    real :: a(1)
    dimension tmp(1),work(1)
    tmin = 9.e28
    do 100 i=1,n
        tmin = min(tmin,abs(a(i)))
    100 END DO
    tmp(1)=tmin
    call gop(tmp,work,'m  ',1)
    glamin=abs(tmp(1))
    return
    end function glamin
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
    function fmdian(a,n,ifok)
!     find the Median of the (global) set A
    use size_m
    DIMENSION A(1)
    DIMENSION WORK1(5),WORK2(5)
    DIMENSION GUES(100)
    LOGICAL :: IFOK

    AMP  =1.5
    AFAC =1.5
    GMIN =GLMIN(A,N)
    GMAX =GLMAX(A,N)
    GMIN0=GLMIN(A,N)
    GMAX0=GLMAX(A,N)
    GUESS=(GMAX+GMIN)/2.0
    EPS  =(GMAX-GMIN)
    IF (EPS == 0.0) THEN
        FMDIAN=GMAX
        return
    ENDIF
    WORK1(1)=N
    CALL GOP(WORK1,WORK2,'+  ',1)
    NTOT=WORK1(1)
    N2 = (NTOT+1)/2
    IF ( .NOT. IFOK) THEN
        WRITE(6,8) NID,N,(A(I),I=1,N)
        WRITE(6,9) NID,NTOT,N2,N,GMIN,GMAX
        8 FORMAT(I5,'N,A:',I5,10(6F10.5,/))
        9 FORMAT(I5,'mnx:',3I6,2F10.5)
    ENDIF

!     This is the trial loop

    ITRY=-1
    10 CONTINUE
    ITRY=ITRY+1
    II=ITRY+1
    IF (II <= 100) GUES(II)=GUESS
!     error check for infinite loop
    IF (ITRY > 2*NTOT) GOTO 9000
    CALL RZERO(WORK1,5)
    NLT=0
    NGT=0
    CLT=GMIN0
    CGT=GMAX0
    DO 100 I=1,N
        AA=A(I)
        IF (AA /= GUESS) THEN
            IF (AA < GUESS) THEN
                NLT=NLT+1
            !              CLT - closest value to GUESS Less Than GUESS
                IF (AA > CLT) CLT=AA
            ENDIF
            IF (AA > GUESS) THEN
                NGT=NGT+1
            !              CGT - closest value to GUESS Greater Than GUESS
                IF (AA < CGT) CGT=AA
            ENDIF
            DUM=1./(EPS+ABS(AA-GUESS))
            WORK1(1)=WORK1(1)+DUM
            WORK1(2)=WORK1(2)+DUM*AA
        ELSE
        !           detected values equaling the guess.
            WORK1(5)=WORK1(5)+1.0
        ENDIF
    100 END DO
!     Invoke vector reduction across processors:
    WORK2(1)=CLT
    CLT=GLMAX(WORK2,1)
    WORK2(1)=CGT
    CGT=GLMIN(WORK2,1)
    WORK1(3)=NLT
    WORK1(4)=NGT
    CALL GOP(WORK1,WORK2,'+  ',5)
    NLT=WORK1(3)
    NGT=WORK1(4)
    IF ( .NOT. IFOK) THEN
        WRITE(6,101) NID,GUESS,CLT,CGT
        WRITE(6,102) NID,(WORK1(I),I=1,5)
        101 FORMAT(I5,'Glg:',3F12.5)
        102 FORMAT(I5,'WORK1:',5F12.5)
    ENDIF

!     Done?

    IF (NLT > N2 .OR. NGT > N2) THEN
    !        we're not done.....
        IF (NGT > NLT) THEN
        !           guess is too low
            GMIN=CGT
            G2=CGT+MAX(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2 > GMAX) G2=0.5*(GUESS+GMAX)
            EPS=AFAC*ABS(G2-GUESS)
        !           see that we move at least as far as the next closest value.
            GUESS=MAX(G2,CGT)
            GOTO 10
        ELSE IF (NLT > NGT) THEN
        !           guess is too high
            GMAX=CLT
            G2=CLT+MIN(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2 < GMIN) G2=0.5*(GUESS+GMIN)
            EPS=AFAC*ABS(G2-GUESS)
        !           see that we move at least as far as the next closest value.
            GUESS=MIN(G2,CLT)
            GOTO 10
        ENDIF
    ELSE
    
    !        we're done....
        IF (WORK1(5) /= 0) THEN
        !           the median is (usually) one of the values
            FMDIAN=GUESS
            IF (WORK1(5) == 1.0) THEN
                IF (MOD(NTOT,2) == 0) THEN
                    IF (NGT > NLT) THEN
                        FMDIAN=0.5*(GUESS+CGT)
                    ELSE
                        FMDIAN=0.5*(GUESS+CLT)
                    ENDIF
                ELSE
                    IF (NGT == NLT) THEN
                        FMDIAN=GUESS
                    ELSE IF(NGT > NLT) THEN
                        FMDIAN=CGT
                    ELSE
                        FMDIAN=CLT
                    ENDIF
                ENDIF
            ENDIF
        ELSE
            IF (MOD(NTOT,2) == 0) THEN
                IF (NGT == NLT) THEN
                    FMDIAN=0.5*(CLT+CGT)
                ELSE IF(NGT > NLT) THEN
                    FMDIAN=0.5*(GUESS+CGT)
                ELSE
                    FMDIAN=0.5*(GUESS+CLT)
                ENDIF
            ELSE
                IF (NGT == NLT) THEN
                    FMDIAN=GUESS
                ELSE IF(NGT > NLT) THEN
                    FMDIAN=CGT
                ELSE
                    FMDIAN=CLT
                ENDIF
            ENDIF
        ENDIF
    
    ENDIF
    IF ( .NOT. IFOK) WRITE(6,*) NID,'FMDIAN2',FMDIAN,(A(I),I=1,N)
    return

!     Error handling

    9000 CONTINUE
    WRITE(6,11) NTOT,GMIN0,GMAX0,GUESS
    11 FORMAT('ABORTING IN FMDIAN: N,AMIN,AMAX:',I6,3G14.6)
    DO 13 I1=1,N,5
        IN=I1+5
        IN=MIN(IN,N)
        WRITE(6,12) NID,(A(I),I=I1,IN)
        12 FORMAT(I4,' FMA:',5G14.6)
    13 END DO
    DO 15 I1=1,ITRY,5
        IN=I1+5
        IN=MIN(IN,ITRY)
        WRITE(6,14) NID,(GUES(I),I=I1,IN)
        14 FORMAT(I4,' FMG:',5G14.6)
    15 END DO
    call exitt
    end function fmdian

!========================================================================
!     Double precision matrix and vector routines
!========================================================================

!-----------------------------------------------------------------------
    subroutine dcadd(a,const,n)
    real*8 :: A(1),CONST

    DO 100 I=1,N
        A(I)=A(I)+CONST
    100 END DO
    return
    end subroutine dcadd
!-----------------------------------------------------------------------
    subroutine dsub2(a,b,n)
    real*8 :: A(1), B(1)

    DO 100 I=1,N
        A(I)=A(I)-B(I)
    100 END DO
    return
    end subroutine dsub2

!-----------------------------------------------------------------------
    subroutine dadd2(a,b,n)
    real*8 :: A(1), B(1)

    DO 100 I=1,N
        A(I)=A(I)+B(I)
    100 END DO
    return
    end subroutine dadd2
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
    subroutine sorts(xout,xin,work,n)
    real :: xout(1),xin(1),work(1)
    call copy(xout,xin,n)
    call sort(xout,work,n)
    return
    end subroutine sorts

!-----------------------------------------------------------------------
    function ivlsum(a,n)
    INTEGER :: A(1)
    INTEGER :: TSUM
    if (n == 0) then
        ivlsum = 0
        return
    endif
    TSUM=A(1)
    DO 100 I=2,N
        TSUM=TSUM+A(I)
    100 END DO
    IVLSUM=TSUM
    return
    end function ivlsum
!-----------------------------------------------------------------------
    subroutine icadd(a,c,n)
    INTEGER :: A(1),C
    DO 100 I = 1, N
        A(I) = A(I) + C
    100 END DO
    return
    end subroutine icadd
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
    subroutine iswap_ip(x,p,n)
    integer :: x(1),xstart
    integer :: p(1)

!     In-place permutation: x' = x(p)

    do k=1,n
        if (p(k) > 0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
                next    = p(last)
                if (next < 0) then
                    write(6,*) 'Hey! iswap_ip problem.',j,k,n,next
                    call exitt
                elseif (next == loop_start) then
                    x(last) = xstart
                    p(last) = -p(last)
                    goto 10
                else
                    x(last) = x(next)
                    p(last) = -p(last)
                    last    = next
                endif
            enddo
            10 continue
        endif
    enddo

    do k=1,n
        p(k) = -p(k)
    enddo
    return
    end subroutine iswap_ip
!-----------------------------------------------------------------------
    subroutine iswapt_ip(x,p,n)
    integer :: x(1),t1,t2
    integer :: p(1)

!     In-place permutation: x'(p) = x


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
!-----------------------------------------------------------------------
    subroutine swap_ip(x,p,n)
    real ::    x(1),xstart
    integer :: p(1)

!     In-place permutation: x' = x(p)

    do k=1,n
        if (p(k) > 0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
                next    = p(last)
                if (next < 0) then
                    write(6,*) 'Hey! swap_ip problem.',j,k,n,next
                    call exitt
                elseif (next == loop_start) then
                    x(last) = xstart
                    p(last) = -p(last)
                    goto 10
                else
                    x(last) = x(next)
                    p(last) = -p(last)
                    last    = next
                endif
            enddo
            10 continue
        endif
    enddo

    do k=1,n
        p(k) = -p(k)
    enddo
    return
    end subroutine swap_ip
!-----------------------------------------------------------------------
    subroutine swapt_ip(x,p,n)
    real ::    x(1),t1,t2
    integer :: p(1)

!     In-place permutation: x'(p) = x


    do k=1,n
        if (p(k) > 0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
                if (next < 0) then
                    write(6,*) 'Hey! swapt_ip problem.',j,k,n,next
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
    end subroutine swapt_ip
!-----------------------------------------------------------------------
    subroutine glvadd(x,w,n)
    real :: x(1),w(1)
    call gop(x,w,'+  ',n)
    return
    end subroutine glvadd
!-----------------------------------------------------------------------
    subroutine add3s12(x,y,z,c1,c2,n)
    real :: x(1),y(1),z(1),c1,c2
    do i=1,n
        x(i) = c1*y(i)+c2*z(i)
    enddo
    return
    end subroutine add3s12
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
    subroutine admcol3(a,b,c,d,n)
    REAL :: A(1),B(1),C(1),D

    DO 100 I=1,N
        A(I)=A(I)+B(I)*C(I)*D
    100 END DO
    return
    end subroutine admcol3
!-----------------------------------------------------------------------
    subroutine add2col2(a,b,c,n)
    real :: a(1),b(1),c(1)

    do i=1,n
        a(i) = a(i) + b(i)*c(i)
    enddo
    return
    end subroutine add2col2
!-----------------------------------------------------------------------
    subroutine add2sxy(x,a,y,b,n)
    real :: x(1),y(1)

    do i=1,n
        x(i) = a*x(i) + b*y(i)
    enddo

    return
    end subroutine add2sxy
!-----------------------------------------------------------------------
    subroutine col2s2(x,y,s,n)
    real :: x(n),y(n)

    do i=1,n
        x(i)=s*x(i)*y(i)
    enddo

    return
    end subroutine col2s2
!-----------------------------------------------------------------------
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
