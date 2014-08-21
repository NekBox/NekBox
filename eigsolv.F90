!**********************************************************************

!     ROUTINES FOR ESITMATING AND CALCULATING EIGENVALUES
!     USED IN NEKTON

!**********************************************************************

    SUBROUTINE ESTEIG
!--------------------------------------------------------------

!     Estimate eigenvalues

!-------------------------------------------------------------
    use size_m
    use eigen
    use geom
    use input
    INCLUDE 'TSTEP'

    NTOT1=NX1*NY1*NZ1*NELFLD(IFIELD)
    XMIN = GLMIN(XM1,NTOT1)
    XMAX = GLMAX(XM1,NTOT1)
    YMIN = GLMIN(YM1,NTOT1)
    YMAX = GLMAX(YM1,NTOT1)
    IF (IF3D) THEN
        ZMIN = GLMIN(ZM1,NTOT1)
        ZMAX = GLMAX(ZM1,NTOT1)
    ELSE
        ZMIN = 0.0
        ZMAX = 0.0
    ENDIF

    XX = XMAX - XMIN
    YY = YMAX - YMIN
    ZZ = ZMAX - ZMIN
    RXY = XX/YY
    RYX = YY/XX
    RMIN = RXY
    IF (RYX < RMIN) RMIN = RYX
    IF (NDIM == 3) THEN
        RXZ = XX/ZZ
        RZX = ZZ/XX
        RYZ = YY/ZZ
        RZY = ZZ/YY
        IF (RXZ < RMIN) RMIN = RXZ
        IF (RZX < RMIN) RMIN = RZX
        IF (RYZ < RMIN) RMIN = RYZ
        IF (RZY < RMIN) RMIN = RZY
    ENDIF

    XX2    = 1./XX**2
    YY2    = 1./YY**2
    XYZMIN = XX2
    XYZMAX = XX2+YY2
    IF (YY2 < XYZMIN) XYZMIN = YY2
    IF (NDIM == 3) THEN
        ZZ2 = 1./ZZ**2
        XYZMAX = XYZMAX+ZZ2
        IF (ZZ2 < XYZMIN) XYZMIN = ZZ2
    ENDIF

    one    = 1.
    PI     = 4.*ATAN(one)
    RATIO  = XYZMIN/XYZMAX
    EIGAE  = PI*PI*XYZMIN
    EIGGE  = EIGGA
    IF (NDIM == 2) EIGAA = PI*PI*(XX2+YY2)/2.
    IF (NDIM == 3) EIGAA = PI*PI*(XX2+YY2+ZZ2)/3.
    IF (IFAXIS)      EIGAA = .25*PI*PI*YY2
    EIGAS  = 0.25*RATIO
    EIGGS  = 2.0

    IF (NID == 0 .AND. ISTEP <= 0) THEN
        WRITE (6,*) ' '
        WRITE (6,*) 'Estimated eigenvalues'
        WRITE (6,*) 'EIGAA = ',EIGAA
        WRITE (6,*) 'EIGGA = ',EIGGA
        IF (IFFLOW) THEN
            WRITE (6,*) 'EIGAE = ',EIGAE
            WRITE (6,*) 'EIGAS = ',EIGAS
            WRITE (6,*) 'EIGGE = ',EIGGE
            WRITE (6,*) 'EIGGS = ',EIGGS
        ENDIF
        WRITE (6,*) ' '
    ENDIF

    RETURN
    END SUBROUTINE ESTEIG

    SUBROUTINE EIGENV
!-------------------------------------------------------------------------

!     Compute the following eigenvalues:
!     EIGAA  = minimum eigenvalue of the matrix A  (=Laplacian)
!     EIGAE  = minimum eigenvalue of the matrix E  (=DB-1DT)
!     EIGAS  = minimum eigenvalue of the matrix S  (=DA-1DT)
!     EIGAST = minimum eigenvalue of the matrix St (=D(A+B/dt)-1DT
!     EIGGA  = maximum eigenvalue of the matrix A
!     EIGGS  = maximum eigenvalue of the matrix S
!     EIGGE  = maximum eigenvalue of the matrix E
!     EIGGST = maximum eigenvalue of the matrix St

!     Method : Power method/Inverse iteration & Rayleigh quotient wo shift

!-------------------------------------------------------------------------
    use size_m
    use eigen
    use input
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'

    COMMON /SCRVH/ H1 (LX1,LY1,LZ1,LELT) &
    ,             H2 (LX1,LY1,LZ1,LELT)
    COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)

    NTOT1  = NX1*NY1*NZ1*NELV

    IF (IFAA) THEN
        NTOT1  = NX1*NY1*NZ1*NELV
        CALL RONE    (H1,NTOT1)
        CALL RZERO   (H2,NTOT1)
        CALL ALPHAM1 (EIGAA1,V1MASK,VMULT,H1,H2,1)
        CALL ALPHAM1 (EIGAA2,V2MASK,VMULT,H1,H2,2)
        EIGAA = MIN  (EIGAA1,EIGAA2)
        IF (NDIM == 3) THEN
            CALL ALPHAM1 (EIGAA3,V3MASK,VMULT,H1,H2,3)
            EIGAA = MIN  (EIGAA,EIGAA3)
        ENDIF
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAA = ',EIGAA
    ENDIF

    IF (IFAS) THEN
        INLOC = 0
        CALL RONE    (H1,NTOT1)
        CALL RZERO   (H2,NTOT1)
        CALL RZERO   (H2INV,NTOT1)
        CALL ALPHAM2  (EIGAS,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAS = ',EIGAS
    ENDIF

    IF (IFAE) THEN
        INLOC = 1
        CALL RZERO   (H1,NTOT1)
        CALL RONE    (H2,NTOT1)
        CALL RONE    (H2INV,NTOT1)
        CALL ALPHAM2  (EIGAE,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAE = ',EIGAE
    ENDIF

    IF (IFAST) THEN
        INLOC = -1
        CALL SETHLM  (H1,H2,INLOC)
        CALL INVERS2 (H2INV,H2,NTOT1)
        CALL ALPHAM2  (EIGAST,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAST = ',EIGAST
    ENDIF

    IF (IFGS) THEN
        INLOC = 0
        CALL RONE    (H1,NTOT1)
        CALL RZERO   (H2,NTOT1)
        CALL RZERO   (H2INV,NTOT1)
        CALL GAMMAM2 (EIGGS,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGS = ',EIGGS
    ENDIF

    IF (IFGE) THEN
        INLOC = 1
        CALL RZERO   (H1,NTOT1)
        CALL RONE    (H2,NTOT1)
        CALL RONE    (H2INV,NTOT1)
        CALL GAMMAM2 (EIGGE,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGE = ',EIGGE
    ENDIF

    IF (IFGST) THEN
        INLOC = -1
        CALL SETHLM  (H1,H2,INLOC)
        CALL INVERS2 (H2INV,H2,NTOT1)
        CALL GAMMAM2 (EIGGST,H1,H2,H2INV,INLOC)
        IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGST = ',EIGGST
    ENDIF

    IF (IFGA) THEN
        NTOT1  = NX1*NY1*NZ1*NELV
        CALL RONE    (H1,NTOT1)
        CALL RZERO   (H2,NTOT1)
        IF ( .NOT. IFSTRS) THEN
            CALL GAMMAM1 (EIGGA1,V1MASK,VMULT,H1,H2,1)
            CALL GAMMAM1 (EIGGA2,V2MASK,VMULT,H1,H2,2)
            EIGGA3 = 0.
            IF (NDIM == 3) &
            CALL GAMMAM1 (EIGGA3,V3MASK,VMULT,H1,H2,3)
            EIGGA = MAX  (EIGGA1,EIGGA2,EIGGA3)
        ELSE
!max            CALL GAMMASF (H1,H2)
        ENDIF
    ENDIF

    RETURN
    END SUBROUTINE EIGENV

    SUBROUTINE ALPHAM1 (ALPHA,MASK,MULT,H1,H2,ISD)
!---------------------------------------------------------------------------

!     Compute minimum eigenvalue, ALPHA, of the discrete Helmholtz operator

!---------------------------------------------------------------------------
    use size_m
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'

    REAL ::            MASK (LX1,LY1,LZ1,1)
    REAL ::            MULT (LX1,LY1,LZ1,1)
    REAL ::            H1   (LX1,LY1,LZ1,1)
    REAL ::            H2   (LX1,LY1,LZ1,1)
    COMMON /SCREV/  X1   (LX1,LY1,LZ1,LELT) &
    ,              Y1   (LX1,LY1,LZ1,LELT)
    CHARACTER NAME*4

    IF (IMESH == 1) NEL  = NELV
    IF (IMESH == 2) NEL  = NELT
    IF (ISD  == 1) NAME = 'EVVX'
    IF (ISD  == 2) NAME = 'EVVX'
    IF (ISD  == 3) NAME = 'EVVX'

    NXYZ1  = NX1*NY1*NZ1
    NTOT1  = NXYZ1*NEL
    EVNEW  = 0.
    CALL STARTX1 (X1,Y1,MASK,MULT,NEL)

    DO 1000 ITER=1,NMXE
        CALL AXHELM (Y1,X1,H1,H2,IMESH,ISD)
        CALL COL2   (Y1,MASK,NTOT1)
        CALL DSSUM  (Y1,NX1,NY1,NZ1)
        RQ     = GLSC3 (X1,Y1,MULT,NTOT1)
        EVOLD  = EVNEW
        EVNEW  = RQ
        write (6,*) 'alphaa = ',rq
        CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
        IF (CRIT < TOLEV)                  GOTO 2000
        CALL COL2    (X1,BM1,NTOT1)
        CALL HMHOLTZ ('NOMG',Y1,X1,H1,H2,MASK,MULT, &
        IMESH,TOLHE,NMXH,ISD)
        CALL COL3    (X1,BM1,Y1,NTOT1)
        CALL DSSUM   (X1,NX1,NY1,NZ1)
        YY = GLSC3  (X1,Y1,MULT,NTOT1)
        YNORM = 1./SQRT(YY)
        CALL CMULT   (Y1,YNORM,NTOT1)
        CALL COPY    (X1,Y1,NTOT1)
    1000 END DO
    2000 CONTINUE

    ALPHA = RQ
    RETURN
    END SUBROUTINE ALPHAM1

    SUBROUTINE GAMMAM1 (GAMMA,MASK,MULT,H1,H2,ISD)
!---------------------------------------------------------------------------

!     Compute maximum eigenvalue of the discrete Helmholtz operator

!---------------------------------------------------------------------------
    use size_m
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'

    REAL ::            MASK (LX1,LY1,LZ1,1)
    REAL ::            MULT (LX1,LY1,LZ1,1)
    REAL ::            H1   (LX1,LY1,LZ1,1)
    REAL ::            H2   (LX1,LY1,LZ1,1)
    COMMON /SCREV/  X1   (LX1,LY1,LZ1,LELT) &
    ,              Y1   (LX1,LY1,LZ1,LELT)

    IF (IMESH == 1) NEL = NELV
    IF (IMESH == 2) NEL = NELT
    NXYZ1  = NX1*NY1*NZ1
    NTOT1  = NXYZ1*NEL
    EVNEW  = 0.
!     pff (2/15/96)
    if (isd == 1) CALL STARTX1 (X1,Y1,MASK,MULT,NEL)

    DO 1000 ITER=1,NMXE
        CALL AXHELM (Y1,X1,H1,H2,IMESH,ISD)
        CALL COL2   (Y1,MASK,NTOT1)
        CALL DSSUM  (Y1,NX1,NY1,NZ1)
        RQ     = GLSC3 (X1,Y1,MULT,NTOT1)
        EVOLD  = EVNEW
        EVNEW  = RQ
        CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
    
    ! HMT removed
    
    !         if (nid.eq.0) then
    !            write(6,*) iter,' eig_max A:',evnew,crit,tolev
    !         endif
        IF (CRIT < TOLEV)                  GOTO 2000
        CALL COL3 (X1,BINVM1,Y1,NTOT1)
        XX     = GLSC3 (X1,Y1,MULT,NTOT1)
        XNORM  = 1./SQRT(XX)
        CALL CMULT (X1,XNORM,NTOT1)
    1000 END DO
    2000 CONTINUE

    GAMMA = RQ
    RETURN
    END SUBROUTINE GAMMAM1

    SUBROUTINE ALPHAM2 (ALPHA,H1,H2,H2INV,INLOC)
!----------------------------------------------------------------------

!     Compute minimum eigenvalue, ALPHA, of one of the matrices
!     defined on the pressure mesh:
!     INLOC =  0  : DA-1DT
!     INLOC =  1  : DB-1DT
!     INLOC = -1  : D(A+B/DT)-1DT

!----------------------------------------------------------------------
    use size_m
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'

    REAL ::           H1   (LX1,LY1,LZ1,1)
    REAL ::           H2   (LX1,LY1,LZ1,1)
    REAL ::           H2INV(LX1,LY1,LZ1,1)
    COMMON /SCREV/ X2   (LX2,LY2,LZ2,LELV) &
    ,             Y2   (LX2,LY2,LZ2,LELV)

    NTOT2  = NX2*NY2*NZ2*NELV
    EVNEW  = 0.
    CALL STARTX2 (X2,Y2)

    DO 1000 ITER=1,NMXE
        CALL CDABDTP (Y2,X2,H1,H2,H2INV,INLOC)
        RQ = GLSC2  (X2,Y2,NTOT2)
        EVOLD  = EVNEW
        EVNEW  = RQ
    !        write (6,*) 'new eigenvalue ************* eigas = ',evnew
        CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
        IF (CRIT < TOLEV)                GOTO 2000
        CALL COL2   (X2,BM2,NTOT2)
        CALL UZAWA  (X2,H1,H2,H2INV,INLOC,ICG)
        CALL COL3   (Y2,BM2,X2,NTOT2)
        XX = GLSC2 (X2,Y2,NTOT2)
        XNORM = 1./SQRT(XX)
        CALL CMULT  (X2,XNORM,NTOT2)
    1000 END DO
    2000 CONTINUE

    ALPHA = RQ
    RETURN
    END SUBROUTINE ALPHAM2

    SUBROUTINE GAMMAM2 (GAMMA,H1,H2,H2INV,INLOC)
!-------------------------------------------------------------------

!     Compute maximum eigenvalue, GAMMA, of one of the matrices
!     defined on the pressure mesh:
!     INLOC =  0  : DA-1DT
!     INLOC =  1  : DB-1DT
!     INLOC = -1  : D(A+B/DT)-1DT

!-------------------------------------------------------------------
    use size_m
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'

    REAL ::           H1    (LX1,LY1,LZ1,1)
    REAL ::           H2    (LX1,LY1,LZ1,1)
    REAL ::           H2INV (LX1,LY1,LZ1,1)
    COMMON /SCREV/ X2 (LX2,LY2,LZ2,LELV) &
    ,             Y2 (LX2,LY2,LZ2,LELV)

    NTOT2  = NX2*NY2*NZ2*NELV
    EVNEW  = 0.
    CALL STARTX2 (X2,Y2)

    DO 1000 ITER=1,NMXE
        CALL CDABDTP (Y2,X2,H1,H2,H2INV,INLOC)
        RQ = GLSC2  (X2,Y2,NTOT2)
        EVOLD  = EVNEW
        EVNEW  = RQ
        CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
        IF (CRIT < TOLEV)                GOTO 2000
        CALL INVCOL3 (X2,Y2,BM2,NTOT2)
        XX     = GLSC2  (Y2,X2,NTOT2)
        XNORM  = 1./SQRT(XX)
        CALL CMULT   (X2,XNORM,NTOT2)
    1000 END DO
    2000 CONTINUE

    GAMMA = RQ
    RETURN
    END SUBROUTINE GAMMAM2
!-----------------------------------------------------------------------
    SUBROUTINE STARTX1 (X1,Y1,MASK,MULT,NEL)

!     Compute startvector for finding an eigenvalue on mesh 1.
!     Normalization: XT*B*X = 1

    use size_m
    INCLUDE 'MASS'

    REAL :: X1   (LX1,LY1,LZ1,1)
    REAL :: Y1   (LX1,LY1,LZ1,1)
    REAL :: MASK (LX1,LY1,LZ1,1)
    REAL :: MULT (LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NEL
    CALL COPY       (X1,BM1,NTOT1)


    call rand_fld_h1(y1)            ! pff 3/21/12
    small = 0.001*glamax(x1,ntot1)
    call add2s2(x1,y1,small,ntot1)


    CALL COL2       (X1,MASK,NTOT1)
    CALL COL3       (Y1,BM1,X1,NTOT1)
    CALL DSSUM      (Y1,NX1,NY1,NZ1)
    XX     = GLSC3 (X1,Y1,MULT,NTOT1)
    XNORM  = 1./SQRT(XX)
    CALL CMULT      (X1,XNORM,NTOT1)

    RETURN
    END SUBROUTINE STARTX1
!-----------------------------------------------------------------------
    SUBROUTINE STARTX2 (X2,Y2)
!------------------------------------------------------------------

!     Compute startvector for finding an eigenvalue on mesh 2.

!------------------------------------------------------------------
    use size_m
    INCLUDE 'MASS'

    REAL :: X2 (LX2,LY2,LZ2,LELV)
    REAL :: Y2 (LX2,LY2,LZ2,LELV)

    NXYZ2  = NX2*NY2*NZ2
    NTOT2  = NXYZ2*NELV
    ICONST = 0
    IF ((NDIM == 2) .AND. (NXYZ2 == 4)) ICONST = 1
    IF ((NDIM == 3) .AND. (NXYZ2 == 8)) ICONST = 1

    IF (ICONST == 1) THEN
        DO 1000 IEL=1,NELV
            DO 1000 K=1,NZ2
                DO 1000 J=1,NY2
                    DO 1000 I=1,NX2
                        X2(I,J,K,IEL) = I*J*K
        1000 END DO
    ELSE
        CALL COPY (X2,BM2,NTOT2)
    ENDIF

    call ortho (x2)
    CALL COL3 (Y2,BM2,X2,NTOT2)
    XX     = GLSC2 (X2,Y2,NTOT2)
    XNORM  = 1./SQRT(XX)
    CALL CMULT (X2,XNORM,NTOT2)

    RETURN
    END SUBROUTINE STARTX2
