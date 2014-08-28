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
    use tstep

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

!-------------------------------------------------------------------------
!> \brief Compute the following eigenvalues:.
!!  EIGAA  = minimum eigenvalue of the matrix A  (=Laplacian)
!!  EIGAE  = minimum eigenvalue of the matrix E  (=DB-1DT)
!!  EIGAS  = minimum eigenvalue of the matrix S  (=DA-1DT)
!!  EIGAST = minimum eigenvalue of the matrix St (=D(A+B/dt)-1DT
!!  EIGGA  = maximum eigenvalue of the matrix A
!!  EIGGS  = maximum eigenvalue of the matrix S
!!  EIGGE  = maximum eigenvalue of the matrix E
!!  EIGGST = maximum eigenvalue of the matrix St
!!  Method : Power method/Inverse iteration & Rayleigh quotient wo shift
!-------------------------------------------------------------------------
SUBROUTINE EIGENV
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, nelv, ndim
  use eigen, only : ifaa, ifas, ifae, ifast, ifgs, ifge, ifgst, ifga, eigga
  use input, only : ifstrs
  use soln, only : vmult, v1mask, v2mask, v3mask
  implicit none

  real(DP) :: H1(LX1,LY1,LZ1,LELT), H2(LX1,LY1,LZ1,LELT)
!    C!OMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
  integer :: ntot1
  real(DP) :: eigga1, eigga2, eigga3
  
  NTOT1  = NX1*NY1*NZ1*NELV

  IF (IFAA) THEN
     write(*,*) "Oops: IFAA"
#if 0
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
#endif
  ENDIF

  IF (IFAS) THEN
     write(*,*) "Oops: IFAA"
#if 0
      INLOC = 0
      CALL RONE    (H1,NTOT1)
      CALL RZERO   (H2,NTOT1)
      CALL RZERO   (H2INV,NTOT1)
      CALL ALPHAM2  (EIGAS,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAS = ',EIGAS
#endif
  ENDIF

  IF (IFAE) THEN
     write(*,*) "Oops: IFAE"
#if 0
      INLOC = 1
      CALL RZERO   (H1,NTOT1)
      CALL RONE    (H2,NTOT1)
      CALL RONE    (H2INV,NTOT1)
      CALL ALPHAM2  (EIGAE,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAE = ',EIGAE
#endif
  ENDIF

  IF (IFAST) THEN
     write(*,*) "Oops: IFAST"
#if 0
      INLOC = -1
      CALL SETHLM  (H1,H2,INLOC)
      CALL INVERS2 (H2INV,H2,NTOT1)
      CALL ALPHAM2  (EIGAST,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGAST = ',EIGAST
#endif
  ENDIF

  IF (IFGS) THEN
     write(*,*) "Oops: IFGS"
#if 0
      INLOC = 0
      CALL RONE    (H1,NTOT1)
      CALL RZERO   (H2,NTOT1)
      CALL RZERO   (H2INV,NTOT1)
      CALL GAMMAM2 (EIGGS,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGS = ',EIGGS
#endif
  ENDIF

  IF (IFGE) THEN
     write(*,*) "Oops: IFGE"
#if 0
      INLOC = 1
      CALL RZERO   (H1,NTOT1)
      CALL RONE    (H2,NTOT1)
      CALL RONE    (H2INV,NTOT1)
      CALL GAMMAM2 (EIGGE,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGE = ',EIGGE
#endif
  ENDIF

  IF (IFGST) THEN
     write(*,*) "Oops: IFGST"
#if 0
      INLOC = -1
      CALL SETHLM  (H1,H2,INLOC)
      CALL INVERS2 (H2INV,H2,NTOT1)
      CALL GAMMAM2 (EIGGST,H1,H2,H2INV,INLOC)
      IF (NID == 0 .AND. ISTEP <= 0) WRITE (6,*) 'EIGGST = ',EIGGST
#endif
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

!---------------------------------------------------------------------------
!> \brief Compute maximum eigenvalue of the discrete Helmholtz operator
!---------------------------------------------------------------------------
SUBROUTINE GAMMAM1 (GAMMA,MASK,MULT,H1,H2,ISD)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelt, nelv
  use mass, only : binvm1
  use tstep, only : imesh, nmxe, tolev
  implicit none

  real(DP) :: gamma
  REAL(DP) ::            MASK (LX1,LY1,LZ1,1)
  REAL(DP) ::            MULT (LX1,LY1,LZ1,1)
  REAL(DP) ::            H1   (LX1,LY1,LZ1,1)
  REAL(DP) ::            H2   (LX1,LY1,LZ1,1)
  integer :: isd

  real(DP) :: X1(LX1,LY1,LZ1,LELT), Y1(LX1,LY1,LZ1,LELT)

  integer :: nel, nxyz1, ntot1
  integer :: iter
  real(DP) :: evnew, rq, evold, crit, xx, xnorm
  real(DP), external :: glsc3

  IF (IMESH == 1) NEL = NELV
  IF (IMESH == 2) NEL = NELT
  NXYZ1  = NX1*NY1*NZ1
  NTOT1  = NXYZ1*NEL
  EVNEW  = 0.
!   pff (2/15/96)
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
!-----------------------------------------------------------------------
    SUBROUTINE STARTX1 (X1,Y1,MASK,MULT,NEL)

!     Compute startvector for finding an eigenvalue on mesh 1.
!     Normalization: XT*B*X = 1

    use size_m
    use mass

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
    use mass

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
