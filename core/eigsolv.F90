!**********************************************************************
!> \file eigsolv.F90
!! \brief routines for estimating and calculating eigenvalues 
!**********************************************************************

!--------------------------------------------------------------
!> \brief Estimate eigenvalues
!-------------------------------------------------------------
SUBROUTINE ESTEIG
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, ndim, nid
  use eigen, only : eigae, eigge, eigga, eigaa, eigas, eiggs
  use geom, only : xm1, ym1, zm1
  use input, only : if3d, ifaxis, ifflow
  use tstep, only : nelfld, ifield, pi, istep
  implicit none

  integer :: ntot1
  real(DP) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(DP) :: xx, yy, zz, rxy, ryx, rxz, rzx, ryz, rzy, rmin
  real(DP) :: xx2, yy2, zz2, xyzmin, xyzmax
  real(DP) :: one, ratio
  real(DP), external :: glmin, glmax

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
!> \brief Compute startvector for finding an eigenvalue on mesh 1.
!! Normalization: XT*B*X = 1
SUBROUTINE STARTX1 (X1,Y1,MASK,MULT,NEL)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1
  use mass, only : bm1
  implicit none

  REAL(DP) :: X1   (LX1,LY1,LZ1,1)
  REAL(DP) :: Y1   (LX1,LY1,LZ1,1)
  REAL(DP) :: MASK (LX1,LY1,LZ1,1)
  REAL(DP) :: MULT (LX1,LY1,LZ1,1)
  integer :: nel

  integer :: ntot1
  real(DP) :: small, xx, xnorm
  real(DP), external :: glamax, glsc3

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
