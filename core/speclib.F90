!==============================================================================
!> \file  speclib.F90
!! \brief LIBRARY ROUTINES FOR SPECTRAL METHODS
!! \date  March 1989
!! \author Einar Malvin Ronquit
!!
!!     ABBRIVIATIONS:
!!
!!     M   - Set of mesh points
!!     Z   - Set of collocation/quadrature points
!!     W   - Set of quadrature weights
!!     H   - Lagrangian interpolant
!!     D   - Derivative operator
!!     I   - Interpolation operator
!!     GL  - Gauss Legendre
!!     GLL - Gauss-Lobatto Legendre
!!     GJ  - Gauss Jacobi
!!     GLJ - Gauss-Lobatto Jacobi
!!
!!     MAIN ROUTINES:
!!
!!     Points and weights:
!!     ZWGL      Compute Gauss Legendre points and weights
!!     ZWGLL     Compute Gauss-Lobatto Legendre points and weights
!!     ZWGJ      Compute Gauss Jacobi points and weights (general)
!!     ZWGLJ     Compute Gauss-Lobatto Jacobi points and weights (general)
!!
!!     Lagrangian interpolants:
!!     HGL       Compute Gauss Legendre Lagrangian interpolant
!!     HGLL      Compute Gauss-Lobatto Legendre Lagrangian interpolant
!!     HGJ       Compute Gauss Jacobi Lagrangian interpolant (general)
!!     HGLJ      Compute Gauss-Lobatto Jacobi Lagrangian interpolant (general)
!!
!!     Derivative operators:
!!     DGLL      Compute Gauss-Lobatto Legendre derivative matrix
!!     DGLLGL    Compute derivative matrix for a staggered mesh (GLL->GL)
!!     DGJ       Compute Gauss Jacobi derivative matrix (general)
!!     DGLJ      Compute Gauss-Lobatto Jacobi derivative matrix (general)
!!     DGLJGJ    Compute derivative matrix for a staggered mesh (GLJ->GJ) (general)
!!
!!     Interpolation operators:
!!     IGLM      Compute interpolation operator GL  -> M
!!     IGLLM     Compute interpolation operator GLL -> M
!!     IGJM      Compute interpolation operator GJ  -> M  (general)
!!     IGLJM     Compute interpolation operator GLJ -> M  (general)
!!
!!     Other:
!!     PNLEG     Compute Legendre polynomial of degree N
!!     PNDLEG    Compute derivative of Legendre polynomial of degree N
!!
!!     Comments:
!!     Note that many of the above routines exist in both single and
!!     double precision. If the name of the single precision routine is
!!     SUB, the double precision version is called SUBD. In most cases
!!     all the "low-level" arithmetic is done in double precision, even
!!     for the single precsion versions.
!!
!!     Useful references:
!! [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
!!     Providence, Rhode Island, 1939.
!! [2] Abramowitz & Stegun: Handbook of Mathematical Functions,
!!     Dover, New York, 1972.
!! [3] Canuto, Hussaini, Quarteroni & Zang: Spectral Methods in Fluid
!!     Dynamics, Springer-Verlag, 1988.
!!==============================================================================

module speclib
  use kinds, only : DP

  integer, PARAMETER :: nmax = 84
  public zwgl, zwgll, igllm, dgll, dgllgl, zwglj
  private

contains

!--------------------------------------------------------------------
    SUBROUTINE ZWGL (Z,W,NP)
!--------------------------------------------------------------------
!> \brief Generate NP Gauss Legendre points (Z) and weights (W)
!!        associated with Jacobi polynomial P(N)(alpha=0,beta=0).
!! The polynomial degree N=NP-1.
!! Z and W are in single precision, but all the arithmetic
!! operations are done in double precision.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: NP
    REAL(DP), intent(out) :: Z(NP), W(NP)
    real(DP) :: alpha, beta
    ALPHA = 0._dp
    BETA  = 0._dp
    CALL ZWGJ (Z,W,NP,ALPHA,BETA)
    RETURN
    END SUBROUTINE ZWGL

    SUBROUTINE ZWGLL (Z,W,NP)
!--------------------------------------------------------------------
!     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
!     The polynomial degree N=NP-1.
!     Z and W are in single precision, but all the arithmetic
!     operations are done in double precision.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: NP
    REAL(DP), intent(out) :: Z(NP),W(NP)
    real(DP) :: alpha, beta
    ALPHA = 0._dp
    BETA  = 0._dp
    CALL ZWGLJ (Z,W,NP,ALPHA,BETA)
    RETURN
    END SUBROUTINE ZWGLL

    SUBROUTINE ZWGJ (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------

!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.

!--------------------------------------------------------------------
    implicit none
    integer, parameter :: NZD = NMAX

    integer,  intent(in)  :: NP
    REAL(DP), intent(out) :: Z(NP), W(NP)
    real(DP), intent(in)  :: ALPHA, BETA
    REAL(DP) :: ZD(NZD),WD(NZD),ALPHAD,BETAD
    integer :: npmax, i

    NPMAX = NZD
    IF (NP > NPMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in ZWGJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here NP=',NP
        call exitt
    ENDIF
    ALPHAD = ALPHA
    BETAD  = BETA
    CALL ZWGJD (ZD,WD,NP,ALPHAD,BETAD)
    DO 100 I=1,NP
        Z(I) = ZD(I)
        W(I) = WD(I)
    100 END DO
    RETURN
    END SUBROUTINE ZWGJ

    SUBROUTINE ZWGJD (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in)   :: NP
    REAL(DP), intent(out) :: Z(NP),W(NP)
    real(DP), intent(in)  :: ALPHA,BETA

    integer :: n, np1, np2, i
    real(DP) :: dn, one, two, apb, dnp1, dnp2, fac1, fac2, fac3, fnorm
    real(DP) :: rcoef, P, PD, PM1, PDM1, PM2, PDM2

    N     = NP-1
    DN    = ((N))
    ONE   = 1._dp
    TWO   = 2._dp
    APB   = ALPHA+BETA

    IF (NP <= 0) THEN
        WRITE (6,*) 'ZWGJD: Minimum number of Gauss points is 1',np
        call exitt
    ENDIF
    IF ((ALPHA <= -ONE) .OR. (BETA <= -ONE)) THEN
        WRITE (6,*) 'ZWGJD: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    IF (NP == 1) THEN
        Z(1) = (BETA-ALPHA)/(APB+TWO)
        W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO) &
        * TWO**(APB+ONE)
        RETURN
    ENDIF

    CALL JACG (Z,NP,ALPHA,BETA)

    NP1   = N+1
    NP2   = N+2
    DNP1  = ((NP1))
    DNP2  = ((NP2))
    FAC1  = DNP1+ALPHA+BETA+ONE
    FAC2  = FAC1+DNP1
    FAC3  = FAC2+ONE
    FNORM = PNORMJ(NP1,ALPHA,BETA)
    RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
    DO 100 I=1,NP
        CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z(I))
        W(I) = -RCOEF/(P*PDM1)
    100 END DO
    RETURN
    END SUBROUTINE ZWGJD

    SUBROUTINE ZWGLJ (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------

!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.

!--------------------------------------------------------------------
    implicit none
    integer, PARAMETER :: NZD = NMAX
    integer,  intent(in)  :: NP
    REAL(DP), intent(out) :: Z(NP),W(NP)
    real(DP), intent(in)  :: ALPHA,BETA
    REAL(DP) ::  ZD(NZD),WD(NZD),ALPHAD,BETAD

    integer :: npmax, i

    NPMAX = NZD
    IF (NP > NPMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in ZWGLJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here NP=',NP
        call exitt
    ENDIF
    ALPHAD = ALPHA
    BETAD  = BETA
    CALL ZWGLJD (ZD,WD,NP,ALPHAD,BETAD)
    DO 100 I=1,NP
        Z(I) = ZD(I)
        W(I) = WD(I)
    100 END DO
    RETURN
    END SUBROUTINE ZWGLJ

    SUBROUTINE ZWGLJD (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: np
    REAL(DP), intent(out) ::  Z(NP),W(NP)
    REAL(DP), intent(in) ::  ALPHA,BETA

    integer :: n, nm1, i
    real(DP) :: one, two, alpg, betg, p, pd, pm1, pdm1, pm2, pdm2
    N     = NP-1
    NM1   = N-1
    ONE   = 1._dp
    TWO   = 2._dp

    IF (NP <= 1) THEN
        WRITE (6,*) 'ZWGLJD: Minimum number of Gauss-Lobatto points is 2'
        WRITE (6,*) 'ZWGLJD: alpha,beta:',alpha,beta,np
        call exitt
    ENDIF
    IF ((ALPHA <= -ONE) .OR. (BETA <= -ONE)) THEN
        WRITE (6,*) 'ZWGLJD: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    IF (NM1 > 0) THEN
        ALPG  = ALPHA+ONE
        BETG  = BETA+ONE
        CALL ZWGJD (Z(2),W(2),NM1,ALPG,BETG)
    ENDIF
    Z(1)  = -ONE
    Z(NP) =  ONE
    DO 100  I=2,NP-1
        W(I) = W(I)/(ONE-Z(I)**2)
    100 END DO
    CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(1))
    W(1)  = ENDW1 (N,ALPHA,BETA)/(TWO*PD)
    CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(NP))
    W(NP) = ENDW2 (N,ALPHA,BETA)/(TWO*PD)

    RETURN
    END SUBROUTINE ZWGLJD

    REAL(DP) FUNCTION ENDW1 (N,ALPHA,BETA)
    implicit none
    integer, intent(in) :: n
    real(DP), intent(in) :: alpha, beta

    integer :: i
    real(DP) :: zero, one, two, three, four, apb
    real(DP) :: f1, f2, fint1, fint2, f3, a1, a2, a3, abn, abnn, di
    ZERO  = 0._dp
    ONE   = 1._dp
    TWO   = 2._dp
    THREE = 3._dp
    FOUR  = 4._dp
    F3    = -1._dp
    APB   = ALPHA+BETA
    IF (N == 0) THEN
        ENDW1 = ZERO
        RETURN
    ENDIF
    F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
    F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
    IF (N == 1) THEN
        ENDW1 = F1
        RETURN
    ENDIF
    FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
    FINT1 = FINT1*TWO**(APB+TWO)
    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
    FINT2 = FINT2*TWO**(APB+THREE)
    F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2) &
    * (APB+THREE)/FOUR
    IF (N == 2) THEN
        ENDW1 = F2
        RETURN
    ENDIF
    DO 100 I=3,N
        DI   = ((I-1))
        ABN  = ALPHA+BETA+DI
        ABNN = ABN+DI
        A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
        A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
        A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
        F3   =  -(A2*F2+A1*F1)/A3
        F1   = F2
        F2   = F3
    100 END DO
    ENDW1  = F3
    RETURN
    END FUNCTION

    REAL(DP) FUNCTION ENDW2 (N,ALPHA,BETA)
    implicit none
    integer, intent(in) :: n
    real(DP), intent(in) :: alpha, beta

    integer :: i
    real(DP) :: zero, one, two, three, four, apb
    real(DP) :: f1, f2, fint1, fint2, f3, a1, a2, a3, abn, abnn, di

    ZERO  = 0._dp
    ONE   = 1._dp
    TWO   = 2._dp
    THREE = 3._dp
    FOUR  = 4._dp
    APB   = ALPHA+BETA
    F3 = -1._dp
    IF (N == 0) THEN
        ENDW2 = ZERO
        RETURN
    ENDIF
    F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
    F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
    IF (N == 1) THEN
        ENDW2 = F1
        RETURN
    ENDIF
    FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
    FINT1 = FINT1*TWO**(APB+TWO)
    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
    FINT2 = FINT2*TWO**(APB+THREE)
    F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2) &
    * (APB+THREE)/FOUR
    IF (N == 2) THEN
        ENDW2 = F2
        RETURN
    ENDIF
    DO 100 I=3,N
        DI   = ((I-1))
        ABN  = ALPHA+BETA+DI
        ABNN = ABN+DI
        A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
        A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
        A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
        F3   =  -(A2*F2+A1*F1)/A3
        F1   = F2
        F2   = F3
    100 END DO
    ENDW2  = F3
    RETURN
    END FUNCTION

    REAL(DP) FUNCTION GAMMAF (X)
    implicit none
    REAL(DP), intent(in) :: X

    real(DP) :: zero, half, one, two, four, pi
    ZERO = 0.0_dp
    HALF = 0.5_dp
    ONE  = 1.0_dp
    TWO  = 2.0_dp
    FOUR = 4.0_dp
    PI   = FOUR*ATAN(ONE)
    GAMMAF = ONE
    IF (X == -HALF) GAMMAF = -TWO*SQRT(PI)
    IF (X == HALF) GAMMAF =  SQRT(PI)
    IF (X == ONE ) GAMMAF =  ONE
    IF (X == TWO ) GAMMAF =  ONE
    IF (X == 1.5_dp  ) GAMMAF =  SQRT(PI)/2._dp
    IF (X == 2.5_dp) GAMMAF =  1.5_dp*SQRT(PI)/2._dp
    IF (X == 3.5_dp) GAMMAF =  0.5_dp*(2.5_dp*(1.5_dp*SQRT(PI)))
    IF (X == 3._dp ) GAMMAF =  2._dp
    IF (X == 4._dp ) GAMMAF = 6._dp
    IF (X == 5._dp ) GAMMAF = 24._dp
    IF (X == 6._dp ) GAMMAF = 120._dp
    RETURN
    END FUNCTION

    REAL(DP)  FUNCTION PNORMJ (N,ALPHA,BETA)
    implicit none
    integer, intent(in) :: n
    REAL(DP), intent(in) ::  ALPHA,BETA

    real(DP) :: one, two, dn, const, prod, dindx, frac
    integer :: i

    ONE   = 1._dp
    TWO   = 2._dp
    DN    = ((N))
    CONST = ALPHA+BETA+ONE
    IF (N <= 1) THEN
        PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
        PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
        PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
        RETURN
    ENDIF
    PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
    PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
    PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
    PROD  = PROD*(ONE+BETA)*(TWO+BETA)
    DO 100 I=3,N
        DINDX = ((I))
        FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
        PROD  = PROD*FRAC
    100 END DO
    PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
    RETURN
    END FUNCTION


    SUBROUTINE JACG (XJAC,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!     Compute NP Gauss points XJAC, which are the zeros of the
!     Jacobi polynomial J(NP) with parameters ALPHA and BETA.
!     ALPHA and BETA determines the specific type of Gauss points.
!     Examples:
!     ALPHA = BETA =  0.0  ->  Legendre points
!     ALPHA = BETA = -0.5  ->  Chebyshev points
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: np
    REAL(DP), intent(out) ::  XJAC(NP)
    real(DP), intent(in) :: alpha, beta

    real(DP) :: eps, one, dth, x, xlast, x1, x2, swap, xmin, delx
    real(DP) :: p, pd, pm1, pdm1, pm2, pdm2, recsum
    integer :: n, kstop, j, k, i, jmin, jm
    DATA KSTOP /10/

    kstop = 10
    eps = 1.e-12_dp

    N   = NP-1
    one = 1._dp
    DTH = 4._dp*ATAN(one)/(2._dp*((N))+2._dp)
    JMIN = -1
    DO 40 J=1,NP
        IF (J == 1) THEN
            X = COS((2._dp*(((J))-1._dp)+1._dp)*DTH)
        ELSE
            X1 = COS((2._dp*(((J))-1._dp)+1._dp)*DTH)
            X2 = XLAST
            X  = (X1+X2)/2._dp
        ENDIF
        DO 30 K=1,KSTOP
            CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X)
            RECSUM = 0._dp
            JM = J-1
            DO 29 I=1,JM
                RECSUM = RECSUM+1./(X-XJAC(NP-I+1))
            29 END DO
            DELX = -P/(PD-RECSUM*P)
            X    = X+DELX
            IF (ABS(DELX) < EPS) GOTO 31
        30 END DO
        31 CONTINUE
        XJAC(NP-J+1) = X
        XLAST        = X
    40 END DO
    DO 200 I=1,NP
        XMIN = 2._dp
        DO 100 J=I,NP
            IF (XJAC(J) < XMIN) THEN
                XMIN = XJAC(J)
                JMIN = J
            ENDIF
        100 END DO
        IF (JMIN /= I) THEN
            SWAP = XJAC(I)
            XJAC(I) = XJAC(JMIN)
            XJAC(JMIN) = SWAP
        ENDIF
    200 END DO
    RETURN
    END SUBROUTINE JACG

    SUBROUTINE JACOBF (POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2, &
    N,ALP,BET,X)
!--------------------------------------------------------------------

!     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
!     of degree N at X.

!--------------------------------------------------------------------
    implicit none
    real(DP), intent(out) :: poly, pder, polym1, pderm1, polym2, pderm2
    real(DP), intent(in) :: alp, bet, x
    integer, intent(in) :: n

    real(DP) :: apb, polyl, pderl, dk, a1, a2, a3, b3, a4, polyn, pdern
    real(DP) :: psave, pdsave
    integer :: k

    APB  = ALP+BET
    POLY = 1._dp
    PDER = 0._dp
    IF (N == 0) RETURN
    POLYL = POLY
    PDERL = PDER
    POLY  = (ALP-BET+(APB+2._dp)*X)/2._dp
    PDER  = (APB+2._dp)/2._dp
    PSAVE = 0._dp
    PDSAVE = 0._dp
    IF (N == 1) RETURN
    DO 20 K=2,N
        DK = ((K))
        A1 = 2._dp*DK*(DK+APB)*(2._dp*DK+APB-2._dp)
        A2 = (2._dp*DK+APB-1._dp)*(ALP**2-BET**2)
        B3 = (2._dp*DK+APB-2._dp)
        A3 = B3*(B3+1._dp)*(B3+2._dp)
        A4 = 2._dp*(DK+ALP-1._dp)*(DK+BET-1._dp)*(2._dp*DK+APB)
        POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
        PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
        PSAVE  = POLYL
        PDSAVE = PDERL
        POLYL  = POLY
        POLY   = POLYN
        PDERL  = PDER
        PDER   = PDERN
    20 END DO
    POLYM1 = POLYL
    PDERM1 = PDERL
    POLYM2 = PSAVE
    PDERM2 = PDSAVE
    RETURN
    END SUBROUTINE JACOBF

    REAL(DP) FUNCTION HGJ (II,Z,ZGJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant HGJ through
!     the NP Gauss Jacobi points ZGJ at the point Z.
!     Single precision version.
!---------------------------------------------------------------------
    implicit none
    integer, PARAMETER :: NZD = NMAX
    integer, intent(in) :: np, ii
    REAL(DP), intent(in) :: Z,ZGJ(np),ALPHA,BETA

    REAL(DP) ::  ZD,ZGJD(NZD),ALPHAD,BETAD
    integer :: i, npmax

    NPMAX = NZD
    IF (NP > NPMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in HGJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here NP=',NP
        call exitt
    ENDIF
    ZD = Z
    DO 100 I=1,NP
        ZGJD(I) = ZGJ(I)
    100 END DO
    ALPHAD = ALPHA
    BETAD  = BETA
    HGJ    = HGJD (II,ZD,ZGJD,NP,ALPHAD,BETAD)
    RETURN
    END FUNCTION HGJ

    REAL(DP) FUNCTION HGJD (II,Z,ZGJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant HGJD through
!     the NZ Gauss-Lobatto Jacobi points ZGJ at the point Z.
!     Double precision version.
!---------------------------------------------------------------------
    implicit none
    integer, intent(in) :: ii, np
    REAL(DP), intent(in) :: Z,ZGJ(np),ALPHA,BETA

    real(DP) :: eps, one, zi, dz
    real(DP) :: pzi, pdzi, pm1, pdm1, pm2, pdm2, pz, pdz

    EPS = 1.e-5_dp
    ONE = 1._dp
    ZI  = ZGJ(II)
    DZ  = Z-ZI
    IF (ABS(DZ) < EPS) THEN
        HGJD = ONE
        RETURN
    ENDIF
    CALL JACOBF (PZI,PDZI,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,ZI)
    CALL JACOBF (PZ,PDZ,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,Z)
    HGJD  = PZ/(PDZI*(Z-ZI))
    RETURN
    END FUNCTION

    REAL(DP) FUNCTION HGLJ (II,Z,ZGLJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant HGLJ through
!     the NZ Gauss-Lobatto Jacobi points ZGLJ at the point Z.
!     Single precision version.
!---------------------------------------------------------------------
    implicit none
    integer, PARAMETER :: NZD = NMAX
    integer,  intent(in) :: ii, NP
    REAL(DP), intent(in) :: Z,ZGLJ(NP)
    REAL(DP), intent(in) :: ALPHA,BETA

    REAL(DP) :: ZD,ZGLJD(NZD),ALPHAD,BETAD
    integer :: npmax, i

    NPMAX = NZD
    IF (NP > NPMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in HGLJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here NP=',NP
        call exitt
    ENDIF
    ZD = Z
    DO 100 I=1,NP
        ZGLJD(I) = ZGLJ(I)
    100 END DO
    ALPHAD = ALPHA
    BETAD  = BETA
    HGLJ   = HGLJD (II,ZD,ZGLJD,NP,ALPHAD,BETAD)
    RETURN
    END FUNCTION HGLJ

    REAL(DP) FUNCTION HGLJD (I,Z,ZGLJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant HGLJD through
!     the NZ Gauss-Lobatto Jacobi points ZJACL at the point Z.
!     Double precision version.
!---------------------------------------------------------------------
    implicit none
    integer,  intent(in) :: i, np
    REAL(DP), intent(in) ::  Z,ZGLJ(NP),ALPHA,BETA

    real(DP) :: eps, one, dz, zi, dn, eigval, const
    real(DP) :: p, pd, pi, pdi, pm1, pdm1, pm2, pdm2
    integer :: n

    EPS = 1.e-5_dp
    ONE = 1._dp
    ZI  = ZGLJ(I)
    DZ  = Z-ZI
    IF (ABS(DZ) < EPS) THEN
        HGLJD = ONE
        RETURN
    ENDIF
    N      = NP-1
    DN     = ((N))
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
    CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZI)
    CONST  = EIGVAL*PI+ALPHA*(ONE+ZI)*PDI-BETA*(ONE-ZI)*PDI
    CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z)
    HGLJD  = (ONE-Z**2)*PD/(CONST*(Z-ZI))
    RETURN
    END FUNCTION

    SUBROUTINE DGJ (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!-----------------------------------------------------------------
    implicit none
    integer, PARAMETER :: NZDD = NMAX
    integer :: nz, nzd
    REAL(DP), intent(out) :: D(NZD,NZD),DT(NZD,NZD)
    real(DP), intent(in) :: Z(nz),ALPHA,BETA

    REAL(DP) ::  DD(NZDD,NZDD),DTD(NZDD,NZDD),ZD(NZDD),ALPHAD,BETAD
    integer :: i, j

    IF (NZ <= 0) THEN
        WRITE (6,*) 'DGJ: Minimum number of Gauss points is 1'
        call exitt
    ENDIF
    IF (NZ > NMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in DGJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here Nz=',Nz
        call exitt
    ENDIF
    IF ((ALPHA <= -1.) .OR. (BETA <= -1.)) THEN
        WRITE (6,*) 'DGJ: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF
    ALPHAD = ALPHA
    BETAD  = BETA
    DO 100 I=1,NZ
        ZD(I) = Z(I)
    100 END DO
    CALL DGJD (DD,DTD,ZD,NZ,NZDD,ALPHAD,BETAD)
    DO 200 I=1,NZ
        DO 200 J=1,NZ
            D(I,J)  = DD(I,J)
            DT(I,J) = DTD(I,J)
    200 END DO
    RETURN
    END SUBROUTINE DGJ

    SUBROUTINE DGJD (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!-----------------------------------------------------------------
    implicit none
    integer,  intent(in)  :: nz, nzd
    REAL(DP), intent(out) :: D(NZD,NZD),DT(NZD,NZD)
    real(DP), intent(in)  :: Z(NZ),ALPHA,BETA

    integer :: n, i, j
    real(DP) :: dn, one, two
    real(DP) :: pi, pdi, pm1, pdm1, pm2, pdm2, pj, pdj

    N    = NZ-1
    DN   = ((N))
    ONE  = 1._dp
    TWO  = 2._dp

    IF (NZ <= 1) THEN
        WRITE (6,*) 'DGJD: Minimum number of Gauss-Lobatto points is 2'
        call exitt
    ENDIF
    IF ((ALPHA <= -ONE) .OR. (BETA <= -ONE)) THEN
        WRITE (6,*) 'DGJD: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    DO 200 I=1,NZ
        DO 200 J=1,NZ
            CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z(I))
            CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z(J))
            IF (I /= J) D(I,J) = PDI/(PDJ*(Z(I)-Z(J)))
            IF (I == J) D(I,J) = ((ALPHA+BETA+TWO)*Z(I)+ALPHA-BETA)/ &
            (TWO*(ONE-Z(I)**2))
            DT(J,I) = D(I,J)
    200 END DO
    RETURN
    END SUBROUTINE DGJD

    SUBROUTINE DGLJ (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!-----------------------------------------------------------------
    implicit none
    integer, PARAMETER :: NZDD = NMAX
    integer, intent(in) :: nz, nzd
    REAL(DP), intent(out) :: D(NZD,NZD),DT(NZD,NZD)
    REAL(DP), intent(in) :: Z(NZ),ALPHA,BETA

    REAL(DP) :: DD(NZDD,NZDD),DTD(NZDD,NZDD),ZD(NZDD),ALPHAD,BETAD
    integer :: i, j

    IF (NZ <= 1) THEN
        WRITE (6,*) 'DGLJ: Minimum number of Gauss-Lobatto points is 2'
        call exitt
    ENDIF
    IF (NZ > NMAX) THEN
        WRITE (6,*) 'Too large polynomial degree in DGLJ'
        WRITE (6,*) 'Maximum polynomial degree is',NMAX
        WRITE (6,*) 'Here NZ=',NZ
        call exitt
    ENDIF
    IF ((ALPHA <= -1.) .OR. (BETA <= -1.)) THEN
        WRITE (6,*) 'DGLJ: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF
    ALPHAD = ALPHA
    BETAD  = BETA
    DO 100 I=1,NZ
        ZD(I) = Z(I)
    100 END DO
    CALL DGLJD (DD,DTD,ZD,NZ,NZDD,ALPHAD,BETAD)
    DO 200 I=1,NZ
        DO 200 J=1,NZ
            D(I,J)  = DD(I,J)
            DT(I,J) = DTD(I,J)
    200 END DO
    RETURN
    END SUBROUTINE DGLJ

    SUBROUTINE DGLJD (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!-----------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz, nzd
    REAL(DP), intent(out) ::  D(NZD,NZD),DT(NZD,NZD)
    REAL(DP), intent(in) ::  Z(NZ),ALPHA,BETA

    integer :: n, i, j
    real(DP) :: dn, one, two, eigval
    real(DP) :: pi, pdi, pm1, pdm1, pm2, pdm2, pj, pdj, ci, cj
    N    = NZ-1
    DN   = ((N))
    ONE  = 1._dp
    TWO  = 2._dp
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

    IF (NZ <= 1) THEN
        WRITE (6,*) 'DGLJD: Minimum number of Gauss-Lobatto points is 2'
        call exitt
    ENDIF
    IF ((ALPHA <= -ONE) .OR. (BETA <= -ONE)) THEN
        WRITE (6,*) 'DGLJD: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    DO 200 I=1,NZ
        DO 200 J=1,NZ
            CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(I))
            CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(J))
            CI = EIGVAL*PI-(BETA*(ONE-Z(I))-ALPHA*(ONE+Z(I)))*PDI
            CJ = EIGVAL*PJ-(BETA*(ONE-Z(J))-ALPHA*(ONE+Z(J)))*PDJ
            IF (I /= J) D(I,J) = CI/(CJ*(Z(I)-Z(J)))
            IF ((I == J) .AND. (I /= 1) .AND. (I /= NZ)) &
            D(I,J) = (ALPHA*(ONE+Z(I))-BETA*(ONE-Z(I)))/ &
            (TWO*(ONE-Z(I)**2))
            IF ((I == J) .AND. (I == 1)) &
            D(I,J) =  (EIGVAL+ALPHA)/(TWO*(BETA+TWO))
            IF ((I == J) .AND. (I == NZ)) &
            D(I,J) = -(EIGVAL+BETA)/(TWO*(ALPHA+TWO))
            DT(J,I) = D(I,J)
    200 END DO
    RETURN
    END SUBROUTINE DGLJD

    SUBROUTINE DGLL (D,DT,Z,NZ,NZD)
!-----------------------------------------------------------------
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Legendre points Z.
!     Note: D and DT are square matrices.
!-----------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz, nzd 
    REAL(DP), intent(out) :: D(NZD,NZD),DT(NZD,NZD)
    REAL(DP), intent(in) :: Z(NZ)

    integer :: n, i, j
    real(DP) :: fn, d0
    N  = NZ-1
    IF (NZ > NMAX) THEN
        WRITE (6,*) 'Subroutine DGLL'
        WRITE (6,*) 'Maximum polynomial degree =',NMAX
        WRITE (6,*) 'Polynomial degree         =',NZ
    ENDIF
    IF (NZ == 1) THEN
        D(1,1) = 0._dp
        RETURN
    ENDIF
    FN = (N)
    d0 = FN*(FN+1._dp)/4._dp
    DO 200 I=1,NZ
        DO 200 J=1,NZ
            D(I,J) = 0._dp
            IF  (I /= J) D(I,J) = PNLEG(Z(I),N)/ &
            (PNLEG(Z(J),N)*(Z(I)-Z(J)))
            IF ((I == J) .AND. (I == 1))  D(I,J) = -d0
            IF ((I == J) .AND. (I == NZ)) D(I,J) =  d0
            DT(J,I) = D(I,J)
    200 END DO
    RETURN
    END SUBROUTINE DGLL

    REAL(DP) FUNCTION HGLL (I,Z,ZGLL,NZ)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant L through
!     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
!---------------------------------------------------------------------
    implicit none
    integer, intent(in) :: i, nz
    REAL(DP), intent(in) :: z, ZGLL(nz)

    integer :: n
    real(DP) :: eps, dz, alfan

    EPS = 1.E-5_dp
    DZ = Z - ZGLL(I)
    IF (ABS(DZ) < EPS) THEN
        HGLL = 1._dp
        RETURN
    ENDIF
    N = NZ - 1
    ALFAN = (N)*((N)+1._dp)
    HGLL = - (1._dp-Z*Z)*PNDLEG(Z,N)/ &
    (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))
    RETURN
    END FUNCTION HGLL

    REAL(DP) FUNCTION HGL (I,Z,ZGL,NZ)
!---------------------------------------------------------------------
!     Compute the value of the Lagrangian interpolant HGL through
!     the NZ Gauss Legendre points ZGL at the point Z.
!---------------------------------------------------------------------
    implicit none
    integer, intent(in) :: i, nz
    REAL(DP), intent(in) :: z, ZGL(nz)

    real(DP) :: eps, dz
    integer :: n

    EPS = 1.E-5_dp
    DZ = Z - ZGL(I)
    IF (ABS(DZ) < EPS) THEN
        HGL = 1._dp
        RETURN
    ENDIF
    N = NZ-1
    HGL = PNLEG(Z,NZ)/(PNDLEG(ZGL(I),NZ)*(Z-ZGL(I)))
    RETURN
    END FUNCTION HGL

    REAL(DP) FUNCTION PNLEG (Z,N)
!---------------------------------------------------------------------
!     Compute the value of the Nth order Legendre polynomial at Z.
!     (Simpler than JACOBF)
!     Based on the recursion formula for the Legendre polynomials.
!---------------------------------------------------------------------
    implicit none
    real(DP), intent(in) :: z
    integer,  intent(in) :: n

    real(DP) :: p1, p2, p3, fk
    integer :: k

    P1   = 1._dp
    IF (N == 0) THEN
        PNLEG = P1
        RETURN
    ENDIF
    P2   = Z
    P3   = P2
    DO 10 K = 1, N-1
        FK  = (K)
        P3  = ((2._dp*FK+1._dp)*Z*P2 - FK*P1)/(FK+1._dp)
        P1  = P2
        P2  = P3
    10 END DO
    PNLEG = P3
    if (n == 0) pnleg = 1._dp
    RETURN
    END FUNCTION PNLEG

    REAL(DP) FUNCTION PNDLEG (Z,N)
!----------------------------------------------------------------------
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     (Simpler than JACOBF)
!     Based on the recursion formula for the Legendre polynomials.
!----------------------------------------------------------------------
    implicit none
    real(DP), intent(in) :: z
    integer,  intent(in) :: n

    real(DP) :: p1, p2, p3, p1d, p2d, p3d
    integer :: k, fk
    P1   = 1._dp
    P2   = Z
    P1D  = 0._dp
    P2D  = 1._dp
    P3D  = 1._dp
    DO 10 K = 1, N-1
        FK  = (K)
        P3  = ((2._dp*FK+1._dp)*Z*P2 - FK*P1)/(FK+1._dp)
        P3D = ((2._dp*FK+1._dp)*P2 + (2._dp*FK+1._dp)*Z*P2D - FK*P1D)/(FK+1._dp)
        P1  = P2
        P2  = P3
        P1D = P2D
        P2D = P3D
    10 END DO
    PNDLEG = P3D
    IF (N == 0) pndleg = 0._dp
    RETURN
    END FUNCTION PNDLEG

    SUBROUTINE DGLLGL (D,DT,ZM1,ZM2,IM12,NZM1,NZM2,ND1,ND2)
!-----------------------------------------------------------------------
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Legendre mesh (M1), and evaluate its
!     derivative on a Guass Legendre mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLLGL).
!     Note: D and DT are rectangular matrices.
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nd1, nd2, nzm1, nzm2
    REAL(DP), intent(out) :: D(ND2,ND1), DT(ND1,ND2)
    REAL(DP), intent(in)  :: ZM1(ND1), ZM2(ND2), IM12(ND2,ND1)

    integer :: ip, jq, nm1
    real(DP) :: zp, zq, eps

    IF (NZM1 == 1) THEN
        D (1,1) = 0._dp
        DT(1,1) = 0._dp
        RETURN
    ENDIF
    EPS = 1.E-6_dp
    NM1 = NZM1-1
    DO 10 IP = 1, NZM2
        DO 10 JQ = 1, NZM1
            ZP = ZM2(IP)
            ZQ = ZM1(JQ)
            IF ((ABS(ZP) < EPS) .AND. (ABS(ZQ) < EPS)) THEN
                D(IP,JQ) = 0._dp
            ELSE
                D(IP,JQ) = (PNLEG(ZP,NM1)/PNLEG(ZQ,NM1) &
                -IM12(IP,JQ))/(ZP-ZQ)
            ENDIF
            DT(JQ,IP) = D(IP,JQ)
    10 END DO
    RETURN
    END SUBROUTINE DGLLGL

    SUBROUTINE DGLJGJ (D,DT,ZGL,ZG,IGLG,NPGL,NPG,ND1,ND2,ALPHA,BETA)
!-----------------------------------------------------------------------
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
!     derivative on a Guass Jacobi mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLJGJ).
!     Note: D and DT are rectangular matrices.
!     Single precision version.
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nd1, nd2, npgl, npg
    REAL(DP), intent(out) :: D(ND2,ND1), DT(ND1,ND2)
    REAL(DP), intent(in)  :: ZGL(ND1), ZG(ND2), IGLG(ND2,ND1)
    real(DP), intent(in)  :: alpha, beta

    integer, PARAMETER :: NDD = NMAX
    REAL(DP) ::  DD(NDD,NDD), DTD(NDD,NDD)
    REAL(DP) ::  ZGD(NDD), ZGLD(NDD), IGLGD(NDD,NDD)
    REAL(DP) ::  ALPHAD, BETAD
    integer :: i, j

    IF (NPGL <= 1) THEN
        WRITE(6,*) 'DGLJGJ: Minimum number of Gauss-Lobatto points is 2'
        call exitt
    ENDIF
    IF (NPGL > NMAX) THEN
        WRITE(6,*) 'Polynomial degree too high in DGLJGJ'
        WRITE(6,*) 'Maximum polynomial degree is',NMAX
        WRITE(6,*) 'Here NPGL=',NPGL
        call exitt
    ENDIF
    IF ((ALPHA <= -1._dp) .OR. (BETA <= -1._dp)) THEN
        WRITE(6,*) 'DGLJGJ: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    ALPHAD = ALPHA
    BETAD  = BETA
    DO 100 I=1,NPG
        ZGD(I) = ZG(I)
        DO 100 J=1,NPGL
            IGLGD(I,J) = IGLG(I,J)
    100 END DO
    DO 200 I=1,NPGL
        ZGLD(I) = ZGL(I)
    200 END DO
    CALL DGLJGJD (DD,DTD,ZGLD,ZGD,IGLGD,NPGL,NPG,NDD,NDD,ALPHAD,BETAD)
    DO 300 I=1,NPG
        DO 300 J=1,NPGL
            D(I,J)  = DD(I,J)
            DT(J,I) = DTD(J,I)
    300 END DO
    RETURN
    END SUBROUTINE DGLJGJ

    SUBROUTINE DGLJGJD (D,DT,ZGL,ZG,IGLG,NPGL,NPG,ND1,ND2,ALPHA,BETA)
!-----------------------------------------------------------------------
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
!     derivative on a Guass Jacobi mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLJGJ).
!     Note: D and DT are rectangular matrices.
!     Double precision version.
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nd1, nd2, npgl, npg
    REAL(DP), intent(out) ::  D(ND2,ND1), DT(ND1,ND2)
    REAL(DP), intent(in) ::  ZGL(ND1), ZG(ND2), IGLG(ND2,ND1), ALPHA, BETA

    real(DP) :: eps, one, two, faci, facj, const, dz, dn, eigval
    real(DP) :: pj, pdj, pm1, pdm1, pm2, pdm2, pi, pdi
    integer :: i, j, ngl
    IF (NPGL <= 1) THEN
        WRITE(6,*) 'DGLJGJD: Minimum number of Gauss-Lobatto points is 2'
        call exitt
    ENDIF
    IF ((ALPHA <= -1.) .OR. (BETA <= -1.)) THEN
        WRITE(6,*) 'DGLJGJD: Alpha and Beta must be greater than -1'
        call exitt
    ENDIF

    EPS    = 1.e-6_dp
    ONE    = 1._dp
    TWO    = 2._dp
    NGL    = NPGL-1
    DN     = ((NGL))
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

    DO 100 I=1,NPG
        DO 100 J=1,NPGL
            DZ = ABS(ZG(I)-ZGL(J))
            IF (DZ < EPS) THEN
                D(I,J) = (ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I)))/ &
                (TWO*(ONE-ZG(I)**2))
            ELSE
                CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,NGL,ALPHA,BETA,ZG(I))
                CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,NGL,ALPHA,BETA,ZGL(J))
                FACI   = ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I))
                FACJ   = ALPHA*(ONE+ZGL(J))-BETA*(ONE-ZGL(J))
                CONST  = EIGVAL*PJ+FACJ*PDJ
                D(I,J) = ((EIGVAL*PI+FACI*PDI)*(ZG(I)-ZGL(J)) &
                -(ONE-ZG(I)**2)*PDI)/(CONST*(ZG(I)-ZGL(J))**2)
            ENDIF
            DT(J,I) = D(I,J)
    100 END DO
    RETURN
    END SUBROUTINE DGLJGJD

    SUBROUTINE IGLM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2)
!----------------------------------------------------------------------
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss Legendre points.
!     Z2 : NZ2 points on mesh M.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz1, nz2, nd1, nd2
    REAL(DP), intent(out) :: I12(ND2,ND1),IT12(ND1,ND2)
    REAL(DP), intent(in) :: Z1(ND1),Z2(ND2)

    integer :: i, j
    real(DP) :: zi

    IF (NZ1 == 1) THEN
        I12 (1,1) = 1._dp
        IT12(1,1) = 1._dp
        RETURN
    ENDIF
    DO 10 I=1,NZ2
        ZI = Z2(I)
        DO 10 J=1,NZ1
            I12 (I,J) = HGL(J,ZI,Z1,NZ1)
            IT12(J,I) = I12(I,J)
    10 END DO
    RETURN
    END SUBROUTINE IGLM

    SUBROUTINE IGLLM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2)
!----------------------------------------------------------------------
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss-Lobatto Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss-Lobatto Legendre points.
!     Z2 : NZ2 points on mesh M.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz1, nz2, nd1, nd2
    REAL(DP), intent(out) :: I12(ND2,ND1),IT12(ND1,ND2)
    REAL(DP), intent(in) :: Z1(ND1),Z2(ND2)

    integer :: i, j
    real(DP) :: zi

    IF (NZ1 == 1) THEN
        I12 (1,1) = 1._dp
        IT12(1,1) = 1._dp
        RETURN
    ENDIF
    DO 10 I=1,NZ2
        ZI = Z2(I)
        DO 10 J=1,NZ1
            I12 (I,J) = HGLL(J,ZI,Z1,NZ1)
            IT12(J,I) = I12(I,J)
    10 END DO
    RETURN
    END SUBROUTINE IGLLM

    SUBROUTINE IGJM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2,ALPHA,BETA)
!----------------------------------------------------------------------
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss Jacobi mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss Jacobi points.
!     Z2 : NZ2 points on mesh M.
!     Single precision version.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz1, nz2, nd1, nd2
    REAL(DP), intent(out) :: I12(ND2,ND1),IT12(ND1,ND2)
    REAL(DP), intent(in) :: Z1(ND1),Z2(ND2), ALPHA, BETA

    integer :: i, j
    real(DP) :: zi
    IF (NZ1 == 1) THEN
        I12 (1,1) = 1._dp
        IT12(1,1) = 1._dp
        RETURN
    ENDIF
    DO 10 I=1,NZ2
        ZI = Z2(I)
        DO 10 J=1,NZ1
            I12 (I,J) = HGJ(J,ZI,Z1,NZ1,ALPHA,BETA)
            IT12(J,I) = I12(I,J)
    10 END DO
    RETURN
    END SUBROUTINE IGJM

    SUBROUTINE IGLJM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2,ALPHA,BETA)
!----------------------------------------------------------------------
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss-Lobatto Jacobi mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss-Lobatto Jacobi points.
!     Z2 : NZ2 points on mesh M.
!     Single precision version.
!--------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nz1, nz2, nd1, nd2
    REAL(DP), intent(out) :: I12(ND2,ND1),IT12(ND1,ND2)
    REAL(DP), intent(in) :: Z1(ND1),Z2(ND2), alpha, beta

    integer :: i, j
    real(DP) :: zi
    IF (NZ1 == 1) THEN
        I12 (1,1) = 1._dp
        IT12(1,1) = 1._dp
        RETURN
    ENDIF
    DO 10 I=1,NZ2
        ZI = Z2(I)
        DO 10 J=1,NZ1
            I12 (I,J) = HGLJ(J,ZI,Z1,NZ1,ALPHA,BETA)
            IT12(J,I) = I12(I,J)
    10 END DO
    RETURN
    END SUBROUTINE IGLJM

end module speclib
