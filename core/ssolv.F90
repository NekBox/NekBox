!-------------------------------------------------------------------
!> \brief  Set tolerances for velocity solver
!-------------------------------------------------------------------
SUBROUTINE SETTOLV
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use eigen, only : eigaa, eigas
  use input, only : iftran, ifnav
  use soln, only : vdiff, vtrans
  use tstep, only : ifield, nelfld, istep, tolabs, dt, vnrml8, vnrmsm 
  use tstep, only : tolps, tolrel, tolhv, tolhr, tolhs, tolpdf
  implicit none

  REAL(DP) :: avvisc, avdens, vnorm, factor
  integer :: ntot
  real(DP), external :: glmin, glmax

  NTOT   = NX1*NY1*NZ1*NELFLD(IFIELD)
  AVVISC = GLMIN(VDIFF(1,1,1,1,IFIELD),NTOT)
  AVDENS = GLMAX(VTRANS(1,1,1,1,IFIELD),NTOT)

  vnorm = 0.
  IF (IFTRAN) THEN
      IF (ISTEP == 1)  VNORM = VNRML8
      IF (ISTEP > 1)  VNORM = VNRMSM
      IF (VNORM == 0.) VNORM = TOLABS
      FACTOR = 1.+(AVDENS/(EIGAA*AVVISC*DT))
  ELSE
      VNORM = VNRML8
      IF (VNORM == 0.) VNORM = TOLABS
      FACTOR = 1.
  ENDIF

  TOLPS  = TOLREL*VNORM * SQRT(EIGAS)/(4.*FACTOR)
  TOLHV  = TOLREL*VNORM * SQRT(EIGAA)*AVVISC/2.
  TOLHV  = TOLHV/3.
  IF ( .NOT. IFTRAN .AND. .NOT. IFNAV) TOLHV = TOLHV/10.
  TOLHR  = TOLHV
  TOLHS  = TOLHV

!   Non-zero default pressure tolerance
!   NOTE: This tolerance may change due to precision problems.
!         See subroutine CHKSSV

  IF (TOLPDF /= 0.) TOLPS = TOLPDF

  RETURN
END SUBROUTINE SETTOLV

!-------------------------------------------------------------------
!> \brief Set tolerances for temerature/passive scalar solver
!-------------------------------------------------------------------
SUBROUTINE SETTOLT
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use eigen, only : eigaa
  use input, only : iftran
  use soln, only : vdiff
  use tstep, only : ifield, nelfld, istep, tnrml8, tnrmsm, tolabs, tolht, tolrel
  implicit none

  REAL(DP) :: avcond, tnorm
  integer :: ntot
  real(DP), external :: glmin

  NTOT   = NX1*NY1*NZ1*NELFLD(IFIELD)
  AVCOND = GLMIN (VDIFF(1,1,1,1,IFIELD),NTOT)

  tnorm = 0.
  IF (IFTRAN) THEN
      IF (ISTEP == 1)  TNORM = TNRML8(IFIELD-1)
      IF (ISTEP > 1)  TNORM = TNRMSM(IFIELD-1)
      IF (TNORM == 0.) TNORM = TOLABS
  ELSE
      TNORM = TNRML8(IFIELD-1)
      IF (TNORM == 0.) TNORM = TOLABS
  ENDIF

  TOLHT(IFIELD) = TOLREL*TNORM * SQRT(EIGAA)*AVCOND

  RETURN
END SUBROUTINE SETTOLT

!-----------------------------------------------------------------------
!> \brief If characteristics, need number of sub-timesteps (DT/DS).
!! Current sub-timeintegration scheme: RK4.
!! If not characteristics, i.e. standard semi-implicit scheme,
!! check user-defined Courant number.
!----------------------------------------------------------------------
SUBROUTINE SETCHAR
  use kinds, only : DP
  use size_m, only : nid
  use input, only : ifchar
  use tstep, only : ctarg, ntaubd
  implicit none

  integer :: ict
  real(DP) :: rict, dct

  IF (IFCHAR) THEN
      ICT    = INT(CTARG)
      RICT   = (ICT)
      DCT    = CTARG-RICT
      IF (DCT == 0.) NTAUBD = ICT
      IF (DCT > 0.) NTAUBD = ICT+1
  !        if (param(78).ne.0) then
  !             ntaupf=int(param(78))
  !             if (nid.eq.0) write(6,*) ' new ntaubd:',ntaubd,ntaupf
  !             ntaubd=max(ntaubd,ntaupf)
  !        endif
  ELSE
      NTAUBD = 0
      IF (CTARG > 0.5) THEN
          IF (NID == 0) &
          WRITE (6,*) 'Reset the target Courant number to .5'
          CTARG = 0.5
      ENDIF
  ENDIF

  RETURN
END SUBROUTINE SETCHAR
