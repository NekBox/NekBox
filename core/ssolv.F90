    SUBROUTINE SETTOLV
!-------------------------------------------------------------------
!     Set tolerances for velocity solver
!-------------------------------------------------------------------
    use size_m
    use eigen
    use input
    use mass
    use soln
    use tstep
    REAL :: LENGTH

    NTOT   = NX1*NY1*NZ1*NELFLD(IFIELD)
    AVVISC = GLMIN(VDIFF(1,1,1,1,IFIELD),NTOT)
    AVDENS = GLMAX(VTRANS(1,1,1,1,IFIELD),NTOT)

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

!     Non-zero default pressure tolerance
!     NOTE: This tolerance may change due to precision problems.
!           See subroutine CHKSSV

    IF (TOLPDF /= 0.) TOLPS = TOLPDF

    RETURN
    END SUBROUTINE SETTOLV

    SUBROUTINE SETTOLT
!-------------------------------------------------------------------

!     Set tolerances for temerature/passive scalar solver

!-------------------------------------------------------------------
    use size_m
    use eigen
    use input
    use mass
    use soln
    use tstep
    REAL :: LENGTH

    NTOT   = NX1*NY1*NZ1*NELFLD(IFIELD)
    AVCOND = GLMIN (VDIFF(1,1,1,1,IFIELD),NTOT)

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

    SUBROUTINE SETCHAR
!-----------------------------------------------------------------------

!     If characteristics, need number of sub-timesteps (DT/DS).
!     Current sub-timeintegration scheme: RK4.
!     If not characteristics, i.e. standard semi-implicit scheme,
!     check user-defined Courant number.

!----------------------------------------------------------------------
    use size_m
    use input
    use tstep

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
