    SUBROUTINE SSTEST (ISSS)
!------------------------------------------------------------------------

!     Test if Steady State Solver should be activated.

!------------------------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    ISSS = 0
    IADV = 0
    DO 100 IFIELD=1,NFIELD
        IF (IFADVC(IFIELD)) IADV = 1
    100 END DO
    IF ( .NOT. IFTRAN .AND. (IADV == 1)) ISSS = 1
    IF (ISSS == 1 .AND. NFIELD > 4 .AND. NID == 0) THEN
        WRITE (6,*) ' '
        WRITE (6,*) 'Trying to activate the steady state solver'
        WRITE (6,*) 'using NFIELD =',NFIELD
        WRITE (6,*) 'Maximum number of fields is 4'
        call exitt
    ENDIF
    RETURN
    END SUBROUTINE SSTEST

    SUBROUTINE SSINIT (KMAX)
!-----------------------------------------------------------

!     Initialize steady state solver

!-----------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'EIGEN'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'

    IFTRAN    = .TRUE. 
    IFCHAR    = .TRUE. 
    CTARG     = 3.
    IF (IFMODEL) CTARG = 2.
    IF (NX1 >= 10) THEN
        CTARG = 5.
        IF (IFMODEL) CTARG = 3.
    ENDIF
    CALL SETPROP
    TAUMIN = 1.E20
    MFIELD = 1
    IF ( .NOT. IFFLOW) MFIELD=2
    DO 10 IFIELD=MFIELD,NFIELD
        DIFFUS = AVDIFF(IFIELD)/AVTRAN(IFIELD)
        TAUSS(IFIELD)  = 1./(EIGAA*DIFFUS)
        TXNEXT(IFIELD) = TAUSS(IFIELD)
        IF (TAUSS(IFIELD) < TAUMIN) TAUMIN = TAUSS(IFIELD)
    10 END DO

    NBDINP    = 1.
    TIME      = 0.
    DT        = 0.
    DTINIT    = TAUMIN/5.
    NSTEPS    = 10000
    IOSTEP    = 10000

    IFMODP = .TRUE. 
    IF (IFNATC)  IFMODP = .FALSE. 
    IF (IFMODEL) IFMODP = .FALSE. 
    IFSKIP = .TRUE. 
    NSSKIP = 1

    PRELAX = 1.E-1
    IF (IFSPLIT) PRELAX = 1.E-4

    IFSSVT = .FALSE. 
    IFEXVT = .FALSE. 

    KMAX   = 5

    CALL SETCHAR


    RETURN
    END SUBROUTINE SSINIT

    SUBROUTINE CHKEXT (IFACCX,Z,S)
!------------------------------------------------------------------

!     Accept extrapolation?

!------------------------------------------------------------------
    use size_m
    INCLUDE 'TSTEP'
    INCLUDE 'INPUT'
    INCLUDE 'STEADY'
    LOGICAL :: IFACCX
    real*8 :: Z(1),S(1)
    REAL :: H1NRM1 (LDIMT1), H1NRM2(LDIMT1)

    CALL RZERO (H1NRM1,NFIELD)
    IF (IFFLOW) THEN
        IFIELD = 1
        CALL UNORM
        H1NRM1(IFIELD) = VNRMH1
    ENDIF
    DO 10 IFIELD=2,NFIELD
        CALL UNORM
        H1NRM1(IFIELD) = TNRMH1(IFIELD-1)
    10 END DO

    CALL MKVEC (Z)
    CALL MKARR (S)

    CALL RZERO (H1NRM2,NFIELD)
    IF (IFFLOW) THEN
        IFIELD = 1
        CALL UNORM
        H1NRM2(IFIELD) = VNRMH1
    ENDIF
    DO 20 IFIELD=2,NFIELD
        CALL UNORM
        H1NRM2(IFIELD) = TNRMH1(IFIELD-1)
    20 END DO

    XLIM   = .2
    IFACCX = .TRUE. 
    RDMAX  = 0.
    RDLIM  = .5*TOLREL+1.E-4
    IF (IFFLOW) THEN
        IFIELD = 1
        RDIFF = ABS((H1NRM2(IFIELD)-H1NRM1(IFIELD))/H1NRM1(IFIELD))
        IF (RDIFF > RDMAX) RDMAX = RDIFF
        IF (NID == 0) WRITE (6,*) ' ifield, rdiff ',ifield,rdiff
        IF (RDIFF > XLIM) IFACCX = .FALSE. 
    ENDIF
    DO 100 IFIELD=2,NFIELD
        RDIFF = ABS((H1NRM2(IFIELD)-H1NRM1(IFIELD))/H1NRM1(IFIELD))
        IF (RDIFF > RDMAX) RDMAX = RDIFF
        IF (NID == 0) WRITE (6,*) ' ifield, rdiff ',ifield,rdiff
        IF (RDIFF > XLIM) IFACCX = .FALSE. 
    100 END DO

    IF ( .NOT. IFACCX) THEN
        IF (NID == 0) THEN
            WRITE (6,*) ' '
            write (6,*) 'Extrapolation attempt rejected'
            write (6,*) ' '
        ENDIF
        CALL MKARR (Z)
    ELSE
        IF (NID == 0) THEN
            write (6,*)  ' '
            write (6,*) 'Extrapolation accepted'
            write (6,*) ' '
        ENDIF
        IF (RDMAX < RDLIM) IFSSVT = .TRUE. 
        CALL FILLLAG
    ENDIF

    RETURN
    END SUBROUTINE CHKEXT

    SUBROUTINE FILLLAG
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    NBDINP = 3
    IF (IFFLOW) THEN
        CALL LAGVEL
        CALL LAGVEL
    ENDIF
    IF (IFHEAT) THEN
        DO 100 IFIELD=2,NFIELD
            CALL LAGSCAL
            CALL LAGSCAL
        100 END DO
    ENDIF
    RETURN
    END SUBROUTINE FILLLAG

    SUBROUTINE GONSTEP (N,ITEST)
!----------------------------------------------------------------

!     Do N steps; return if steady state

!----------------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'STEADY'
    EXTERNAL GOSTEP

    DO 1000 JSTEP=1,N
        IF (ITEST == 0 .AND. IFSSVT) GOTO 1001
        IF (ITEST == 1 .AND. (IFSSVT .OR. IFEXVT)) GOTO 1001
        CALL GOSTEP
    1000 END DO
    1001 CONTINUE

    RETURN
    END SUBROUTINE GONSTEP

    SUBROUTINE GO1STEP (X,Y,NVEC)
!----------------------------------------------------------------

!     Advance one (or more) time step(s)

!----------------------------------------------------------------
    use size_m
    INCLUDE 'TOTAL'
    real*8 :: X(1), Y(1)

    CALL MKARR (X)
    IF ( .NOT. IFSKIP) NJSTEP=1
    IF (     IFSKIP) NJSTEP=NSSKIP

    DO 9000 JSTEP=1,NJSTEP
    
        ISTEP = ISTEP+1
        CALL SETTIME
        CALL SETPROP
        IF (IFMODP) CALL MODPROP
        CALL SETSOLV
        CALL COMMENT
        DO 100 IGEOM=1,2
            IF (IFGEOM) THEN
                CALL GENGEOM (IGEOM)
                CALL GENEIG  (IGEOM)
            ENDIF
            IF (IFFLOW) CALL FLUID (IGEOM)
            IF (IFHEAT) CALL HEAT  (IGEOM)
            IF (IFMVBD) CALL MESHV (IGEOM)
        100 END DO
        CALL PREPOST( .FALSE. )
        CALL USERCHK
    
    9000 END DO

    IF (ISTEP > 1) CALL CHKSSVT
    CALL MKVEC (Y)

    RETURN
    END SUBROUTINE GO1STEP

    SUBROUTINE GOSTEP
!----------------------------------------------------------------

!     Advance one (or more) time step(s)

!----------------------------------------------------------------
    use size_m
    INCLUDE 'TOTAL'

    IF ( .NOT. IFSKIP) NJSTEP=1
    IF (     IFSKIP) NJSTEP=NSSKIP

    DO 9000 JSTEP=1,NJSTEP
    
        ISTEP = ISTEP+1
        CALL SETTIME
        CALL SETPROP
        IF (IFMODP) CALL MODPROP
        CALL SETSOLV
        CALL COMMENT
        DO 100 IGEOM=1,2
            IF (IFGEOM) THEN
                CALL GENGEOM (IGEOM)
                CALL GENEIG  (IGEOM)
            ENDIF
            IF (IFFLOW) CALL FLUID (IGEOM)
            IF (IFHEAT) CALL HEAT  (IGEOM)
            IF (IFMVBD) CALL MESHV (IGEOM)
        100 END DO
        CALL PREPOST( .FALSE. )
        CALL USERCHK
    
    9000 END DO

    IF (ISTEP > 1) CALL CHKSSVT

    RETURN
    END SUBROUTINE GOSTEP

    SUBROUTINE MODPROP
!------------------------------------------------------------------

!     Modify the properties

!------------------------------------------------------------------
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'
    INCLUDE 'INPUT'

    MFIELD=1
    IF ( .NOT. IFFLOW) MFIELD=2
    DO 100 IFIELD=MFIELD,NFIELD
        NTOT  = NX1*NY1*NZ1*NELFLD(IFIELD)
        TAU   = .02*TAUSS(IFIELD)
        DECAY = 1.+99.*EXP(-TIME/TAU)
        CALL CMULT (VDIFF(1,1,1,1,IFIELD),DECAY,NTOT)
    !         if (nid.eq.0)
    !     $   write (6,*) '.......... diff = ',IFIELD,vdiff(1,1,1,1,IFIELD)
    100 END DO

    RETURN
    END SUBROUTINE MODPROP

    SUBROUTINE MKVEC (X)
!-------------------------------------------------------------

!     Fill up the vector X with VX, VY, ....

!-------------------------------------------------------------
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    real*8 :: X(1)

    NTOTV = NX1*NY1*NZ1*NELV

    IF (IFFLOW) THEN
        DO 100 I=1,NTOTV
            X(I) = VX(I,1,1,1)
        100 END DO
        DO 200 I=1,NTOTV
            X(I+NTOTV) = VY(I,1,1,1)
        200 END DO
        IF (NDIM == 3) THEN
            IOFF = 2*NTOTV
            DO 300 I=1,NTOTV
                X(I+IOFF) = VZ(I,1,1,1)
            300 END DO
        ENDIF
    ENDIF

    IF (IFHEAT) THEN
        IOFF = NDIM*NTOTV
        DO 401 IFIELD=2,NFIELD
            NTOT = NX1*NY1*NZ1*NELFLD(IFIELD)
            DO 400 I=1,NTOT
                X(I+IOFF) = T(I,1,1,1,IFIELD-1)
            400 END DO
            IOFF = IOFF+NTOT
        401 END DO
    ENDIF

    RETURN
    END SUBROUTINE MKVEC

    SUBROUTINE MKARR (X)
!------------------------------------------------------------------

!     Split the vector X into VX, VY, .....

!------------------------------------------------------------------
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'INPUT'
    real*8 :: X(1)

    NTOTV = NX1*NY1*NZ1*NELV

    IF (IFFLOW) THEN
        DO 10 I=1,NTOTV
            VX(I,1,1,1) = X(I)
        10 END DO
        DO 20 I=1,NTOTV
            VY(I,1,1,1) = X(I+NTOTV)
        20 END DO
        IF (NDIM == 3) THEN
            IOFF = 2*NTOTV
            DO 30 I=1,NTOTV
                VZ(I,1,1,1) = X(I+IOFF)
            30 END DO
        ENDIF
    ENDIF

    IF (IFHEAT) THEN
        IOFF = NDIM*NTOTV
        DO 41 IFIELD=2,NFIELD
            NTOT = NX1*NY1*NZ1*NELFLD(IFIELD)
            DO 40 I=1,NTOT
                T(I,1,1,1,IFIELD-1) = X(I+IOFF)
            40 END DO
            IOFF = IOFF+NTOT
        41 END DO
    ENDIF

    RETURN
    END SUBROUTINE MKARR

    SUBROUTINE SSPARAM (KMAX,L)
!------------------------------------------------------------------------------

!     Set steady state parameters

!------------------------------------------------------------------------------
    use size_m
    INCLUDE 'TSTEP'
    INCLUDE 'INPUT'
    INCLUDE 'STEADY'

    IF (L == 0) THEN
        CALL SSINIT (KMAX)
    ELSEIF (L == 1) THEN
        ISTEP  = 0
    
        PRELAX = 1.E-2
        IF (IFSPLIT) PRELAX = 1.E-5
    
        CTARG  = 1.
        IF (IFSPLIT) CTARG  = 1.
        IF (NX1 >= 10) THEN
            CTARG  = 2.
            IF (IFSPLIT) CTARG = 2.
        ENDIF
    
        KMAX   = 5
        NBDINP = 3
        IF (IFMODEL) NBDINP = 2
        NSSKIP = 2
        IFSKIP = .TRUE. 
        IFMODP = .FALSE. 
    
    ELSEIF (L == 2) THEN
    
        PRELAX = 1.E-3
        IF (IFSPLIT) PRELAX = 1.E-5
    
        CTARG = 1.
        IF (IFSPLIT) CTARG = 1.
        IF (NX1 >= 10) THEN
            CTARG = 2.
            IF (IFSPLIT) CTARG = 2.
        ENDIF
    
        KMAX   = 5
        NBDINP = 3
        IF (IFMODEL) NBDINP = 2
        NSSKIP = 2
        IFSKIP = .TRUE. 
        IFMODP = .FALSE. 
    
    ELSE
    ENDIF
    CALL SETCHAR
    RETURN
    END SUBROUTINE SSPARAM

    SUBROUTINE CHKSSVT
!-----------------------------------------------------------------------

!     Check for global steady state (velocity and temp/passive scalar)

!-----------------------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'

    IF (IFFLOW) THEN
        IMESH  = 1
        IFIELD = 1
        CALL CHKSSV
    ENDIF

    IMESH = 2
    DO 100 IFIELD=2,NFIELD
        CALL CHKSST
    100 END DO

    IFSSVT = .TRUE. 
    IFEXVT = .TRUE. 
    MFIELD = 1
    IF ( .NOT. IFFLOW) MFIELD=2
    DO 200 IFIELD=MFIELD,NFIELD
        IF( .NOT. IFSTST(IFIELD)) IFSSVT = .FALSE. 
        IF( .NOT. IFEXTR(IFIELD)) IFEXVT = .FALSE. 
    200 END DO
    IF (IFNATC) THEN
        IF (IFSTST(2) .AND. ( .NOT. IFSTST(1))) IFSTST(2) = .FALSE. 
    ENDIF
    RETURN
    END SUBROUTINE CHKSSVT

    SUBROUTINE CHKSSV
!--------------------------------------------------------------------

!     Check steady state for velocity

!--------------------------------------------------------------------
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'MASS'
    INCLUDE 'INPUT'
    INCLUDE 'EIGEN'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'
    COMMON /CTOLPR/ DIVEX
    COMMON /CPRINT/ IFPRINT
    LOGICAL ::         IFPRINT

    COMMON /SCRSS2/ DV1 (LX1,LY1,LZ1,LELV) &
    ,              DV2 (LX1,LY1,LZ1,LELV) &
    ,              DV3 (LX1,LY1,LZ1,LELV)
    COMMON /SCRUZ/  W1  (LX1,LY1,LZ1,LELV) &
    ,              W2  (LX1,LY1,LZ1,LELV) &
    ,              W3  (LX1,LY1,LZ1,LELV) &
    ,              BDIVV(LX2,LY2,LZ2,LELV)
    COMMON /SCRMG/  T1  (LX1,LY1,LZ1,LELV) &
    ,              T2  (LX1,LY1,LZ1,LELV) &
    ,              T3  (LX1,LY1,LZ1,LELV) &
    ,              DIVV(LX2,LY2,LZ2,LELV)
    COMMON /SCRVH/  H1  (LX1,LY1,LZ1,LELV) &
    ,              H2  (LX1,LY1,LZ1,LELV)

    CALL OPSUB3 (DV1,DV2,DV3,VX,VY,VZ,VXLAG,VYLAG,VZLAG)
    CALL NORMVC (DVNNH1,DVNNSM,DVNNL2,DVNNL8,DV1,DV2,DV3)
    INTYPE = -1
    CALL SETHLM (H1,H2,INTYPE)
    CALL OPHX   (W1,W2,W3,DV1,DV2,DV3,H1,H2)
    CALL OPDSSUM(W1,W2,W3)
    CALL OPMASK (W1,W2,W3)
    CALL OPCOLV3(T1,T2,T3,W1,W2,W3,BINVM1)
    CALL OPHX   (W1,W2,W3,T1,T2,T3,H1,H2)
    CALL OPCOL2 (W1,W2,W3,DV1,DV2,DV3)
    NTOT1  = NX1*NY1*NZ1*NELV
    USNRM1 = GLSUM(W1,NTOT1)
    USNRM2 = GLSUM(W2,NTOT1)
    USNRM3 = 0.
    IF (NDIM == 3) USNRM3 = GLSUM(W3,NTOT1)
    USNORM = SQRT( (USNRM1+USNRM2+USNRM3)/VOLVM1 )

    NTOT2 = NX2*NY2*NZ2*NELV
    CALL OPDIV (BDIVV,VX,VY,VZ)
    CALL COL3 (DIVV,BDIVV,BM2INV,NTOT2)
    DNORM = SQRT(GLSC2(DIVV,BDIVV,NTOT2)/VOLVM2)

    TOLOLD = TOLPS
    CALL SETTOLV
    TOLHV3 = TOLHV*(NDIM)
    IF (IFSTRS) TOLHV3 = TOLHV
    IF (NID == 0 .AND. IFPRINT) THEN
        WRITE (6,*) 'USNORM, TOLHV',USNORM,TOLHV3
        WRITE (6,*) 'DNORM, TOLPS',DNORM,TOLPS
    ENDIF
    IF (DNORM > (1.1*DIVEX) .AND. DIVEX > 0. &
     .AND. TOLPDF == 0.) TOLPDF = 5.*DNORM
    USREL = USNORM/TOLHV3
    DREL  = DNORM/TOLPS

    IF (TOLREL > 0.) THEN
        EXFAC = .3/TOLREL
    ELSE
        WRITE (6,*) 'WARNING: TOLREL=0. Please modify *.rea'
        call exitt
    ENDIF
    IFEXTR(IFIELD) = .FALSE. 
!      IF ((USREL.LT.EXFAC).OR.(TIME.GT.TXNEXT(IFIELD)))
!     $                                       IFEXTR(IFIELD) = .TRUE.
    IF (USREL < EXFAC)                    IFEXTR(IFIELD) = .TRUE. 
    if (nid == 0 .AND. ifprint) &
    WRITE (6,*) 'Tau, Txnext ',IFIELD,tauss(ifield),txnext(ifield)

    IFSTST(IFIELD) = .FALSE. 
    USLIM = 2.*TOLHV3
    DLIM  = 2.*TOLPS
    IF (USNORM < USLIM .AND. DNORM < DLIM .AND. .NOT. IFSPLIT) &
    IFSTST(IFIELD) = .TRUE. 
    IF (USNORM < USLIM .AND. IFSPLIT) &
    IFSTST(IFIELD) = .TRUE. 

    RETURN
    END SUBROUTINE CHKSSV

    SUBROUTINE CHKSST
!----------------------------------------------------------------------

!     Check for steady state for temperature/passive scalar

!----------------------------------------------------------------------
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'
    COMMON /SCRUZ/  DELTAT (LX1,LY1,LZ1,LELT) &
    ,              WA     (LX1,LY1,LZ1,LELT) &
    ,              WB     (LX1,LY1,LZ1,LELT)
    COMMON /SCRVH/  H1     (LX1,LY1,LZ1,LELT) &
    ,              H2     (LX1,LY1,LZ1,LELT)
    COMMON /CPRINT/ IFPRINT
    LOGICAL ::         IFPRINT

    NTOT = NX1*NY1*NZ1*NELT
    CALL SUB3 (DELTAT(1,1,1,1),T(1,1,1,1,IFIELD-1), &
    TLAG(1,1,1,1,1,IFIELD-1),NTOT)
    CALL NORMSC (DVNNH1,DVNNSM,DVNNL2,DVNNL8,DELTAT,IMESH)
    INTYPE = -1
    ISD    = 1
    CALL SETHLM (H1,H2,INTYPE)
    CALL AXHELM (WA,DELTAT,H1,H2,IMESH,ISD)
    CALL DSSUM  (WA,NX1,NY1,NZ1)
    CALL COL2   (WA,TMASK(1,1,1,1,IFIELD-1),NTOT)
    CALL COL3   (WB,WA,BINTM1,NTOT)
    CALL AXHELM (WA,WB,H1,H2,IMESH,ISD)
    CALL COL2   (WA,DELTAT,NTOT)
    USNORM = SQRT(GLSUM(WA,NTOT)/VOLTM1)

    CALL SETTOLT
    IF (NID == 0 .AND. IFPRINT) &
    WRITE (6,*) 'USNORM, TOLHT',USNORM,TOLHT(IFIELD)
    USREL = USNORM/TOLHT(IFIELD)

    IF (TOLREL > 0.) THEN
        EXFAC = .3/TOLREL
    ELSE
        WRITE (6,*) 'WARNING: TOLREL=0. Please modify *.rea'
        call exitt
    ENDIF
    IFEXTR(IFIELD) = .FALSE. 
!      IF ((USREL.LT.EXFAC).OR.(TIME.GT.TXNEXT(IFIELD)))
!     $                     IFEXTR(IFIELD) = .TRUE.
    IF (USREL < EXFAC)  IFEXTR(IFIELD) = .TRUE. 
    IF(NID == 0 .AND. IFPRINT) &
    WRITE (6,*) 'Tau, Txnext ',IFIELD,tauss(ifield),txnext(ifield)

    IFSTST(IFIELD) = .FALSE. 
    USLIM = 2.*TOLHT(IFIELD)
    IF (USNORM < USLIM) IFSTST(IFIELD) = .TRUE. 

    RETURN
    END SUBROUTINE CHKSST

    SUBROUTINE SSNORMD (DV1,DV2,DV3)
    use size_m
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'
    REAL :: DV1(1),DV2(1),DV3(1)
    CALL NORMVC (DVDFH1,DVDFSM,DVDFL2,DVDFL8,DV1,DV2,DV3)
    RETURN
    END SUBROUTINE SSNORMD

    SUBROUTINE SSNORMP (DV1,DV2,DV3)
    use size_m
    INCLUDE 'TSTEP'
    INCLUDE 'STEADY'
    REAL :: DV1(1),DV2(1),DV3(1)
    CALL NORMVC (DVPRH1,DVPRSM,DVPRL2,DVPRL8,DV1,DV2,DV3)
    RETURN
    END SUBROUTINE SSNORMP

    SUBROUTINE SETTOLV
!-------------------------------------------------------------------

!     Set tolerances for velocity solver

!-------------------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'EIGEN'
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'
    INCLUDE 'SOLN'
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
    INCLUDE 'INPUT'
    INCLUDE 'EIGEN'
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'
    INCLUDE 'SOLN'
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

    SUBROUTINE CHKTOLP (TOLMIN)
    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'MASS'
    INCLUDE 'TSTEP'
    COMMON /SCRMG/ DIVFLD (LX2,LY2,LZ2,LELV) &
    ,             WORK   (LX2,LY2,LZ2,LELV)
    NTOT2 = NX2*NY2*NZ2*NELV
    CALL OPDIV   (DIVFLD,VX,VY,VZ)
    CALL COL3    (WORK,DIVFLD,BM2INV,NTOT2)
    CALL COL2    (WORK,DIVFLD,NTOT2)
    DIVV  = SQRT(GLSUM(WORK,NTOT2)/VOLVM2)

    IFIELD = 1
    CALL SETTOLV
    TOLMIN = DIVV/100.
    IF (TOLMIN < TOLPS) TOLMIN = TOLPS
    RETURN
    END SUBROUTINE CHKTOLP

    SUBROUTINE SETCHAR
!-----------------------------------------------------------------------

!     If characteristics, need number of sub-timesteps (DT/DS).
!     Current sub-timeintegration scheme: RK4.
!     If not characteristics, i.e. standard semi-implicit scheme,
!     check user-defined Courant number.

!----------------------------------------------------------------------
    use size_m
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'

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

    SUBROUTINE PROJECT
!--------------------------------------------------------------------

!     Project current solution onto the closest incompressible field

!--------------------------------------------------------------------
    use size_m
    INCLUDE 'TOTAL'
    COMMON /SCRNS/ W1    (LX1,LY1,LZ1,LELV) &
    ,             W2    (LX1,LY1,LZ1,LELV) &
    ,             W3    (LX1,LY1,LZ1,LELV) &
    ,             DV1   (LX1,LY1,LZ1,LELV) &
    ,             DV2   (LX1,LY1,LZ1,LELV) &
    ,             DV3   (LX1,LY1,LZ1,LELV) &
    ,             RESPR (LX2,LY2,LZ2,LELV)
    COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV) &
    ,             H2    (LX1,LY1,LZ1,LELV)

    IF (NID == 0) WRITE(6,5)
    5 FORMAT(/,'  Project',/)

    NTOT1  = NX1*NY1*NZ1*NELV
    NTOT2  = NX2*NY2*NZ2*NELV
    INTYPE = 1
    CALL RZERO   (H1,NTOT1)
    CALL RONE    (H2,NTOT1)
    CALL OPDIV   (RESPR,VX,VY,VZ)
    CALL CHSIGN  (RESPR,NTOT2)
    CALL ORTHO   (RESPR)
    CALL UZAWA   (RESPR,H1,H2,INTYPE,ICG)
    CALL OPGRADT (W1,W2,W3,RESPR)
    CALL OPBINV  (DV1,DV2,DV3,W1,W2,W3,H2)
    CALL OPADD2  (VX,VY,VZ,DV1,DV2,DV3)
    RETURN
    END SUBROUTINE PROJECT
