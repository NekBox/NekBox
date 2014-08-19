    SUBROUTINE SETTURB
!------------------------------------------------------------------------

!     Set parameters for k-e turbulence models

!------------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'

    IFHEAT = .TRUE. 
    IFADVC(NFIELD-1) = .TRUE. 
    IFADVC(NFIELD)   = .TRUE. 
    IFTMSH(NFIELD-1) = .FALSE. 
    IFTMSH(NFIELD)   = .FALSE. 

    IFLDTK  = NFIELD - 2
    IFLDTE  = NFIELD - 1
    IFLDK   = NFIELD - 1
    IFLDE   = NFIELD

    RETURN
    END SUBROUTINE SETTURB
    SUBROUTINE PRETMIC

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    COMMON /COMLMX/ ISTPCL

    NTOT1 = NX1*NY1*NZ1*NELV
    NTOTF = NX1*NZ1*6*NELV
    ZWINI = 0.09*TLMAX
    UWINI = 0.0
    TKINI = 1.0
    TEINI = 1.0
    TLIMUL= 1.0
    ISTPCL= -1

    CALL CFILL (UWALL,UWINI,NTOTF)
    CALL CFILL (ZWALL,ZWINI,NTOTF)
    CALL RZERO (TURBL,NTOT1)

    IF (IFMODEL .AND. IFKEPS) THEN
        CALL CFILL (T(1,1,1,1,IFLDTK),TKINI,NTOT1)
        CALL CFILL (T(1,1,1,1,IFLDTE),TEINI,NTOT1)
    ENDIF

    RETURN
    END SUBROUTINE PRETMIC
    SUBROUTINE POSTMIC

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    COMMON /COMLMX/ ISTPCL

    IF ( .NOT. IFKEPS .AND. ISTPCL /= 0) CALL COMLSQ

    RETURN
    END SUBROUTINE POSTMIC
    SUBROUTINE CBCTURB
!---------------------------------------------------------------------

!     Set b.c. for k and e if k-e turbulence model.
!     Legal types of b.c.: Dirichlet, Neumann, Robin.

!---------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TURBO'
    CHARACTER CB*3

    IFLD   = 1
    NFACES = 2*NDIM

!     Set b.c. for k and e based on user specified velocity b.c.

    DO 100 IEL=1,NELV
        DO 100 IFACE=1,NFACES
            CB  = CBC (IFACE,IEL,IFLD)
            BC1 = BC(1,IFACE,IEL,IFLD)
            BC2 = BC(2,IFACE,IEL,IFLD)
            BC3 = BC(3,IFACE,IEL,IFLD)
        
            IF (CB == 'E  ' .OR. CB == 'P  ' .OR. cb == 'p  ') THEN
                CBC(  IFACE,IEL,IFLDK) = CB
                BC (1,IFACE,IEL,IFLDK) = BC1
                BC (2,IFACE,IEL,IFLDK) = BC2
                BC (3,IFACE,IEL,IFLDK) = BC3
                CBC(  IFACE,IEL,IFLDE) = CB
                BC (1,IFACE,IEL,IFLDE) = BC1
                BC (2,IFACE,IEL,IFLDE) = BC2
                BC (3,IFACE,IEL,IFLDE) = BC3
            ELSEIF (CB == 'V  ' .OR. CB == 'VL ') THEN
                CBC(IFACE,IEL,IFLDK) = 'KD '
                CBC(IFACE,IEL,IFLDE) = 'ED '
            ELSEIF (CB == 'v  ' .OR. CB == 'vl ') THEN
                CBC(IFACE,IEL,IFLDK) = 'kd '
                CBC(IFACE,IEL,IFLDE) = 'ed '
            ELSEIF (CB == 'W  ') THEN
                CBC(IFACE,IEL,IFLDK) = 'KW '
                CBC(IFACE,IEL,IFLDE) = 'EW '
            ELSEIF (CB == 'WS ') THEN
                CBC(IFACE,IEL,IFLDK) = 'KWS'
                CBC(IFACE,IEL,IFLDE) = 'EWS'
            ELSE
                CBC(IFACE,IEL,IFLDK) = 'KN '
                CBC(IFACE,IEL,IFLDE) = 'EN '
            ENDIF
        
    100 END DO

    RETURN
    END SUBROUTINE CBCTURB
    SUBROUTINE WHATFLD (IFTURB)
!--------------------------------------------------------------------

!     Find out if the current field is part of a turbulence model.

!--------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    LOGICAL ::  IFTURB, IFKFLD, IFEFLD

    IFTURB = .FALSE. 
    IF (IFMODEL) THEN
        CALL TURBFLD (IFKFLD,IFEFLD)
        IF (IFKFLD .OR. IFEFLD) IFTURB = .TRUE. 
    ENDIF

    RETURN
    END SUBROUTINE WHATFLD
    SUBROUTINE TURBFLD (IFKFLD,IFEFLD)
!------------------------------------------------------------------

!     Check if field is k or e. Set appropriate logicals.

!------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    LOGICAL ::  IFKFLD, IFEFLD

    IFKFLD = .FALSE. 
    IFEFLD = .FALSE. 
    IF (IFKEPS .AND. IFIELD == IFLDK) IFKFLD = .TRUE. 
    IF (IFKEPS .AND. IFIELD == IFLDE) IFEFLD = .TRUE. 

    RETURN
    END SUBROUTINE TURBFLD
    SUBROUTINE TVISCOS
!------------------------------------------------------------------------

!     Turbulence model active. Set effective viscosity.
!     Current turbulence models: k-e and algebraic.

!------------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TURBO'

    IF (IFKEPS) THEN
        CALL TVISCKE
    ELSE
        CALL TVISCA
    ENDIF

    RETURN
    END SUBROUTINE TVISCOS
    SUBROUTINE TVISCKE
!------------------------------------------------------------------------

!     Set effective viscosity (k-e model).

!------------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'

    IFLD  = 1
    NTOT1 = NX1*NY1*NZ1*NELV

    CALL CMULT2 (VTURB,T(1,1,1,1,IFLDTK),CMU,NTOT1)
    CALL COL2   (VTURB,T(1,1,1,1,IFLDTK),NTOT1)
    CALL INVCHK2(VTURB,T(1,1,1,1,IFLDTE),NTOT1)
    CALL ADDCOL3(VDIFF(1,1,1,1,IFLD),VTRANS(1,1,1,1,IFLD),VTURB,NTOT1)

    RETURN
    END SUBROUTINE TVISCKE
    SUBROUTINE TVISCA
!------------------------------------------------------------------------

!     Set effective viscosity (algebraic model).

!------------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'MASS'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    COMMON /SCRVH/ TA(LX1,LY1,LZ1,LELV)
    COMMON /COMLMX/ ISTPCL

    NTOT1  = NX1*NY1*NZ1*NELV

    IF (ISTPCL /= ISTEP) THEN
        CALL COMLSQ
        IF (ISTEP == 0) THEN
            CALL INIPHI (TA)
        ELSE
            CALL COMPHI (TA)
            CALL VSQRT  (TA,NTOT1)
        ENDIF
        CALL COL2   (TA,BM1,NTOT1)
        CALL DSSUM  (TA,NX1,NY1,NZ1)
        CALL COL2   (TA,BINVM1,NTOT1)
        CALL COL3   (VTURB,TA,TURBL,NTOT1)
        CALL COL2   (VTURB,TURBL,NTOT1)
        ISTPCL = ISTEP
    ENDIF

!     Sum eddy vscosity

    IF (IFIELD == 1) THEN
        CALL ADDCOL3 (VDIFF(1,1,1,1,IFIELD),VTRANS(1,1,1,1,IFIELD), &
        VTURB,NTOT1)
    ELSE
        CALL CMULT2  (TA,VTURB,STI,NTOT1)
        CALL ADDCOL3 (VDIFF(1,1,1,1,IFIELD),VTRANS(1,1,1,1,IFIELD), &
        TA,NTOT1)
    ENDIF

    RETURN
    END SUBROUTINE TVISCA
    SUBROUTINE TPROPK
!-----------------------------------------------------------------------

!     Set "viscosity" and "density" for k transport equation

!-----------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    COMMON /SCRCH/ TA(LX1,LY1,LZ1,LELV),TB(LX1,LY1,LZ1,LELV)

    NTOT1   = NX1*NY1*NZ1*NELV
    VISKIN  = PARAM(2)/PARAM(1)

    CALL CFILL  (TA,VISKIN,NTOT1)
    CALL CMULT2 (TB,VTURB,SKI,NTOT1)
    CALL ADD3   (VDIFF(1,1,1,1,IFLDK),TA,TB,NTOT1)
    CALL RONE   (VTRANS(1,1,1,1,IFLDK),NTOT1)

    RETURN
    END SUBROUTINE TPROPK
    SUBROUTINE TPROPE
!-----------------------------------------------------------------------

!     Set "viscosity" and "density" for epsilon transport equation

!-----------------------------------------------------------------------
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    COMMON /SCRCH/ TA(LX1,LY1,LZ1,LELV),TB(LX1,LY1,LZ1,LELV)

    NTOT1   = NX1*NY1*NZ1*NELV
    VISKIN  = PARAM(2)/PARAM(1)

    CALL CFILL  (TA,VISKIN,NTOT1)
    CALL CMULT2 (TB,VTURB,SEI,NTOT1)
    CALL ADD3   (VDIFF(1,1,1,1,IFLDE),TA,TB,NTOT1)
    CALL RONE   (VTRANS(1,1,1,1,IFLDE),NTOT1)

    RETURN
    END SUBROUTINE TPROPE
    SUBROUTINE MAKETQ
!-------------------------------------------------------------------

!      Set source term for k and epsilon transport equations

!-------------------------------------------------------------------
    LOGICAL :: IFKFLD,IFEFLD

    CALL TURBFLD (IFKFLD,IFEFLD)
    IF (IFKFLD) CALL SETQK
    IF (IFEFLD) CALL SETQE

    RETURN
    END SUBROUTINE MAKETQ
    SUBROUTINE SETQK

    INCLUDE 'SIZE'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'

    NTOT1 = NX1*NY1*NZ1*NELV

    CALL RZERO (BQ(1,1,1,1,IFLDTK),NTOT1)
    IF (ISTEP > 0) CALL TURBQK

    RETURN
    END SUBROUTINE SETQK
    SUBROUTINE SETQE

    INCLUDE 'SIZE'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'

    NTOT1 = NX1*NY1*NZ1*NELV

    CALL RZERO (BQ(1,1,1,1,IFLDTE),NTOT1)
    IF (ISTEP > 0) CALL TURBQE

    RETURN
    END SUBROUTINE SETQE
    SUBROUTINE TURBQK

    INCLUDE 'SIZE'
    INCLUDE 'MASS'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    COMMON /SCRCH/ TA(LX1,LY1,LZ1,LELV)

    NTOT1  = NX1*NY1*NZ1*NELV

    CALL COMPHI  (TA)
    CALL COL2    (TA,VTURB,NTOT1)
    CALL SUB2    (TA,T(1,1,1,1,IFLDTE),NTOT1)
    CALL ADDCOL3 (BQ(1,1,1,1,IFLDTK),BM1,TA,NTOT1)

    RETURN
    END SUBROUTINE TURBQK
    SUBROUTINE TURBQE

    INCLUDE 'SIZE'
    INCLUDE 'MASS'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    COMMON /SCRUZ/ TA(LX1,LY1,LZ1,LELV) &
    , TB(LX1,LY1,LZ1,LELV)

    NTOT1  = NX1*NY1*NZ1*NELV

    CALL COMPHI  (TA)
    CALL COL2    (TA,VTURB,NTOT1)
    CALL CMULT   (TA,CE1,NTOT1)
    CALL CMULT2  (TB,T(1,1,1,1,IFLDTE),CE2,NTOT1)
    CALL SUB2    (TA,TB,NTOT1)
    CALL COL2    (TA,T(1,1,1,1,IFLDTE),NTOT1)
    CALL INVCHK2 (TA,T(1,1,1,1,IFLDTK),NTOT1)
    CALL ADDCOL3 (BQ(1,1,1,1,IFLDTE),BM1,TA,NTOT1)

    RETURN
    END SUBROUTINE TURBQE
    SUBROUTINE TURBWBC (TMP,TMA,SMU)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'

    DIMENSION TMP(LX1,LY1,LZ1,1) &
    , TMA(LX1,LY1,LZ1,1) &
    , SMU(LX1,LY1,LZ1,1)
    CHARACTER CB*3

    NFACE = 2*NDIM
    NTOT1 = NX1*NY1*NZ1*NELV
    CALL RZERO (TMA,NTOT1)
    CALL RZERO (SMU,NTOT1)

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB=CBC(IFC,IEL,IFIELD)
            IF (CB /= 'KWS' .AND. CB /= 'EWS') GOTO 100
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFC)
            DO 120 IZ=KZ1,KZ2
                DO 120 IY=KY1,KY2
                    DO 120 IX=KX1,KX2
                        SMU(IX,IY,IZ,IEL)=SMU(IX,IY,IZ,IEL) + 1.0
            120 END DO
            IF (CB == 'KWS') CALL FACEWSK (TMA(1,1,1,IEL),IEL,IFC)
            IF (CB == 'EWS') CALL FACEWSE (TMA(1,1,1,IEL),IEL,IFC)
    100 END DO

    CALL DSSUM   (TMA,NX1,NY1,NZ1)
    CALL DSSUM   (SMU,NX1,NY1,NZ1)
    CALL INVCHK2 (TMA,SMU,NTOT1)

    DO 200 IEL=1,NELV
        DO 200 IFC=1,NFACE
            CB=CBC(IFC,IEL,IFIELD)
            IF (CB /= 'KWS' .AND. CB /= 'EWS') GOTO 200
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFC)
            DO 220 IZ=KZ1,KZ2
                DO 220 IY=KY1,KY2
                    DO 220 IX=KX1,KX2
                        TMP(IX,IY,IZ,IEL)=TMA(IX,IY,IZ,IEL)
            220 END DO
    200 END DO

    RETURN
    END SUBROUTINE TURBWBC
    SUBROUTINE FACEWSK (S,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'TURBO'
    DIMENSION S(LX1,LY1,LZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
    I = 0

    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I   = I + 1
            ZZ  = ZWALL(I,1,IFC,IEL)
            USS = UWALL(I,1,IFC,IEL)**2
            IF (ZZ < 1.E-20) USS = 0.0
            S(J1,J2,1) = S(J1,J2,1) + CMI*USS
    100 END DO

    RETURN
    END SUBROUTINE FACEWSK
    SUBROUTINE FACEWSE (S,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'TURBO'
    DIMENSION S(LX1,LY1,LZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
    I = 0

    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I   = I + 1
            ZZ  = ZWALL(I,1,IFC,IEL)
            IF (ZZ > 1.E-20) THEN
                USS = ABS( UWALL(I,1,IFC,IEL) )
                USS = USS**3 / ZZ
            ELSE
                USS = 0.0
            ENDIF
            S(J1,J2,1) = S(J1,J2,1) + VKI*USS
    100 END DO

    RETURN
    END SUBROUTINE FACEWSE
    SUBROUTINE SETTMC

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TURBO'

    NTOT1 = NX1*NY1*NZ1*NELV

    CMU = 0.09
    CMT = 0.09
    SGK = 1.00
    SGE = 1.30
    CE1 = 1.44
    CE2 = 1.92
    VKC = 0.40
    BTA = 0.13
    TGA = 0.0385
    SGT = 0.90

    BETA1 =  7.0
    BETA2 = 14.0

    CMI = 1.0 / SQRT(CMU)
    SKI = 1.0 / SGK
    SEI = 1.0 / SGE
    VKI = 1.0 / VKC
    BTI = 1.0 / BTA
    STI = 1.0 / SGT

    ZPLDAT =  30.0
    ZPUDAT = 100.0
    ZPVDAT =   5.0

    RETURN
    END SUBROUTINE SETTMC
    SUBROUTINE COMPHI (PHI)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    COMMON /CTMP0/ EXZ(LX1,LY1,LZ1,LELT) &
    , EYZ(LX1,LY1,LZ1,LELT)
    COMMON /CTMP1/ EXX(LX1,LY1,LZ1,LELT) &
    , EXY(LX1,LY1,LZ1,LELT) &
    , EYY(LX1,LY1,LZ1,LELT) &
    , EZZ(LX1,LY1,LZ1,LELT)

    DIMENSION PHI(LX1,LY1,LZ1,1)

    NTOT1  = NX1*NY1*NZ1*NELV
    MATMOD = 0
    CALL STNRATE (VX,VY,VZ,NELV,MATMOD)

    CALL RZERO   (PHI,NTOT1)
    CALL ADDCOL3 (PHI,EXX,EXX,NTOT1)
    CALL ADDCOL3 (PHI,EYY,EYY,NTOT1)
    IF (IFAXIS .OR. NDIM == 3) CALL ADDCOL3 (PHI,EZZ,EZZ,NTOT1)

    CONST = 2.0
    CALL CMULT   (PHI,CONST,NTOT1)
    CALL ADDCOL3 (PHI,EXY,EXY,NTOT1)
    IF (NDIM == 3) THEN
        CALL ADDCOL3 (PHI,EXZ,EXZ,NTOT1)
        CALL ADDCOL3 (PHI,EYZ,EYZ,NTOT1)
    ENDIF

    RETURN
    END SUBROUTINE COMPHI
    SUBROUTINE INIPHI (PHI)
    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    COMMON /SCRMG/ TA(LX1,LY1,LZ1,LELV) &
    ,  TB(LX1,LY1,LZ1,LELV) &
    ,  TC(LX1,LY1,LZ1,LELV) &
    ,  TD(LX1,LY1,LZ1,LELV)

    DIMENSION PHI(LX1,LY1,LZ1,1)

    DELTA  = 1.0E-9
    X      = 1.0 + DELTA
    Y      = 1.0
    DIFF   = ABS(X - Y)
    IF (DIFF == 0.0) EPS = 1.0E-6
    IF (DIFF > 0.0) EPS = 1.0E-13

    NTOT1  = NX1*NY1*NZ1*NELV

    CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
    CALL OPDOT   (TA,VX,VY,VZ,VX,VY,VZ,NTOT1)
    UMAX = GLMAX(TA,NTOT1)
    UMIN = GLMIN(TA,NTOT1)
    UDEL = SQRT(UMAX) - SQRT(UMIN)

    IF (UDEL > EPS) THEN
        UDEL = UDEL / TLMAX
    ELSE
        CALL NEKUF (TA,TB,TC)
        CALL OPDOT (TD,TA,TB,TC,TA,TB,TC,NTOT1)
        FMAX = GLMAX(TD,NTOT1)
        FMAX = SQRT(FMAX)
        UDEL = SQRT( 2.0 * FMAX ) / TLMAX
    ENDIF

    CALL CFILL (PHI,UDEL,NTOT1)

    RETURN
    END SUBROUTINE INIPHI
    SUBROUTINE TWALLUZ (IGEOM)

    INCLUDE 'SIZE'
    INCLUDE 'GEOM'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    COMMON /SCRMG/ VISKIN(LX1,LY1,LZ1,LELV) &
    , XWLL(LX1,LZ1) &
    , YWLL(LX1,LZ1) &
    , ZWLL(LX1,LZ1)
    COMMON /SCRUZ/ V1(LX1,LY1,LZ1,LELV) &
    , V2(LX1,LY1,LZ1,LELV) &
    , V3(LX1,LY1,LZ1,LELV)
    CHARACTER CB*3

    IF (IGEOM == 1 .OR. .NOT. IFCWUZ) RETURN

    IFLD  = 1
    NFACE = 2*NDIM
    NTOT1 = NX1*NY1*NZ1*NELV

    CALL RZERO3  (V1,V2,V3,NTOT1)
    CALL BCTWALL (V1,V2,V3)
    CALL COPY    (VISKIN,VDIFF(1,1,1,1,IFLD),NTOT1)
    CALL INVCOL2 (VISKIN,VTRANS(1,1,1,1,IFLD),NTOT1)
    CALL SUB2    (VISKIN,VTURB,NTOT1)
    VISMIN = GLMIN(VISKIN,NTOT1)

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFLD)
            IF (CB /= 'WS ' .AND. CB /= 'ws '   .AND. &
            CB /= 'WSL' .AND. CB /= 'wsl' ) GOTO 100
            CALL FACEXV (XWLL,YWLL,ZWLL,XM1(1,1,1,IEL),YM1(1,1,1,IEL), &
            ZM1(1,1,1,IEL),IFC,0)
            CALL COMWUZ (XWLL,YWLL,ZWLL,V1(1,1,1,IEL),V2(1,1,1,IEL), &
            V3(1,1,1,IEL),VISKIN(1,1,1,IEL),VISMIN,IEL,IFC)
    100 END DO

    RETURN
    END SUBROUTINE TWALLUZ
    SUBROUTINE TWALLSH

    INCLUDE 'SIZE'
    INCLUDE 'GEOM'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    COMMON /SCRSF/ TRX(LX1,LY1,LZ1) &
    , TRY(LX1,LY1,LZ1) &
    , TRZ(LX1,LY1,LZ1)
    CHARACTER CB*3

    IF ( ISTEP <= 1 .OR. .NOT. IFSWALL ) RETURN

    IFLD    = 1
    NFACE   = 2*NDIM
    NXYZ1   = NX1*NY1*NZ1

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFLD)
            IF (CB /= 'WS ' .AND. CB /= 'ws '   .AND. &
            CB /= 'WSL' .AND. CB /= 'wsl' ) GOTO 100
            CALL RZERO3  (TRX,TRY,TRZ,NXYZ1)
            CALL FACEWS  (TRX,TRY,TRZ,IEL,IFC)
            CALL FACCVS  (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
            CALL ADD2    (BFX(1,1,1,IEL),TRX,NXYZ1)
            CALL ADD2    (BFY(1,1,1,IEL),TRY,NXYZ1)
            IF (NDIM == 3) CALL ADD2 (BFZ(1,1,1,IEL),TRZ,NXYZ1)
    100 END DO

    RETURN
    END SUBROUTINE TWALLSH
    SUBROUTINE FACEWS (TRX,TRY,TRZ,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),TRZ(LX1,LY1,LZ1)

    IFLD = 1
    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    IF (NDIM == 2) THEN
        I = 0
        DO 100 J2=JS2,JF2,JSKIP2
            DO 100 J1=JS1,JF1,JSKIP1
                I    = I + 1
                WSH  = VTRANS(J1,J2,1,IEL,IFLD) * UWALL(I,1,IFC,IEL)**2
                TRX(J1,J2,1) = TWX(I,1,IFC,IEL) * WSH
                TRY(J1,J2,1) = TWY(I,1,IFC,IEL) * WSH
        100 END DO
    ELSE
        I = 0
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                I    = I + 1
                WSH  = VTRANS(J1,J2,1,IEL,IFLD) * UWALL(I,1,IFC,IEL)**2
                TRX(J1,J2,1) = TWX(I,1,IFC,IEL) * WSH
                TRY(J1,J2,1) = TWY(I,1,IFC,IEL) * WSH
                TRZ(J1,J2,1) = TWZ(I,1,IFC,IEL) * WSH
        200 END DO
    ENDIF

    RETURN
    END SUBROUTINE FACEWS
    SUBROUTINE COMWUZ (XWLL,YWLL,ZWLL,V1,V2,V3,VISKIN,VISMIN,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    INCLUDE 'TSTEP'
    COMMON /SCRCH/ UTW(LX1),ZNW(LX1),USTR(LX1),ZSTR(LX1)

    DIMENSION XWLL(LX1,LZ1),YWLL(LX1,LZ1),ZWLL(LX1,LZ1) &
    , VISKIN(LX1,LY1,LZ1) &
    , V1(LX1,LY1,LZ1) &
    , V2(LX1,LY1,LZ1) &
    , V3(LX1,LY1,LZ1) &
    , ISHDAT(6)

    DATA ISHDAT / 1,-1,-1, 1, 1,-1/

    NXM1  = NX1 - 1
    NXM2  = NX1 - 2
    ISCH  = ISHDAT(IFC)

    DELTA = 1.0E-9
    X     = 1.0 + DELTA
    Y     = 1.0
    DIFF  = ABS(X - Y)
    IF (DIFF == 0.0) EPSA = 1.0E-5
    IF (DIFF > 0.0) EPSA = 1.0E-12

    EPSV = EPSA * VNRML8
    IF (EPSV < EPSA) EPSV = EPSA
    ZLOG = 0.50 * TLMAX

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    DO 500 JWY=1,NZ1
        DO 500 JWX=1,NX1
            IWX    = JS1 + (JWX-1)*JSKIP1
            IWY    = JS2 + (JWY-1)*JSKIP2
            CALL NORIND (JS3,JF3,JSKIP3,IWX,IWY,IFC,ISCH)
            VTAN1  = V1(IWX,IWY,1)
            VTAN2  = V2(IWX,IWY,1)
            VTAN3  = V3(IWX,IWY,1)
            CALL CWREF  (XWLL,YWLL,ZWLL,UTW,ZNW,VTAN1,VTAN2,VTAN3,IEL,IFC, &
            JWX,JWY,JS3,JSKIP3)
            IZER0   = 0
            IULOG   = 0
            UST1    = 5.0 * VISMIN / ZNW(2)
            USTR(1) = 0.0
            ZSTR(1) = 0.0
            UL      = UTW(1)
            ZL      = 0.0
            DO 510 JWZ=2,NX1
                J3 = JS3 + (JWZ-1)*JSKIP3
                VS = VISKIN(J3,1,1)
                ZW = ZNW(JWZ)
                UW = UTW(JWZ)
                IF (UW < EPSV) THEN
                    IZER0 = IZER0 + 1
                    IF (IZER0 == NXM2) THEN
                        UST = 0.0
                        ZW  = 0.1*TLMAX
                        GOTO 550
                    ENDIF
                    GOTO 510
                ENDIF
                IF (ZW > ZLOG) GOTO 520
                CALL COMUFR (UST,UW,ZW,VS,UST1,IEL,IFC,JWX,JWY)
                ZPW       = UST*ZW / VS
                USTR(JWZ) = UST
                ZSTR(JWZ) = ZPW
                IF (ZPW < ZPLDAT) THEN
                    UL     = UW
                    ZL     = ZW
                    UST1   = UST
                    IULOG  = IULOG + 1
                ELSEIF (ZPW > ZPUDAT) THEN
                    UU     = UW
                    ZU     = ZW
                    UST1   = UST
                    NWZ    = JWZ
                    GOTO 540
                ELSE
                    NWZ    = JWZ
                    GOTO 550
                ENDIF
            510 END DO
        
            520 CONTINUE
            NWZ = JWZ - 1
            ZPMAX = -1.0E20
            DO 530 JJ=2,NWZ
                IF (ZPMAX < ZSTR(JJ)) THEN
                    ZPMAX = ZSTR(JJ)
                    JMAX  = JJ
                ENDIF
            530 END DO
            UST = USTR(JMAX)
            ZW  = ZNW (JMAX)
        
        !         IF (IULOG .EQ. NXM1) THEN
        !            WRITE (6,*) ' '
        !            WRITE (6,*) ' *** WARNING MESSAGE ***'
        !            WRITE (6,'(1X,A4,I5,5X,A6,I5,5X,A9,I3,A1,I3)')
        !     $      'EL =',IEL,'FACE =',IFC,'WALL PT =',JWX,',',JWY
        !            WRITE (6,*) ' ELEMENT HEIGHT IS BELOW THE LOG LAYER.'
        !         ENDIF
        
            GOTO 550
        
            540 CONTINUE
        
        !         IF (NWZ.EQ.2) THEN
        !            WRITE (6,*) ' '
        !            WRITE (6,*) ' *** WARNING MESSAGE ***'
        !            WRITE (6,'(1X,A4,I5,5X,A6,I5,5X,A9,I3,A1,I3)')
        !     $      'EL =',IEL,'FACE =',IFC,'WALL PT =',JWX,',',JWY
        !            WRITE (6,*) ' LOG LAYER IS WITHIN ONE GRIDSPACE OF WALL.'
        !            WRITE (6,*) ' MESH REFINEMENT MAY BE REQUIRED.'
        !         ENDIF
        
            CALL SUBGRID (UTW,ZNW,UST,ZW,ZL,ZU,UL,UU,VS,UST1,IEL,IFC, &
            JWX,JWY,NWZ)
        
            550 CONTINUE
            UWALL(JWX,JWY,IFC,IEL) = UST
            ZWALL(JWX,JWY,IFC,IEL) = ZW
            IF (IZER0 == NXM2) THEN
                TWX(JWX,JWY,IFC,IEL) = 0.0
                TWY(JWX,JWY,IFC,IEL) = 0.0
            ELSE
                IWZ = JS3 + (NWZ - 1)*JSKIP3
                CALL COMTDIR (VTAN1,VTAN2,VTAN3,JWX,JWY,IWZ,IEL,IFC)
            ENDIF
        
    500 END DO

    RETURN
    END SUBROUTINE COMWUZ
    SUBROUTINE NORIND (JS3,JF3,JSKIP3,IWX,IWY,IFC,ISCH)

    INCLUDE 'SIZE'

    IF (IFC == 1 .OR. IFC == 3) THEN
        JSKIP3 = NX1
    ELSEIF (IFC == 2 .OR. IFC == 4) THEN
        JSKIP3 = 1
    ELSE
        JSKIP3 = NX1**2
    ENDIF
    IF (ISCH == -1) JSKIP3 = -JSKIP3

    JS3 = IWX + NX1 * (IWY - 1)
    JF3 = JS3 + JSKIP3 * (NX1 - 1)

    RETURN
    END SUBROUTINE NORIND
    SUBROUTINE CWREF (XWLL,YWLL,ZWLL,UTW,ZNW,VTAN1,VTAN2,VTAN3, &
    IEL,IFC,JWX,JWY,JS3,JSKIP3)

    INCLUDE 'SIZE'
    INCLUDE 'GEOM'
    INCLUDE 'SOLN'
    DIMENSION XWLL(LX1,LZ1),YWLL(LX1,LZ1),ZWLL(LX1,LZ1) &
    , UTW(LX1),ZNW(LX1)

    IF (NDIM == 2) THEN
        DO 100 II=1,NX1
            J3   =  JS3 + (II-1)*JSKIP3
            VRX  =  VX(J3,1,1,IEL) - VTAN1
            VRY  =  VY(J3,1,1,IEL) - VTAN2
            VDOT =  VRX**2 + VRY**2
            VNOR =  VRX * UNX(JWX,JWY,IFC,IEL) + &
            VRY * UNY(JWX,JWY,IFC,IEL)
            DIS2 =( XM1(J3,1,1,IEL) - XWLL(JWX,JWY) )**2 + &
            ( YM1(J3,1,1,IEL) - YWLL(JWX,JWY) )**2
            UTW(II) = SQRT(VDOT - VNOR**2)
            ZNW(II) = SQRT(DIS2)
        100 END DO
    ELSE
        DO 200 II=1,NX1
            J3   =  JS3 + (II-1)*JSKIP3
            VRX  =  VX(J3,1,1,IEL) - VTAN1
            VRY  =  VY(J3,1,1,IEL) - VTAN2
            VRZ  =  VZ(J3,1,1,IEL) - VTAN3
            VDOT =  VRX**2 + VRY**2 + VRZ**2
            VNOR =  VRX * UNX(JWX,JWY,IFC,IEL) + &
            VRY * UNY(JWX,JWY,IFC,IEL) + &
            VRZ * UNZ(JWX,JWY,IFC,IEL)
            DIS2 =( XM1(J3,1,1,IEL) - XWLL(JWX,JWY) )**2 + &
            ( YM1(J3,1,1,IEL) - YWLL(JWX,JWY) )**2 + &
            ( ZM1(J3,1,1,IEL) - ZWLL(JWX,JWY) )**2
            UTW(II) = SQRT(VDOT - VNOR**2)
            ZNW(II) = SQRT(DIS2)
        200 END DO
    ENDIF

    RETURN
    END SUBROUTINE CWREF
    SUBROUTINE COMTDIR (VTAN1,VTAN2,VTAN3,JWX,JWY,IWZ,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'GEOM'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'

    VRX = VX(IWZ,1,1,IEL) - VTAN1
    VRY = VY(IWZ,1,1,IEL) - VTAN2
    VRZ = VZ(IWZ,1,1,IEL) - VTAN3
    UN1 = UNX(JWX,JWY,IFC,IEL)
    UN2 = UNY(JWX,JWY,IFC,IEL)
    UN3 = UNZ(JWX,JWY,IFC,IEL)

    IF (NDIM == 2) THEN
        BINZ = VRY*UN1 -  VRX*UN2
        RBSQ = ABS( BINZ )
        TWX(JWX,JWY,IFC,IEL) =  UN2*BINZ / RBSQ
        TWY(JWX,JWY,IFC,IEL) = -UN1*BINZ / RBSQ
    ELSE
        BINX = VRZ*UN2 - VRY*UN3
        BINY = VRX*UN3 - VRZ*UN1
        BINZ = VRY*UN1 - VRX*UN2
        RBSQ = SQRT(BINX**2 + BINY**2 + BINZ**2)
        TWX(JWX,JWY,IFC,IEL) = (UN2*BINZ - UN3*BINY) / RBSQ
        TWY(JWX,JWY,IFC,IEL) = (UN3*BINX - UN1*BINZ) / RBSQ
        TWZ(JWX,JWY,IFC,IEL) = (UN1*BINY - UN2*BINX) / RBSQ
    ENDIF

    RETURN
    END SUBROUTINE COMTDIR
    SUBROUTINE SUBGRID (UTW,ZNW,UW,ZW,ZL,ZU,UL,UU,AKVIS,UST1,IEL,IFC, &
    JWX,JWY,JWZ)

    INCLUDE 'SIZE'
    INCLUDE 'TSTEP'
    INCLUDE 'TURBO'
    INCLUDE 'WZ'
    DIMENSION UTW(LX1),ZNW(LX1)

    IC     = 0
    NTMAX  = 100
    JWN    = JWZ - 1
    RL     = ZGM1(JWN,1)
    RU     = ZGM1(JWZ,1)
    EPSABS = 1.0e-20

    100 IC = IC + 1
    IF (IC > NTMAX) THEN
        ICODE = 1
        GOTO 9000
    ENDIF

    ZM = 0.5*(ZL + ZU)
    IF (ZM < EPSABS) THEN
        ICODE = 2
        GOTO 9000
    ENDIF
    RM = 0.5*(RL + RU)
    CALL GETVAR (UTW,UM,RM)
    CALL COMUFR (UST,UM,ZM,AKVIS,UST1,IEL,IFC,JWX,JWY)
    ZPM =  UST * ZM / AKVIS

    IF (ZPM < ZPLDAT) THEN
        ZL = ZM
        RL = RM
        UL = UM
        GOTO 100
    ELSEIF (ZPM > ZPUDAT) THEN
        ZU = ZM
        RU = RM
        UU = UM
        GOTO 100
    ELSE
        UW = UST
        ZW = ZM
    ENDIF

    RETURN

    9000 CONTINUE
    WRITE (6,*) '   '
    WRITE (6,*) ' *** ABNORMAL TERMINATION ***'
    WRITE (6,'(1X,A4,I5,5X,A6,I5,5X,A9,I3,A1,I3)') &
    'EL =',IEL,'FACE =',IFC,'WALL PT =',JWX,',',JWY
    WRITE (6,*) ' SUBGRID INTERPOLATION OF FRICTIONAL VELOCITY FAILED'
    IF (ICODE == 1) &
    WRITE (6,*) ' MAX NUMBER OF SUBDIVISION (',NTMAX,')  IS REACHED.'
    IF (ICODE == 2) &
    WRITE (6,*) ' DISTANCE FROM WALL VANISHES ( =',ZM,').'
    CALL EMERXIT
    call exitt

    END SUBROUTINE SUBGRID
    SUBROUTINE COMUFR (UST,U,Z,AKVIS,UST1,IEL,IFC,JWX,JWY)

    INCLUDE 'SIZE'
    INCLUDE 'TURBO'

    IC     = 0
    NTMAX  = 1000
    EPS    = 1.E-05
    EPSABS = 1.E-20

    AU = ABS(U)
    AZ = ABS(Z)
    US = ABS(UST1)

    100 IC = IC + 1
    IF (IC > NTMAX) THEN
        ICODE = 1
        GOTO 9000
    ENDIF

    TL = LOG( BTI*AZ*US/AKVIS )
    FF = US*TL - VKC*AU
    FP = 1.0 + TL
    IF (FP < EPSABS) THEN
        ICODE = 2
        GOTO 9000
    ENDIF
    DU = -FF/FP
    DE = ABS( DU/US )
    IF (DE < EPS) GOTO 500
    US = US + DU
    GOTO 100

    500 CONTINUE
    IF (US < -EPSABS) THEN
        ICODE = 3
        GOTO 9000
    ENDIF

    UST = US

    RETURN

    9000 CONTINUE
    WRITE (6,*) '   '
    WRITE (6,*) ' *** ABNORMAL TERMINATION ***'
    WRITE (6,'(1X,A4,I5,5X,A6,I5,5X,A9,I3,A1,I3)') &
    'EL =',IEL,'FACE =',IFC,'WALL PT =',JWX,',',JWY
    WRITE (6,*) ' ITERATION OF FRICTIONAL VELOCITY FAILED.'
    IF (ICODE == 1) &
    WRITE (6,*) ' MAX NUMBER OF ITERATION =',NTMAX,'  IS REACHED.'
    IF (ICODE == 2) &
    WRITE (6,*) ' DIVERGENT NEWTON-RAPHSON DERIVATIVE ( =',FP,' ).'
    IF (ICODE == 3) &
    WRITE (6,*) ' NEGATIVE FRICTIONAL VELOCITY ( =',US,' ).'
    CALL EMERXIT
    call exitt

    END SUBROUTINE COMUFR
    SUBROUTINE COMLSQ

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'MASS'
    INCLUDE 'TURBO'
    INCLUDE 'TSTEP'

    NTOT1 = NX1*NY1*NZ1*NELV
    TLFAC = PARAM(49)

    IF (IFCWUZ) THEN
        CALL LSCALE
        CL = GLMAX (TURBL,NTOT1)
        IF (ISTEP == 0) TLMAX = SQRT( CL )
        TLIMUL = TLFAC / SQRT( CL )
        CALL CMULT (TURBL,TLIMUL,NTOT1)
    ELSE
        TLMAX = VOLVM1**( 1./NDIM )
        TLENG = TLFAC*TLMAX
        CALL CFILL (TURBL,TLENG,NTOT1)
    ENDIF

    RETURN
    END SUBROUTINE COMLSQ
    SUBROUTINE LSCALE

    INCLUDE 'SIZE'
    INCLUDE 'EIGEN'
    INCLUDE 'MASS'
    INCLUDE 'SOLN'
    INCLUDE 'TURBO'
    INCLUDE 'TSTEP'
    INCLUDE 'INPUT'
    COMMON /SCRNS/ RHS(LX1,LY1,LZ1,LELV) &
    , H1 (LX1,LY1,LZ1,LELV) &
    , H2 (LX1,LY1,LZ1,LELV) &
    , SMS(LX1,LY1,LZ1,LELV) &
    , TA (LX1,LY1,LZ1,LELV)

    CHARACTER NAME*4

    IMESH = 1
    NTOT1 = NX1*NY1*NZ1*NELV
    TOLLS = TOLREL / SQRT(EIGAA)
    TTINV = 1.0 / TLIMUL
    MAXIT = 1000
    ISD   = 1
    NAME  = 'NOMG'

    CONST = 2.0
    IF (IFAXIS) CONST = 4.0
    CALL CFILL   (RHS,CONST,NTOT1)
    CALL COL2    (RHS,BM1,NTOT1)
    CALL RONE    (H1, NTOT1)
    CALL RZERO   (H2, NTOT1)
    CALL TLMASK  (SMS)
    CALL CMULT   (TURBL,TTINV,NTOT1)
    CALL BCDIRTL (TURBL,SMS,TA)
    CALL AXHELM  (TA,TURBL,H1,H2,IMESH,ISD)
    CALL SUB2    (RHS,TA,NTOT1)
    CALL HMHOLTZ (NAME,TA,RHS,H1,H2,SMS,VMULT, &
    IMESH,TOLLS,MAXIT,ISD)
    CALL ADD2    (TURBL,TA,NTOT1)

    RETURN
    END SUBROUTINE LSCALE
    SUBROUTINE TLMASK (SMASK)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'

    DIMENSION SMASK(LX1,LY1,LZ1,LELV)
    CHARACTER CB*3

    IFLD  = 1
    NFACE = 2*NDIM
    NTOT1 = NX1*NY1*NZ1*NELV
    CALL RONE (SMASK,NTOT1)

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFLD)
            IF (CB /= 'WS ' .AND. CB /= 'ws '   .AND. &
            CB /= 'WSL' .AND. CB /= 'wsl' ) GOTO 100
            CALL FACEV (SMASK,IEL,IFC,0.0,NX1,NY1,NZ1)
    100 END DO

    CALL DSOP (SMASK,'MUL',NX1,NX1,NZ1)

    RETURN
    END SUBROUTINE TLMASK
    SUBROUTINE BCDIRTL (TLS,SMS,TMP)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TSTEP'

    DIMENSION TLS(LX1,LY1,LZ1,LELV) &
    , SMS(LX1,LY1,LZ1,LELV) &
    , TMP(LX1,LY1,LZ1,LELV)
    CHARACTER CB*3

    IFLD  = 1
    NFACE = 2*NDIM
    NTOT1 = NX1*NY1*NZ1*NELV
    CALL RZERO (TMP,NTOT1)

    DO 1000 ISWEEP=1,2
        DO 100 IEL=1,NELV
            DO 100 IFC=1,NFACE
                CB = CBC(IFC,IEL,IFLD)
                IF (CB /= 'WS ' .AND. CB /= 'ws '   .AND. &
                CB /= 'WSL' .AND. CB /= 'wsl' ) GOTO 100
                IF (ISTEP == 0) THEN
                    CALL FACEV  (TMP,IEL,IFC,0.0,NX1,NY1,NZ1)
                ELSE
                    CALL FACEWL (TMP(1,1,1,IEL),IEL,IFC)
                ENDIF
        100 END DO
        IF (ISWEEP == 1) CALL DSOP (TMP,'MXA',NX1,NY1,NZ1)
        IF (ISWEEP == 2) CALL DSOP (TMP,'MNA',NX1,NY1,NZ1)
    1000 END DO

    CALL COL2 (TLS,SMS,NTOT1)
    CALL ADD2 (TLS,TMP,NTOT1)

    RETURN
    END SUBROUTINE BCDIRTL
    SUBROUTINE FACEWL (S,IEL,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'INPUT'
    INCLUDE 'TURBO'
    DIMENSION S(LX1,LY1,LZ1)

    TLFAC = PARAM(49)
    ZLMAX = 0.09 * TLMAX
    ZFAC  = 2.00 * TLFAC / TLIMUL

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
    I = 0

    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I   = I + 1
            ZL  = MIN( ZWALL(I,1,IFC,IEL),ZLMAX )
            S(J1,J2,1) = ZL * ZFAC
    100 END DO

    RETURN
    END SUBROUTINE FACEWL
    SUBROUTINE GETVAR (V,VP,RP)

    INCLUDE 'SIZE'
    INCLUDE 'WZ'
    DIMENSION V(LX1)

    VP = 0.0

    DO 100 I=1,NX1
        HI = HGLL(I,RP,ZGM1,NX1)
        VP = VP + HI*V(I)
    100 END DO

    RETURN
    END SUBROUTINE GETVAR
    SUBROUTINE INVCHK2 (A,B,N)
    DIMENSION A(1),B(1)
    DO 100 I=1,N
        IF (B(I) > 1.e-20) THEN
            A(I)=A(I)/B(I)
        ELSE
            A(I)=0.0
        ENDIF
    100 END DO
    RETURN
    END SUBROUTINE INVCHK2
    SUBROUTINE FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    INCLUDE 'SIZE'
    INCLUDE 'TOPOL'

    CALL DSSET (NX1,NY1,NZ1)
    IFACE  = EFACE1(IFC)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    RETURN
    END SUBROUTINE FACIND2
