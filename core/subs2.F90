    SUBROUTINE SETAXW1 (IFAXWG)

    use size_m
    use wz_m

    LOGICAL :: IFAXWG

    IF (IFAXWG) THEN
        CALL COPY (W3M1,W2AM1,NX1*NY1)
    ELSE
        CALL COPY (W3M1,W2CM1,NX1*NY1)
    ENDIF

    RETURN
    END SUBROUTINE SETAXW1
    SUBROUTINE SETAXW2 (IFAXWG)

    use size_m
    use wz_m

    LOGICAL :: IFAXWG

    IF (IFAXWG) THEN
        CALL COPY (W3M2,W2AM2,NX2*NY2)
    ELSE
        CALL COPY (W3M2,W2CM2,NX2*NY2)
    ENDIF

    RETURN
    END SUBROUTINE SETAXW2
    SUBROUTINE STNRINV

!     Calculate 2nd and 3rd strain-rate invariants

    use size_m
    use soln
    use tstep
    common /screv/ ei2(lx1,ly1,lz1,lelt) &
    , ei3(lx1,ly1,lz1,lelt)
    common /ctmp1/ exx(lx1,ly1,lz1,lelt) &
    , exy(lx1,ly1,lz1,lelt) &
    , eyy(lx1,ly1,lz1,lelt) &
    , ezz(lx1,ly1,lz1,lelt)
    common /ctmp0/ exz(lx1,ly1,lz1,lelt) &
    , eyz(lx1,ly1,lz1,lelt)

    NTOT1  = NX1*NY1*NZ1*NELV
    CALL RZERO (EI2,NTOT1)
    CALL RZERO (EI3,NTOT1)
    IF (ISTEP == 0) RETURN

    MATMOD = 0
    CALL STNRATE (VX,VY,VZ,NELV,MATMOD)

    IF (NDIM == 2) THEN
        CALL COL3    (EI2,EXX,EYY,NTOT1)
        CALL SUBCOL3 (EI2,EXY,EXY,NTOT1)
        CALL RZERO   (EI3,NTOT1)
    ELSE
        CONST = 2.0
        CALL COL4    (EI3,EXX,EYY,EZZ,NTOT1)
        CALL COL4    (EI2,EXY,EXZ,EYZ,NTOT1)
        CALL ADD2S2  (EI3,EI2,CONST,NTOT1)
        CALL SUBCOL4 (EI3,EXX,EYZ,EYZ,NTOT1)
        CALL SUBCOL4 (EI3,EYY,EXZ,EXZ,NTOT1)
        CALL SUBCOL4 (EI3,EZZ,EXY,EXY,NTOT1)
        CALL COL3    (EI2,EXX,EYY,NTOT1)
        CALL ADDCOL3 (EI2,EXX,EZZ,NTOT1)
        CALL ADDCOL3 (EI2,EYY,EZZ,NTOT1)
        CALL SUBCOL3 (EI2,EXY,EXY,NTOT1)
        CALL SUBCOL3 (EI2,EXZ,EXZ,NTOT1)
        CALL SUBCOL3 (EI2,EYZ,EYZ,NTOT1)
    ENDIF

    RETURN
    END SUBROUTINE STNRINV
    SUBROUTINE OPDOT (DP,A1,A2,A3,B1,B2,B3,N)

    use size_m

    DIMENSION DP(LX1,LY1,LZ1,1) &
    , A1(LX1,LY1,LZ1,1) &
    , A2(LX1,LY1,LZ1,1) &
    , A3(LX1,LY1,LZ1,1) &
    , B1(LX1,LY1,LZ1,1) &
    , B2(LX1,LY1,LZ1,1) &
    , B3(LX1,LY1,LZ1,1)

    IF (NDIM == 2) THEN
        CALL VDOT2 (DP,A1,A2,B1,B2,N)
    ELSE
        CALL VDOT3 (DP,A1,A2,A3,B1,B2,B3,N)
    ENDIF

    RETURN
    END SUBROUTINE OPDOT
    SUBROUTINE OPADDS (A1,A2,A3,B1,B2,B3,CONST,N,ISC)

    use size_m

    DIMENSION A1(LX1,LY1,LZ1,1) &
    , A2(LX1,LY1,LZ1,1) &
    , A3(LX1,LY1,LZ1,1) &
    , B1(LX1,LY1,LZ1,1) &
    , B2(LX1,LY1,LZ1,1) &
    , B3(LX1,LY1,LZ1,1)

    IF (ISC == 1) THEN
        CALL ADD2S1 (A1,B1,CONST,N)
        CALL ADD2S1 (A2,B2,CONST,N)
        IF (NDIM == 3) CALL ADD2S1 (A3,B3,CONST,N)
    ELSEIF (ISC == 2) THEN
        CALL ADD2S2 (A1,B1,CONST,N)
        CALL ADD2S2 (A2,B2,CONST,N)
        IF (NDIM == 3) CALL ADD2S2 (A3,B3,CONST,N)
    ENDIF

    RETURN
    END SUBROUTINE OPADDS
    SUBROUTINE FACEXS (A,B,IFACE1,IOP)

!     IOP = 0
!     Extract scalar A from B on face IFACE1.

!     IOP = 1
!     Extract scalar B from A on face IFACE1.

!     A has the (NX,NY,NFACE) data structure
!     B has the (NX,NY,NZ)    data structure
!     IFACE1 is in the preprocessor notation
!     IFACE  is the dssum notation.

    use size_m
    use topol

    DIMENSION A(LX1,LY1),B(LX1,LY1,LZ1)

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    IF (IOP == 0) THEN
        DO 100 J2=JS2,JF2,JSKIP2
            DO 100 J1=JS1,JF1,JSKIP1
                I = I+1
                A(I,1) = B(J1,J2,1)
        100 END DO
    ELSE
        DO 150 J2=JS2,JF2,JSKIP2
            DO 150 J1=JS1,JF1,JSKIP1
                I = I+1
                B(J1,J2,1) = A(I,1)
        150 END DO
    ENDIF

    RETURN
    END SUBROUTINE FACEXS
    SUBROUTINE FACEXV (A1,A2,A3,B1,B2,B3,IFACE1,IOP)

!     IOP = 0
!     Extract vector (A1,A2,A3) from (B1,B2,B3) on face IFACE1.

!     IOP = 1
!     Extract vector (B1,B2,B3) from (A1,A2,A3) on face IFACE1.

!     A1, A2, A3 have the (NX,NY,NFACE) data structure
!     B1, B2, B3 have the (NX,NY,NZ)    data structure
!     IFACE1 is in the preprocessor notation
!     IFACE  is the dssum notation.

    use size_m
    use topol

    DIMENSION A1(LX1,LY1),A2(LX1,LY1),A3(LX1,LY1), &
    B1(LX1,LY1,LZ1),B2(LX1,LY1,LZ1),B3(LX1,LY1,LZ1)

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)
    I = 0

    IF (IOP == 0) THEN
        DO 100 J2=JS2,JF2,JSKIP2
            DO 100 J1=JS1,JF1,JSKIP1
                I = I+1
                A1(I,1) = B1(J1,J2,1)
                A2(I,1) = B2(J1,J2,1)
                A3(I,1) = B3(J1,J2,1)
        100 END DO
    ELSE
        DO 150 J2=JS2,JF2,JSKIP2
            DO 150 J1=JS1,JF1,JSKIP1
                I = I+1
                B1(J1,J2,1) = A1(I,1)
                B2(J1,J2,1) = A2(I,1)
                B3(J1,J2,1) = A3(I,1)
        150 END DO
    ENDIF

    RETURN
    END SUBROUTINE FACEXV
    SUBROUTINE FACSUB2 (A1,A2,A3,B1,B2,B3,IFACE1)

!     Subtract B1,B2,B3 from A1,A2,A3 on surface IFACE1 of element IE.

!     A1, A2, A3 have the (NX,NY,NZ)    data structure
!     B1, B2, B3 have the (NX,NY,NFACE) data structure
!     IFACE1 is in the preprocessor notation
!     IFACE  is the dssum notation.

    use size_m
    use topol

    DIMENSION A1(LX1,LY1,LZ1),A2(LX1,LY1,LZ1),A3(LX1,LY1,LZ1), &
    B1(LX1,LY1),B2(LX1,LY1),B3(LX1,LY1)

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A1(J1,J2,1) = A1(J1,J2,1) - B1(I,1)
            A2(J1,J2,1) = A2(J1,J2,1) - B2(I,1)
            A3(J1,J2,1) = A3(J1,J2,1) - B3(I,1)
    100 END DO

    RETURN
    END SUBROUTINE FACSUB2

#if 0
    SUBROUTINE GAMMASF (H1,H2)
!-----------------------------------------------------------------------

!     Compute lagrest eigenvalue of coupled Helmholtz operator

!-----------------------------------------------------------------------
    use size_m
    use eigen
    use input
    use mass
    use mvgeom
    use soln
    use tstep
    use wz_m
    common /scrmg/ ae1(lx1,ly1,lz1,lelv) &
    , ae2(lx1,ly1,lz1,lelv) &
    , ae3(lx1,ly1,lz1,lelv)
    common /scruz/ e1(lx1,ly1,lz1,lelv) &
    , e2(lx1,ly1,lz1,lelv) &
    , e3(lx1,ly1,lz1,lelv)

    DIMENSION H1(LX1,LY1,LZ1,1),H2(LX1,LY1,LZ1,1)

    NTOT1  = NX1*NY1*NZ1*NELV
    IMESH  = 1
    MATMOD = 0

    IF (ISTEP == 0) THEN
        EIGGA = 0.0
        CALL STX1SF
    ELSE
        CALL COPY (E1,EV1,NTOT1)
        CALL COPY (E2,EV2,NTOT1)
        IF (NDIM == 3) CALL COPY (E3,EV3,NTOT1)
    ENDIF

    EVNEW = EIGGA

    DO 1000 ITER=1,NMXE
    
        CALL AXHMSF  (AE1,AE2,AE3,E1,E2,E3,H1,H2,MATMOD)
        CALL RMASK   (AE1,AE2,AE3,NELV)
        CALL OPDSSUM (AE1,AE2,AE3)
    
        EVOLD = EVNEW
        EVNEW = GLSC3(E1,AE1,VMULT,NTOT1) + GLSC3(E2,AE2,VMULT,NTOT1)
        IF (NDIM == 3) EVNEW = EVNEW + GLSC3(E3,AE3,VMULT,NTOT1)
        CRIT = ABS( (EVNEW - EVOLD)/EVNEW )
        IF ( CRIT < TOLEV ) GOTO 2000
    
        CALL COL3 (E1,BINVM1,AE1,NTOT1)
        CALL COL3 (E2,BINVM1,AE2,NTOT1)
        IF (NDIM == 3) CALL COL3 (E3,BINVM1,AE3,NTOT1)
        XX = GLSC3(E1,AE1,VMULT,NTOT1) + GLSC3(E2,AE2,VMULT,NTOT1)
        IF (NDIM == 3) XX = XX + GLSC3(E3,AE3,VMULT,NTOT1)
        IF (XX < 0.0) GO TO 9000
    
        XNORM=1./SQRT( XX )
        CALL CMULT (E1,XNORM,NTOT1)
        CALL CMULT (E2,XNORM,NTOT1)
        IF (NDIM == 3) CALL CMULT (E3,XNORM,NTOT1)
    
    1000 END DO

!     Save eigenvalue for all cases.
!     Save eigenvectors for deforming geometries.

    2000 EIGGA = EVNEW
    IF (IFMVBD) THEN
        CALL COPY (EV1,E1,NTOT1)
        CALL COPY (EV2,E2,NTOT1)
        IF (NDIM == 3) CALL COPY (EV3,E3,NTOT1)
    ENDIF

    RETURN

    9000 CONTINUE
    IF (NID == 0) &
    WRITE ( 6,*) ' Non +ve def. A-operator detected during eigenvalue &
    computation : tran(x)Ax =',XX
    CALL EMERXIT
    call exitt

    END SUBROUTINE GAMMASF
#endif


    SUBROUTINE CMULT2 (A,B,CONST,N)
    DIMENSION A(1),B(1)
    DO 100 I=1,N
        A(I)=B(I)*CONST
    100 END DO
    RETURN
    END SUBROUTINE CMULT2
    SUBROUTINE ADD3S (A,B,C,CONST,N)
    DIMENSION A(1),B(1),C(1)
    DO 100 I=1,N
        A(I)=B(I)+CONST*C(I)
    100 END DO
    RETURN
    END SUBROUTINE ADD3S
    SUBROUTINE EMERXIT

    use size_m
    use input
    use parallel
    use tstep

!     Try to hang in there on the first few time steps (pff 8/92)
    IF (IFTRAN .AND. ISTEP < 9) RETURN

    LASTEP = 1
    CALL PREPOST( .TRUE. ,'   ')

    IF (NP == 1) THEN
        WRITE (6,*) '       '
        WRITE (6,*) &
        ' Emergency exit:',ISTEP,'   time =',TIME
        WRITE (6,*) &
        ' Latest solution and data are dumped for post-processing.'
        WRITE (6,*) ' *** STOP ***'
    ELSE
        WRITE (6,*) '       '
        WRITE (6,*) NID, &
        ' Emergency exit:',ISTEP,'   time =',TIME
        WRITE (6,*) &
        ' Latest solution and data are dumped for post-processing.'
        WRITE (6,*) ' *** STOP ***'
    ENDIF

    call runstat

    call exitt
    RETURN
    END SUBROUTINE EMERXIT
    SUBROUTINE FACCVS (A1,A2,A3,B,IFACE1)

!     Collocate scalar B with vector A, components A1,A2,A3,
!     on the surface IFACE1 of an element.

!         A1,A2,A3 have the (NX,NY,NZ) data structure
!         B has the (NX,NY,IFACE) data structure
!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.

    use size_m
    use topol
    DIMENSION A1(LX1,LY1,LZ1),A2(LX1,LY1,LZ1),A3(LX1,LY1,LZ1), &
    B(LX1,LY1)

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)
    I = 0

    IF (NDIM == 2) THEN
        DO 100 J2=JS2,JF2,JSKIP2
            DO 100 J1=JS1,JF1,JSKIP1
                I = I+1
                A1(J1,J2,1) = A1(J1,J2,1)*B(I,1)
                A2(J1,J2,1) = A2(J1,J2,1)*B(I,1)
        100 END DO
    ELSE
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                I = I+1
                A1(J1,J2,1) = A1(J1,J2,1)*B(I,1)
                A2(J1,J2,1) = A2(J1,J2,1)*B(I,1)
                A3(J1,J2,1) = A3(J1,J2,1)*B(I,1)
        200 END DO
    ENDIF

    RETURN
    END SUBROUTINE FACCVS
    SUBROUTINE STX1SF
!------------------------------------------------------------------

!     Compute startvector for finding an eigenvalue on mesh 1.
!     Normalization: XT*B*X = 1

!------------------------------------------------------------------
    use size_m
    use mass
    use soln
    common /scrmg/ ae1(lx1,ly1,lz1,lelv) &
    , ae2(lx1,ly1,lz1,lelv) &
    , ae3(lx1,ly1,lz1,lelv)
    common /scruz/ e1(lx1,ly1,lz1,lelv) &
    , e2(lx1,ly1,lz1,lelv) &
    , e3(lx1,ly1,lz1,lelv)

    NTOT1 = NX1*NY1*NZ1*NELV
    CALL RZERO3 (E1 ,E2 ,E3 ,NTOT1)
    CALL RZERO3 (AE1,AE2,AE3,NTOT1)

    CALL COPY  (E1,BM1,NTOT1)
    CALL COPY  (E2,BM1,NTOT1)
    IF (NDIM == 3) CALL COPY  (E3,BM1,NTOT1)

    CALL RMASK (E1,E2,E3,NELV)
    CALL COL3  (AE1,BM1,E1,NTOT1)
    CALL COL3  (AE2,BM1,E2,NTOT1)
    IF (NDIM == 3) CALL COL3  (AE3,BM1,E3,NTOT1)

    CALL OPDSSUM (AE1,AE2,AE3)

    XX = GLSC3 (E1,AE1,VMULT,NTOT1) + GLSC3 (E2,AE2,VMULT,NTOT1)
    IF (NDIM == 3) XX = XX + GLSC3 (E3,AE3,VMULT,NTOT1)
    XNORM  = 1./SQRT(XX)
    CALL CMULT (E1,XNORM,NTOT1)
    CALL CMULT (E2,XNORM,NTOT1)
    IF (NDIM == 3) CALL CMULT (E3,XNORM,NTOT1)

!     call exitti   ('quit in stx1sf$,',nel)


    RETURN
    END SUBROUTINE STX1SF

    SUBROUTINE SOLVEL

    use size_m
    use geom
    use soln
    use tstep

    DO 100 IEL=1,NELV
        DO 100 K=1,NZ1
            DO 100 J=1,NY1
                DO 100 I=1,NX1
                    CALL VSOLN (VX (I,J,K,IEL),VY (I,J,K,IEL),VZ (I,J,K,IEL), &
                    XM1(I,J,K,IEL),YM1(I,J,K,IEL),ZM1(I,J,K,IEL),PI)
    100 END DO

    RETURN
    END SUBROUTINE SOLVEL
    SUBROUTINE VSOLN (UX,UY,UZ,X,Y,Z,PI)

!       URR=(1.-0.75/SQRT(X**2+Y**2)+0.0625/(SQRT(X**2+Y**2)**3))*(X/
!     $     SQRT(X**2+Y**2))
!       UTETA=-(1.-0.375/SQRT(X**2+Y**2)-0.03125/(SQRT(X**2+Y**2)**3))*
!     $       (Y/SQRT(X**2+Y**2))
!       UX=URR*(X/SQRT(X**2+Y**2))-UTETA*(Y/SQRT(X**2+Y**2))
!       UY=URR*(Y/SQRT(X**2+Y**2))+UTETA*(X/SQRT(X**2+Y**2))

!       UX = 2.*COS( PI*X )
!       UY = PI*Y*SIN( PI*X )

    UX = 0.0
    UY = 0.0

    RETURN
    END SUBROUTINE VSOLN
    SUBROUTINE SOLPRES

    use size_m
    use geom
    use soln
    use tstep

    DO 100 IEL=1,NELV
        DO 100 K=1,NZ2
            DO 100 J=1,NY2
                DO 100 I=1,NX2
                    CALL PRSOLN (PR (I,J,K,IEL),XM2(I,J,K,IEL),YM2(I,J,K,IEL), &
                    ZM2(I,J,K,IEL),PI)
    100 END DO

    RETURN
    END SUBROUTINE SOLPRES
    SUBROUTINE PRSOLN (P,X,Y,Z,PI)

!      R  = SQRT( X**2 + Y**2 )
!      CS = X/R
!      P  = -0.75 * CS / R**2

!      P  = -SIN( PI*X )*COS( PI*Y )

    P = 0.0

    RETURN
    END SUBROUTINE PRSOLN
    SUBROUTINE STORE
!-----------------------------------------------------------------------

!     Store the results

!-----------------------------------------------------------------------
    use size_m
    INCLUDE 'TOTAL'
    LOGICAL :: IFPREL(LELT)

    NZ1I   =  1
    NZ1J   = NZ1
    NZ1INC =  1
    NZ2I   =  1
    NZ2J   = NZ2
    NZ2INC =  1
    NFACE  = 2*NDIM

    CALL LFALSE (IFPREL,NELTOT)

    IFPREL(1)= .TRUE. 
    IFPREL(2)= .TRUE. 

    IF (IFFLOW) THEN
        WRITE (21,*) ' '
        WRITE (21,*) 'FLUID FLOW FIELD'
        DO 9001 IEL = 1,NELV
            IF ( .NOT. IFPREL(IEL) ) GOTO 9001
            WRITE (21,*) 'ELEMENT NUMBER ',IEL
            DO 101 IPL=NZ1I,NZ1J,NZ1INC
                CALL OUTM1 (VX,'VX        ',NZ1,IEL,IPL)
            101 END DO
            DO 102 IPL=NZ1I,NZ1J,NZ1INC
                CALL OUTM1 (VY,'VY        ',NZ1,IEL,IPL)
            102 END DO
        !           DO 103 IPL=NZ1I,NZ1J,NZ1INC
        !              CALL OUTM1 (VZ,'VZ        ',NZ1,IEL,IPL)
        ! 103       CONTINUE
        
        !           DO 310 IPL=NZ2I,NZ2J,NZ2INC
        !              CALL OUTM2 (PR,'PR        ',NZ2,IEL,IPL)
        ! 310       CONTINUE
        
            IF (IFMODEL .AND. IFSWALL) THEN
                DO 410 IFC=1,NFACE
                    CALL OUTF1 (UWALL,'U*-(1)    ',IEL,IFC)
                410 END DO
                DO 430 IFC=1,NFACE
                    CALL OUTF1 (ZWALL,'Z0        ',IEL,IFC)
                430 END DO
            ENDIF
            IF (IFMODEL .AND. .NOT. IFKEPS) THEN
                DO 450 IPL=NZ1I,NZ1J,NZ1INC
                    CALL OUTM1 (TURBL,'TURBLENGTH',NZ1,IEL,IPL)
                450 END DO
            ENDIF
        
        9001 END DO
    ENDIF

    IF (IFHEAT) THEN
        DO 500 IFIELD=2,NFIELD
            DO 9002 IEL = 1,NELT
                IF ( .NOT. IFPREL(IEL) ) GOTO 9002
                WRITE (21,*) ' '
                WRITE (21,*) 'FIELD   NUMBER ',IFIELD
                WRITE (21,*) 'ELEMENT NUMBER ',IEL
                DO 601 IPL=NZ1I,NZ1J,NZ1INC
                    CALL OUTM1 (T(1,1,1,1,IFIELD-1),'T         ',NZ1,IEL,IPL)
                601 END DO
            9002 END DO
        500 END DO
    ENDIF

    RETURN
    END SUBROUTINE STORE
    SUBROUTINE PRINTEL (TA,A,IEL)

    use size_m
    DIMENSION TA(LX1,LY1,LZ1,LELT)
    CHARACTER A*10

    NZ1I   =  1
    NZ1J   = NZ1
    NZ1INC =  1

    WRITE (21,*) 'ELEMENT NUMBER ',IEL
    DO 101 IPL=NZ1I,NZ1J,NZ1INC
        CALL OUTM1 (TA,A,NZ1,IEL,IPL)
    101 END DO

    RETURN
    END SUBROUTINE PRINTEL
    SUBROUTINE PRINTV (TA,A,NEL)
!-----------------------------------------------------------------------

!     Store the results

!-----------------------------------------------------------------------
    use size_m
    DIMENSION TA(LX1,LY1,LZ1,LELT)
    CHARACTER A*10

    NZ1I   =  1
    NZ1J   = NZ1
    NZ1INC =  1

    DO 9001 IEL = 1,NEL
        WRITE (21,*) 'ELEMENT NUMBER ',IEL
        DO 101 IPL=NZ1I,NZ1J,NZ1INC
            CALL OUTM1 (TA,A,NZ1,IEL,IPL)
        101 END DO
    9001 END DO

    RETURN
    END SUBROUTINE PRINTV
    SUBROUTINE OUTF1 (X,TXT,IEL,IFC)
    use size_m
    DIMENSION X(LX1,LZ1,6,LELT)
    CHARACTER(10) :: TXT

    NFACE = 2*NDIM
    NZI   = NZ1
    NZJ   =  1
    NZINC = -1
    NXI   =  1
    NXJ   = NX1
    NXINC =  1

    WRITE(21,106) TXT,IFC,NFACE
    DO 100 J=NZI,NZJ,NZINC
        WRITE(21,105) (X(I,J,IFC,IEL),I=NXI,NXJ,NXINC)
    100 END DO

    105 FORMAT(5E15.6)
    106 FORMAT(///,5X,'     ^              ',/, &
    &            5X,'   S |              ',/, &
    &            5X,'     |              ',A10,/, &
    &            5X,'     +---->         ','Plane = ',I2,'/',I2,/, &
    &            5X,'       R            ',/)

    RETURN
    END SUBROUTINE OUTF1
    SUBROUTINE OUTM1 (X,TXT,NP,IEL,IP)
    use size_m
    DIMENSION X(LX1,LY1,LZ1,LELT)
    CHARACTER(10) :: TXT

    NYI   = NY1
    NYJ   =  1
    NYINC = -1
    NXI   =  1
    NXJ   = NX1
    NXINC =  1

    WRITE(6,106) TXT,IP,NP
    DO 100 J=NYI,NYJ,NYINC
        WRITE(6,105) (X(I,J,IP,IEL),I=NXI,NXJ,NXINC)
    100 END DO

! 105 FORMAT(1p8e10.3)
    105 FORMAT(8f10.3)
    106 FORMAT(///,5X,'     ^              ',/, &
    &            5X,'   Y |              ',/, &
    &            5X,'     |              ',A10,/, &
    &            5X,'     +---->         ','Plane = ',I2,'/',I2,/, &
    &            5X,'       X            ',/)

    RETURN
    END SUBROUTINE OUTM1
    SUBROUTINE OUTM2 (X,TXT,NP,IEL,IP)
    use size_m
    DIMENSION X(LX2,LY2,LZ2,LELV)
    CHARACTER(10) :: TXT

    NYI   = NY2
    NYJ   =  1
    NYINC = -1
    NXI   =  1
    NXJ   = NX2
    NXINC =  1

    WRITE(21,106) TXT,IP,NP
    DO 100 J=NYI,NYJ,NYINC
        WRITE(21,105) (X(I,J,IP,IEL),I=NXI,NXJ,NXINC)
    100 END DO

    105 FORMAT(5E15.6)
    106 FORMAT(///,5X,'     ^              ',/, &
    &            5X,'   Y |              ',/, &
    &            5X,'     |              ',A10,/, &
    &            5X,'     +---->         ','Plane = ',I2,'/',I2,/, &
    &            5X,'       X            ',/)

    RETURN
    END SUBROUTINE OUTM2
    SUBROUTINE STSMASK (C1MASK,C2MASK,C3MASK)

    use size_m
    use geom
    use input
    use tstep
    common /screv/ hfmask(lx1,lz1,6,lelt) &
    , hvmask(lx1,ly1,lz1,lelt)

    DIMENSION C1MASK(LX1,LY1,LZ1,1) &
    , C2MASK(LX1,LY1,LZ1,1) &
    , C3MASK(LX1,LY1,LZ1,1)
    INTEGER ::   IMDATA
    SAVE      IMDATA
    DATA      IMDATA /0/

    IFLD = IFIELD
    NEL  = NELFLD(IFIELD)

    IF (IMDATA == 0) THEN
        CALL SETCDAT
        IMDATA=1
    ENDIF

    IF (IFLD == 1) CALL SKIPCNR (NEL)
    CALL SETHMSK (HVMASK,HFMASK,IFLD,NEL)
    CALL SETMLOG (HVMASK,HFMASK,IFLD,NEL)
    CALL SETMASK (C1MASK,C2MASK,C3MASK,HVMASK,NEL)
    IF (IFLMSF(IFLD)) CALL SETCSYS (HVMASK,HFMASK,NEL)
    IF (IFLD == 0)    CALL FIXWMSK (C2MASK,C3MASK,HVMASK,HFMASK,NEL)

    if (ifaxis .AND. ifld == 1)  call fixmska (c1mask,c2mask,c3mask)

    RETURN
    END SUBROUTINE STSMASK
    SUBROUTINE UPDMSYS (IFLD)

    use size_m
    use geom
    use tstep
    common /screv/ hfmask(lx1,lz1,6,lelt) &
    , hvmask(lx1,ly1,lz1,lelt)

    IF ( .NOT. IFLMSF(IFLD)) RETURN

    NEL  = NELFLD(IFLD)
    CALL SETHMSK (HVMASK,HFMASK,IFLD,NEL)
    CALL SETCSYS (HVMASK,HFMASK,NEL)

    RETURN
    END SUBROUTINE UPDMSYS
    SUBROUTINE SETHMSK (HVMASK,HFMASK,IFLD,NEL)

    use size_m
    use input
    use tstep

    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)
    CHARACTER CB*3

    NTOT1 = NX1*NY1*NZ1*NEL
    NXZ1  = NX1*NZ1
    NTOTF = NXZ1*6*NEL
    NFACE = 2*NDIM
    CONST = 5.0
    CALL CFILL (HVMASK,CONST,NTOT1)
    CALL CFILL (HFMASK,CONST,NTOTF)

    IF (IFLD == 1) THEN
    
        DO 110 IEL=1,NEL
            DO 110 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'ON ' .OR. CB == 'on ' .OR. &
                CB == 'MM ' .OR. CB == 'mm ' ) THEN
                    CALL FACEV (HVMASK,IEL,IFC,3.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),3.0,NXZ1)
                ENDIF
        110 END DO
    
        DO 120 IEL=1,NEL
            DO 120 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'SYM' .OR. CB == 'A  ' .OR. CB == 'WS ' .OR. &
                CB == 'ws ' .OR. CB == 'WSL' .OR. CB == 'wsl' .OR. &
                CB == 'SH ' .OR. CB == 'sh ' .OR. CB == 'SHL' .OR. &
                CB == 'shl')                                  THEN
                    CALL FACEV (HVMASK,IEL,IFC,2.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),2.0,NXZ1)
                ENDIF
        120 END DO
    
        DO 130 IEL=1,NEL
            DO 130 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'MF ' .OR. CB == 'V  ' .OR. CB == 'v  ' .OR. &
                CB == 'VL ' .OR. CB == 'vl ' .OR. CB(1:2) == 'mv') THEN
                    CALL FACEV (HVMASK,IEL,IFC,1.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),1.0,NXZ1)
                ENDIF
        130 END DO
    
        DO 140 IEL=1,NEL
            DO 140 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'W  ') THEN
                    CALL FACEV (HVMASK,IEL,IFC,0.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),0.0,NXZ1)
                ENDIF
        140 END DO
    
    ELSE
    
        DO 210 IEL=1,NEL
            DO 210 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'SYM') THEN
                    CALL FACEV (HVMASK,IEL,IFC,2.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),2.0,NXZ1)
                ENDIF
        210 END DO
    
    !     write(6,*) 'MASK this is ifield:',ifield
        DO 220 IEL=1,NEL
            DO 220 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB(1:1) == 'M' .OR. CB(1:1) == 'm') THEN
                !            CALL FACEV (HVMASK,IEL,IFC,1.0,NX1,NY1,NZ1)
                !            CALL CFILL (HFMASK(1,1,IFC,IEL),1.0,NXZ1)
                    CALL FACEV (HVMASK,IEL,IFC,2.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),2.0,NXZ1)
                ENDIF
        220 END DO
    
        DO 230 IEL=1,NEL
            DO 230 IFC=1,NFACE
                CB=CBC(IFC,IEL,IFLD)
                IF (CB == 'FIX') THEN
                    CALL FACEV (HVMASK,IEL,IFC,0.0,NX1,NY1,NZ1)
                    CALL CFILL (HFMASK(1,1,IFC,IEL),0.0,NXZ1)
                ENDIF
        230 END DO
    
    ENDIF

    CALL DSOP (HVMASK,'MNA',NX1,NY1,NZ1)

    RETURN
    END SUBROUTINE SETHMSK
    SUBROUTINE SKIPCNR (NEL)

    use size_m
    use geom
    use input
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    NFACE = 2*NDIM
    NCRFC = NFACE - 2
    NMXCR = 8*NEL
    CALL LFALSE (IFNSKP,NMXCR)

    DO 100 IEL=1,NEL
        DO 100 IFC=1,NFACE
            IF (CDOF(IFC,IEL) == '1') THEN
                ICR=MCRFC(1,IFC)
                IFNSKP(ICR,IEL)= .TRUE. 
            ELSEIF (CDOF(IFC,IEL) == '2') THEN
                ICR=MCRFC(2,IFC)
                IFNSKP(ICR,IEL)= .TRUE. 
            ELSEIF (CDOF(IFC,IEL) == '3') THEN
                ICR=MCRFC(3,IFC)
                IFNSKP(ICR,IEL)= .TRUE. 
            ELSEIF (CDOF(IFC,IEL) == '4') THEN
                ICR=MCRFC(4,IFC)
                IFNSKP(ICR,IEL)= .TRUE. 
            ENDIF
            IF (CDOF(IFC,IEL) == '*') THEN
                DO 160 ICRFC=1,NCRFC
                    ICR=MCRFC(ICRFC,IFC)
                    IFNSKP(ICR,IEL)= .TRUE. 
                160 END DO
            ENDIF
    100 END DO

    RETURN
    END SUBROUTINE SKIPCNR
    SUBROUTINE SETMASK (C1MASK,C2MASK,C3MASK,HVMASK,NEL)

    use size_m

    DIMENSION HVMASK (LX1,LY1,LZ1,1) &
    , C1MASK(LX1,LY1,LZ1,1) &
    , C2MASK(LX1,LY1,LZ1,1) &
    , C3MASK(LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NEL
    CALL RZERO3  (C1MASK,C2MASK,C3MASK,NTOT1)

    DO 100 IEL=1,NEL
        DO 100 IZ=1,NZ1
            DO 100 IY=1,NY1
                DO 100 IX=1,NX1
                    HMV=ABS( HVMASK(IX,IY,IZ,IEL) )
                    IF (HMV > 2.9) THEN
                        C1MASK(IX,IY,IZ,IEL) = 1.0
                    ENDIF
                    IF ((HMV > 1.9 .AND. HMV < 2.1) .OR. HMV > 4.9) THEN
                        C2MASK(IX,IY,IZ,IEL) = 1.0
                    ENDIF
    100 END DO

    IF (NDIM == 3) CALL COPY (C3MASK,C2MASK,NTOT1)

    RETURN
    END SUBROUTINE SETMASK
    SUBROUTINE SETMLOG (HVMASK,HFMASK,IFLD,NEL)

    use size_m
    use geom
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NFACE = 2*NDIM
    NEDGE = 12
    NCRNR = 2**NDIM
    NTOTF = NFACE*NEL
    NTOTS = NEDGE*NEL
    NTOTC = NCRNR*NEL
    EPSA  = 1.E-6

    IFLMSF(IFLD) = .FALSE. 
    IFLMSE(IFLD) = .FALSE. 
    IFLMSC(IFLD) = .FALSE. 
    CALL LFALSE (IFMSFC(1,1,IFLD),NTOTF)
    CALL LFALSE (IFMSEG(1,1,IFLD),NTOTS)
    CALL LFALSE (IFMSCR(1,1,IFLD),NTOTC)

    DO 100 IEL=1,NEL
        DO 100 IFC=1,NFACE
            HMF = ABS( HFMASK(1,1,IFC,IEL) )
            IF (HMF > 1.9  .AND. HMF < 3.1 ) THEN
                IFLMSF(IFLD)         = .TRUE. 
                IFMSFC(IFC,IEL,IFLD) = .TRUE. 
            ENDIF
    100 END DO
    CALL GLLOG(IFLMSF(IFLD), .TRUE. )

    IF (NDIM == 3) THEN
        DO 200 IEL=1,NEL
            DO 200 ISD=1,NEDGE
                IX  = MIDRST(1,ISD)
                IY  = MIDRST(2,ISD)
                IZ  = MIDRST(3,ISD)
                HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
                IF (HMV < 1.9  .OR.  HMV > 3.1) GOTO 200
                IDIFF = 0
                DO 220 II=1,2
                    IFC = MFCEG(II,ISD)
                    HMF = ABS( HFMASK(1,1,IFC,IEL) )
                    IF (ABS(HMV - HMF) > EPSA) IDIFF=IDIFF + 1
                220 END DO
                IF (IDIFF == 2) THEN
                    IFLMSE(IFLD)         = .TRUE. 
                    IFMSEG(ISD,IEL,IFLD) = .TRUE. 
                ENDIF
        200 END DO
        CALL GLLOG(IFLMSE(IFLD), .TRUE. )
    ENDIF

    DO 300 IEL=1,NEL
        DO 300 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            IZ  = MCRRST(3,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR.  HMV > 3.1) GOTO 300
            IDIFF = 0
            DO 330 II=1,NDIM
                IFC = MFCCR(II,ICR)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) > EPSA) IDIFF=IDIFF + 1
            330 END DO
            IF (NDIM == 3) THEN
                DO 360 II=1,NDIM
                    ISD = MEGCR(II,ICR)
                    IXS = MIDRST(1,ISD)
                    IYS = MIDRST(2,ISD)
                    IZS = MIDRST(3,ISD)
                    HMS = ABS( HVMASK(IXS,IYS,IZS,IEL) )
                    IF (ABS(HMV - HMS) > EPSA) IDIFF=IDIFF + 1
                360 END DO
            ENDIF
            IF ( (NDIM == 2 .AND. IDIFF == 2)   .OR. &
            (NDIM == 3 .AND. IDIFF == 6) ) THEN
                IFLMSC(IFLD)         = .TRUE. 
                IFMSCR(ICR,IEL,IFLD) = .TRUE. 
            ENDIF
    300 END DO
    CALL GLLOG(IFLMSC(IFLD), .TRUE. )

    RETURN
    END SUBROUTINE SETMLOG
    SUBROUTINE SETCSYS (HVMASK,HFMASK,NEL)

    use size_m
    use geom
    use tstep
    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NFACE  = 2*NDIM
    NTOT1  = NX1*NY1*NZ1*NEL

    CALL RZERO3  (VNX,VNY,VNZ,NTOT1)
    CALL RZERO3  (V1X,V1Y,V1Z,NTOT1)
    CALL RZERO3  (V2X,V2Y,V2Z,NTOT1)

    DO 10 IEL=1,NEL
        DO 10 IFC=1,NFACE
            HMF = ABS( HFMASK(1,1,IFC,IEL) )
            IF (HMF > 1.9  .AND.  HMF < 3.1) &
            CALL FACEXV(UNX(1,1,IFC,IEL),UNY(1,1,IFC,IEL),UNZ(1,1,IFC,IEL), &
            VNX(1,1,1,IEL),VNY(1,1,1,IEL),VNZ(1,1,1,IEL),IFC,1)
    10 END DO

    IF (NDIM == 2) THEN
        CALL COMAVN2 (HVMASK,HFMASK,NEL)
    ELSE
        CALL COMAVN3 (HVMASK,HFMASK,NEL)
    ENDIF

    RETURN
    END SUBROUTINE SETCSYS
    SUBROUTINE COMAVN2 (HVMASK,HFMASK,NEL)

    use size_m
    use geom
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NTOT1  = NX1*NY1*NZ1*NEL
    NFACE  = 2*NDIM
    NCRNR  = 2**NDIM
    EPSA   = 1.0E-06
    CALL RZERO (VNZ,NTOT1)

    IZ = 1
    DO 100 IEL=1,NEL
        DO 100 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 100
            VNX(IX,IY,IZ,IEL) = 0.0
            VNY(IX,IY,IZ,IEL) = 0.0
    100 END DO

    DO 200 IEL=1,NEL
        DO 200 ICR=1,NCRNR
            IF (IFNSKP(ICR,IEL)) GOTO 200
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 200
            DO 220 II=1,2
                IFC = MFCCR(II,ICR)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    IR = MCRRST(II,ICR)
                    VNX(IX,IY,IZ,IEL)=VNX(IX,IY,IZ,IEL) + UNX(IR,IZ,IFC,IEL)
                    VNY(IX,IY,IZ,IEL)=VNY(IX,IY,IZ,IEL) + UNY(IR,IZ,IFC,IEL)
                ENDIF
            220 END DO
    200 END DO

    CALL DSSUM   (VNX,NX1,NY1,NZ1)
    CALL DSSUM   (VNY,NX1,NY1,NZ1)
    CALL UNITVEC (VNX,VNY,VNZ,NTOT1)

    CALL COPY   (V1Y,VNX,NTOT1)
    CALL COPY   (V1X,VNY,NTOT1)
    CALL CHSIGN (V1X,NTOT1)

    RETURN
    END SUBROUTINE COMAVN2
    SUBROUTINE COMAVN3 (HVMASK,HFMASK,NEL)

    use size_m
    use geom
    common /scrcg/  vnmag(lx1,ly1,lz1,lelt)
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NTOT1  = NX1*NY1*NZ1*NEL
    NFACE  = 2*NDIM
    NCRNR  = 2**NDIM
    NEDGE  = 12
    NMID   =(NX1 + 1)/2
    NXM1   = NX1 - 1
    EPSA   = 1.0E-06
    EPSN   = 1.0E-03

    DO 100 IEL=1,NEL
        DO 100 ISD=1,NEDGE
            IX  = MIDRST(1,ISD)
            IY  = MIDRST(2,ISD)
            IZ  = MIDRST(3,ISD)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 100
            CALL EDGINDV (LV1,LV2,LVSKIP,ISD)
            DO 120 I=2,NXM1
                LV = LV1 + (I-1)*LVSKIP
                VNX(LV,1,1,IEL) = 0.0
                VNY(LV,1,1,IEL) = 0.0
                VNZ(LV,1,1,IEL) = 0.0
            120 END DO
    100 END DO

    DO 150 IEL=1,NEL
        DO 150 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            IZ  = MCRRST(3,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 150
            VNX(IX,IY,IZ,IEL) = 0.0
            VNY(IX,IY,IZ,IEL) = 0.0
            VNZ(IX,IY,IZ,IEL) = 0.0
    150 END DO

!     (1) All Edges

    DO 200 IEL=1,NEL
        DO 200 ISD=1,NEDGE
            IX  = MIDRST(1,ISD)
            IY  = MIDRST(2,ISD)
            IZ  = MIDRST(3,ISD)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 200
            DO 220 II=1,2
                IFC = MFCEG(II,ISD)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    CALL EDGINDV (LV1,LV2,LVSKIP,ISD)
                    CALL EDGINDF (LF1,LF2,LFSKIP,ISD,II)
                    DO 240 I=2,NXM1
                        LV = LV1 + (I-1)*LVSKIP
                        LF = LF1 + (I-1)*LFSKIP
                        VNX(LV,1,1,IEL)=VNX(LV,1,1,IEL)+UNX(LF,1,IFC,IEL)
                        VNY(LV,1,1,IEL)=VNY(LV,1,1,IEL)+UNY(LF,1,IFC,IEL)
                        VNZ(LV,1,1,IEL)=VNZ(LV,1,1,IEL)+UNZ(LF,1,IFC,IEL)
                    240 END DO
                ENDIF
            220 END DO
    200 END DO

!     (2) All Corners

    DO 300 IEL=1,NEL
        DO 300 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            IZ  = MCRRST(3,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 3.1) GOTO 300
            DO 320 II=1,3
                IFC = MFCCR(II,ICR)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    IRA = NTCRF(1,II)
                    ISA = NTCRF(2,II)
                    IR  = MCRRST(IRA,ICR)
                    IS  = MCRRST(ISA,ICR)
                    VNX(IX,IY,IZ,IEL)=VNX(IX,IY,IZ,IEL)+UNX(IR,IS,IFC,IEL)
                    VNY(IX,IY,IZ,IEL)=VNY(IX,IY,IZ,IEL)+UNY(IR,IS,IFC,IEL)
                    VNZ(IX,IY,IZ,IEL)=VNZ(IX,IY,IZ,IEL)+UNZ(IR,IS,IFC,IEL)
                ENDIF
            320 END DO
    300 END DO

    CALL DSSUM   (VNX,NX1,NY1,NZ1)
    CALL DSSUM   (VNY,NX1,NY1,NZ1)
    CALL DSSUM   (VNZ,NX1,NY1,NZ1)
    CALL UNITVEC (VNX,VNY,VNZ,NTOT1)
    CALL VDOT3   (VNMAG,VNX,VNY,VNZ,VNX,VNY,VNZ,NTOT1)

    DO 500 IEL=1,NEL
        DO 500 IFC=1,NFACE
            HMF = ABS( HFMASK(1,1,IFC,IEL) )
            IF (HMF < 1.9  .OR. HMF > 3.1) GOTO 500
            CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
            DO 520 J2=JS2,JF2,JSKIP2
                DO 520 J1=JS1,JF1,JSKIP1
                    IF (VNMAG(J1,J2,1,IEL) < EPSA) GOTO 520
                    VNZDIF = ABS(VNZ(J1,J2,1,IEL)) - 1.0
                    IF (ABS(VNZDIF) < EPSN) THEN
                        V1X(J1,J2,1,IEL) = 1.0
                        V1Y(J1,J2,1,IEL) = 0.0
                        V1Z(J1,J2,1,IEL) = 0.0
                    ELSE
                        SSN = SQRT(VNX(J1,J2,1,IEL)**2 + VNY(J1,J2,1,IEL)**2)
                        V1X(J1,J2,1,IEL) = -VNY(J1,J2,1,IEL) / SSN
                        V1Y(J1,J2,1,IEL) =  VNX(J1,J2,1,IEL) / SSN
                        V1Z(J1,J2,1,IEL) =  0.0
                    ENDIF
            520 END DO
    500 END DO

    CALL DSSUM   (V1X,NX1,NY1,NZ1)
    CALL DSSUM   (V1Y,NX1,NY1,NZ1)
    CALL DSSUM   (V1Z,NX1,NY1,NZ1)
    CALL UNITVEC (V1X,V1Y,V1Z,NTOT1)

    CALL VCROSS (V2X,V2Y,V2Z,VNX,VNY,VNZ,V1X,V1Y,V1Z,NTOT1)

    RETURN
    END SUBROUTINE COMAVN3
    SUBROUTINE FIXWMSK (W2MASK,W3MASK,HVMASK,HFMASK,NEL)

    use size_m
    DIMENSION HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    DIMENSION W2MASK(LX1,LY1,LZ1,1) &
    , W3MASK(LX1,LY1,LZ1,1)

    IF (NDIM == 2) THEN
        CALL FXWMS2 (W2MASK,HVMASK,HFMASK,NEL)
    ELSE
        CALL FXWMS3 (W2MASK,W3MASK,HVMASK,HFMASK,NEL)
    ENDIF

    CALL DSOP(W2MASK,'MUL',NX1,NY1,NZ1)
    IF (NDIM == 3) CALL DSOP(W3MASK,'MUL',NX1,NY1,NZ1)

    RETURN
    END SUBROUTINE FIXWMSK
    SUBROUTINE FXWMS2 (W2MASK,HVMASK,HFMASK,NEL)

    use size_m
    use geom
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    DIMENSION W2MASK(LX1,LY1,LZ1,1) &
    , HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NCRNR  = 2**NDIM
    EPSA   = 1.0E-06

    IZ = 1
    DO 100 IEL=1,NEL
        DO 100 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR.  HMV > 2.1) GOTO 100
            DO 120 II=1,2
                IFC = MFCCR(II,ICR)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    IR  = MCRRST(II,ICR)
                    DOT = VNX(IX,IY,IZ,IEL)*UNX(IR,IZ,IFC,IEL) + &
                    VNY(IX,IY,IZ,IEL)*UNY(IR,IZ,IFC,IEL)
                    IF (DOT < 0.99) THEN
                        W2MASK(IX,IY,IZ,IEL) = 0.0
                        GOTO 100
                    ENDIF
                ENDIF
            120 END DO
    100 END DO

    RETURN
    END SUBROUTINE FXWMS2
    SUBROUTINE FXWMS3 (W2MASK,W3MASK,HVMASK,HFMASK,NEL)

    use size_m
    use geom
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    DIMENSION W2MASK(LX1,LY1,LZ1,1) &
    , W3MASK(LX1,LY1,LZ1,1) &
    , HVMASK(LX1,LY1,LZ1,1) &
    , HFMASK(LX1,LZ1,6,1)

    NCRNR = 2**NDIM
    NEDGE = 12
    NMID  = (NX1 + 1)/2
    EPSA  = 1.0E-06

    DO 100 IEL=1,NEL
        DO 100 ISD=1,12
            IX  = MIDRST(1,ISD)
            IY  = MIDRST(2,ISD)
            IZ  = MIDRST(3,ISD)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR.  HMV > 2.1) GOTO 100
            DO 120 II=1,2
                IFC = MFCEG(II,ISD)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    CALL EDGINDF (LF1,LF2,LFSKIP,ISD,II)
                    LF  = LF1 + (NMID-1)*LFSKIP
                    DOT = VNX(IX,IY,IZ,IEL)*UNX(LF,1,IFC,IEL) + &
                    VNY(IX,IY,IZ,IEL)*UNY(LF,1,IFC,IEL) + &
                    VNZ(IX,IY,IZ,IEL)*UNZ(LF,1,IFC,IEL)
                    IF (DOT < 0.99) THEN
                        CALL EDGINDV (LV1,LV2,LVSKIP,ISD)
                        DO 140 LV=LV1,LV2,LVSKIP
                            W3MASK(LV,1,1,IEL) = 0.0
                        140 END DO
                    ENDIF
                ENDIF
            120 END DO
    100 END DO

!     All Corners

    DO 300 IEL=1,NEL
        DO 300 ICR=1,NCRNR
            IX  = MCRRST(1,ICR)
            IY  = MCRRST(2,ICR)
            IZ  = MCRRST(3,ICR)
            HMV = ABS( HVMASK(IX,IY,IZ,IEL) )
            IF (HMV < 1.9  .OR. HMV > 2.1) GOTO 300
            DO 320 II=1,3
                IFC = MFCCR(II,ICR)
                HMF = ABS( HFMASK(1,1,IFC,IEL) )
                IF (ABS(HMV - HMF) < EPSA) THEN
                    IRA = NTCRF(1,II)
                    ISA = NTCRF(2,II)
                    IR  = MCRRST(IRA,ICR)
                    IS  = MCRRST(ISA,ICR)
                    DOT = VNX(IX,IY,IZ,IEL)*UNX(IR,IS,IFC,IEL) + &
                    VNY(IX,IY,IZ,IEL)*UNY(IR,IS,IFC,IEL) + &
                    VNZ(IX,IY,IZ,IEL)*UNZ(IR,IS,IFC,IEL)
                ENDIF
                IF (DOT < 0.99) THEN
                    W2MASK(IX,IY,IZ,IEL) = 0.0
                    GOTO 300
                ENDIF
            320 END DO
    300 END DO

    RETURN
    END SUBROUTINE FXWMS3
    SUBROUTINE SETCDAT

    use size_m
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    NMID = (NX1 +1)/2

!     Corners on faces

    MCRFC(1,1) = 1
    MCRFC(2,1) = 2
    MCRFC(3,1) = 6
    MCRFC(4,1) = 5
    MCRFC(1,2) = 2
    MCRFC(2,2) = 3
    MCRFC(3,2) = 7
    MCRFC(4,2) = 6
    MCRFC(1,3) = 3
    MCRFC(2,3) = 4
    MCRFC(3,3) = 8
    MCRFC(4,3) = 7
    MCRFC(1,4) = 4
    MCRFC(2,4) = 1
    MCRFC(3,4) = 5
    MCRFC(4,4) = 8
    MCRFC(1,5) = 4
    MCRFC(2,5) = 3
    MCRFC(3,5) = 2
    MCRFC(4,5) = 1
    MCRFC(1,6) = 5
    MCRFC(2,6) = 6
    MCRFC(3,6) = 7
    MCRFC(4,6) = 8

!     Faces at corners

    MFCCR(1,1) = 4
    MFCCR(2,1) = 1
    MFCCR(3,1) = 5
    MFCCR(1,2) = 1
    MFCCR(2,2) = 2
    MFCCR(3,2) = 5
    MFCCR(1,3) = 2
    MFCCR(2,3) = 3
    MFCCR(3,3) = 5
    MFCCR(1,4) = 3
    MFCCR(2,4) = 4
    MFCCR(3,4) = 5
    MFCCR(1,5) = 4
    MFCCR(2,5) = 1
    MFCCR(3,5) = 6
    MFCCR(1,6) = 1
    MFCCR(2,6) = 2
    MFCCR(3,6) = 6
    MFCCR(1,7) = 2
    MFCCR(2,7) = 3
    MFCCR(3,7) = 6
    MFCCR(1,8) = 3
    MFCCR(2,8) = 4
    MFCCR(3,8) = 6

!     Edges at corners

    MEGCR(1,1) = 4
    MEGCR(2,1) = 1
    MEGCR(3,1) = 9
    MEGCR(1,2) = 1
    MEGCR(2,2) = 2
    MEGCR(3,2) = 10
    MEGCR(1,3) = 2
    MEGCR(2,3) = 3
    MEGCR(3,3) = 11
    MEGCR(1,4) = 3
    MEGCR(2,4) = 4
    MEGCR(3,4) = 12
    MEGCR(1,5) = 8
    MEGCR(2,5) = 5
    MEGCR(3,5) = 9
    MEGCR(1,6) = 5
    MEGCR(2,6) = 6
    MEGCR(3,6) = 10
    MEGCR(1,7) = 6
    MEGCR(2,7) = 7
    MEGCR(3,7) = 11
    MEGCR(1,8) = 7
    MEGCR(2,8) = 8
    MEGCR(3,8) = 12

!     Faces on edges

    MFCEG(1,1)  = 1
    MFCEG(2,1)  = 5
    MFCEG(1,2)  = 2
    MFCEG(2,2)  = 5
    MFCEG(1,3)  = 3
    MFCEG(2,3)  = 5
    MFCEG(1,4)  = 4
    MFCEG(2,4)  = 5
    MFCEG(1,5)  = 1
    MFCEG(2,5)  = 6
    MFCEG(1,6)  = 2
    MFCEG(2,6)  = 6
    MFCEG(1,7)  = 3
    MFCEG(2,7)  = 6
    MFCEG(1,8)  = 4
    MFCEG(2,8)  = 6
    MFCEG(1,9)  = 4
    MFCEG(2,9)  = 1
    MFCEG(1,10) = 1
    MFCEG(2,10) = 2
    MFCEG(1,11) = 2
    MFCEG(2,11) = 3
    MFCEG(1,12) = 3
    MFCEG(2,12) = 4

!     Corners at edges

    MCREG(1,1)  = 1
    MCREG(2,1)  = 2
    MCREG(1,2)  = 2
    MCREG(2,2)  = 3
    MCREG(1,3)  = 4
    MCREG(2,3)  = 3
    MCREG(1,4)  = 1
    MCREG(2,4)  = 4
    MCREG(1,5)  = 5
    MCREG(2,5)  = 6
    MCREG(1,6)  = 6
    MCREG(2,6)  = 7
    MCREG(1,7)  = 8
    MCREG(2,7)  = 7
    MCREG(1,8)  = 5
    MCREG(2,8)  = 8
    MCREG(1,9)  = 1
    MCREG(2,9)  = 5
    MCREG(1,10) = 2
    MCREG(2,10) = 6
    MCREG(1,11) = 3
    MCREG(2,11) = 7
    MCREG(1,12) = 4
    MCREG(2,12) = 8

!     Corner indices (Vol array)

    MCRRST(1,1) = 1
    MCRRST(2,1) = 1
    MCRRST(3,1) = 1
    MCRRST(1,2) = NX1
    MCRRST(2,2) = 1
    MCRRST(3,2) = 1
    MCRRST(1,3) = NX1
    MCRRST(2,3) = NX1
    MCRRST(3,3) = 1
    MCRRST(1,4) = 1
    MCRRST(2,4) = NX1
    MCRRST(3,4) = 1
    MCRRST(1,5) = 1
    MCRRST(2,5) = 1
    MCRRST(3,5) = NX1
    MCRRST(1,6) = NX1
    MCRRST(2,6) = 1
    MCRRST(3,6) = NX1
    MCRRST(1,7) = NX1
    MCRRST(2,7) = NX1
    MCRRST(3,7) = NX1
    MCRRST(1,8) = 1
    MCRRST(2,8) = NX1
    MCRRST(3,8) = NX1

!     Mid-edge indcies (Vol array)

    MIDRST(1,1)  = NMID
    MIDRST(1,2)  = NX1
    MIDRST(1,3)  = NMID
    MIDRST(1,4)  = 1
    MIDRST(1,5)  = NMID
    MIDRST(1,6)  = NX1
    MIDRST(1,7)  = NMID
    MIDRST(1,8)  = 1
    MIDRST(1,9)  = 1
    MIDRST(1,10) = NX1
    MIDRST(1,11) = NX1
    MIDRST(1,12) = 1
    MIDRST(2,1)  = 1
    MIDRST(2,2)  = NMID
    MIDRST(2,3)  = NX1
    MIDRST(2,4)  = NMID
    MIDRST(2,5)  = 1
    MIDRST(2,6)  = NMID
    MIDRST(2,7)  = NX1
    MIDRST(2,8)  = NMID
    MIDRST(2,9)  = 1
    MIDRST(2,10) = 1
    MIDRST(2,11) = NX1
    MIDRST(2,12) = NX1
    MIDRST(3,1)  = 1
    MIDRST(3,2)  = 1
    MIDRST(3,3)  = 1
    MIDRST(3,4)  = 1
    MIDRST(3,5)  = NX1
    MIDRST(3,6)  = NX1
    MIDRST(3,7)  = NX1
    MIDRST(3,8)  = NX1
    MIDRST(3,9)  = NMID
    MIDRST(3,10) = NMID
    MIDRST(3,11) = NMID
    MIDRST(3,12) = NMID

!     1-D corners indices (Vol array)

    MCRIND(1) = 1
    MCRIND(2) = NX1
    MCRIND(3) = NX1**2
    MCRIND(7) = NX1**3
    MCRIND(4) = MCRIND(3) - NX1 + 1
    MCRIND(5) = MCRIND(7) - MCRIND(3) + 1
    MCRIND(6) = MCRIND(5) + NX1 - 1
    MCRIND(8) = MCRIND(7) - NX1 + 1

!     1-D  edge indices (Face array)

    MEDIND(1,1) = 1
    MEDIND(2,1) = NX1
    MEDIND(1,2) = NX1**2 - NX1 + 1
    MEDIND(2,2) = NX1**2
    MEDIND(1,3) = 1
    MEDIND(2,3) = MEDIND(1,2)
    MEDIND(1,4) = NX1
    MEDIND(2,4) = NX1**2

!     1-D edge index type (Face array)

    NTEFC(1,1)  = 1
    NTEFC(2,1)  = 1
    NTEFC(1,2)  = 1
    NTEFC(2,2)  = 4
    NTEFC(1,3)  = 1
    NTEFC(2,3)  = 2
    NTEFC(1,4)  = 1
    NTEFC(2,4)  = 3
    NTEFC(1,5)  = 2
    NTEFC(2,5)  = 1
    NTEFC(1,6)  = 2
    NTEFC(2,6)  = 4
    NTEFC(1,7)  = 2
    NTEFC(2,7)  = 2
    NTEFC(1,8)  = 2
    NTEFC(2,8)  = 3
    NTEFC(1,9)  = 3
    NTEFC(2,9)  = 3
    NTEFC(1,10) = 4
    NTEFC(2,10) = 3
    NTEFC(1,11) = 4
    NTEFC(2,11) = 4
    NTEFC(1,12) = 3
    NTEFC(2,12) = 4

!     Corner index address on face in MCRRST

    NTCRF(1,1) = 1
    NTCRF(2,1) = 3
    NTCRF(1,2) = 2
    NTCRF(2,2) = 3
    NTCRF(1,3) = 1
    NTCRF(2,3) = 2

    RETURN
    END SUBROUTINE SETCDAT
    SUBROUTINE EDGINDF (LF1,LF2,LFSKIP,ISD,IFCN)

    use size_m
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    ITYP = NTEFC(IFCN,ISD)

    LF1 = MEDIND(1,ITYP)
    LF2 = MEDIND(2,ITYP)

    LFSKIP = 1
    IF (ITYP >= 3) LFSKIP = NX1

    RETURN
    END SUBROUTINE EDGINDF
    SUBROUTINE EDGINDV (LV1,LV2,LVSKIP,ISD)

    use size_m
    common /indxfc/ mcrfc(4,6) &
    , mfccr(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)

    IODD = ISD - ISD/2*2
    ICR1 = MCREG(1,ISD)
    ICR2 = MCREG(2,ISD)

    LV1  = MCRIND(ICR1)
    LV2  = MCRIND(ICR2)

    IF (ISD >= 9) THEN
        LVSKIP = NX1**2
    ELSE
        IF (IODD == 0) THEN
            LVSKIP = NX1
        ELSE
            LVSKIP = 1
        ENDIF
    ENDIF

    RETURN
    END SUBROUTINE EDGINDV
    SUBROUTINE SETCDOF

    use size_m
    use input

    NFACE = 2*NDIM

    DO 100 IEL=1,NELT
        DO 100 IFC=1,NFACE
            CDOF(IFC,IEL)=CBC(IFC,IEL,0)(1:1)
    100 END DO

    RETURN
    END SUBROUTINE SETCDOF
    SUBROUTINE AMASK (VB1,VB2,VB3,V1,V2,V3,NEL)

    use size_m
    use geom
    use input
    use soln
    use tstep
    common /scrsf/ a1mask(lx1,ly1,lz1,lelt) &
    , a2mask(lx1,ly1,lz1,lelt) &
    , a3mask(lx1,ly1,lz1,lelt)
    common /ctmp0/ wa(lx1,ly1,lz1,lelt)

    DIMENSION VB1(LX1,LY1,LZ1,1) &
    , VB2(LX1,LY1,LZ1,1) &
    , VB3(LX1,LY1,LZ1,1) &
    , V1(LX1,LY1,LZ1,1) &
    , V2(LX1,LY1,LZ1,1) &
    , V3(LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NEL
    CALL RONE (WA,NTOT1)
    CALL COPY (VB1,V1,NTOT1)
    CALL COPY (VB2,V2,NTOT1)
    IF (NDIM == 3) CALL COPY (VB3,V3,NTOT1)

    IF (IFIELD == 1) THEN
        CALL SUB3  (A1MASK,WA,V1MASK,NTOT1)
        CALL SUB3  (A2MASK,WA,V2MASK,NTOT1)
        IF (NDIM == 3) CALL SUB3 (A3MASK,WA,V3MASK,NTOT1)
    ELSEIF (IFIELD == ifldmhd) THEN
        CALL SUB3  (A1MASK,WA,B1MASK,NTOT1)
        CALL SUB3  (A2MASK,WA,B2MASK,NTOT1)
        IF (NDIM == 3) CALL SUB3 (A3MASK,WA,B3MASK,NTOT1)
    ELSE
        CALL SUB3  (A1MASK,WA,W1MASK,NTOT1)
        CALL SUB3  (A2MASK,WA,W2MASK,NTOT1)
        IF (NDIM == 3) CALL SUB3 (A3MASK,WA,W3MASK,NTOT1)
    ENDIF

    CALL QMASK (VB1,VB2,VB3,A1MASK,A2MASK,A3MASK,NEL)

    RETURN
    END SUBROUTINE AMASK
    SUBROUTINE RMASK (R1,R2,R3,NEL)

    use size_m
    use input
    use mvgeom
    use soln
    use tstep

    DIMENSION R1  (LX1,LY1,LZ1,1) &
    , R2  (LX1,LY1,LZ1,1) &
    , R3  (LX1,LY1,LZ1,1)

    if (ifsplit) then
        call opmask(r1,r2,r3)
        return
    endif

!     call outfldro (v1mask,'v1mask rmk',0)
!     call outfldro (v2mask,'v2mask rmk',1)

    IF (IFIELD == 1) THEN
        CALL QMASK (R1,R2,R3,V1MASK,V2MASK,V3MASK,NEL)
    ELSEIF (ifield == ifldmhd) then
        CALL QMASK (R1,R2,R3,B1MASK,B2MASK,B3MASK,NEL)
    ELSE
        CALL QMASK (R1,R2,R3,W1MASK,W2MASK,W3MASK,NEL)
    ENDIF

!     call outfldro (r1,'r1   rmask',0)
!     call outfldro (r2,'r2   rmask',1)
!     call exitti   ('quit in rmask$,',nel)

    RETURN
    END SUBROUTINE RMASK
    SUBROUTINE QMASK (R1,R2,R3,R1MASK,R2MASK,R3MASK,NEL)

    use size_m
    use geom
    use tstep
    common /ctmp1/ s1(lx1,ly1,lz1,lelt) &
    , s2(lx1,ly1,lz1,lelt) &
    , s3(lx1,ly1,lz1,lelt)

    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , R3(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1) &
    , R3MASK(LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NEL

!     (0) Collocate Volume Mask

    CALL COPY  (S1,R1,NTOT1)
    CALL COPY  (S2,R2,NTOT1)
    CALL COL2  (R1,R1MASK,NTOT1)
    CALL COL2  (R2,R2MASK,NTOT1)
    IF (NDIM == 3) THEN
        CALL COPY (S3,R3,NTOT1)
        CALL COL2 (R3,R3MASK,NTOT1)
    ENDIF

!     (1) Face Mask

!     write(6,*) 'iflmsf:',iflmsf(ifield),ifield
    IF (IFLMSF(IFIELD)) THEN
        IF (NDIM == 2) THEN
            CALL FCMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
        ELSE
            CALL FCMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
        ENDIF
    ENDIF

!     (2) Edge Mask  (3-D only)

    IF (NDIM == 3 .AND. IFLMSE(IFIELD)) &
    CALL EGMASK (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)

!     (3) Corner Mask

    IF (IFLMSC(IFIELD)) THEN
        IF (NDIM == 2) THEN
            CALL CRMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
        ELSE
            CALL CRMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
        ENDIF
    ENDIF

    RETURN
    END SUBROUTINE QMASK
    SUBROUTINE FCMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)

    use size_m
    use geom
    use tstep
    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , S1(LX1,LY1,LZ1,1) &
    , S2(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1)

    NFACE = 2*NDIM

    DO 100 IEL=1,NEL
        DO 100 IFC=1,NFACE
            IF ( .NOT. IFMSFC(IFC,IEL,IFIELD)) GO TO 100
            CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
            DO 120 J2=JS2,JF2,JSKIP2
                DO 120 J1=JS1,JF1,JSKIP1
                    RNOR = ( S1(J1,J2,1,IEL)*VNX(J1,J2,1,IEL) + &
                    S2(J1,J2,1,IEL)*VNY(J1,J2,1,IEL) ) * &
                    R1MASK(J1,J2,1,IEL)
                    RTN1 = ( S1(J1,J2,1,IEL)*V1X(J1,J2,1,IEL) + &
                    S2(J1,J2,1,IEL)*V1Y(J1,J2,1,IEL) ) * &
                    R2MASK(J1,J2,1,IEL)
                    R1(J1,J2,1,IEL) = RNOR*VNX(J1,J2,1,IEL) + &
                    RTN1*V1X(J1,J2,1,IEL)
                    R2(J1,J2,1,IEL) = RNOR*VNY(J1,J2,1,IEL) + &
                    RTN1*V1Y(J1,J2,1,IEL)
            120 END DO
    100 END DO

    RETURN
    END SUBROUTINE FCMSK2
    SUBROUTINE FCMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)

    use size_m
    use geom
    use tstep
    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , R3(LX1,LY1,LZ1,1) &
    , S1(LX1,LY1,LZ1,1) &
    , S2(LX1,LY1,LZ1,1) &
    , S3(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1) &
    , R3MASK(LX1,LY1,LZ1,1)

    NFACE = 2*NDIM

    DO 100 IEL=1,NEL
        DO 100 IFC=1,NFACE
            IF ( .NOT. IFMSFC(IFC,IEL,IFIELD)) GO TO 100
            CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
            DO 120 J2=JS2,JF2,JSKIP2
                DO 120 J1=JS1,JF1,JSKIP1
                    RNOR = ( S1(J1,J2,1,IEL)*VNX(J1,J2,1,IEL) + &
                    S2(J1,J2,1,IEL)*VNY(J1,J2,1,IEL) + &
                    S3(J1,J2,1,IEL)*VNZ(J1,J2,1,IEL) ) * &
                    R1MASK(J1,J2,1,IEL)
                    RTN1 = ( S1(J1,J2,1,IEL)*V1X(J1,J2,1,IEL) + &
                    S2(J1,J2,1,IEL)*V1Y(J1,J2,1,IEL) + &
                    S3(J1,J2,1,IEL)*V1Z(J1,J2,1,IEL) ) * &
                    R2MASK(J1,J2,1,IEL)
                    RTN2 = ( S1(J1,J2,1,IEL)*V2X(J1,J2,1,IEL) + &
                    S2(J1,J2,1,IEL)*V2Y(J1,J2,1,IEL) + &
                    S3(J1,J2,1,IEL)*V2Z(J1,J2,1,IEL) ) * &
                    R3MASK(J1,J2,1,IEL)
                    R1(J1,J2,1,IEL) = RNOR*VNX(J1,J2,1,IEL) + &
                    RTN1*V1X(J1,J2,1,IEL) + &
                    RTN2*V2X(J1,J2,1,IEL)
                    R2(J1,J2,1,IEL) = RNOR*VNY(J1,J2,1,IEL) + &
                    RTN1*V1Y(J1,J2,1,IEL) + &
                    RTN2*V2Y(J1,J2,1,IEL)
                    R3(J1,J2,1,IEL) = RNOR*VNZ(J1,J2,1,IEL) + &
                    RTN1*V1Z(J1,J2,1,IEL) + &
                    RTN2*V2Z(J1,J2,1,IEL)
            120 END DO
    100 END DO

    RETURN
    END SUBROUTINE FCMSK3
    SUBROUTINE EGMASK (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)

    use size_m
    use geom
    use tstep
    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , R3(LX1,LY1,LZ1,1) &
    , S1(LX1,LY1,LZ1,1) &
    , S2(LX1,LY1,LZ1,1) &
    , S3(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1) &
    , R3MASK(LX1,LY1,LZ1,1)

    NEDGE = 12

    DO 100 IEL=1,NEL
        DO 100 ISD=1,NEDGE
            IF ( .NOT. IFMSEG(ISD,IEL,IFIELD)) GOTO 100
            CALL EDGINDV (LV1,LV2,LVSKIP,ISD)
            DO 120 LV=LV1,LV2,LVSKIP
                RNOR = ( S1(LV,1,1,IEL)*VNX(LV,1,1,IEL) + &
                S2(LV,1,1,IEL)*VNY(LV,1,1,IEL) + &
                S3(LV,1,1,IEL)*VNZ(LV,1,1,IEL) ) * &
                R1MASK(LV,1,1,IEL)
                RTN1 = ( S1(LV,1,1,IEL)*V1X(LV,1,1,IEL) + &
                S2(LV,1,1,IEL)*V1Y(LV,1,1,IEL) + &
                S3(LV,1,1,IEL)*V1Z(LV,1,1,IEL) ) * &
                R2MASK(LV,1,1,IEL)
                RTN2 = ( S1(LV,1,1,IEL)*V2X(LV,1,1,IEL) + &
                S2(LV,1,1,IEL)*V2Y(LV,1,1,IEL) + &
                S3(LV,1,1,IEL)*V2Z(LV,1,1,IEL) ) * &
                R3MASK(LV,1,1,IEL)
                R1(LV,1,1,IEL) = RNOR*VNX(LV,1,1,IEL) + &
                RTN1*V1X(LV,1,1,IEL) + &
                RTN2*V2X(LV,1,1,IEL)
                R2(LV,1,1,IEL) = RNOR*VNY(LV,1,1,IEL) + &
                RTN1*V1Y(LV,1,1,IEL) + &
                RTN2*V2Y(LV,1,1,IEL)
                R3(LV,1,1,IEL) = RNOR*VNZ(LV,1,1,IEL) + &
                RTN1*V1Z(LV,1,1,IEL) + &
                RTN2*V2Z(LV,1,1,IEL)
            120 END DO
    100 END DO

    RETURN
    END SUBROUTINE EGMASK
    SUBROUTINE CRMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)

    use size_m
    use geom
    use tstep
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)
    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , S1(LX1,LY1,LZ1,1) &
    , S2(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1)

    NCRNR = 2**NDIM

    DO 100 IEL=1,NEL
        DO 100 ICR=1,NCRNR
            IF ( .NOT. IFMSCR(ICR,IEL,IFIELD)) GO TO 100
            IX = MCRRST(1,ICR)
            IY = MCRRST(2,ICR)
            IZ = 1
            RNOR = ( S1(IX,IY,IZ,IEL)*VNX(IX,IY,IZ,IEL) + &
            S2(IX,IY,IZ,IEL)*VNY(IX,IY,IZ,IEL) ) * &
            R1MASK(IX,IY,IZ,IEL)
            RTN1 = ( S1(IX,IY,IZ,IEL)*V1X(IX,IY,IZ,IEL) + &
            S2(IX,IY,IZ,IEL)*V1Y(IX,IY,IZ,IEL) ) * &
            R2MASK(IX,IY,IZ,IEL)
            R1(IX,IY,IZ,IEL) = RNOR*VNX(IX,IY,IZ,IEL) + &
            RTN1*V1X(IX,IY,IZ,IEL)
            R2(IX,IY,IZ,IEL) = RNOR*VNY(IX,IY,IZ,IEL) + &
            RTN1*V1Y(IX,IY,IZ,IEL)
    100 END DO

    RETURN
    END SUBROUTINE CRMSK2
    SUBROUTINE CRMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)

    use size_m
    use geom
    use tstep
    common /indxfc/ mcrfc(4,6) &
    , MFCCR(3,8) &
    , MEGCR(3,8) &
    , MFCEG(2,12) &
    , MCREG(2,12) &
    , MCRRST(3,8) &
    , MIDRST(3,12) &
    , MCRIND(8) &
    , MEDIND(2,4) &
    , NTEFC(2,12) &
    , NTCRF(2,3)
    DIMENSION R1(LX1,LY1,LZ1,1) &
    , R2(LX1,LY1,LZ1,1) &
    , R3(LX1,LY1,LZ1,1) &
    , S1(LX1,LY1,LZ1,1) &
    , S2(LX1,LY1,LZ1,1) &
    , S3(LX1,LY1,LZ1,1) &
    , R1MASK(LX1,LY1,LZ1,1) &
    , R2MASK(LX1,LY1,LZ1,1) &
    , R3MASK(LX1,LY1,LZ1,1)

    NCRNR = 2**NDIM

    DO 100 IEL=1,NEL
        DO 100 ICR=1,NCRNR
            IF ( .NOT. IFMSCR(ICR,IEL,IFIELD)) GO TO 100
            IX = MCRRST(1,ICR)
            IY = MCRRST(2,ICR)
            IZ = MCRRST(3,ICR)
            RNOR = ( S1(IX,IY,IZ,IEL)*VNX(IX,IY,IZ,IEL) + &
            S2(IX,IY,IZ,IEL)*VNY(IX,IY,IZ,IEL) + &
            S3(IX,IY,IZ,IEL)*VNZ(IX,IY,IZ,IEL) ) * &
            R1MASK(IX,IY,IZ,IEL)
            RTN1 = ( S1(IX,IY,IZ,IEL)*V1X(IX,IY,IZ,IEL) + &
            S2(IX,IY,IZ,IEL)*V1Y(IX,IY,IZ,IEL) + &
            S3(IX,IY,IZ,IEL)*V1Z(IX,IY,IZ,IEL) ) * &
            R2MASK(IX,IY,IZ,IEL)
            RTN2 = ( S1(IX,IY,IZ,IEL)*V2X(IX,IY,IZ,IEL) + &
            S2(IX,IY,IZ,IEL)*V2Y(IX,IY,IZ,IEL) + &
            S3(IX,IY,IZ,IEL)*V2Z(IX,IY,IZ,IEL) ) * &
            R3MASK(IX,IY,IZ,IEL)
            R1(IX,IY,IZ,IEL) = RNOR*VNX(IX,IY,IZ,IEL) + &
            RTN1*V1X(IX,IY,IZ,IEL) + &
            RTN2*V2X(IX,IY,IZ,IEL)
            R2(IX,IY,IZ,IEL) = RNOR*VNY(IX,IY,IZ,IEL) + &
            RTN1*V1Y(IX,IY,IZ,IEL) + &
            RTN2*V2Y(IX,IY,IZ,IEL)
            R3(IX,IY,IZ,IEL) = RNOR*VNZ(IX,IY,IZ,IEL) + &
            RTN1*V1Z(IX,IY,IZ,IEL) + &
            RTN2*V2Z(IX,IY,IZ,IEL)
    100 END DO

    RETURN
    END SUBROUTINE CRMSK3

    subroutine getSnormal(sn,ix,iy,iz,iside,e)

!     calculate surface normal

    use size_m
    use geom
    use topol

    real :: sn(3)
    integer :: e,f

    f = eface1(iside)

    if (1 <= f .AND. f <= 2) then     ! "r face"
        sn(1) = unx(iy,iz,iside,e)
        sn(2) = uny(iy,iz,iside,e)
        sn(3) = unz(iy,iz,iside,e)
    elseif (3 <= f .AND. f <= 4) then ! "s face"
        sn(1) = unx(ix,iz,iside,e)
        sn(2) = uny(ix,iz,iside,e)
        sn(3) = unz(ix,iz,iside,e)
    elseif (5 <= f .AND. f <= 6) then ! "t face"
        sn(1) = unx(ix,iy,iside,e)
        sn(2) = uny(ix,iy,iside,e)
        sn(3) = unz(ix,iy,iside,e)
    endif

    return
    end subroutine getSnormal

    subroutine fixmska (c1mask,c2mask,c3mask)

!     fixes masks for A/SYM face corners

    use size_m
    use input
     
    real ::   c1mask(lx1,ly1,lz1,1) &
    ,c2mask(lx1,ly1,lz1,1) &
    ,c3mask(lx1,ly1,lz1,1)

    common /ctmp0/ im1(lx1,ly1,lz1),im2(lx1,ly1,lz1)
    integer :: e,f,val,im1,im2

    character(3) :: cb

    n = nx1*ny1*nz1

    do e=1,nelv
        call izero (im1,n)
        call izero (im2,n)
        do f=1,nface
            cb  = cbc (f,e,1)
            if (cb == 'SYM')  call iface_e(im1,f,1,nx1,ny1,nz1)
            if (cb == 'A  ')  call iface_e(im2,f,2,nx1,ny1,nz1)
        enddo
        call icol2(im2,im1,n)

        k = 1
        do j=1,ny1,ny1-1
            do i=1,nx1,nx1-1
                if  ( im2(i,j,k) == 2) then  ! corner of SYM & 'A  ' faces
                    c1mask(i,j,k,e) = 0.
                    c2mask(i,j,k,e) = 0.
                endif
            enddo
        enddo
    enddo

    return
    end subroutine fixmska
!-----------------------------------------------------------------------
    subroutine icol2(a,b,n)
    integer :: a(1),b(1)

    do i=1,n
        a(i)=a(i)*b(i)
    enddo

    return
    end subroutine icol2
!-----------------------------------------------------------------------
    subroutine iface_e(a,iface,val,nx,ny,nz)

!     Assign the value VAL to face(IFACE,IE) of array A.
!     IFACE is the input in the pre-processor ordering scheme.

    use size_m
    integer :: a(nx,ny,nz),val
    call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
    do 100 iz=kz1,kz2
        do 100 iy=ky1,ky2
            do 100 ix=kx1,kx2
                a(ix,iy,iz)=val
    100 END DO
    return
    end subroutine iface_e
!-----------------------------------------------------------------------
    function op_vlsc2_wt(b1,b2,b3,x1,x2,x3,wt)
    use size_m
    use input
    use tstep
    real :: b1(1),b2(1),b3(1),x1(1),x2(1),x3(1),wt(1)

    nel = nelfld(ifield)
    n   = nx1*ny1*nz1*nel

    s = 0
    if (if3d) then
        do i=1,n
            s=s+wt(i)*(b1(i)*x1(i)+b2(i)*x2(i)+b3(i)*x3(i))
        enddo
    else
        do i=1,n
            s=s+wt(i)*(b1(i)*x1(i)+b2(i)*x2(i))
        enddo
    endif
    op_vlsc2_wt = s
          
    return
    end function op_vlsc2_wt
!-----------------------------------------------------------------------
    function op_glsc2_wt(b1,b2,b3,x1,x2,x3,wt)
    use size_m
    use input
    use tstep
    real :: b1(1),b2(1),b3(1),x1(1),x2(1),x3(1),wt(1)

    nel = nelfld(ifield)
    n   = nx1*ny1*nz1*nel

    s = 0
    if (if3d) then
        do i=1,n
            s=s+wt(i)*(b1(i)*x1(i)+b2(i)*x2(i)+b3(i)*x3(i))
        enddo
    else
        do i=1,n
            s=s+wt(i)*(b1(i)*x1(i)+b2(i)*x2(i))
        enddo
    endif
    op_glsc2_wt = glsum(s,1)
          
    return
    end function op_glsc2_wt
!-----------------------------------------------------------------------
      SUBROUTINE FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

      use size_m
      use topol

      CALL DSSET (NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)

      RETURN
      END
!-----------------------------------------------------------------------
