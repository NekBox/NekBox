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

    SUBROUTINE CMULT2 (A,B,CONST,N)
    DIMENSION A(1),B(1)
    DO 100 I=1,N
        A(I)=B(I)*CONST
    100 END DO
    RETURN
    END SUBROUTINE CMULT2

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

    SUBROUTINE UPDMSYS (IFLD)

    use size_m
    use geom
    use tstep
    common /screv/ hfmask(lx1,lz1,6,lelt) &
    , hvmask(lx1,ly1,lz1,lelt)

    IF ( .NOT. IFLMSF(IFLD)) RETURN
#if 0
    NEL  = NELFLD(IFLD)
    CALL SETHMSK (HVMASK,HFMASK,IFLD,NEL)
    CALL SETCSYS (HVMASK,HFMASK,NEL)
#endif

    RETURN
    END SUBROUTINE UPDMSYS

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
