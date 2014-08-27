!-----------------------------------------------------------------------
    SUBROUTINE SETLOG

!     Subroutine to initialize logical flags

    use ctimer
    use size_m
    use geom
    use input
    use tstep
    use turbo
    COMMON  /CPRINT/ IFPRINT

    common  /nekcb/ cb
    CHARACTER CB*3
    LOGICAL :: IFALGN,IFNORX,IFNORY,IFNORZ,IFPRINT

    NFACE  = 2*NDIM
    NMXV   = NFACE*NELV
    NMXT   = NFACE*NELT

    IFPRINT = .TRUE. 
    IFVCOR  = .TRUE. 
    IFGEOM  = .FALSE. 
    IFINTQ  = .FALSE. 
    IFSURT  = .FALSE. 
    IFWCNO  = .FALSE. 
    IFSWALL = .FALSE. 
    DO 10 IFIELD=1,NFIELD
        IFNONL(IFIELD) = .FALSE. 
    10 END DO

    CALL LFALSE (IFEPPM,NMXV)
    CALL LFALSE (IFQINP,NMXV)

!max    IF (IFMODEL) CALL SETSHL

    IF (IFMVBD) THEN
        IFGEOM = .TRUE. 
        IF ( IFFLOW .AND. .NOT. IFNAV )       IFWCNO          = .TRUE. 
        IF ( IFMELT .AND. .NOT. IFFLOW )      IFWCNO          = .TRUE. 
    ENDIF

    IF (IFFLOW) THEN

    ! k         call check_cyclic  ! fow now; set in .rea file

        IFIELD = 1
        DO 100 IEL=1,NELV
            DO 100 IFC=1,NFACE
                CB = CBC(IFC,IEL,IFIELD)
                CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
                IF ( .NOT. IFSTRS ) CALL CHKCBC  (CB,IEL,IFC,IFALGN)
                IF  (CB == 'O  ' .OR. CB == 'o  ' .OR. &
                CB == 'ON ' .OR. CB == 'on ' .OR. &
                CB == 'S  ' .OR. CB == 's  ' .OR. &
                CB == 'SL ' .OR. CB == 'sl ' .OR. &
                CB == 'MM ' .OR. CB == 'mm ' .OR. &
                CB == 'MS ' .OR. CB == 'ms ')  THEN
                    IFVCOR          = .FALSE. 
                    IFEPPM(IFC,IEL) = .TRUE. 
                ENDIF
                IF  (CB == 'VL ' .OR. CB == 'vl ' .OR. &
                CB == 'WSL' .OR. CB == 'wsl' .OR. &
                CB == 'SL ' .OR. CB == 'sl ' .OR. &
                CB == 'SHL' .OR. CB == 'shl' .OR. &
                CB == 'MM ' .OR. CB == 'mm ' .OR. &
                CB == 'MS ' .OR. CB == 'ms ' .OR. &
                CB == 'O  ' .OR. CB == 'o  ' .OR. &
                CB == 'ON ' .OR. CB == 'on ')  THEN
                    IFQINP(IFC,IEL) = .TRUE. 
                ENDIF
                IF  (CB == 'MS ' .OR. CB == 'ms ' .OR. &
                CB == 'MM ' .OR. CB == 'mm ' .OR. &
                CB == 'MSI' .OR. CB == 'msi' ) THEN
                    IFSURT          = .TRUE. 
                ENDIF
                IF  (CB == 'WS ' .OR. CB == 'ws ' .OR. &
                CB == 'WSL' .OR. CB == 'wsl') THEN
                    IFSWALL         = .TRUE. 
                    IFCWUZ          = .TRUE. 
                ENDIF
        100 END DO
    ENDIF

    IF (IFHEAT) THEN
    
        DO 250 IFIELD=2,NFIELD
            DO 250 IEL=1,NELFLD(IFIELD)
                DO 250 IFC=1,NFACE
                    CB=CBC(IFC,IEL,IFIELD)
                    IF  (CB == 'r  ' .OR. CB == 'R  ') THEN
                        IFNONL(IFIELD)  = .TRUE. 
                    ENDIF
        250 END DO
    
    ENDIF

!max    if (ifmhd) call set_ifbcor

    IF (NHIS > 0) THEN
        IQ = 0
        DO 300 IH=1,NHIS
            IF ( HCODE(10,IH) == 'I' ) THEN
                IFINTQ = .TRUE. 
                IOBJ   = LOCHIS(1,IH)
                IQ     = IQ + 1
                IF (IOBJ > NOBJ .OR. IOBJ < 0)  THEN
                    WRITE (6,*) &
                    'ERROR : Undefined Object for integral',IQ
                    call exitt
                ENDIF
            ENDIF
        300 END DO
    ENDIF

!     Establish global consistency of LOGICALS amongst all processors.

    CALL GLLOG(IFVCOR , .FALSE. )
    CALL GLLOG(IFSURT , .TRUE. )
    CALL GLLOG(IFSWALL, .TRUE. )
    CALL GLLOG(IFCWUZ , .TRUE. )
    CALL GLLOG(IFWCNO , .TRUE. )
    DO 400 IFIELD=2,NFIELD
        CALL GLLOG(IFNONL(IFIELD), .TRUE. )
    400 END DO

    IF (NID == 0) THEN
        WRITE (6,*) 'IFTRAN   =',IFTRAN
        WRITE (6,*) 'IFFLOW   =',IFFLOW
        WRITE (6,*) 'IFHEAT   =',IFHEAT
        WRITE (6,*) 'IFSPLIT  =',IFSPLIT
        WRITE (6,*) 'IFLOMACH =',IFLOMACH
        WRITE (6,*) 'IFUSERVP =',IFUSERVP
        WRITE (6,*) 'IFUSERMV =',IFUSERMV
        WRITE (6,*) 'IFSTRS   =',IFSTRS
        WRITE (6,*) 'IFCHAR   =',IFCHAR
        WRITE (6,*) 'IFCYCLIC =',IFCYCLIC
        WRITE (6,*) 'IFAXIS   =',IFAXIS
        WRITE (6,*) 'IFMVBD   =',IFMVBD
        WRITE (6,*) 'IFMELT   =',IFMELT
        WRITE (6,*) 'IFMODEL  =',IFMODEL
        WRITE (6,*) 'IFKEPS   =',IFKEPS
        WRITE (6,*) 'IFMOAB   =',IFMOAB
        WRITE (6,*) 'IFNEKNEK =',IFNEKNEK
        WRITE (6,*) 'IFSYNC   =',IFSYNC
        WRITE (6,*) '  '
        WRITE (6,*) 'IFVCOR   =',IFVCOR
        WRITE (6,*) 'IFINTQ   =',IFINTQ
        WRITE (6,*) 'IFCWUZ   =',IFCWUZ
        WRITE (6,*) 'IFSWALL  =',IFSWALL
        WRITE (6,*) 'IFGEOM   =',IFGEOM
        WRITE (6,*) 'IFSURT   =',IFSURT
        WRITE (6,*) 'IFWCNO   =',IFWCNO
        DO 500 IFIELD=1,NFIELD
            WRITE (6,*) '  '
            WRITE (6,*) 'IFTMSH for field',IFIELD,'   = ',IFTMSH(IFIELD)
            WRITE (6,*) 'IFADVC for field',IFIELD,'   = ',IFADVC(IFIELD)
            WRITE (6,*) 'IFNONL for field',IFIELD,'   = ',IFNONL(IFIELD)
        500 END DO
        WRITE (6,*) '  '
        if (param(99) > 0) write(6,*) 'Dealiasing enabled, lxd=', lxd
    ENDIF

    RETURN
    END SUBROUTINE SETLOG

!-----------------------------------------------------------------------
    SUBROUTINE SETRZER
!-------------------------------------------------------------------

!     Check for axisymmetric case.
!     Are some of the elements close to the axis?

!-------------------------------------------------------------------
    use size_m
    use geom
    use input

!     Single or double precision???

    DELTA = 1.E-9
    X     = 1.+DELTA
    Y     = 1.
    DIFF  = ABS(X-Y)
    IF (DIFF == 0.) EPS = 1.E-7
    IF (DIFF > 0.) EPS = 1.E-14
    eps1 = 1.e-6 ! for prenek mesh in real*4

    DO 100 IEL=1,NELT
        IFRZER(IEL) = .FALSE. 
        IF (IFAXIS) THEN
            NVERT = 0
            DO 10 IC=1,4
                IF(ABS(YC(IC,IEL)) < EPS1) THEN
                    NVERT = NVERT+1
                    YC(IC,IEL) = 0.0  ! exactly on the axis
                ENDIF
            10 END DO
        ENDIF
        IEDGE = 1
        IF ((NVERT == 2) .AND. (CCURVE(IEDGE,IEL) == ' ')) &
        IFRZER(IEL) = .TRUE. 
    100 END DO
    RETURN
    END SUBROUTINE SETRZER

!-----------------------------------------------------------------------
    SUBROUTINE CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)

!      Check direction of normal of an element face for
!      alignment with the X, Y, or Z axis.

    use size_m
    use geom

    LOGICAL :: IFALGN,IFNORX,IFNORY,IFNORZ

    SUMX    = 0.0
    SUMY    = 0.0
    SUMZ    = 0.0
    TOLNOR  = 1.0e-3
    IFALGN  = .FALSE. 
    IFNORX  = .FALSE. 
    IFNORY  = .FALSE. 
    IFNORZ  = .FALSE. 

    IF (NDIM == 2) THEN
    
        NCPF = NX1
        DO 100 IX=1,NX1
            SUMX = SUMX + ABS( ABS(UNX(IX,1,IFC,IEL)) - 1.0 )
            SUMY = SUMY + ABS( ABS(UNY(IX,1,IFC,IEL)) - 1.0 )
        100 END DO
        SUMX = SUMX / NCPF
        SUMY = SUMY / NCPF
        IF ( SUMX < TOLNOR ) THEN
            IFNORX  = .TRUE. 
            IFALGN = .TRUE. 
        ENDIF
        IF ( SUMY < TOLNOR ) THEN
            IFNORY  = .TRUE. 
            IFALGN = .TRUE. 
        ENDIF
    
    ELSE
    
        NCPF = NX1*NX1
        DO 200 IX=1,NX1
            DO 200 IY=1,NY1
                SUMX = SUMX + ABS( ABS(UNX(IX,IY,IFC,IEL)) - 1.0 )
                SUMY = SUMY + ABS( ABS(UNY(IX,IY,IFC,IEL)) - 1.0 )
                SUMZ = SUMZ + ABS( ABS(UNZ(IX,IY,IFC,IEL)) - 1.0 )
        200 END DO
        SUMX = SUMX / NCPF
        SUMY = SUMY / NCPF
        SUMZ = SUMZ / NCPF
        IF ( SUMX < TOLNOR ) THEN
            IFNORX  = .TRUE. 
            IFALGN = .TRUE. 
        ENDIF
        IF ( SUMY < TOLNOR ) THEN
            IFNORY  = .TRUE. 
            IFALGN = .TRUE. 
        ENDIF
        IF ( SUMZ < TOLNOR ) THEN
            IFNORZ  = .TRUE. 
            IFALGN = .TRUE. 
        ENDIF
    
    ENDIF

    RETURN
    END SUBROUTINE CHKNORD
!-----------------------------------------------------------------------
    SUBROUTINE CHKAXCB

    use size_m
    use input
    CHARACTER CB*3

    IFLD  = 1
    NFACE = 2*NDIM

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFLD)
            IF  (CB == 'A  ' .AND. IFC /= 1)  GOTO 9000
    100 END DO

    RETURN

    9000 WRITE (6,*) ' Element face on the axis of symmetry must be FACE 1'
    WRITE (6,*) ' Element',IEL,'   face',IFC,'  is on the axis.'
    call exitt

    END SUBROUTINE CHKAXCB
!-----------------------------------------------------------------------
    SUBROUTINE CHKCBC (CB,IEL,IFC,IFALGN)

!     Check for illegal boundary conditions

    CHARACTER CB*3
    LOGICAL :: IFALGN

!     Laplacian formulation only

    IF  (CB == 'SH ' .OR.  CB == 'sh ' .OR. &
    CB == 'SHL' .OR.  CB == 'shl' .OR. &
    CB == 'S  ' .OR.  CB == 's  ' .OR. &
    CB == 'SL ' .OR.  CB == 'sl ' .OR. &
    CB == 'MM ' .OR.  CB == 'mm ' .OR. &
    CB == 'MS ' .OR.  CB == 'ms ' .OR. &
    CB == 'MSI' .OR.  CB == 'msi'    )                GOTO 9001
    IF ( .NOT. IFALGN .AND. &
    (CB == 'ON ' .OR.  CB == 'on ' .OR. CB == 'SYM') ) GOTO 9010
    RETURN

    9001 WRITE (6,*) ' Illegal traction boundary conditions detected for'
    GOTO 9999
    9010 WRITE (6,*) ' Mixed B.C. on a side nonaligned with either the X,Y, &
    or Z axis detected for'
    9999 WRITE (6,*) ' Element',IEL,'   side',IFC,'.'
    WRITE (6,*) ' Selected option only allowed for STRESS FORMULATION'
    WRITE (6,*) ' Execution terminates'
    call exitt
    END SUBROUTINE CHKCBC
!-----------------------------------------------------------------------
    SUBROUTINE BCMASK

!     Zero out masks corresponding to Dirichlet boundary points.

    use size_m
    use input
    use mvgeom
    use soln
    use topol
    use tstep

    common  /nekcb/ cb
    character(3) :: cb
    character(1) :: cb1(3)
    equivalence (cb1,cb)

    logical :: ifalgn,ifnorx,ifnory,ifnorz
    integer :: e,f

    NFACES=2*NDIM
    NXYZ  =NX1*NY1*NZ1


!     Masks for moving mesh

    IF (IFMVBD) THEN
#if 0
        IFIELD = 0
        CALL STSMASK (W1MASK,W2MASK,W3MASK)
        do e=1,nelv
            do f=1,nfaces
                if (cbc(f,e,1) == 'msi' .OR. cbc(f,e,1) == 'msi') then
                    call facev(w1mask,e,f,0.0,nx1,ny1,nz1)
                    call facev(w2mask,e,f,0.0,nx1,ny1,nz1)
                    call facev(w3mask,e,f,0.0,nx1,ny1,nz1)
                endif
            enddo
        enddo
#endif
    ENDIF

!     Masks for flow variables

    IF (IFFLOW) THEN
        IFIELD = 1
        NEL    = NELFLD(IFIELD)
        NTOT   = NXYZ*NEL
    
    !        Pressure mask
    
        CALL RONE(PMASK,NTOT)
        DO 50 IEL=1,NELV
            DO 50 IFACE=1,NFACES
                CB=CBC(IFACE,IEL,IFIELD)
                IF (CB == 'O  ' .OR. CB == 'ON ') &
                CALL FACEV(PMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
        50 END DO
    
    !        Zero out mask at Neumann-Dirichlet interfaces
    
        CALL DSOP(PMASK,'MUL',NX1,NY1,NZ1)
    
    !        Velocity masks
    
    !        write(6,*) 'MASK ifstrs',ifstrs,ifield
    !        call exitt
        IF (IFSTRS) THEN
!max            CALL STSMASK (V1MASK,V2MASK,V3MASK)
        ELSE
        
            CALL RONE(V1MASK,NTOT)
            CALL RONE(V2MASK,NTOT)
            CALL RONE(V3MASK,NTOT)
            CALL RONE( OMASK,NTOT)
        
            DO 100 IEL=1,NELV
                DO 100 IFACE=1,NFACES
                    CB =CBC(IFACE,IEL,IFIELD)
                    CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFACE,IEL)
                
                !            All-Dirichlet boundary conditions
                
                    IF (CB == 'v  ' .OR. CB == 'V  ' .OR. CB == 'vl ' .OR. &
                    CB == 'VL ' .OR. CB == 'W  ') THEN
                        CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        GOTO 100
                    ENDIF
                
                !        Mixed-Dirichlet-Neumann boundary conditions
                
                    IF (CB == 'SYM') THEN
                        IF ( .NOT. IFALGN .OR. IFNORX ) &
                        CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        IF ( IFNORY ) &
                        CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        IF ( IFNORZ ) &
                        CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        GOTO 100
                    ENDIF
                    IF (CB == 'ON ') THEN
                        IF ( IFNORY .OR. IFNORZ ) &
                        CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        IF ( .NOT. IFALGN .OR. IFNORX .OR. IFNORZ ) &
                        CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        IF ( .NOT. IFALGN .OR. IFNORX .OR. IFNORY ) &
                        CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        GOTO 100
                    ENDIF
                    IF (CB == 'A  ') THEN
                        CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                        CALL FACEV ( OMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                    ENDIF
            100 END DO

            CALL DSOP  ( OMASK,'MUL',NX1,NY1,NZ1)
            call opdsop(v1mask,v2mask,v3mask,'MUL') ! no rotation for mul


        ENDIF
    
    ENDIF

!     Masks for passive scalars +
!     k and e if k-e turbulence modem:
!     k = nfield-1
!     e = nfield

    IF (IFHEAT) THEN
    
        DO 1200 IFIELD=2,NFIELD
            IPSCAL = IFIELD-1
            NEL    = NELFLD(IFIELD)
            NTOT   = NXYZ*NEL
            CALL RONE (TMASK(1,1,1,1,IPSCAL),NTOT)
            DO 1100 IEL=1,NEL
                DO 1100 IFACE=1,NFACES
                    CB =CBC(IFACE,IEL,IFIELD)
                
                !           Assign mask values.
                
                    IF  (CB == 'T  ' .OR. CB == 't  ' .OR. &
                    (CB == 'A  ' .AND. IFAZIV)    .OR. &
                    CB == 'MCI' .OR. CB == 'MLI' .OR. &
                    CB == 'KD ' .OR. CB == 'kd ' .OR. &
                    CB == 'ED ' .OR. CB == 'ed ' .OR. &
                    CB == 'KW ' .OR. CB == 'KWS' .OR. CB == 'EWS') &
                    CALL FACEV (TMASK(1,1,1,1,IPSCAL), &
                    IEL,IFACE,0.0,NX1,NY1,NZ1)
            1100 END DO
            CALL DSOP (TMASK(1,1,1,1,IPSCAL),'MUL',NX1,NY1,NZ1)
        1200 END DO
    
    ENDIF

!     Masks for B-field

    if (ifmhd) then
        ifield = ifldmhd
        nel    = nelfld(ifield)
        ntot   = nxyz*nel
    
    !        B-field pressure mask
    
        call rone(bpmask,ntot)
        do iel=1,nelv
            do iface=1,nfaces
                cb=cbc(iface,iel,ifield)
                if (cb == 'O  ' .OR. cb == 'ON ') &
                call facev(bpmask,iel,iface,0.0,nx1,ny1,nz1)
            enddo
        enddo
    
    !        Zero out mask at Neumann-Dirichlet interfaces
    
        call dsop(bpmask,'MUL',nx1,ny1,nz1)
    
    !        B-field masks
    
        if (ifstrs) then
!max            call stsmask (b1mask,b2mask,b3mask)
        else
        
            call rone(b1mask,ntot)
            call rone(b2mask,ntot)
            call rone(b3mask,ntot)
        
            do iel=1,nelv
                do iface=1,nfaces
                    cb =cbc(iface,iel,ifield)
                    call chknord (ifalgn,ifnorx,ifnory,ifnorz,iface,iel)
                
                    if (cb == 'v  ' .OR. cb == 'V  ' .OR. cb == 'vl ' .OR. &
                    cb == 'VL ' .OR. cb == 'W  ') then
                    
                    !               All-Dirichlet boundary conditions
                    
                        call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                        call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                        call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                    
                    elseif (cb == 'SYM') then
                    
                    !               Mixed-Dirichlet-Neumann boundary conditions
                    
                        if ( .NOT. ifalgn .OR. ifnorx ) &
                        call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( ifnory ) &
                        call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( ifnorz ) &
                        call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                    
                    elseif (cb == 'ON ') then
                    
                    !               Mixed-Dirichlet-Neumann boundary conditions
                    
                        if ( ifnory .OR. ifnorz ) &
                        call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( .NOT. ifalgn .OR. ifnorx .OR. ifnorz ) &
                        call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( .NOT. ifalgn .OR. ifnorx .OR. ifnory ) &
                        call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                    
                    elseif (cb == 'A  ') then
                    
                    !               axisymmetric centerline
                    
                        call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                    
                    else
                    
                        if ( cb1(1) == 'd' ) &
                        call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( cb1(2) == 'd' ) &
                        call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                        if ( cb1(3) == 'd' .AND. if3d ) &
                        call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                    
                    endif
                enddo
            enddo
        
            call dsop(b1mask,'MUL',nx1,ny1,nz1)
            call dsop(b2mask,'MUL',nx1,ny1,nz1)
            if (ndim == 3) call dsop(b3mask,'MUL',nx1,ny1,nz1)
        endif
    endif

    RETURN
    END SUBROUTINE BCMASK
!-----------------------------------------------------------------------
    SUBROUTINE BCDIRVC(V1,V2,V3,mask1,mask2,mask3)

!     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
!     Use IFIELD as a guide to which boundary conditions are to be applied.

    use ctimer
    use size_m
    use geom
    use input
    use soln
    use topol
    use tstep
    COMMON /SCRUZ/ TMP1(LX1,LY1,LZ1,LELV) &
    , TMP2(LX1,LY1,LZ1,LELV) &
    , TMP3(LX1,LY1,LZ1,LELV)
    COMMON /SCRMG/ TMQ1(LX1,LY1,LZ1,LELV) &
    , TMQ2(LX1,LY1,LZ1,LELV) &
    , TMQ3(LX1,LY1,LZ1,LELV)

    REAL :: V1(NX1,NY1,NZ1,LELV),V2(NX1,NY1,NZ1,LELV) &
    ,V3(NX1,NY1,NZ1,LELV)
    real :: mask1(nx1,ny1,nz1,lelv),mask2(nx1,ny1,nz1,lelv) &
    ,mask3(nx1,ny1,nz1,lelv)

    common  /nekcb/ cb
    character cb*3
    character(1) :: cb1(3)
    equivalence (cb1,cb)

    logical :: ifonbc

    ifonbc = .FALSE. 

#ifndef NOTIMER
    if (icalld == 0) then
        tusbc=0.0
        nusbc=0
        icalld=icalld+1
    endif
    nusbc=nusbc+1
    etime1=dnekclock()
#endif


    NFACES=2*NDIM
    NXYZ  =NX1*NY1*NZ1
    NEL   =NELFLD(IFIELD)
    NTOT  =NXYZ*NEL

    CALL RZERO(TMP1,NTOT)
    CALL RZERO(TMP2,NTOT)
    IF (IF3D) CALL RZERO(TMP3,NTOT)

!     Velocity boundary conditions

!     write(6,*) 'BCDIRV: ifield',ifield
    DO 2100 ISWEEP=1,2
        DO 2000 IE=1,NEL
            DO 2000 IFACE=1,NFACES
                CB  = CBC(IFACE,IE,IFIELD)
                BC1 = BC(1,IFACE,IE,IFIELD)
                BC2 = BC(2,IFACE,IE,IFIELD)
                BC3 = BC(3,IFACE,IE,IFIELD)

                IF (CB == 'V  ' .OR. CB == 'VL '  .OR. &
                CB == 'WS ' .OR. CB == 'WSL') THEN
                    CALL FACEV (TMP1,IE,IFACE,BC1,NX1,NY1,NZ1)
                    CALL FACEV (TMP2,IE,IFACE,BC2,NX1,NY1,NZ1)
                    IF (IF3D) CALL FACEV (TMP3,IE,IFACE,BC3,NX1,NY1,NZ1)
                    IF ( IFQINP(IFACE,IE) ) &
                    CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                    TMP3(1,1,1,IE),IE,IFACE)
                ENDIF

                IF (CB == 'v  ' .OR. CB == 'vl ' .OR. &
                CB == 'ws ' .OR. CB == 'wsl' .OR. &
                CB == 'mv ' .OR. CB == 'mvn' .OR. &
                cb1(1) == 'd' .OR. cb1(2) == 'd' .OR. cb1(3) == 'd') then

                    call faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie), &
                    tmp3(1,1,1,ie),ie,iface,nx1,ny1,nz1)

                    IF ( IFQINP(IFACE,IE) ) &
                    CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                    TMP3(1,1,1,IE),IE,IFACE)
                ENDIF

                IF (CB == 'ON ' .OR. CB == 'on ') then   ! 5/21/01 pff
                    ifonbc = .TRUE. 
                    CALL FACEIV ('v  ',TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                    TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
                ENDIF

        2000 END DO
        DO 2010 IE=1,NEL
            DO 2010 IFACE=1,NFACES
                IF (CBC(IFACE,IE,IFIELD) == 'W  ') THEN
                    CALL FACEV (TMP1,IE,IFACE,0.0,NX1,NY1,NZ1)
                    CALL FACEV (TMP2,IE,IFACE,0.0,NX1,NY1,NZ1)
                    IF (IF3D) CALL FACEV (TMP3,IE,IFACE,0.0,NX1,NY1,NZ1)
                ENDIF
        2010 END DO
    
    !        Take care of Neumann-Dirichlet shared edges...
    
        if (isweep == 1) then
            call opdsop(tmp1,tmp2,tmp3,'MXA')
        else
            call opdsop(tmp1,tmp2,tmp3,'MNA')
        endif
    2100 END DO

!     Copy temporary array to velocity arrays.

    IF ( .NOT. IFSTRS ) THEN
        CALL COL2(V1,mask1,NTOT)
        CALL COL2(V2,mask2,NTOT)
        IF (IF3D) CALL COL2(V3,mask3,NTOT)
        if (ifonbc) then
            call antimsk1(tmp1,mask1,ntot)
            call antimsk1(tmp2,mask2,ntot)
            if (if3d) call antimsk1(tmp3,mask3,ntot)
        endif
    ELSE
#if 0
        IF (IFMODEL) THEN
            CALL COPY (TMQ1,TMP1,NTOT)
            CALL COPY (TMQ2,TMP2,NTOT)
            IF (NDIM == 3) CALL COPY (TMQ3,TMP3,NTOT)
            CALL AMASK (TMP1,TMP2,TMP3,TMQ1,TMQ2,TMQ3,NELV)
        ENDIF
        CALL RMASK (V1,V2,V3,NELV)
#endif
    ENDIF

    CALL ADD2(V1,TMP1,NTOT)
    CALL ADD2(V2,TMP2,NTOT)
    IF (IF3D) CALL ADD2(V3,TMP3,NTOT)


#ifndef NOTIMER
    tusbc=tusbc+(dnekclock()-etime1)
#endif

    RETURN
    END SUBROUTINE BCDIRVC
!-----------------------------------------------------------------------
!> \brief Apply Dirichlet boundary conditions to surface of scalar, S.
!! Use IFIELD as a guide to which boundary conditions are to be applied.
SUBROUTINE BCDIRSC(S)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, ndim, nfield
  use ctimer, only : icalld, tusbc, nusbc, etime1, dnekclock
  use input, only : cbc, bc
  use soln, only : tmask
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: S(LX1,LY1,LZ1,LELT)
  real(DP) :: tmp(LX1,LY1,LZ1,LELT) 

  CHARACTER CB*3
  common  /nekcb/ cb

  integer :: ifld, nfaces, nxyz, nel, ntot, nvldt, isweep, ie, iface, nfldt
  real(DP) :: BC1, BC2, BC3, BC4, BCK, BCE

#ifndef NOTIMER
  if (icalld == 0) then
      tusbc=0.0
      nusbc=0
      icalld=icalld+1
  endif
  nusbc=nusbc+1
  etime1=dnekclock()
#endif

  IFLD   = 1
  NFACES = 2*NDIM
  NXYZ   = NX1*NY1*NZ1
  NEL    = NELFLD(IFIELD)
  NTOT   = NXYZ*NEL
  NFLDT  = NFIELD - 1

  CALL RZERO(TMP,NTOT)

!     Temperature boundary condition

  DO ISWEEP=1,2
    
#if 0
      IF (IFMODEL .AND. IFKEPS .AND. IFIELD >= NFLDT) &
      CALL TURBWBC (TMP,TMA,SMU)
#endif
    
      DO IE=1,NEL
          DO IFACE=1,NFACES
              CB=CBC(IFACE,IE,IFIELD)
              BC1=BC(1,IFACE,IE,IFIELD)
              BC2=BC(2,IFACE,IE,IFIELD)
              BC3=BC(3,IFACE,IE,IFIELD)
              BC4=BC(4,IFACE,IE,IFIELD)
              BCK=BC(4,IFACE,IE,IFLD)
              BCE=BC(5,IFACE,IE,IFLD)
              IF (CB == 'T  ') CALL FACEV (TMP,IE,IFACE,BC1,NX1,NY1,NZ1)
              IF (CB == 'MCI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
              IF (CB == 'MLI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
              IF (CB == 'KD ') CALL FACEV (TMP,IE,IFACE,BCK,NX1,NY1,NZ1)
              IF (CB == 'ED ') CALL FACEV (TMP,IE,IFACE,BCE,NX1,NY1,NZ1)
              IF (CB == 't  ' .OR. CB == 'kd ' .OR. CB == 'ed ') &
              CALL FACEIS (CB,TMP(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
          END DO
      enddo
  
  !        Take care of Neumann-Dirichlet shared edges...
  
      IF (ISWEEP == 1) CALL DSOP(TMP,'MXA',NX1,NY1,NZ1)
      IF (ISWEEP == 2) CALL DSOP(TMP,'MNA',NX1,NY1,NZ1)
  END DO

!     Copy temporary array to temperature array.

  CALL COL2(S,TMASK(1,1,1,1,IFIELD-1),NTOT)
  CALL ADD2(S,TMP,NTOT)

#ifndef NOTIMER
  tusbc=tusbc+(dnekclock()-etime1)
#endif

  RETURN
END SUBROUTINE BCDIRSC

!-----------------------------------------------------------------------
    SUBROUTINE BCNEUSC(S,ITYPE)

!     Apply Neumann boundary conditions to surface of scalar, S.
!     Use IFIELD as a guide to which boundary conditions are to be applied.

!     If ITYPE = 1, then S is returned as the rhs contribution to the
!                   volumetric flux.

!     If ITYPE =-1, then S is returned as the lhs contribution to the
!                   diagonal of A.


    use ctimer
    use size_m
    use nekuse
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    DIMENSION S(LX1,LY1,LZ1,LELT)
    common  /nekcb/ cb
    CHARACTER CB*3

#ifndef NOTIMER
    if (icalld == 0) then
        tusbc=0.0
        nusbc=0
        icalld=icalld+1
    endif
    nusbc=nusbc+1
    etime1=dnekclock()
#endif

    NFACES=2*NDIM
    NXYZ  =NX1*NY1*NZ1
    NEL   =NELFLD(IFIELD)
    NTOT  =NXYZ*NEL
    CALL RZERO(S,NTOT)

    IF (ITYPE == -1) THEN
    
    !        Compute diagonal contributions to accomodate Robin boundary conditions
    
        DO 1000 IE=1,NEL
            DO 1000 IFACE=1,NFACES
                ieg=lglel(ie)
                CB =CBC(IFACE,IE,IFIELD)
                IF (CB == 'C  ' .OR. CB == 'c  ' .OR. &
                CB == 'R  ' .OR. CB == 'r  ') THEN
                
                    IF (CB == 'C  ') HC   = BC(2,IFACE,IE,IFIELD)
                    IF (CB == 'R  ') THEN
                        TINF = BC(1,IFACE,IE,IFIELD)
                        HRAD = BC(2,IFACE,IE,IFIELD)
                    ENDIF
                    IA=0
                
                ! IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
                
                    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
                    DO 100 IZ=KZ1,KZ2
                        DO 100 IY=KY1,KY2
                            DO 100 IX=KX1,KX2
                                IA = IA + 1
                                TS = T(IX,IY,IZ,IE,IFIELD-1)
                                IF (CB == 'c  ' .OR. CB == 'r  ') THEN
                                    CALL NEKASGN (IX,IY,IZ,IE)
                                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                                ENDIF
                                IF (CB == 'r  ' .OR. CB == 'R  ') &
                                HC = HRAD * (TINF**2 + TS**2) * (TINF + TS)
                                S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) + &
                                HC*AREA(IA,1,IFACE,IE)/BM1(IX,IY,IZ,IE)
                    100 END DO
                ENDIF
        1000 END DO
    ENDIF
    IF (ITYPE == 1) THEN
    
    !        Add passive scalar fluxes to rhs
    
        DO 2000 IE=1,NEL
            DO 2000 IFACE=1,NFACES
                ieg=lglel(ie)
                CB =CBC(IFACE,IE,IFIELD)
                IF (CB == 'F  ' .OR. CB == 'f  ' .OR. &
                CB == 'C  ' .OR. CB == 'c  ' .OR. &
                CB == 'R  ' .OR. CB == 'r  ' ) THEN
                
                    IF (CB == 'F  ') FLUX=BC(1,IFACE,IE,IFIELD)
                    IF (CB == 'C  ') FLUX=BC(1,IFACE,IE,IFIELD) &
                    *BC(2,IFACE,IE,IFIELD)
                    IF (CB == 'R  ') THEN
                        TINF=BC(1,IFACE,IE,IFIELD)
                        HRAD=BC(2,IFACE,IE,IFIELD)
                    ENDIF
                
                !              Add local weighted flux values to rhs, S.
                
                ! IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
                    IA=0
                    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
                    DO 200 IZ=KZ1,KZ2
                        DO 200 IY=KY1,KY2
                            DO 200 IX=KX1,KX2
                                IA = IA + 1
                                TS = T(IX,IY,IZ,IE,IFIELD-1)
                                IF (CB == 'f  ') THEN
                                    CALL NEKASGN (IX,IY,IZ,IE)
                                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                                ENDIF
                                IF (CB == 'c  ') THEN
                                    CALL NEKASGN (IX,IY,IZ,IE)
                                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                                    FLUX = TINF*HC
                                ENDIF
                                IF (CB == 'r  ') THEN
                                    CALL NEKASGN (IX,IY,IZ,IE)
                                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                                ENDIF
                                IF (CB == 'R  ' .OR. CB == 'r  ') &
                                FLUX = HRAD*(TINF**2 + TS**2)*(TINF + TS) * TINF
                            
                            !                 Add computed fluxes to boundary surfaces:
                            
                                S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) &
                                + FLUX*AREA(IA,1,IFACE,IE)
                    200 END DO
                ENDIF
        2000 END DO
    ENDIF

#ifndef NOTIMER
    tusbc=tusbc+(dnekclock()-etime1)
#endif

    RETURN
    END SUBROUTINE BCNEUSC
!-----------------------------------------------------------------------
    SUBROUTINE FACEIS (CB,S,IEL,IFACE,NX,NY,NZ)

!     Assign inflow boundary conditions to face(IE,IFACE)
!     for scalar S.

    use size_m
    use nekuse
    use parallel
    use soln      ! tmask()   11/19/2010
    use tstep     ! ifield    11/19/2010

    DIMENSION S(LX1,LY1,LZ1)
    CHARACTER CB*3

    common  /nekcb/ cb3
    character(3) :: cb3
    cb3 = cb

    ifld1 = ifield-1


!     Passive scalar term

    ieg = lglel(iel)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)

    IF (CB == 't  ') THEN

        DO 100 IZ=KZ1,KZ2                           !  11/19/2010: The tmask() screen
            DO 100 IY=KY1,KY2                           !  added here so users can leave
                DO 100 IX=KX1,KX2                           !  certain points to be Neumann,
                    if (tmask(ix,iy,iz,iel,ifld1) == 0) then !  if desired.
                        CALL NEKASGN (IX,IY,IZ,IEL)
                        CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                        S(IX,IY,IZ) = TEMP
                    endif
        100 END DO
        RETURN
    
    ELSEIF (CB == 'ms ' .OR. CB == 'msi') THEN
    
        DO 200 IZ=KZ1,KZ2
            DO 200 IY=KY1,KY2
                DO 200 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    S(IX,IY,IZ) = SIGMA
        200 END DO
    
    ELSEIF (CB == 'kd ') THEN
    
        DO 300 IZ=KZ1,KZ2
            DO 300 IY=KY1,KY2
                DO 300 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    S(IX,IY,IZ) = TURBK
        300 END DO
    
    ELSEIF (CB == 'ed ') THEN
    
        DO 400 IZ=KZ1,KZ2
            DO 400 IY=KY1,KY2
                DO 400 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    S(IX,IY,IZ) = TURBE
        400 END DO
    
    ENDIF

    RETURN
    END SUBROUTINE FACEIS
!-----------------------------------------------------------------------
    SUBROUTINE FACEIV (CB,V1,V2,V3,IEL,IFACE,NX,NY,NZ)

!     Assign fortran function boundary conditions to
!     face IFACE of element IEL for vector (V1,V2,V3).

    use size_m
    use nekuse
    use parallel

    dimension v1(nx,ny,nz),v2(nx,ny,nz),v3(nx,ny,nz)
    character cb*3

    character(1) :: cb1(3)

    common  /nekcb/ cb3
    character(3) :: cb3
    cb3 = cb

    call chcopy(cb1,cb,3)

    ieg = lglel(iel)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)

    IF (CB == 'v  ' .OR. CB == 'ws ' .OR. CB == 'mv ' .OR. &
    CB == 'mvn') THEN
    
        DO 100 IZ=KZ1,KZ2
            DO 100 IY=KY1,KY2
                DO 100 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = UX
                    V2(IX,IY,IZ) = UY
                    V3(IX,IY,IZ) = UZ
        100 END DO
        RETURN
    
    elseif (cb1(1) == 'd' .OR. cb1(2) == 'd' .OR. cb1(3) == 'd') then
    
        do iz=kz1,kz2
            do iy=ky1,ky2
                do ix=kx1,kx2
                    call nekasgn (ix,iy,iz,iel)
                    call userbc  (ix,iy,iz,iface,ieg)
                    if (cb1(1) == 'd') v1(ix,iy,iz) = ux
                    if (cb1(2) == 'd') v2(ix,iy,iz) = uy
                    if (cb1(3) == 'd') v3(ix,iy,iz) = uz
                enddo
            enddo
        enddo
        return
    
    ELSEIF (CB == 'vl ' .OR. CB == 'wsl') THEN
    
        DO 120 IZ=KZ1,KZ2
            DO 120 IY=KY1,KY2
                DO 120 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = UN
                    V2(IX,IY,IZ) = U1
                    V3(IX,IY,IZ) = U2
        120 END DO
        RETURN
    
    ELSEIF (CB == 's  ' .OR. CB == 'sh ') THEN
    
        DO 200 IZ=KZ1,KZ2
            DO 200 IY=KY1,KY2
                DO 200 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = TRX
                    V2(IX,IY,IZ) = TRY
                    V3(IX,IY,IZ) = TRZ
        200 END DO
        RETURN
    
    ELSEIF (CB == 'sl ' .OR. CB == 'shl') THEN
    
        DO 220 IZ=KZ1,KZ2
            DO 220 IY=KY1,KY2
                DO 220 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = TRN
                    V2(IX,IY,IZ) = TR1
                    V3(IX,IY,IZ) = TR2
        220 END DO
    
    ELSEIF (CB == 'ms ') THEN
    
        DO 240 IZ=KZ1,KZ2
            DO 240 IY=KY1,KY2
                DO 240 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = -PA
                    V2(IX,IY,IZ) = TR1
                    V3(IX,IY,IZ) = TR2
        240 END DO
    
    ELSEIF (CB == 'on ' .OR. CB == 'o  ') THEN
    
        DO 270 IZ=KZ1,KZ2
            DO 270 IY=KY1,KY2
                DO 270 IX=KX1,KX2
                    CALL NEKASGN (IX,IY,IZ,IEL)
                    CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                    V1(IX,IY,IZ) = -PA
                    V2(IX,IY,IZ) = 0.0
                    V3(IX,IY,IZ) = 0.0
        270 END DO
    
    ENDIF

    RETURN
    END SUBROUTINE FACEIV
!-----------------------------------------------------------------------
    SUBROUTINE NEKASGN (IX,IY,IZ,IEL)

!     Assign NEKTON variables for definition (by user) of
!     boundary conditions at collocation point (IX,IY,IZ)
!     of element IEL.

!       X             X-coordinate
!       Y             Y-coordinate
!       Z             Z-coordinate
!       UX            X-velocity
!       UY            Y-velocity
!       UZ            Z-velocity
!       TEMP          Temperature
!       PS1           Passive scalar No. 1
!       PS2           Passive scalar No. 2
!        .             .
!        .             .
!       PS9           Passive scalar No. 9
!       SI2           Strainrate invariant II
!       SI3           Strainrate invariant III

!     Variables to be defined by user for imposition of
!     boundary conditions :

!       SH1           Shear component No. 1
!       SH2           Shear component No. 2
!       TRX           X-traction
!       TRY           Y-traction
!       TRZ           Z-traction
!       SIGMA         Surface-tension coefficient
!       FLUX          Flux
!       HC            Convection heat transfer coefficient
!       HRAD          Radiation  heat transfer coefficient
!       TINF          Temperature at infinity

    use size_m
    use geom
    use input
    use nekuse
    use soln
    use tstep

    common  /nekcb/ cb
    CHARACTER CB*3

    COMMON /SCREV / SII (LX1,LY1,LZ1,LELT) &
    , SIII(LX1,LY1,LZ1,LELT)

    X     = XM1(IX,IY,IZ,IEL)
    Y     = YM1(IX,IY,IZ,IEL)
    Z     = ZM1(IX,IY,IZ,IEL)
    R     = X**2+Y**2
    IF (R > 0.0) R=SQRT(R)
    IF (X /= 0.0 .OR. Y /= 0.0) THETA = ATAN2(Y,X)

    UX    = VX(IX,IY,IZ,IEL)
    UY    = VY(IX,IY,IZ,IEL)
    UZ    = VZ(IX,IY,IZ,IEL)
    TEMP  = T(IX,IY,IZ,IEL,1)
    DO 100 IPS=1,NPSCAL
        PS(IPS) = T(IX,IY,IZ,IEL,IPS+1)
    100 END DO
    SI2   = SII (IX,IY,IZ,IEL)
    SI3   = SIII(IX,IY,IZ,IEL)
    UDIFF = VDIFF (IX,IY,IZ,IEL,IFIELD)
    UTRANS= VTRANS(IX,IY,IZ,IEL,IFIELD)

    cbu   = cb

    RETURN
    END SUBROUTINE NEKASGN
!-----------------------------------------------------------------------
    SUBROUTINE TRSTAX (TRX,TRY,SIGST,IEL,IFC)

!     Compute taction due to surface tension (axisymmetric)

    use size_m
    use dxyz
    use geom
    use topol
    use wz_m
    COMMON /CTMP1/ A1X(LX1),A1Y(LX1),A2X(LX1),A2Y(LX1) &
    , STX(LX1),STY(LX1),XJM1(LX1)
    COMMON /CTMP0/ XFM1(LX1),YFM1(LX1),T1XF(LX1),T1YF(LX1)

    DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),SIGST(LX1,LY1)
    DIMENSION CANG(2),SANG(2)
    DIMENSION IXN(2),IYN(2),IAN(2)
    LOGICAL :: IFGLJ

    IFGLJ = .FALSE. 
    IF ( IFRZER(IEL) .AND. (IFC == 2 .OR. IFC == 4) ) IFGLJ = .TRUE. 
    CALL FACEC2 (XFM1,YFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),IFC)

    IF (IFGLJ) THEN
        CALL MXM (DAM1,NY1,XFM1,NY1,T1XF,1)
        CALL MXM (DAM1,NY1,YFM1,NY1,T1YF,1)
        YS0 = T1YF(1)
    ELSE
        CALL MXM (DXM1,NX1,XFM1,NX1,T1XF,1)
        CALL MXM (DXM1,NX1,YFM1,NX1,T1YF,1)
    ENDIF

    DO 10 IX=1,NX1
        XJM1(IX)=SQRT( T1XF(IX)**2 + T1YF(IX)**2 )
        T1XF(IX)=T1XF(IX) / XJM1(IX)
        T1YF(IX)=T1YF(IX) / XJM1(IX)
    10 END DO

    IF ( IFGLJ ) THEN
        CALL MXM (DAM1,1,T1XF,NY1,T1XS0,1)
        CALL MXM (DAM1,1,UNY(1,1,IFC,IEL),NY1,UNYS0,1)
        DDX    = WAM1(1)*SIGST(1,1)*T1XS0*YS0
        DDY    = WAM1(1)*SIGST(1,1)*T1YF(1)*YS0*2.0
        A2X(1) = WAM1(1)*SIGST(1,1)*XJM1(1)*UNX(1,1,IFC,IEL)*UNYS0
        A2Y(1) = 0.0
        STX(1) = 0.0
        STY(1) = 0.0
        DO 100 IY=2,NY1
            AA = WAM1(IY) * SIGST(IY,1) / (1.0 + ZAM1(IY))
            STX(IY) = T1XF(IY) * AA
            STY(IY) = T1YF(IY) * AA
            AA = AA * XJM1(IY) * UNY(IY,1,IFC,IEL)
            A2X(IY) = UNX(IY,1,IFC,IEL) * AA
            A2Y(IY) = UNY(IY,1,IFC,IEL) * AA
        100 END DO
    ELSE
        DO 200 IX=1,NX1
            AA = SIGST(IX,1) * WXM1(IX)
            STX(IX) = T1XF(IX) * AA
            STY(IX) = T1YF(IX) * AA
            AA = AA * XJM1(IX) * UNY(IX,1,IFC,IEL)
            A2X(IX) = UNX(IX,1,IFC,IEL) * AA
            A2Y(IX) = UNY(IX,1,IFC,IEL) * AA
        200 END DO
    ENDIF

    IF (IFGLJ) THEN
        DO 220 IY=1,NY1
            YSIY = T1YF(IY)*XJM1(IY)
            DTX1 = 0.0
            DTY1 = DATM1(IY,1)*DDY
            DTX2 = YSIY*STX(IY)
            DTY2 = YSIY*STY(IY)
            DTY3 = 0.0
            DO 240 J=2,NY1
                DTYS = DATM1(IY,J)*YFM1(J)
                DTX1 = DTX1 + DTYS*STX(J)
                DTY3 = DTY3 + DTYS*STY(J)
            240 END DO
            A1X(IY) = DTX1 + DTX2
            A1Y(IY) = DTY1 + DTY2 + DTY3
        220 END DO
        A1X(1)  = A1X(1) + DDX
    ELSE
        CALL MXM  (DXTM1,NX1,STX,NX1,A1X,1)
        CALL MXM  (DXTM1,NX1,STY,NX1,A1Y,1)
        CALL COL2 (A1X,YFM1,NX1)
        CALL COL2 (A1Y,YFM1,NX1)
    ENDIF

    CALL DSSET (NX1,NY1,NZ1)
    IFACE  = EFACE1(IFC)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)
    I = 0

    DO 300 J2=JS2,JF2,JSKIP2
        DO 300 J1=JS1,JF1,JSKIP1
            I  = I + 1
            TRX(J1,J2,1) = TRX(J1,J2,1) - A2X(I) - A1X(I)
            TRY(J1,J2,1) = TRY(J1,J2,1) - A2Y(I) - A1Y(I)
    300 END DO

!     Contact angle corrections

    CALL CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)
    DO 500 I=1,2
        IX = IXN(I)
        IY = IYN(I)
        IA = IAN(I)
        AA = SIGST(IA,1)*YM1(IX,IY,1,IEL)
        TRX(IX,IY,1)=TRX(IX,IY,1) + AA*CANG(I)
        TRY(IX,IY,1)=TRY(IX,IY,1) + AA*SANG(I)
    500 END DO

    RETURN
    END SUBROUTINE TRSTAX
!-----------------------------------------------------------------------
    SUBROUTINE CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)

    use size_m
    use geom
    use input
    use soln

    DIMENSION CANG(2),SANG(2)
    DIMENSION IXN(2),IYN(2),IAN(2),ISN(2),NEBPT(4,2)
    CHARACTER CBN*3

    DATA NEBPT /4,1,2,3, 2,3,4,1/
    IFLD = 1
    EPS  = 1.e-6

    DO 100 I=1,2
        IFCN    = NEBPT(IFC,I)
        CBN     = CBC(IFCN,IEL,IFLD)
        IXN(I)  = 1
        IYN(I)  = 1
        IAN(I)  = 1
        ISN(I)  = 1
        CANG(I) = 0.0
        SANG(I) = 0.0
        IF (CBN == 'E  ' .OR. CBN == 'P  ' .OR. cbn == 'p  ' .OR. &
        CBN(1:1) == 'M' .OR. CBN(1:1) == 'm') GOTO 100
        NC = IFC
        IF (I == 2) NC=IFCN
        IF (NC  == 2 .OR. NC  == 3) IXN(I) = NX1
        IF (NC  == 3 .OR. NC  == 4) IYN(I) = NY1
        IF (IFC == 2 .OR. IFC == 3) ISN(I) = NX1
        IF (IFCN == 2 .OR. IFCN == 3) IAN(I) = NX1
        IX = IXN(I)
        IY = IYN(I)
        IA = IAN(I)
        IS = ISN(I)
        IF (CBN(1:1) == 'V'   .OR. CBN(1:1) == 'v'   .OR. &
        CBN     == 'S  ' .OR. CBN     == 's  ' .OR. &
        CBN     == 'SL ' .OR. CBN     == 'sl ' .OR. &
        CBN(1:1) == 'O'   .OR. CBN(1:1) == 'o' ) THEN
            UX=VX(IX,IY,1,IEL)
            UY=VY(IX,IY,1,IEL)
            UM=UX**2 + UY**2
            IF (UM > EPS) THEN
                UNLX=UNX(IS,1,IFCN,IEL)
                UNLY=UNY(IS,1,IFCN,IEL)
                UM=SQRT(UM)
                DOT =UX*UNLX + UY*UNLY
                IF (DOT < 0.0) UM=-UM
                CANG(I)=UX/UM
                SANG(I)=UY/UM
                GOTO 100
            ENDIF
        ENDIF
        CANG(I)=UNX(IS,1,IFCN,IEL)
        SANG(I)=UNY(IS,1,IFCN,IEL)
    100 END DO

    RETURN
    END SUBROUTINE CTANG2D
!-----------------------------------------------------------------------
    SUBROUTINE GLOBROT (R1,R2,R3,IEL,IFC)

!     Rotate vector components R1,R2,R3 at face IFC
!     of element IEL from local to global system.

!     R1, R2, R3 have the (NX,NY,NZ) data structure
!     IFACE1 is in the preprocessor notation
!     IFACE  is the dssum notation.

    use size_m
    use geom
    use topol

    DIMENSION R1(LX1,LY1,LZ1) &
    , R2(LX1,LY1,LZ1) &
    , R3(LX1,LY1,LZ1)

    CALL DSSET (NX1,NY1,NZ1)
    IFACE  = EFACE1(IFC)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)
    I = 0

    IF (NDIM == 2) THEN
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                I = I+1
                RNORL = R1(J1,J2,1)
                RTAN1 = R2(J1,J2,1)
                R1(J1,J2,1) = RNORL*UNX(I,1,IFC,IEL) + &
                RTAN1*T1X(I,1,IFC,IEL)
                R2(J1,J2,1) = RNORL*UNY(I,1,IFC,IEL) + &
                RTAN1*T1Y(I,1,IFC,IEL)
        200 END DO
    ELSE
        DO 300 J2=JS2,JF2,JSKIP2
            DO 300 J1=JS1,JF1,JSKIP1
                I = I+1
                RNORL = R1(J1,J2,1)
                RTAN1 = R2(J1,J2,1)
                RTAN2 = R3(J1,J2,1)
                R1(J1,J2,1) = RNORL*UNX(I,1,IFC,IEL) + &
                RTAN1*T1X(I,1,IFC,IEL) + &
                RTAN2*T2X(I,1,IFC,IEL)
                R2(J1,J2,1) = RNORL*UNY(I,1,IFC,IEL) + &
                RTAN1*T1Y(I,1,IFC,IEL) + &
                RTAN2*T2Y(I,1,IFC,IEL)
                R3(J1,J2,1) = RNORL*UNZ(I,1,IFC,IEL) + &
                RTAN1*T1Z(I,1,IFC,IEL) + &
                RTAN2*T2Z(I,1,IFC,IEL)
        300 END DO
    ENDIF

    RETURN
    END SUBROUTINE GLOBROT
!-----------------------------------------------------------------------
    SUBROUTINE FACEC2 (A1,A2,B1,B2,IFC)

!     2-D Geometry only
!     Extract A1,A2 from B1,B2 on surface IFC.

!     A1, A2 have the (NX1,  1,NFACE) data structure
!     B1, B2 have the (NX1,NY1,    1) data structure

    use size_m

    DIMENSION A1(LX1),A2(LX1),B1(LX1,LY1),B2(LX1,LY1)

    IX=1
    IY=1
    IF (IFC == 1 .OR. IFC == 3) THEN
        IF (IFC == 3) IY = NY1
        DO 10 IX=1,NX1
            A1(IX)=B1(IX,IY)
            A2(IX)=B2(IX,IY)
        10 END DO
    ELSE
        IF (IFC == 2) IX = NX1
        DO 20 IY=1,NY1
            A1(IY)=B1(IX,IY)
            A2(IY)=B2(IX,IY)
        20 END DO
    ENDIF

    RETURN
    END SUBROUTINE FACEC2
!-----------------------------------------------------------------------
    SUBROUTINE LFALSE (IFA,N)
    LOGICAL :: IFA(1)
    DO 100 I=1,N
        IFA(I)= .FALSE. 
    100 END DO
    RETURN
    END SUBROUTINE LFALSE
!-----------------------------------------------------------------------
    SUBROUTINE RZERO3 (A,B,C,N)
    DIMENSION A(1),B(1),C(1)
    DO 100 I=1,N
        A(I)=0.0
        B(I)=0.0
        C(I)=0.0
    100 END DO
    RETURN
    END SUBROUTINE RZERO3
!-----------------------------------------------------------------------
    SUBROUTINE UNITVEC (X,Y,Z,N)
    DIMENSION X(1),Y(1),Z(1)
    DO 100 I=1,N
        XLNGTH = SQRT( X(I)**2 + Y(I)**2 + Z(I)**2 )
        IF (XLNGTH /= 0.0) THEN
            X(I) = X(I)/XLNGTH
            Y(I) = Y(I)/XLNGTH
            Z(I) = Z(I)/XLNGTH
        ENDIF
    100 END DO
    RETURN
    END SUBROUTINE UNITVEC
!-----------------------------------------------------------------------
    SUBROUTINE ANTIMSK1(X,XMASK,N)
!------------------------------------------------------------------

!     Return only Dirichlet boundary values of X

!-------------------------------------------------------------------
    use opctr
    REAL ::  X(1),XMASK(1)

    DO 100 I=1,N
        X(I) = X(I)*(1.-XMASK(I))
    100 END DO
    RETURN
    END SUBROUTINE ANTIMSK1
!-----------------------------------------------------------------------
