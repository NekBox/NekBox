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
            CALL STSMASK (V1MASK,V2MASK,V3MASK)
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
            call stsmask (b1mask,b2mask,b3mask)
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
        IF (IFMODEL) THEN
            CALL COPY (TMQ1,TMP1,NTOT)
            CALL COPY (TMQ2,TMP2,NTOT)
            IF (NDIM == 3) CALL COPY (TMQ3,TMP3,NTOT)
            CALL AMASK (TMP1,TMP2,TMP3,TMQ1,TMQ2,TMQ3,NELV)
        ENDIF
        CALL RMASK (V1,V2,V3,NELV)
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
    SUBROUTINE BCDIRSC(S)

!     Apply Dirichlet boundary conditions to surface of scalar, S.
!     Use IFIELD as a guide to which boundary conditions are to be applied.

    use ctimer
    use size_m
    use input
    use soln
    use topol
    use tstep

    DIMENSION S(LX1,LY1,LZ1,LELT)
    COMMON /SCRSF/ TMP(LX1,LY1,LZ1,LELT) &
    , TMA(LX1,LY1,LZ1,LELT) &
    , SMU(LX1,LY1,LZ1,LELT)
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

    IFLD   = 1
    NFACES = 2*NDIM
    NXYZ   = NX1*NY1*NZ1
    NEL    = NELFLD(IFIELD)
    NTOT   = NXYZ*NEL
    NFLDT  = NFIELD - 1

    CALL RZERO(TMP,NTOT)

!     Temperature boundary condition

    DO 2100 ISWEEP=1,2
    
#if 0
        IF (IFMODEL .AND. IFKEPS .AND. IFIELD >= NFLDT) &
        CALL TURBWBC (TMP,TMA,SMU)
#endif
    
        DO 2010 IE=1,NEL
            DO 2010 IFACE=1,NFACES
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
        2010 END DO
    
    !        Take care of Neumann-Dirichlet shared edges...
    
        IF (ISWEEP == 1) CALL DSOP(TMP,'MXA',NX1,NY1,NZ1)
        IF (ISWEEP == 2) CALL DSOP(TMP,'MNA',NX1,NY1,NZ1)
    2100 END DO

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
    INCLUDE 'TOTAL'

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
    SUBROUTINE BCNEUTR

    use size_m
    use geom
    use input
    use soln
    COMMON /SCRSF/ TRX(LX1,LY1,LZ1) &
    , TRY(LX1,LY1,LZ1) &
    , TRZ(LX1,LY1,LZ1)
    COMMON /CTMP0/ STC(LX1,LY1,LZ1)
    REAL :: SIGST(LX1,LY1)

    LOGICAL :: IFALGN,IFNORX,IFNORY,IFNORZ
    common  /nekcb/ cb
    CHARACTER CB*3

    IFLD  = 1
    NFACE = 2*NDIM
    NXY1  = NX1*NY1
    NXYZ1 = NX1*NY1*NZ1

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
        
            CB  = CBC (IFC,IEL,IFLD)
            BC1 = BC(1,IFC,IEL,IFLD)
            BC2 = BC(2,IFC,IEL,IFLD)
            BC3 = BC(3,IFC,IEL,IFLD)
            BC4 = BC(4,IFC,IEL,IFLD)
            CALL RZERO3 (TRX,TRY,TRZ,NXYZ1)
        
        !        Prescribed tractions and shear tractions
        
            IF (CB == 'S  ' .OR. CB == 'SL ' .OR. &
            CB == 'SH ' .OR. CB == 'SHL' ) THEN
                CALL TRCON (TRX,TRY,TRZ,BC1,BC2,BC3,IEL,IFC)
                IF (IFQINP(IFC,IEL)) CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                GOTO 120
            ENDIF
            IF (CB == 's  ' .OR. CB == 'sl ' .OR. &
            CB == 'sh ' .OR. CB == 'shl' ) THEN
                CALL FACEIV (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
                CALL FACCVS (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
                IF (IFQINP(IFC,IEL)) CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                GOTO 120
            ENDIF
        
        !        Prescribed outflow ambient pressure
        
            IF (CB == 'ON ' .OR. CB == 'O  ') THEN
                BCN = -BC1
                BC2 =  0.0
                BC3 =  0.0
                CALL TRCON   (TRX,TRY,TRZ,BCN,BC2,BC3,IEL,IFC)
                CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                GOTO 120
            ENDIF
            IF (CB == 'on ' .OR. CB == 'o  ') THEN
                CALL FACEIV  (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
                CALL FACCVS  (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
                CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                GOTO 120
            ENDIF
        
        !     Surface-tension
        
            IF (CB == 'MS ' .OR. CB == 'MSI' .OR. &
            CB == 'MM ' .OR. CB == 'mm ' .OR. &
            CB == 'ms ' .OR. CB == 'msi') THEN
                IF (CB == 'MS ' .OR. cb == 'MM ') THEN
                    BCN = -BC1
                    CALL TRCON   (TRX,TRY,TRZ,BCN,BC2,BC3,IEL,IFC)
                    CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                ENDIF
            !            IF (CB.EQ.'ms '.or.cb.eq.'mm ') THEN
                IF (CB == 'ms ' .OR. cb == 'msi') THEN
                    CALL FACEIV  (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
                    CALL FACCVS  (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
                    CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
                ENDIF
                IF (CB(1:1) == 'M') THEN
                    CALL CFILL  (SIGST,BC4,NXY1)
                ELSE
                    CALL FACEIS (CB,STC,IEL,IFC,NX1,NY1,NZ1)
                    CALL FACEXS (SIGST,STC,IFC,0)
                ENDIF
                IF (IFAXIS) THEN
                    CALL TRSTAX (TRX,TRY,SIGST,IEL,IFC)
                ELSEIF (NDIM == 2) THEN
                    CALL TRST2D (TRX,TRY,SIGST,IEL,IFC)
                ELSE
                    CALL TRST3D (TRX,TRY,TRZ,SIGST,IEL,IFC)
                ENDIF
            ENDIF
        
            120 CALL ADD2 (BFX(1,1,1,IEL),TRX,NXYZ1)
            CALL ADD2 (BFY(1,1,1,IEL),TRY,NXYZ1)
            IF (NDIM == 3) CALL ADD2 (BFZ(1,1,1,IEL),TRZ,NXYZ1)
        
    100 END DO

    RETURN
    END SUBROUTINE BCNEUTR
!-----------------------------------------------------------------------
    SUBROUTINE TRCON (TRX,TRY,TRZ,TR1,TR2,TR3,IEL,IFC)

    use size_m
    use geom
    use topol

    DIMENSION TRX(LX1,LY1,LZ1) &
    , TRY(LX1,LY1,LZ1) &
    , TRZ(LX1,LY1,LZ1)

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFC)
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
                I = I + 1
                TRX(J1,J2,1) = TR1*AREA(I,1,IFC,IEL)
                TRY(J1,J2,1) = TR2*AREA(I,1,IFC,IEL)
        100 END DO
    ELSE
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                I = I + 1
                TRX(J1,J2,1) = TR1*AREA(I,1,IFC,IEL)
                TRY(J1,J2,1) = TR2*AREA(I,1,IFC,IEL)
                TRZ(J1,J2,1) = TR3*AREA(I,1,IFC,IEL)
        200 END DO
    ENDIF

    RETURN
    END SUBROUTINE TRCON
!-----------------------------------------------------------------------
    SUBROUTINE TRST2D (TRX,TRY,SIGST,IEL,IFC)

!     Compute taction due to surface tension (2D)

    use size_m
    use dxyz
    use geom
    use topol
    use wz_m
    COMMON /CTMP1/ A1X(LX1),A1Y(LX1),STX(LX1),STY(LX1)

    DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),SIGST(LX1,1)
    DIMENSION CANG(2),SANG(2)
    DIMENSION IXN(2),IYN(2),IAN(2)

    DO 100 IX=1,NX1
        AA = SIGST(IX,1) * WXM1(IX)
        STX(IX) = T1X(IX,1,IFC,IEL) * AA
        STY(IX) = T1Y(IX,1,IFC,IEL) * AA
    100 END DO

    IF (IFC == 3 .OR. IFC == 4) THEN
        CALL CHSIGN (STX,NX1)
        CALL CHSIGN (STY,NX1)
    ENDIF

    IF (IFC == 1 .OR. IFC == 3) THEN
        CALL MXM (DXTM1,NX1,STX,NX1,A1X,1)
        CALL MXM (DXTM1,NX1,STY,NX1,A1Y,1)
    ELSE
        CALL MXM (DYTM1,NY1,STX,NY1,A1X,1)
        CALL MXM (DYTM1,NY1,STY,NY1,A1Y,1)
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

    DO 200 J2=JS2,JF2,JSKIP2
        DO 200 J1=JS1,JF1,JSKIP1
            I = I + 1
            TRX(J1,J2,1) = TRX(J1,J2,1) - A1X(I)
            TRY(J1,J2,1) = TRY(J1,J2,1) - A1Y(I)
    200 END DO

!     Contact angle corrections

    CALL CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)
    DO 500 I=1,2
        IX = IXN(I)
        IY = IYN(I)
        IA = IAN(I)
        TRX(IX,IY,1)=TRX(IX,IY,1) + SIGST(IA,1)*CANG(I)
        TRY(IX,IY,1)=TRY(IX,IY,1) + SIGST(IA,1)*SANG(I)
    500 END DO

    RETURN
    END SUBROUTINE TRST2D
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
    SUBROUTINE TRST3D (TRX,TRY,TRZ,SIGST,IEL,IFC)

!     Compute taction due to surface tension (3D)

    use size_m
    use geom
    use wz_m
    COMMON /CTMP0/  XFM1(LX1,LY1),YFM1(LX1,LY1),ZFM1(LX1,LY1)
    COMMON /CTMP1/  DRM1(LX1,LX1),DRTM1(LX1,LY1) &
    ,  DSM1(LX1,LX1),DSTM1(LX1,LY1) &
    ,  WGS(LX1,LY1)
    COMMON /SCRMG/  XRM1(LX1,LY1),YRM1(LX1,LY1),ZRM1(LX1,LY1) &
    ,  XSM1(LX1,LY1),YSM1(LX1,LY1),ZSM1(LX1,LY1)
    COMMON /SCRUZ/  S1X(LX1,LY1),S1Y(LX1,LY1),S1Z(LX1,LY1) &
    ,  S2X(LX1,LY1),S2Y(LX1,LY1),S2Z(LX1,LY1)
    COMMON /SCRNS/  G1X(LX1,LY1),G1Y(LX1,LY1),G1Z(LX1,LY1) &
    ,  G2X(LX1,LY1),G2Y(LX1,LY1),G2Z(LX1,LY1) &
    ,  GBS(LX1,LY1),GB1L(LX1,LY1),GB2L(LX1,LY1)

    DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),TRZ(LX1,LY1,LZ1)
    DIMENSION SIGST(LX1,LY1)

    NXY1 = NX1*NY1

    CALL RZERO3 (S1X,S1Y,S1Z,NXY1)
    CALL RZERO3 (S2X,S2Y,S2Z,NXY1)
    CALL FACEXV (XFM1,YFM1,ZFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL), &
    ZM1(1,1,1,IEL),IFC,0)
    CALL SETDRS (DRM1,DRTM1,DSM1,DSTM1,IFC)

    CALL MXM (DRM1,NX1, XFM1,NX1,XRM1,NY1)
    CALL MXM (DRM1,NX1, YFM1,NX1,YRM1,NY1)
    CALL MXM (DRM1,NX1, ZFM1,NX1,ZRM1,NY1)
    CALL MXM (XFM1,NX1,DSTM1,NY1,XSM1,NY1)
    CALL MXM (YFM1,NX1,DSTM1,NY1,YSM1,NY1)
    CALL MXM (ZFM1,NX1,DSTM1,NY1,ZSM1,NY1)

    DO 100 IX=1,NX1
        DO 100 IY=1,NY1
            GB1X=XRM1(IX,IY)
            GB1Y=YRM1(IX,IY)
            GB1Z=ZRM1(IX,IY)
            GB2X=XSM1(IX,IY)
            GB2Y=YSM1(IX,IY)
            GB2Z=ZSM1(IX,IY)
            GB11=GB1X*GB1X + GB1Y*GB1Y + GB1Z*GB1Z
            GB12=GB1X*GB2X + GB1Y*GB2Y + GB1Z*GB2Z
            GB22=GB2X*GB2X + GB2Y*GB2Y + GB2Z*GB2Z
            GDET=GB11*GB22 - GB12*GB12
            IF (GDET < 1.E-20) GO TO 9001
            GT11= GB22/GDET
            GT12=-GB12/GDET
            GT22= GB11/GDET
            GB1L(IX,IY)=SQRT(GB11)
            GB2L(IX,IY)=SQRT(GB22)
            GBS (IX,IY)=SQRT(GDET)
            WGS (IX,IY)=WXM1(IX)*WYM1(IY)*SIGST(IX,IY)
            BB = GBS(IX,IY) * WGS(IX,IY)
            G1X(IX,IY) = BB * ( GT11*GB1X + GT12*GB2X )
            G1Y(IX,IY) = BB * ( GT11*GB1Y + GT12*GB2Y )
            G1Z(IX,IY) = BB * ( GT11*GB1Z + GT12*GB2Z )
            G2X(IX,IY) = BB * ( GT12*GB1X + GT22*GB2X )
            G2Y(IX,IY) = BB * ( GT12*GB1Y + GT22*GB2Y )
            G2Z(IX,IY) = BB * ( GT12*GB1Z + GT22*GB2Z )
    100 END DO

    CALL MXM (DRTM1,NX1,G1X,NX1,S1X,NY1)
    CALL MXM (DRTM1,NX1,G1Y,NX1,S1Y,NY1)
    CALL MXM (DRTM1,NX1,G1Z,NX1,S1Z,NY1)

    CALL MXM (G2X,NX1,DSM1,NY1,S2X,NY1)
    CALL MXM (G2Y,NX1,DSM1,NY1,S2Y,NY1)
    CALL MXM (G2Z,NX1,DSM1,NY1,S2Z,NY1)

    CALL ADD2 (S1X,S2X,NXY1)
    CALL ADD2 (S1Y,S2Y,NXY1)
    CALL ADD2 (S1Z,S2Z,NXY1)

!     Contact angle option on hold

!      ICONTAC=INT(BC2)
!      IF (ICONTAC.NE.0) THEN
!         IX=1
!         IY=1
!         IF (ICONTAC.GE.3) IY=NY1
!         IF (ICONTAC.EQ.2 .OR. ICONTAC.EQ.3) IX=NX1
!         ANG = BC3 * PI / 180.00
!         RR  = YM1(IX,IY,IZ,IEL)
!         TRX(IX,IY,IZ)=TRX(IX,IY,IZ) + RR*SIGST*COS( ANG )
!         TRY(IX,IY,IZ)=TRY(IX,IY,IZ) + RR*SIGST*SIN( ANG )
!      ENDIF

    CALL FACSUB2 (TRX,TRY,TRZ,S1X,S1Y,S1Z,IFC)

    RETURN

    9001 WRITE ( 6,*) 'Zero area for Element=',IEL,'    Face=',IFC
    call exitt

    END SUBROUTINE TRST3D
!-----------------------------------------------------------------------
    SUBROUTINE SETDRS (DRM1,DRTM1,DSM1,DSTM1,IFC)

    use size_m
    use dxyz

    DIMENSION DRM1(LX1,LX1),DRTM1(LX1,LX1) &
    , DSM1(LY1,LY1),DSTM1(LY1,LY1)

    NXY1=NX1*NY1

    IF (IFC == 5 .OR. IFC == 6) THEN
        CALL COPY (DRM1 ,DXM1 ,NXY1)
        CALL COPY (DSM1 ,DYM1 ,NXY1)
        CALL COPY (DRTM1,DXTM1,NXY1)
        CALL COPY (DSTM1,DYTM1,NXY1)
    ELSEIF (IFC == 2 .OR. IFC == 4) THEN
        CALL COPY (DRM1 ,DYM1 ,NXY1)
        CALL COPY (DSM1 ,DZM1 ,NXY1)
        CALL COPY (DRTM1,DYTM1,NXY1)
        CALL COPY (DSTM1,DZTM1 ,NXY1)
    ELSE
        CALL COPY (DRM1 ,DZM1 ,NXY1)
        CALL COPY (DSM1 ,DXM1 ,NXY1)
        CALL COPY (DRTM1,DZTM1,NXY1)
        CALL COPY (DSTM1,DXTM1,NXY1)
    ENDIF

    RETURN
    END SUBROUTINE SETDRS
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
#if 0
    SUBROUTINE SETSHL

    use size_m
    use input
    use soln
    use tstep
    COMMON /SCRMG/ V1(LX1,LY1,LZ1,LELV) &
    , V2(LX1,LY1,LZ1,LELV) &
    , V3(LX1,LY1,LZ1,LELV) &
    , VV(LX1,LY1,LZ1,LELV)

    common  /nekcb/ cb
    CHARACTER CB*3

    IFIELD = 1
    NFACE  = 2*NDIM
    NTOT1  = NX1*NY1*NZ1*NELV
    DELTA  = 1.E-9
    X      = 1.+DELTA
    Y      = 1.
    DIFF   = ABS(X-Y)
    IF (DIFF == 0.) EPSA = 1.E-06
    IF (DIFF > 0.) EPSA = 1.E-13

    CALL RZERO3  (V1,V2,V3,NTOT1)
    CALL BCTWALL (V1,V2,V3)
    CALL OPDOT   (VV,V1,V2,V3,V1,V2,V3,NTOT1)
    VDOT  = GLMAX(VV,NTOT1)
    VMAX  = SQRT(VDOT)
    IF (VMAX < EPSA) VMAX = -EPSA

    DO 100 IEL=1,NELV
        DO 100 IFC=1,NFACE
            CB=CBC(IFC,IEL,IFIELD)
            IF (CB /= 'V  ' .AND. CB /= 'v  '  .AND. CB /= 'VL ' .AND. &
            CB /= 'vl ') GOTO 100
            IF (VMAX > 0.0) THEN
                CALL CHKZVN (VMAX,IEL,IFC,IVNORL)
                IF (IVNORL == 1) GOTO 100
            ENDIF
            IF (CB == 'V  ') CBC(IFC,IEL,IFIELD)='WS '
            IF (CB == 'VL ') CBC(IFC,IEL,IFIELD)='WSL'
            IF (CB == 'v  ') CBC(IFC,IEL,IFIELD)='ws '
            IF (CB == 'vl ') CBC(IFC,IEL,IFIELD)='wsl'
    100 END DO

    RETURN
    END SUBROUTINE SETSHL
#endif
!-----------------------------------------------------------------------
    SUBROUTINE CHKZVN (VMAX,IEL,IFC,IVNORL)

    use size_m
    use geom
    use soln
    COMMON /SCRMG/ V1(LX1,LY1,LZ1,LELV) &
    , V2(LX1,LY1,LZ1,LELV) &
    , V3(LX1,LY1,LZ1,LELV) &
    , VV(LX1,LY1,LZ1,LELV)

    NXZ1  = NX1*NZ1
    TOLV  = 0.01*VMAX

    VNOR1 = FACDOT(V1(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
    VNOR2 = FACDOT(V2(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
    VNOR  = VNOR1 + VNOR2
    IF (NDIM == 3) THEN
        VNOR3 = FACDOT(V3(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
        VNOR  = VNOR + VNOR3
    ENDIF
    VNOR = ABS(VNOR) / NXZ1

    IVNORL = 1
    IF (VNOR < TOLV) IVNORL = 0

    RETURN
    END SUBROUTINE CHKZVN
!-----------------------------------------------------------------------
    SUBROUTINE BCTWALL (TMP1,TMP2,TMP3)

!     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3)
!     (No antimask operation is applied).

    use size_m
    use geom
    use input
    use tstep

    DIMENSION TMP1(NX1,NY1,NZ1,1) &
    , TMP2(NX1,NY1,NZ1,1) &
    , TMP3(NX1,NY1,NZ1,1)
    common  /nekcb/ cb
    CHARACTER CB*3

    NFACE = 2*NDIM
    NTOT1 = NX1*NY1*NZ1*NELV

    CALL RZERO (TMP1,NTOT1)
    CALL RZERO (TMP2,NTOT1)
    IF (IF3D) CALL RZERO (TMP3,NTOT1)

    DO 2000 IEL=1,NELV
        DO 2000 IFC=1,NFACE
            CB  = CBC (IFC,IEL,IFIELD)
            BC1 = BC(1,IFC,IEL,IFIELD)
            BC2 = BC(2,IFC,IEL,IFIELD)
            BC3 = BC(3,IFC,IEL,IFIELD)
            IF (CB == 'V  ' .OR. CB == 'VL '  .OR. &
            CB == 'WS ' .OR. CB == 'WSL') THEN
                CALL FACEV (TMP1,IEL,IFC,BC1,NX1,NY1,NZ1)
                CALL FACEV (TMP2,IEL,IFC,BC2,NX1,NY1,NZ1)
                IF (NDIM == 3) CALL FACEV (TMP3,IEL,IFC,BC3,NX1,NY1,NZ1)
                IF (CB == 'VL ' .OR. CB == 'WSL') &
                CALL GLOBROT (TMP1(1,1,1,IEL),TMP2(1,1,1,IEL), &
                TMP3(1,1,1,IEL),IEL,IFC)
            ENDIF
            IF (CB == 'v  ' .OR. CB == 'vl ' .OR. &
            CB == 'ws ' .OR. CB == 'wsl' .OR. &
            CB == 'mv ' .OR. CB == 'mvn') THEN
                CALL FACEIV (CB,TMP1(1,1,1,IEL),TMP2(1,1,1,IEL), &
                TMP3(1,1,1,IEL),IEL,IFC,NX1,NY1,NZ1)
                IF (CB == 'vl ' .OR. CB == 'wsl') &
                CALL GLOBROT (TMP1(1,1,1,IEL),TMP2(1,1,1,IEL), &
                TMP3(1,1,1,IEL),IEL,IFC)
            ENDIF
    2000 END DO

    RETURN
    END SUBROUTINE BCTWALL
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
    subroutine check_cyclic  ! check for cyclic bcs
    use size_m
    include 'TOTAL'

    common /scrmg/ v1(lx1,ly1,lz1,lelt) &
    , v2(lx1,ly1,lz1,lelt) &
    , v3(lx1,ly1,lz1,lelt)

    integer :: e,f

    nface = 2*ndim

    n = nx1*ny1*nz1*nelt
    call rzero(v1,n)
    call rzero(v2,n)
    call rzero(v3,n)

    ifield = 1
    do e=1,nelt   ! possibly U or B field
        do f=1,nface

            if (cbc(f,e,ifield) == 'P  ' .OR. cbc(f,e,ifield) == 'p  ') then
                write(*,*) "Ooopps: cyclic"
!                call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,f)
                k = 0
                do j2=js2,jf2,jskip2
                    do j1=js1,jf1,jskip1
                        k = k+1
                        v1(j1,j2,1,e) = unx(j1,j2,1,e)
                        v2(j1,j2,1,e) = uny(j1,j2,1,e)
                        v3(j1,j2,1,e) = unz(j1,j2,1,e)
                    enddo
                enddo
            endif

        enddo
    enddo

    ifcyclic = .FALSE. 
    call opdssum(v1,v2,v3)

    eps = 1.e-4
    if (ndim == 2) call rzero(v3,n)

    do e=1,nelt   ! Check for turning angle
        do f=1,nface

            if (cbc(f,e,ifield) == 'P  ' .OR. cbc(f,e,ifield) == 'p  ') then

                call facindr(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f) ! restricted indx
                snorm = 0.
                dnorm = 0.
                do k=k0,k1
                    do j=j0,j1
                        do i=i0,i1
                            snorm = abs(v1(i,j,k,e)) &
                            + abs(v2(i,j,k,e)) &
                            + abs(v3(i,j,k,e))
                        enddo
                    enddo
                enddo
                if (snorm > eps) ifcyclic = .TRUE. 

            endif

        enddo
    enddo

    itest = 0
    if (ifcyclic) itest = 1
    itest = iglmax(itest,1)

    if (itest > 0) ifcyclic = .TRUE. 

    return
    end subroutine check_cyclic
!-----------------------------------------------------------------------
