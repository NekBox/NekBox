!-----------------------------------------------------------------------
    subroutine cbcmesh

!     Generate boundary conditions (CBC arrays) for mesh solver

    include 'SIZE'
    include 'GEOM'
    include 'INPUT'
    include 'TSTEP'

    CHARACTER CBM*1,CBF*3,CBT*3,CB*3

    IFLD   = 0
    NFACE  = 2*NDIM
    IFMELT = .FALSE. 
    IF (IFTMSH(IFLD)) IFMELT= .TRUE. 

    DO 100 IEL=1,NELT
        DO 100 IFC=1,NFACE
            CBM = CBC(IFC,IEL,0)
            CBF = CBC(IFC,IEL,1)
            CBT = CBC(IFC,IEL,2)
        
            IF (CBT(1:1) == 'M') THEN
                IFLD = 2
                CB   = CBT
                GOTO 200
            ENDIF
            IF (CBF(1:1) == 'M' .OR. CBF(1:1) == 'm') THEN
                IFLD = 1
                CB   = CBF
            !             IF (CBF.EQ.'mv ' .AND. CBT.EQ.'E  ') THEN
            !                 IFTMSH(0)=.TRUE.
            !             ENDIF
                IF (CBF == 'mv ' .AND. CBM == '+'  ) THEN
                    CB = 'mvn'
                    CBC(IFC,IEL,1) = CB
                ENDIF
                GOTO 200
            ENDIF
            IF (CBF(1:1) == 'V' .OR. CBF(1:1) == 'v' .OR. &
            CBF(1:1) == 'W' ) THEN
                IFLD = 1
                CB   = 'FIX'
                IF (IFMELT .OR. CBM == '+') CB='SYM'
                GOTO 200
            ENDIF
            IF (CBT == 'T  ' .OR. CBT == 't  ') THEN
                IFLD = 2
                CB   = 'FIX'
                IF (CBM == '+') CB='SYM'
                GOTO 200
            ENDIF
            IF (CBF == 'P  ' .OR. CBF == 'E  ') THEN
                IFLD = 1
                CB   = CBF
                IF (CBM == '-') CB='FIX'
                IF (CBM == '+') CB='SYM'
                GOTO 200
            ENDIF
            IF (CBT == 'P  ' .OR. CBT == 'E  ') THEN
                IFLD = 2
                CB   = CBT
                IF (CBM == '-') CB='FIX'
                IF (CBM == '+') CB='SYM'
                GOTO 200
            ENDIF
            IFLD = 1
            IF (CBF == '   ') IFLD = 2
            CB   = 'SYM'
            IF (CBM == '-') CB = 'FIX'
        
            200 CBC(IFC,IEL,0) = CB
            DO 250 I=1,5
                BC(I,IFC,IEL,0)=BC(I,IFC,IEL,IFLD)
            250 END DO
    100 END DO

    return
    end subroutine cbcmesh
!-----------------------------------------------------------------------
    subroutine admeshv

    include 'SIZE'
    include 'SOLN'
    include 'TSTEP'

    COMMON /SCRUZ/ FM1(LX1,LY1,LZ1,LELT) &
    , FM2(LX1,LY1,LZ1,LELT) &
    , FM3(LX1,LY1,LZ1,LELT) &
    , PHI(LX1,LY1,LZ1,LELT)

    NTOT1=NX1*NY1*NZ1*NELV

    CALL RZERO (FM1,NTOT1)
    CALL RZERO (FM2,NTOT1)
    CALL RZERO (FM3,NTOT1)

    CALL DIVWS (FM1,VX,PHI,NELV,1)
    CALL DIVWS (FM2,VY,PHI,NELV,2)
    CALL ADD2  (BFX,FM1,NTOT1)
    CALL ADD2  (BFY,FM2,NTOT1)
    IF (NDIM == 3) THEN
        CALL DIVWS (FM3,VZ,PHI,NELV,3)
        CALL ADD2  (BFZ,FM3,NTOT1)
    ENDIF

    return
    end subroutine admeshv
!-----------------------------------------------------------------------
    subroutine admesht

    include 'SIZE'
    include 'SOLN'
    include 'TSTEP'

    COMMON /SCRUZ/ FMT(LX1,LY1,LZ1,LELT) &
    , PHI(LX1,LY1,LZ1,LELT)

    IFLD = 0
    NEL  = NELFLD(IFLD)
    NTOT1= NX1*NY1*NZ1*NEL

    CALL RZERO   (FMT,NTOT1)
    CALL DIVWS   (FMT,T(1,1,1,1,IFIELD-1),PHI,NEL,1)
    CALL ADDCOL3 (BQ(1,1,1,1,IFIELD-1),FMT,VTRANS(1,1,1,1,IFIELD), &
    NTOT1)

    return
    end subroutine admesht
!-----------------------------------------------------------------------
    subroutine divws (fms,sfv,phi,nel,idir)

    include 'SIZE'
    include 'GEOM'
    include 'MASS'
    include 'MVGEOM'
    include 'WZ'
    include 'INPUT'

    COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT) &
    , PHS(LX1,LY1,LZ1,LELT) &
    , PHT(LX1,LY1,LZ1,LELT)

    DIMENSION FMS(LX1,LY1,LZ1,1) &
    , SFV(LX1,LY1,LZ1,1) &
    , PHI(LX1,LY1,LZ1,1)

    NXYZ1 = NX1*NY1*NZ1
    NTOT1 = NXYZ1*NEL

    CALL COL3    (PHI,SFV,WX,NTOT1)
    CALL URST    (PHI,PHR,PHS,PHT,NEL)
    CALL ADDCOL3 (FMS,RXM1,PHR,NTOT1)
    CALL ADDCOL3 (FMS,SXM1,PHS,NTOT1)
    IF (NDIM == 3) CALL ADDCOL3 (FMS,TXM1,PHT,NTOT1)

    CALL COL3    (PHI,SFV,WY,NTOT1)
    CALL URST    (PHI,PHR,PHS,PHT,NEL)
    CALL ADDCOL3 (FMS,RYM1,PHR,NTOT1)
    CALL ADDCOL3 (FMS,SYM1,PHS,NTOT1)
    IF (NDIM == 3) CALL ADDCOL3 (FMS,TYM1,PHT,NTOT1)

    IF (NDIM == 3) THEN
        CALL COL3    (PHI,SFV,WZ,NTOT1)
        CALL URST    (PHI,PHR,PHS,PHT,NEL)
        CALL ADDCOL3 (FMS,RZM1,PHR,NTOT1)
        CALL ADDCOL3 (FMS,SZM1,PHS,NTOT1)
        CALL ADDCOL3 (FMS,TZM1,PHT,NTOT1)
    ENDIF

    CALL COL2    (FMS,BM1,NTOT1)
    CALL INVCOL2 (FMS,JACM1,NTOT1)

    IF (IFAXIS) CALL AXIFMS (FMS,SFV,PHI,NEL,IDIR)

    return
    end subroutine divws
!-----------------------------------------------------------------------
    subroutine axifms (fms,sfv,phi,nel,idir)

    include 'SIZE'
    include 'DXYZ'
    include 'GEOM'
    include 'INPUT'
    include 'MASS'
    include 'MVGEOM'
    include 'WZ'
    COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT) &
    , PHS(LX1,LY1,LZ1,LELT) &
    , PHT(LX1,LY1,LZ1,LELT)

    DIMENSION FMS(LX1,LY1,LZ1,1) &
    , PHI(LX1,LY1,LZ1,1) &
    , SFV(LX1,LY1,LZ1,1) &
    , WYS(LX1)
    EQUIVALENCE (WYS(1),PHT(1,1,1,1))

    NXYZ1 = NX1*NY1*NZ1
    NTOT1 = NXYZ1*NEL
    CALL COL3 (PHI,SFV,WY,NTOT1)

    DO 100 IEL=1,NEL
        IF ( IFRZER(IEL) ) THEN
            IF (IDIR == 1) THEN
                CALL MXM (WY(1,1,1,IEL),NX1,DATM1,NY1,WYS,1)
                DO 220 IX=1,NX1
                    FMS(IX,1,1,IEL)= FMS(IX,1,1,IEL) + WXM1(IX)*WAM1(1)* &
                    WYS(IX)*SFV(IX,1,1,IEL)*JACM1(IX,1,1,IEL)
                220 END DO
            ENDIF
            DO 320 IX=1,NX1
                DO 320 IY=2,NY1
                    FMS(IX,IY,1,IEL)=FMS(IX,IY,1,IEL) + PHI(IX,IY,1,IEL) * &
                    BM1(IX,IY,1,IEL) / YM1(IX,IY,1,IEL)
            320 END DO
        ELSE
            CALL ADDCOL4 (FMS(1,1,1,IEL),PHI(1,1,1,IEL),JACM1(1,1,1,IEL), &
            W2CM1,NXYZ1)
        ENDIF
    100 END DO

    return
    end subroutine axifms
!-----------------------------------------------------------------------
    subroutine updcoor
!-----------------------------------------------------------------------

!     Subroutine to update geometry for moving boundary problems

!-----------------------------------------------------------------------
    include 'SIZE'
    include 'TSTEP'

    IFIELD = 0
    NEL    = NELFLD(IFIELD)

!     Update collocation points coordinates

    CALL UPDXYZ (NEL)

!     Shift lagged mesh velocity

    CALL LAGMSHV (NEL)

    return
    end subroutine updcoor
!-----------------------------------------------------------------------
    subroutine mvbdry (nel)

!     Evaluate mesh velocities at all moving boundaries

    include 'SIZE'
    include 'GEOM'
    include 'INPUT'
    include 'MVGEOM'
    include 'SOLN'
    include 'TSTEP'

    common /scrsf/ wvx(lx1*ly1*lz1,lelt) &
    , wvy(lx1*ly1*lz1,lelt) &
    , wvz(lx1*ly1*lz1,lelt)
    common /scrch/ wtx(lx1*ly1*lz1,lelt) &
    , wty(lx1*ly1*lz1,lelt)
    common /scrmg/ wtz(lx1*ly1*lz1,lelt) &
    , rnx(lx1*ly1*lz1,lelt) &
    , rny(lx1*ly1*lz1,lelt) &
    , rnz(lx1*ly1*lz1,lelt)
    common /scruz/ dsa(lx1*ly1*lz1,lelt) &
    , qni(lx1*ly1*lz1,lelt) &
    , smt(lx1*ly1*lz1,lelt) &
    , ta (lx1*ly1*lz1,lelt)

    logical :: ifalgn,ifnorx,ifnory,ifnorz,ifdsmv,ifregw
    character cb*3
    integer :: e,f

    ifield = 0
    nxyz1  = nx1*ny1*nz1
    n      = nx1*ny1*nz1*nel
    nface  = 2*ndim
    call rzero3  (rnx,rny,rnz,n)

    do 100 e=1,nel
        do 100 f=1,nface
            cb = cbc(f,e,ifield)
            if (cb == 'ms ' .OR. cb == 'MS ' .OR. &
            cb == 'msi' .OR. cb == 'MSI' .OR. &
            cb == 'mm ' .OR. cb == 'MM ' .OR. &
            cb == 'mv ' .OR. cb == 'mvn' .OR. &
            cb == 'MLI') then
                call facexv (unx(1,1,f,e),uny(1,1,f,e), &
                unz(1,1,f,e),rnx(1,e), &
                rny(1,e),rnz(1,e),f,1)
            endif
    100 END DO

    call opdssum (rnx,rny,rnz)
    call unitvec (rnx,rny,rnz,n)

    call rzero3 (wvx,wvy,wvz,n)
    call rzero3 (wtx,wty,wtz,n)
    do 1000 isweep=1,2

        ifregw = .FALSE. 
        ifdsmv = .FALSE. 
        call rzero  (dsa,n)
        call rzero  (ta,n)

        if (ifflow) then
            ifield = 1
            call rzero  (smt,n)
            do 210 e=1,nelv
                do 210 f=1,nface
                    cb = cbc(f,e,ifield)
                    if (cb == 'mv ' .OR. cb == 'mvn'  .OR. &
                    cb == 'mm ' .OR. cb == 'MM '  .OR. &
                    cb == 'msi' .OR. cb == 'MSI'  .OR. &
                    cb == 'ms ' .OR. cb == 'MS ') then
                        ifregw = .TRUE. 
                        call facec3 (wvx(1,e),wvy(1,e),wvz(1,e), &
                        vx (1,1,1,e),vy (1,1,1,e),vz (1,1,1,e),f)
                        if (cb /= 'mv ') &
                        call norcmp2(wvx(1,e),wvy(1,e),wvz(1,e),e,f)
                    endif
                !          if (cb.eq.'msi' .or. cb.eq.'MSI') then
                !             ifdsmv = .true.
                !             call facsmt (smt(1,e),f)
                !          endif
            210 END DO

            call dsavg(wvx)
            call dsavg(wvy)
            if (if3d) call dsavg(wvz)

            if (istep == 0) call opcopy(wx,wy,wz,wvx,wvy,wvz)
        !         if (istep.eq.4) then
        !            call opcopy(wx,wy,wz,wvx,wvy,wvz)
        !            call outpost(vx,vy,vz,wtx,wtx,'   ')
        !            call outpost(wx,wy,wz,wtx,wtx,'   ')
        !            write(6,*) 'quit1',istep
        !            stop
        !         endif

            iregw = 0    ! Global handshake on logicals ifregw, ifdsmv
            if (ifregw) iregw = 1       ! pff  9/11/07
            iregw = iglmax(iregw,1)
            if (iregw == 1) ifregw = .TRUE. 

            idsmv = 0
            if (ifdsmv) idsmv = 1
            idsmv = iglmax(idsmv,1)
            if (idsmv == 1) ifdsmv = .TRUE. 

            ifdsmv = .FALSE. 
            ifregw = .TRUE. 

        !         if (ifdsmv) then
        !           call dssum (smt,nx1,ny1,nz1)
        !           do 215 e=1,nelv
        !           do 215 f=1,nface
        !              cb = cbc(f,e,ifield)
        !              if (cb.eq.'msi' .or. cb.eq.'MSI') then
        !               call facec3 (wtx(1,e),wty(1,e),wtz(1,e),
        !    $                       vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),f)
        !               call facemv (wtx(1,e),wty(1,e),wtz(1,e),
        !    $                       rnx(1,e),rny(1,e),rnz(1,e),
        !    $                       smt(1,e),f)
        !              endif
        ! 215       continue
        !           call opdssum(wtx,wty,wtz)
        !         endif

        endif
    
        if (ifmelt .AND. istep > 0) then
            ifield = 2
            call rzero (smt,n)
            call cqnet (qni,ta,nel)
            do 220 e=1,nelt
                do 220 f=1,nface
                    cb   = cbc(f,e,ifield)
                    if (cb == 'MLI') call facsmt (smt(1,e),f)
                    if (cb == 'MLI' .OR. cb == 'MCI') then
                        call facexs (area(1,1,f,e),ta,f,1)
                        call add2   (dsa(1,e),ta,nxyz1)
                    endif
            220 END DO
            call dssum (smt,nx1,ny1,nz1)
            call dssum (dsa,nx1,ny1,nz1)
            do 280 e=1,nelt
                do 280 f=1,nface
                    cb = cbc(f,e,ifield)
                    if (cb == 'MLI') then
                        rhola = -0.5 * bc(5,f,e,ifield)
                        call facemt (wtx(1,e),wty(1,e),wtz(1,e), &
                        rnx(1,e),rny(1,e),rnz(1,e), &
                        qni(1,e),dsa(1,e),smt(1,e), &
                        rhola,f)
                    endif
            280 END DO
            call opdssum (wtx,wty,wtz)
        endif

        ifield = 0
        do 330 e=1,nel
            do 330 f=1,nface
                cb = cbc(f,e,ifield)
                if (cb == 'SYM') then
                    call chknord (ifalgn,ifnorx,ifnory,ifnorz,f,e)
                    if (ifregw) then
                        if (ifnorx) call facev (wvx,e,f,0.0,nx1,ny1,nz1)
                        if (ifnory) call facev (wvy,e,f,0.0,nx1,ny1,nz1)
                        if (ifnorz) call facev (wvz,e,f,0.0,nx1,ny1,nz1)
                        if ( .NOT. ifalgn) call faczqn (wvx(1,e),wvy(1,e), &
                        wvz(1,e),f,e)
                    endif
                    if (ifdsmv .OR. ifmelt) then
                        if (ifnorx) call facev (wtx,e,f,0.0,nx1,ny1,nz1)
                        if (ifnory) call facev (wty,e,f,0.0,nx1,ny1,nz1)
                        if (ifnorz) call facev (wtz,e,f,0.0,nx1,ny1,nz1)
                        if ( .NOT. ifalgn) call faczqn (wtx(1,e),wty(1,e), &
                        wtz(1,e),f,e)
                    endif
                endif
        330 END DO

        do 350 e=1,nel
            do 350 f=1,nface
                cb = cbc(f,e,ifield)
                if (cb == 'FIX') then
                    if (ifregw) then
                        call facev (wvx,e,f,0.0,nx1,ny1,nz1)
                        call facev (wvy,e,f,0.0,nx1,ny1,nz1)
                        if (ndim == 3) call facev (wvz,e,f,0.0,nx1,ny1,nz1)
                    endif
                    if (ifdsmv .OR. ifmelt) then
                        call facev (wtx,e,f,0.0,nx1,ny1,nz1)
                        call facev (wty,e,f,0.0,nx1,ny1,nz1)
                        if (ndim == 3) call facev (wtz,e,f,0.0,nx1,ny1,nz1)
                    endif
                endif
        350 END DO

        if (isweep == 1) then
            if (ifregw) then
                call dsop (wvx,'MXA',nx1,ny1,nz1)
                call dsop (wvy,'MXA',nx1,ny1,nz1)
                if (ndim == 3) call dsop (wvz,'MXA',nx1,ny1,nz1)
            endif
            if (ifdsmv .OR. ifmelt) then
                call dsop (wtx,'MXA',nx1,ny1,nz1)
                call dsop (wty,'MXA',nx1,ny1,nz1)
                if (ndim == 3) call dsop (wtz,'MXA',nx1,ny1,nz1)
            endif
        else
            if (ifregw) then
                call dsop (wvx,'MNA',nx1,ny1,nz1)
                call dsop (wvy,'MNA',nx1,ny1,nz1)
                if (ndim == 3) call dsop (wvz,'MNA',nx1,ny1,nz1)
            endif
            if (ifdsmv .OR. ifmelt) then
                call dsop (wtx,'MNA',nx1,ny1,nz1)
                call dsop (wty,'MNA',nx1,ny1,nz1)
                if (ndim == 3) call dsop (wtz,'MNA',nx1,ny1,nz1)
            endif
        endif

    1000 END DO

    call rmask (wx,wy,wz,nel)
!     if (istep.eq.2) then
!        call outpost(wx,wy,wz,wtx,wtx,'   ')
!        write(6,*) 'quit2'
!        stop
!     endif

    if (ifregw) then
        call add2  (wx,wvx,n)
        call add2  (wy,wvy,n)
        if (ndim == 3) call add2  (wz,wvz,n)
    endif
    if (ifdsmv .OR. ifmelt) then
        call add2  (wx,wtx,n)
        call add2  (wy,wty,n)
        if (ndim == 3) call add2  (wz,wtz,n)
    endif
     
!     if (istep.gt.4) then
!        call outpost(wx,wy,wz,wtx,wtx,'   ')
!        write(6,*) 'quit3'
!        stop
!     endif
     
    return
    end subroutine mvbdry
!-----------------------------------------------------------------------
    subroutine norcmp2(wvx,wvy,wvz,e,f)
    include 'SIZE'
    include 'GEOM'
    include 'INPUT'


    real :: wvx(lx1,ly1,lz1),wvy(lx1,ly1,lz1),wvz(lx1,ly1,lz1)

    integer :: e,f

    common /scruz/ r1(lx1,ly1,lz1),r2(lx1,ly1,lz1),r3(lx1,ly1,lz1)

    call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)

    l=0
    do k=k0,k1
        do j=j0,j1
            do i=i0,i1
                l=l+1
                scale=wvx(i,j,k)*unx(l,1,f,e) &
                +wvy(i,j,k)*uny(l,1,f,e) &
                +wvz(i,j,k)*unz(l,1,f,e)
                wvx(i,j,k) = scale*unx(l,1,f,e)
                wvy(i,j,k) = scale*uny(l,1,f,e)
                wvz(i,j,k) = scale*unz(l,1,f,e)
            enddo
        enddo
    enddo

!     wvxm = vlamax(wvx,nx1*ny1)
!     write(6,*) f,e,wvxm,' w-max'

    return
    end subroutine norcmp2
!-----------------------------------------------------------------------
    subroutine norcmp (wt1,wt2,wt3,rnx,rny,rnz,ifc)

    include 'SIZE'
    COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)

    DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1) &
    , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1)

    NXYZ1 = NX1*NY1*NZ1

    CALL COPY (R1,WT1,NXYZ1)
    CALL COPY (R2,WT2,NXYZ1)
    IF (NDIM == 3) CALL COPY (R3,WT3,NXYZ1)
    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    IF (NDIM == 2) THEN
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                WN          = R1(J1,J2,1)*RNX(J1,J2,1) + &
                R2(J1,J2,1)*RNY(J1,J2,1)
                WT1(J1,J2,1) = WN *RNX(J1,J2,1)
                WT2(J1,J2,1) = WN *RNY(J1,J2,1)
        200 END DO
    ELSE
        DO 300 J2=JS2,JF2,JSKIP2
            DO 300 J1=JS1,JF1,JSKIP1
                WN          = R1(J1,J2,1)*RNX(J1,J2,1) + &
                R2(J1,J2,1)*RNY(J1,J2,1) + &
                R3(J1,J2,1)*RNZ(J1,J2,1)
                WT1(J1,J2,1) = WN *RNX(J1,J2,1)
                WT2(J1,J2,1) = WN *RNY(J1,J2,1)
                WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
        300 END DO
    ENDIF

    return
    end subroutine norcmp
!-----------------------------------------------------------------------
    subroutine facemv (wt1,wt2,wt3,rnx,rny,rnz,smt,ifc)

    include 'SIZE'
    COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)

    DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1) &
    , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1) &
    , SMT(LX1,LY1,LZ1)

    NXYZ1 = NX1*NY1*NZ1

    CALL COPY (R1,WT1,NXYZ1)
    CALL COPY (R2,WT2,NXYZ1)
    IF (NDIM == 3) CALL COPY (R3,WT3,NXYZ1)
    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    IF (NDIM == 2) THEN
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + &
                R2(J1,J2,1)*RNY(J1,J2,1) ) / SMT(J1,J2,1)
                WT1(J1,J2,1) = WN *RNX(J1,J2,1)
                WT2(J1,J2,1) = WN *RNY(J1,J2,1)
        200 END DO
    ELSE
        DO 300 J2=JS2,JF2,JSKIP2
            DO 300 J1=JS1,JF1,JSKIP1
                WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + &
                R2(J1,J2,1)*RNY(J1,J2,1) + &
                R3(J1,J2,1)*RNZ(J1,J2,1) ) / SMT(J1,J2,1)
                WT1(J1,J2,1) = WN *RNX(J1,J2,1)
                WT2(J1,J2,1) = WN *RNY(J1,J2,1)
                WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
        300 END DO
    ENDIF

    return
    end subroutine facemv
!-----------------------------------------------------------------------
    subroutine faczqn (wt1,wt2,wt3,ifc,iel)

    include 'SIZE'
    include 'GEOM'
    include 'TOPOL'
    COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)

    DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)

    NXYZ1 = NX1*NY1*NZ1
    CALL COPY (R1,WT1,NXYZ1)
    CALL COPY (R2,WT2,NXYZ1)
    IF (NDIM == 3) CALL COPY (R3,WT3,NXYZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
    I = 0

    IF (NDIM == 2) THEN
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                I = I+1
                W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + &
                R2(J1,J2,1)*T1Y(I,1,IFC,IEL)
                WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
                WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
        200 END DO
    ELSE
        DO 300 J2=JS2,JF2,JSKIP2
            DO 300 J1=JS1,JF1,JSKIP1
                I = I+1
                W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + &
                R2(J1,J2,1)*T1Y(I,1,IFC,IEL) + &
                R3(J1,J2,1)*T1Z(I,1,IFC,IEL)
                WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
                WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
                WT3(J1,J2,1) = W1 *T1Z(I,1,IFC,IEL)
        300 END DO
    ENDIF

    return
    end subroutine faczqn
!-----------------------------------------------------------------------
    subroutine facsmt (smt,ifc)

    include 'SIZE'
    DIMENSION SMT(LX1,LY1,LZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            SMT(J1,J2,1)=SMT(J1,J2,1) + 1.0
    100 END DO

    return
    end subroutine facsmt
!-----------------------------------------------------------------------
    subroutine cqnet (qni,ta,nel)

    include 'SIZE'
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'
    COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT) &
    , H2(LX1,LY1,LZ1,LELT)

    DIMENSION QNI(LX1,LY1,LZ1,1) &
    , TA (LX1,LY1,LZ1,1)

    INTLOC = -1
    IMSHL  =  2
    NTOT1  = NX1*NY1*NZ1*NEL

    CALL SETHLM (H1,H2,INTLOC)
    CALL AXHELM (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMSHL,1)
    CALL SUB3   (QNI,TA,BQ(1,1,1,1,IFIELD-1),NTOT1)
    CALL DSSUM  (QNI,NX1,NY1,NZ1)

    return
    end subroutine cqnet
!-----------------------------------------------------------------------
    subroutine facemt (w1,w2,w3,rnx,rny,rnz,qni,dsa,smt,rhola,ifc)

    include 'SIZE'
    include 'GEOM'

    DIMENSION  W1 (LX1,LY1,LZ1) &
    ,  W2 (LX1,LY1,LZ1) &
    ,  W3 (LX1,LY1,LZ1) &
    ,  RNX(LX1,LY1,LZ1) &
    ,  RNY(LX1,LY1,LZ1) &
    ,  RNZ(LX1,LY1,LZ1) &
    ,  QNI(LX1,LY1,LZ1) &
    ,  DSA(LX1,LY1,LZ1) &
    ,  SMT(LX1,LY1,LZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    IF (NDIM == 2) THEN
        DO 200 J2=JS2,JF2,JSKIP2
            DO 200 J1=JS1,JF1,JSKIP1
                AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
                W1(J1,J2,1) = RNX(J1,J2,1) * AA
                W2(J1,J2,1) = RNY(J1,J2,1) * AA
        200 END DO
    ELSE
        DO 300 J2=JS2,JF2,JSKIP2
            DO 300 J1=JS1,JF1,JSKIP1
                AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
                W1(J1,J2,1) = RNX(J1,J2,1) * AA
                W2(J1,J2,1) = RNY(J1,J2,1) * AA
                W3(J1,J2,1) = RNZ(J1,J2,1) * AA
        300 END DO
    ENDIF

    return
    end subroutine facemt
!-----------------------------------------------------------------------
    subroutine elasolv (nel)

!     Elastostatic solver for mesh deformation

    include 'SIZE'
    include 'GEOM'
    include 'INPUT'
    include 'MVGEOM'
    include 'SOLN'
    include 'TSTEP'

    COMMON /SCRNS/ DW1  (LX1,LY1,LZ1,LELT) &
    , DW2  (LX1,LY1,LZ1,LELT) &
    , DW3  (LX1,LY1,LZ1,LELT) &
    , AW1  (LX1,LY1,LZ1,LELT) &
    , AW2  (LX1,LY1,LZ1,LELT) &
    , AW3  (LX1,LY1,LZ1,LELT)
    COMMON /SCRVH/ H1   (LX1,LY1,LZ1,LELT) &
    , H2   (LX1,LY1,LZ1,LELT)
    common /scruz/ prt  (lx1,ly1,lz1,lelt)
    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV

!     Set up parameters for mesh solver

!     return
!     if (istep.gt.1) then
!        ifldx = ifield
!        ifield = 1
!        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
!        ifield = ifldx
!        return
!     endif

!     call quickmv

    if (ifusermv) return  ! Compute wx,wy,wz in userchk.

    NTOT1  = NX1*NY1*NZ1*NEL
    MAXIT  = 1000
    MATMOD = -1
    IFH2   = .FALSE. 
    IFSOLV = .TRUE. 
    IMSOLV = 0
    VNU    = 0.0
    VNU    = param(47)
    if (vnu == 0) VNU    = 0.4
    vnu    = max(0.00,vnu)
    vnu    = min(0.499,vnu)

!     Set up elastic material constants

    CE = 1./(1. + VNU)
    C2 = VNU * CE / (1. - 2.*VNU)
    C3 = 0.5 * CE
    CALL CFILL (H1,C2,NTOT1)
    CALL CFILL (H2,C3,NTOT1)

!     Solve for interior mesh velocities

    CALL MESHTOL (AW1,TOLMSH,NEL,IMSOLV)
    IF (IMSOLV == 1) return


    CALL AXHMSF  (AW1,AW2,AW3,WX,WY,WZ,H1,H2,MATMOD)

!     if (istep.eq.2) then
!        call outpost(wx,wy,wz,h1,h2,'   ')
!        call outpost(aw1,aw2,aw3,h1,h2,'   ')
!        write(6,*) 'quit elas 1'
!        stop
!     endif

    CALL CHSIGN  (AW1,NTOT1)
    CALL CHSIGN  (AW2,NTOT1)
    IF (NDIM == 3) CALL CHSIGN (AW3,NTOT1)
    CALL HMHZSF  ('NOMG',DW1,DW2,DW3,AW1,AW2,AW3,H1,H2, &
    W1MASK,W2MASK,W3MASK,WMULT,TOLMSH, &
    MAXIT,MATMOD)

!     Update mesh velocities

    CALL ADD2 (WX,DW1,NTOT1)
    CALL ADD2 (WY,DW2,NTOT1)
    IF (NDIM == 3) CALL ADD2 (WZ,DW3,NTOT1)

    ifldt = ifield
    ifield=1
    if (ifheat) ifield=2
    call dsavg(wx)
    call dsavg(wy)
    call dsavg(wz)
    ifield = ifldt

!     if (istep.gt.1) then
!        ifldx = ifield
!        ifield = 1
!        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
!        ifield = ifldx
!        return
!     endif

    return
    end subroutine elasolv
!-----------------------------------------------------------------------
    subroutine meshtol (ta,tolmsh,nel,imsolv)

    include 'SIZE'
    include 'EIGEN'
    include 'MVGEOM'
    include 'TSTEP'
    DIMENSION TA(LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NEL
    TOLAB = TOLREL

    DELTA  = 1.0E-9
    X      = 1.0 + DELTA
    Y      = 1.0
    DIFF   = ABS(X - Y)
    IF (DIFF == 0.0) EPS = 1.0E-05
    IF (DIFF > 0.0) EPS = 1.0E-12

    CALL OPDOT  (TA,WX,WY,WZ,WX,WY,WZ,NTOT1)

    WDOT = GLMAX(TA,NTOT1)
    WMAX = SQRT(WDOT)
    IF (WMAX < EPS) THEN
        IMSOLV = 1
        return
    ENDIF

    TOLMSH = TOLAB * WMAX * SQRT(EIGAA)

    return
    end subroutine meshtol
!-----------------------------------------------------------------------
    subroutine updxyz (nel)

    include 'SIZE'
    include 'TSTEP'
    include 'MVGEOM'
    include 'GEOM'
    COMMON /SCRSF/ UX(LX1,LY1,LZ1,LELT) &
    , UY(LX1,LY1,LZ1,LELT) &
    , UZ(LX1,LY1,LZ1,LELT)
    DIMENSION ABM(3)

    NTOT1 = NX1*NY1*NZ1*NEL

    DO 10 I=1,NBD
        ABM(I) = DT*ABMSH(I)
    10 END DO

    IF (ISTEP == 0) THEN
        CALL COPY (UX,WX,NTOT1)
        CALL COPY (UY,WY,NTOT1)
        IF (NDIM == 3) CALL COPY (UZ,WZ,NTOT1)
    ELSE
        CALL CMULT2 (UX,WX,ABM(1),NTOT1)
        CALL CMULT2 (UY,WY,ABM(1),NTOT1)
        IF (NDIM == 3) CALL CMULT2 (UZ,WZ,ABM(1),NTOT1)
        DO 100 ILAG=2,NBD
            CALL ADD2S2 (UX,WXLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            CALL ADD2S2 (UY,WYLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            IF (NDIM == 3) &
            CALL ADD2S2 (UZ,WZLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
        100 END DO
    ENDIF

    CALL ADD2 (XM1,UX,NTOT1)
    CALL ADD2 (YM1,UY,NTOT1)
    IF (NDIM == 3) CALL ADD2 (ZM1,UZ,NTOT1)

    return
    end subroutine updxyz
!-----------------------------------------------------------------------
    subroutine lagmshv (nel)
!-----------------------------------------------------------------------

!     Keep old mesh velocity

!-----------------------------------------------------------------------
    include 'SIZE'
    include 'INPUT'
    include 'MVGEOM'
    include 'TSTEP'

    NTOT1 = NX1*NY1*NZ1*NEL

    DO 100 ILAG=NBDINP-1,2,-1
        CALL COPY (WXLAG(1,1,1,1,ILAG),WXLAG(1,1,1,1,ILAG-1),NTOT1)
        CALL COPY (WYLAG(1,1,1,1,ILAG),WYLAG(1,1,1,1,ILAG-1),NTOT1)
        IF (NDIM == 3) &
        CALL COPY (WZLAG(1,1,1,1,ILAG),WZLAG(1,1,1,1,ILAG-1),NTOT1)
    100 END DO

    CALL COPY (WXLAG(1,1,1,1,1),WX,NTOT1)
    CALL COPY (WYLAG(1,1,1,1,1),WY,NTOT1)
    IF (NDIM == 3) &
    CALL COPY (WZLAG(1,1,1,1,1),WZ,NTOT1)

    return
    end subroutine lagmshv
!-----------------------------------------------------------------------
    subroutine facec3 (a1,a2,a3,b1,b2,b3,ifc)

!     Copy the face (IFC) of B1,B2,B3 to A1,A2,A3.
!     IFACE is the input in the pre-processor ordering scheme.

    include 'SIZE'
    DIMENSION A1(LX1,LY1,LZ1) &
    , A2(LX1,LY1,LZ1) &
    , A3(LX1,LY1,LZ1) &
    , B1(LX1,LY1,LZ1) &
    , B2(LX1,LY1,LZ1) &
    , B3(LX1,LY1,LZ1)

    CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)

    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            A1(J1,J2,1)=B1(J1,J2,1)
            A2(J1,J2,1)=B2(J1,J2,1)
            A3(J1,J2,1)=B3(J1,J2,1)
    100 END DO
    return
    end subroutine facec3
!-----------------------------------------------------------------------
    subroutine ptbgeom
!-----------------------------------------------------------------------

!     Subroutine to impose perturbation to geometry before solution
!     for moving boundary problems

!-----------------------------------------------------------------------
    include 'SIZE'
    include 'GEOM'
    include 'MVGEOM'
    include 'SOLN'
    include 'TSTEP'
    include 'INPUT'
    COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT) &
    ,             YM3 (LX3,LY3,LZ3,LELT) &
    ,             ZM3 (LX3,LY3,LZ3,LELT)

    IF (ISTEP == 0) return
    IFIELD = 0
    NEL    = NELFLD(IFIELD)
    NTOT1  = NX1*NY1*NZ1*NEL
    IMESH  = 1
    IF ( IFTMSH(IFIELD) ) IMESH = 2

    CALL IBDGEOM (NEL)
    CALL ELASOLV (NEL)
    CALL UPDXYZ  (NEL)
    CALL GEOM1 (XM3,YM3,ZM3)
    CALL GEOM2
    CALL UPDMSYS (0)
    CALL VOLUME
    CALL SETINVM
    CALL LAGMASS

    CALL RZERO (WX,NTOT1)
    CALL RZERO (WY,NTOT1)
    IF (NDIM == 3) CALL RZERO (WZ,NTOT1)

    return
    end subroutine ptbgeom
!-----------------------------------------------------------------------
    subroutine ibdgeom (nel)

!     Routine to evaluate mesh velocities at all moving boundaries

    include 'SIZE'
    include 'GEOM'
    include 'INPUT'
    include 'PARALLEL'
    include 'MVGEOM'
    include 'TSTEP'

    CHARACTER CB*1

    NFACE  = 2*NDIM
    NTOT1  = NX1*NY1*NZ1*NEL

    CALL RZERO (WX,NTOT1)
    CALL RZERO (WY,NTOT1)
    CALL RZERO (WZ,NTOT1)

    DO 1000 ISWEEP=1,2
    
        IFLD = 0
        DO 110 IEL=1,NEL
            DO 110 IFC=1,NFACE
                ieg = lglel(iel)
                CB  = CBC(IFC,IEL,IFLD)
                IF (CB == 'M' .OR. CB == 'm') THEN
                    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFC)
                    DO 140 IZ=KZ1,KZ2
                        DO 140 IY=KY1,KY2
                            DO 140 IX=KX1,KX2
                                CALL INIGEOM (WX(IX,IY,IZ,IEL),WY(IX,IY,IZ,IEL), &
                                WZ(IX,IY,IZ,IEL),XM1(IX,IY,IZ,IEL), &
                                YM1(IX,IY,IZ,IEL),ZM1(IX,IY,IZ,IEL), &
                                IFC,IEG)
                    140 END DO
                ENDIF
        110 END DO
    
        IF (ISWEEP == 1) THEN
            CALL DSOP (WX,'MXA',NX1,NY1,NZ1)
            CALL DSOP (WY,'MXA',NX1,NY1,NZ1)
            IF (NDIM == 3) CALL DSOP (WZ,'MXA',NX1,NY1,NZ1)
        ELSE
            CALL DSOP (WX,'MNA',NX1,NY1,NZ1)
            CALL DSOP (WY,'MNA',NX1,NY1,NZ1)
            IF (NDIM == 3) CALL DSOP (WZ,'MNA',NX1,NY1,NZ1)
        ENDIF
    
    1000 END DO

    return
    end subroutine ibdgeom
!-----------------------------------------------------------------------
    subroutine inigeom (ux,uy,uz,x,y,z,iside,iel)

    include 'SIZE'
    include 'TSTEP'

    UX  =  0.0
    UY  =  0.0
    UZ  =  0.0

    return
    end subroutine inigeom
!-----------------------------------------------------------------------
    subroutine quickmv
    include 'SIZE'
    include 'TOTAL'
    include 'ZPER'

    if (if3d) then
        call quickmv3d
    else
        call quickmv2d
    endif
    return
    end subroutine quickmv
!-----------------------------------------------------------------------
    subroutine quickmv2d
    include 'SIZE'
    include 'TOTAL'
    include 'ZPER'

    integer :: e,ex,ey,ez,eg
    common /surfa/ zsurf(lx1,lz1,lelx,lely) &
    , wsurf(lx1,lz1,lelx,lely)
    real :: nxs,nys,nzs

    icount = 0
    do ex=1,nelx
        do ix=1,nx1
            zsurf(ix,1,ex,1) = -1.e20
            wsurf(ix,1,ex,1) = -1.e20
            ey=nely
            eg  = ex + nelx*(ey-1)
            mid = gllnid(eg)
            e   = gllel (eg)
            if (mid == nid) then
                zsurf(ix,1,ex,1) = ym1(ix,ny1,1,e)
                vxs              = vx (ix,ny1,1,e)
                vys              = vy (ix,ny1,1,e)
                nxs              = unx(ix,1,3,e)          ! Face 3 is on top in 2D
                nys              = uny(ix,1,3,e)
                gamma_s          = (vxs*nxs + vys*nys)/(nys)
                wsurf(ix,1,ex,1) = gamma_s           ! vertical component of wsurf
            endif
            zsurf(ix,1,ex,1) = glmax(zsurf(ix,1,ex,1),1)
            wsurf(ix,1,ex,1) = glmax(wsurf(ix,1,ex,1),1)
            icount = icount+1

        !        write(6,6) ex,e,ix,xm1(ix,ny1,1,e),ym1(ix,ny1,1,e)
        !    $   ,vxs,vys,nxs,nys,gamma_s,wsurf(ix,1,ex,1),zsurf(ix,1,ex,1)
        !   6        format(3i3,1p9e12.4,' srf')

        enddo
    enddo
    zmin = glmin(ym1,nx1*ny1*nz1*nelv)

    do ex=1,nelx
        do ix=1,nx1
            do ey=1,nely
                eg  = ex + nelx*(ey-1)
                mid = gllnid(eg)
                e   = gllel (eg)
                if (mid == nid) then
                    do iy=1,ny1
                        wy (ix,iy,1,e) = wsurf(ix,1,ex,1) &
                        * (ym1(ix,iy,1,e)-zmin)/(zsurf(ix,1,ex,1)-zmin)
                    enddo
                endif
            enddo
        enddo
    enddo

    n = nx1*ny1*nz1*nelv
    call rzero(wx,n)
    call dsavg(wy)

!     call opcopy(vx,vy,vz,wx,wy,wz)
!     call outpost(vx,vy,vz,pr,t,'   ')
!     call exitt

    return
    end subroutine quickmv2d
!-----------------------------------------------------------------------
    subroutine quickmv3d
    include 'SIZE'
    include 'TOTAL'
    include 'ZPER'

    integer :: e,ex,ey,ez,eg
    common /surfa/ zsurf(lx1,lz1,lelx,lely) &
    , wsurf(lx1,lz1,lelx,lely)
    real :: nxs,nys,nzs

    icount = 0
    do ey=1,nely
        do ex=1,nelx
            do iy=1,ny1
                do ix=1,nx1
                    zsurf(ix,iy,ex,ey) = -1.e20
                    wsurf(ix,iy,ex,ey) = -1.e20
                    ez  = nelz
                    eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
                    mid = gllnid(eg)
                    e   = gllel (eg)
                    if (mid == nid) then
                        zsurf(ix,iy,ex,ey) = zm1(ix,iy,nz1,e)
                        vxs                = vx (ix,iy,nz1,e)
                        vys                = vy (ix,iy,nz1,e)
                        vzs                = vz (ix,iy,nz1,e)
                        nxs                = unx(ix,iy,6,e)     ! Face 6 is on top in 3D
                        nys                = uny(ix,iy,6,e)
                        nzs                = unz(ix,iy,6,e)
                        gamma_s            = (vxs*nxs+vys*nys+vzs*nzs)/(nzs)
                        wsurf(ix,iy,ex,ey) = gamma_s  ! vertical component of wsurf
                    endif
                    zsurf(ix,iy,ex,ey) = glmax(zsurf(ix,iy,ex,ey),1)
                    wsurf(ix,iy,ex,ey) = glmax(wsurf(ix,iy,ex,ey),1)
                    icount = icount+1
                enddo
            enddo
        enddo
    enddo

    n = nx1*ny1*nz1*nelv
    zmin = glmin(zm1,n)

    do ey=1,nely
        do ex=1,nelx
            do iy=1,ny1
                do ix=1,nx1
                    do ez=1,nelz
                        eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
                        mid = gllnid(eg)
                        e   = gllel (eg)
                        if (mid == nid) then
                            do iz=1,nz1
                                wz (ix,iy,iz,e) = wsurf(ix,iy,ex,ey) &
                                * (zm1(ix,iy,iz,e)-zmin)/(zsurf(ix,iy,ex,ey)-zmin)
                            enddo
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo

    n = nx1*ny1*nz1*nelv
    call rzero(wx,n)
    call rzero(wy,n)
    call dsavg(wz)

    return
    end subroutine quickmv3d
!-----------------------------------------------------------------------
