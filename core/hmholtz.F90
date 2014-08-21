!=======================================================================
    subroutine hmholtz(name,u,rhs,h1,h2,mask,mult,imsh,tli,maxit,isd)
    use ctimer
    use size_m
    use fdmh1
    use input
    use mass
    use soln
    use tstep

    CHARACTER      NAME*4
    REAL ::           U    (LX1,LY1,LZ1,1)
    REAL ::           RHS  (LX1,LY1,LZ1,1)
    REAL ::           H1   (LX1,LY1,LZ1,1)
    REAL ::           H2   (LX1,LY1,LZ1,1)
    REAL ::           MASK (LX1,LY1,LZ1,1)
    REAL ::           MULT (LX1,LY1,LZ1,1)

    logical :: iffdm
    character(3) :: nam3

    tol = abs(tli)

    if (icalld == 0) thmhz=0.0

    iffdm = .FALSE. 
!     iffdm = .true.
    if (ifsplit) iffdm = .TRUE. 

    if (icalld == 0 .AND. iffdm) call set_fdm_prec_h1A

    icalld=icalld+1
    nhmhz=icalld
    etime1=dnekclock()
    ntot = nx1*ny1*nz1*nelfld(ifield)
    if (imsh == 1) ntot = nx1*ny1*nz1*nelv
    if (imsh == 2) ntot = nx1*ny1*nz1*nelt

!     Determine which field is being computed for FDM based preconditioner bc's

    call chcopy(nam3,name,3)

    kfldfdm = -1
!     if (nam3.eq.'TEM' ) kfldfdm =  0
!     if (name.eq.'TEM1') kfldfdm =  0  ! hardcode for temp only, for mpaul
!     if (name.eq.'VELX') kfldfdm =  1
!     if (name.eq.'VELY') kfldfdm =  2
!     if (name.eq.'VELZ') kfldfdm =  3
    if (name == 'PRES') kfldfdm =  ndim+1
!     if (.not.iffdm) kfldfdm=-1

    call dssum   (rhs,nx1,ny1,nz1)
    call col2    (rhs,mask,ntot)
    if (nid == 0 .AND. istep <= 10) &
    write(6,*) param(22),' p22 ',istep,imsh
    if (param(22) == 0 .OR. istep <= 10) &
    call chktcg1 (tol,rhs,h1,h2,mask,mult,imsh,isd)

    if (tli < 0) tol=tli ! caller-specified relative tolerance

    if (imsh == 1) call cggo &
    (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1,name)
    if (imsh == 2) call cggo &
    (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,bintm1,name)


    thmhz=thmhz+(dnekclock()-etime1)
    return
    end subroutine hmholtz

!=======================================================================
    subroutine axhelm (au,u,helm1,helm2,imesh,isd)
!------------------------------------------------------------------

!     Compute the (Helmholtz) matrix-vector product,
!     AU = helm1*[A]u + helm2*[B]u, for NEL elements.

!------------------------------------------------------------------
    use ctimer
    use size_m
    use dxyz
    use geom
    use input
    use mass
    use parallel
    use wz_m

    COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV

    REAL ::           AU    (LX1,LY1,LZ1,1) &
    ,             U     (LX1,LY1,LZ1,1) &
    ,             HELM1 (LX1,LY1,LZ1,1) &
    ,             HELM2 (LX1,LY1,LZ1,1)
    COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1) &
    ,             DUDS  (LX1,LY1,LZ1) &
    ,             DUDT  (LX1,LY1,LZ1) &
    ,             TMP1  (LX1,LY1,LZ1) &
    ,             TMP2  (LX1,LY1,LZ1) &
    ,             TMP3  (LX1,LY1,LZ1)

    REAL ::           TM1   (LX1,LY1,LZ1)
    REAL ::           TM2   (LX1,LY1,LZ1)
    REAL ::           TM3   (LX1,LY1,LZ1)
    REAL ::           DUAX  (LX1)
    REAL ::           YSM1  (LX1)
    EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

    integer :: e

    nel=nelt
    if (imesh == 1) nel=nelv

    NXY=NX1*NY1
    NYZ=NY1*NZ1
    NXZ=NX1*NZ1
    NXYZ=NX1*NY1*NZ1
    NTOT=NXYZ*NEL

    if (icalld == 0) taxhm=0.0
    icalld=icalld+1
    naxhm=icalld
    etime1=dnekclock()

    IF ( .NOT. IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
    CALL RZERO (AU,NTOT)

    do 100 e=1,nel
    
        if (ifaxis) call setaxdy ( ifrzer(e) )
    
        IF (NDIM == 2) THEN
            write(*,*) "Whoops! axhelm"
#if 0
        
        !       2-d case ...............
        
            if (iffast(e)) then
            
            !          Fast 2-d mode: constant properties and undeformed element
            
                h1 = helm1(1,1,1,e)
                call mxm   (wddx,nx1,u(1,1,1,e),nx1,tm1,nyz)
                call mxm   (u(1,1,1,e),nx1,wddyt,ny1,tm2,ny1)
                call col2  (tm1,g4m1(1,1,1,e),nxyz)
                call col2  (tm2,g5m1(1,1,1,e),nxyz)
                call add3  (au(1,1,1,e),tm1,tm2,nxyz)
                call cmult (au(1,1,1,e),h1,nxyz)
            
            else
            
            !          General case, speed-up for undeformed elements
            
                call mxm  (dxm1,nx1,u(1,1,1,e),nx1,dudr,nyz)
                call mxm  (u(1,1,1,e),nx1,dytm1,ny1,duds,ny1)
                call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
                call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
                if (ifdfrm(e)) then
                    call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
                    call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
                endif
                call col2 (tmp1,helm1(1,1,1,e),nxyz)
                call col2 (tmp2,helm1(1,1,1,e),nxyz)
                call mxm  (dxtm1,nx1,tmp1,nx1,tm1,nyz)
                call mxm  (tmp2,nx1,dym1,ny1,tm2,ny1)
                call add2 (au(1,1,1,e),tm1,nxyz)
                call add2 (au(1,1,1,e),tm2,nxyz)

            endif
        
#endif
        else
        
        !       3-d case ...............
        
            if (iffast(e)) then
            
            !          Fast 3-d mode: constant properties and undeformed element
            
                h1 = helm1(1,1,1,e)

                call mxm   (wddx,nx1,u(1,1,1,e),nx1,tm1,nyz)
                do 5 iz=1,nz1
                    call mxm   (u(1,1,iz,e),nx1,wddyt,ny1,tm2(1,1,iz),ny1)
                5 END DO
                call mxm   (u(1,1,1,e),nxy,wddzt,nz1,tm3,nz1)
#if 0
                call col2  (tm1,g4m1(1,1,1,e),nxyz)
                call col2  (tm2,g5m1(1,1,1,e),nxyz)
                call col2  (tm3,g6m1(1,1,1,e),nxyz)

                call add3  (au(1,1,1,e),tm1,tm2,nxyz)
                call add2  (au(1,1,1,e),tm3,nxyz)

                call cmult (au(1,1,1,e),h1,nxyz)
#else
                 
                au(:,:,:,e) = &
                h1 * ( &
                au(:,:,:,e) + tm1*g4m1(:,:,:,e) &
                + tm2*g5m1(:,:,:,e) + tm3*g6m1(:,:,:,e) &
                )
#endif

            
            else
                write(*,*) "Woops! axhelm"
#if 0
            
            !          General case, speed-up for undeformed elements
            
                call mxm(dxm1,nx1,u(1,1,1,e),nx1,dudr,nyz)
                do 10 iz=1,nz1
                    call mxm(u(1,1,iz,e),nx1,dytm1,ny1,duds(1,1,iz),ny1)
                10 END DO
                call mxm     (u(1,1,1,e),nxy,dztm1,nz1,dudt,nz1)
                call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
                call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
                call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
                if (ifdfrm(e)) then
                    call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
                    call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
                    call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
                    call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
                    call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
                    call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
                endif
                call col2 (tmp1,helm1(1,1,1,e),nxyz)
                call col2 (tmp2,helm1(1,1,1,e),nxyz)
                call col2 (tmp3,helm1(1,1,1,e),nxyz)
                call mxm  (dxtm1,nx1,tmp1,nx1,tm1,nyz)
                do 20 iz=1,nz1
                    call mxm(tmp2(1,1,iz),nx1,dym1,ny1,tm2(1,1,iz),ny1)
                20 END DO
                call mxm  (tmp3,nxy,dzm1,nz1,tm3,nz1)
                call add2 (au(1,1,1,e),tm1,nxyz)
                call add2 (au(1,1,1,e),tm2,nxyz)
                call add2 (au(1,1,1,e),tm3,nxyz)
            
#endif
            endif
        
        endif
    
    100 END DO

    call addcol4 (au,helm2,bm1,u,ntot)

!     If axisymmetric, add a diagonal term in the radial direction (ISD=2)

    if (ifaxis .AND. (isd == 2)) then
        write(*,*) "Whoops! axhelm 3"
#if 0
        do 200 e=1,nel
        
            if (ifrzer(e)) then
                call mxm(u  (1,1,1,e),nx1,datm1,ny1,duax,1)
                call mxm(ym1(1,1,1,e),nx1,datm1,ny1,ysm1,1)
            endif
        
            do 190 j=1,ny1
                do 190 i=1,nx1
                !               if (ym1(i,j,1,e).ne.0.) then
                    if (ifrzer(e)) then
                        term1 = 0.0
                        if(j /= 1) &
                        term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                        term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i) &
                        *jacm1(i,1,1,e)/ysm1(i)
                    else
                        term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                        term2 = 0.
                    endif
                    au(i,j,1,e) = au(i,j,1,e) &
                    + helm1(i,j,1,e)*(term1+term2)
                !               endif
            190 END DO
        200 END DO
#endif
    endif

    taxhm=taxhm+(dnekclock()-etime1)
    return
    end subroutine axhelm

!=======================================================================
    subroutine setfast (helm1,helm2,imesh)
!-------------------------------------------------------------------

!     Set logicals for fast evaluation of A*x

!-------------------------------------------------------------------
    use size_m
    use input
    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV
    REAL :: HELM1(NX1,NY1,NZ1,1), HELM2(NX1,NY1,NZ1,1)

    IF (IMESH == 1) NEL=NELV
    IF (IMESH == 2) NEL=NELT
    NXYZ = NX1*NY1*NZ1
    NTOT = NXYZ*NEL

    DELTA = 1.E-9
    X    = 1.+DELTA
    Y    = 1.
    DIFF = ABS(X-Y)
    IF (DIFF == 0.0) EPSM = 1.E-6
    IF (DIFF > 0.0) EPSM = 1.E-13

    DO 100 ie=1,NEL
        IFFAST(ie) = .FALSE. 
        IF (IFDFRM(ie) .OR. IFAXIS .OR. IFMODEL ) THEN
            IFFAST(ie) = .FALSE. 
        !            write(*,*) IFDFRM(ie), IFAXIS, IFMODEL
        ELSE
            H1MIN  = VLMIN(HELM1(1,1,1,ie),NXYZ)
            H1MAX  = VLMAX(HELM1(1,1,1,ie),NXYZ)
            den    = abs(h1max)+abs(h1min)
            if (den > 0) then
                TESTH1 = ABS((H1MAX-H1MIN)/(H1MAX+H1MIN))
                IF (TESTH1 < EPSM) IFFAST(ie) = .TRUE. 
            else
                iffast(ie) = .TRUE. 
            endif
        ENDIF
    100 END DO

    IFH2   = .FALSE. 
    TESTH2 =  VLAMAX(HELM2,NTOT)
    IF (TESTH2 > 0.) IFH2 = .TRUE. 
    return
    end subroutine setfast

!=======================================================================
    subroutine sfastax
!----------------------------------------------------------------------

!     For undeformed elements, set up appropriate elemental matrices
!     and geometric factors for fast evaluation of Ax.

!----------------------------------------------------------------------
    use size_m
    use dxyz
    use geom
    use wz_m
    COMMON /FASTAX/ WDDX(LX1,LY1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV
    LOGICAL :: IFIRST
    SAVE    IFIRST
    DATA    IFIRST / .TRUE. /

    NXX=NX1*NX1
    IF (IFIRST) THEN
        CALL RZERO(WDDX,NXX)
        DO 100 I=1,NX1
            DO 100 J=1,NX1
                DO 100 IP=1,NX1
                    WDDX(I,J) = WDDX(I,J) + WXM1(IP)*DXM1(IP,I)*DXM1(IP,J)
        100 END DO
        NYY=NY1*NY1
        CALL RZERO(WDDYT,NYY)
        DO 200 I=1,NY1
            DO 200 J=1,NY1
                DO 200 IP=1,NY1
                    WDDYT(J,I) = WDDYT(J,I) + WYM1(IP)*DYM1(IP,I)*DYM1(IP,J)
        200 END DO
        NZZ=NZ1*NZ1
        CALL RZERO(WDDZT,NZZ)
        DO 300 I=1,NZ1
            DO 300 J=1,NZ1
                DO 300 IP=1,NZ1
                    WDDZT(J,I) = WDDZT(J,I) + WZM1(IP)*DZM1(IP,I)*DZM1(IP,J)
        300 END DO
        IFIRST= .FALSE. 
    ENDIF

    IF (NDIM == 3) THEN
        DO 1001 IE=1,NELT
            IF ( .NOT. IFDFRM(IE)) THEN
                DO 1000 IZ=1,NZ1
                    DO 1000 IY=1,NY1
                        DO 1000 IX=1,NX1
                            G4M1(IX,IY,IZ,IE)=G1M1(IX,IY,IZ,IE)/WXM1(IX)
                            G5M1(IX,IY,IZ,IE)=G2M1(IX,IY,IZ,IE)/WYM1(IY)
                            G6M1(IX,IY,IZ,IE)=G3M1(IX,IY,IZ,IE)/WZM1(IZ)
                1000 END DO
            ENDIF
        1001 END DO
    ELSE
        DO 2001 IE=1,NELT
            IF ( .NOT. IFDFRM(IE)) THEN
                DO 2000 IY=1,NY1
                    DO 2000 IX=1,NX1
                        G4M1(IX,IY,1,IE)=G1M1(IX,IY,1,IE)/WXM1(IX)
                        G5M1(IX,IY,1,IE)=G2M1(IX,IY,1,IE)/WYM1(IY)
                2000 END DO
            ENDIF
        2001 END DO
    ENDIF
    return
    end subroutine sfastax

!=======================================================================
    subroutine setprec (dpcm1,helm1,helm2,imsh,isd)
!-------------------------------------------------------------------

!     Generate diagonal preconditioner for the Helmholtz operator.

!-------------------------------------------------------------------
    use size_m
    use dxyz
    use geom
    use input
    use mass
    use tstep
    use wz_m
    REAL ::            DPCM1 (LX1,LY1,LZ1,1)
    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV
    REAL ::            HELM1(NX1,NY1,NZ1,1), HELM2(NX1,NY1,NZ1,1)
    REAL :: YSM1(LY1)

    nel=nelt
    if (imsh == 1) nel=nelv

    ntot = nel*nx1*ny1*nz1

!     The following lines provide a convenient debugging option
!     call rone(dpcm1,ntot)
!     if (ifield.eq.1) call copy(dpcm1,binvm1,ntot)
!     if (ifield.eq.2) call copy(dpcm1,bintm1,ntot)
!     return

    CALL RZERO(DPCM1,NTOT)
    DO 1000 IE=1,NEL

        IF (IFAXIS) CALL SETAXDY ( IFRZER(IE) )

        DO 320 IQ=1,NX1
            DO 320 IZ=1,NZ1
                DO 320 IY=1,NY1
                    DO 320 IX=1,NX1
                        DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                        G1M1(IQ,IY,IZ,IE) * DXTM1(IX,IQ)**2
        320 END DO
        DO 340 IQ=1,NY1
            DO 340 IZ=1,NZ1
                DO 340 IY=1,NY1
                    DO 340 IX=1,NX1
                        DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                        G2M1(IX,IQ,IZ,IE) * DYTM1(IY,IQ)**2
        340 END DO

        IF (NDIM == 3) THEN
            DO 360 IQ=1,NZ1
                DO 360 IZ=1,NZ1
                    DO 360 IY=1,NY1
                        DO 360 IX=1,NX1
                            DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                            G3M1(IX,IY,IQ,IE) * DZTM1(IZ,IQ)**2
            360 END DO
        
        !       Add cross terms if element is deformed.
        
            IF (IFDFRM(IE)) THEN
                write(*,*) "Whoops!"
#if 0
                DO 600 IY=1,NY1
                    DO 600 IZ=1,NZ1
                        DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE) &
                        + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY) &
                        + G5M1(1,IY,IZ,IE) * DXTM1(1,1)*DZTM1(IZ,IZ)
                        DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE) &
                        + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY) &
                        + G5M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DZTM1(IZ,IZ)
                600 END DO
                DO 700 IX=1,NX1
                    DO 700 IZ=1,NZ1
                        DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE) &
                        + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX) &
                        + G6M1(IX,1,IZ,IE) * DYTM1(1,1)*DZTM1(IZ,IZ)
                        DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE) &
                        + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX) &
                        + G6M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DZTM1(IZ,IZ)
                700 END DO
                DO 800 IX=1,NX1
                    DO 800 IY=1,NY1
                        DPCM1(IX,IY,1,IE) = DPCM1(IX,IY,1,IE) &
                        + G5M1(IX,IY,1,IE) * DZTM1(1,1)*DXTM1(IX,IX) &
                        + G6M1(IX,IY,1,IE) * DZTM1(1,1)*DYTM1(IY,IY)
                        DPCM1(IX,IY,NZ1,IE) = DPCM1(IX,IY,NZ1,IE) &
                        + G5M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DXTM1(IX,IX) &
                        + G6M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DYTM1(IY,IY)
                800 END DO
#endif
            ENDIF
        ELSE
        
            IF (IFDFRM(IE)) THEN
                write(*,*) "Whoops!"
#if 0
                IZ=1
                DO 602 IY=1,NY1
                    DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE) &
                    + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY)
                    DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE) &
                    + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY)
                602 END DO
                DO 702 IX=1,NX1
                    DO 702 IZ=1,NZ1
                        DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE) &
                        + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX)
                        DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE) &
                        + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX)
                702 END DO
#endif
            ENDIF
        ENDIF
    1000 END DO

    CALL COL2    (DPCM1,HELM1,NTOT)
    CALL ADDCOL3 (DPCM1,HELM2,BM1,NTOT)

!     If axisymmetric, add a diagonal term in the radial direction (ISD=2)

    IF (IFAXIS .AND. (ISD == 2)) THEN
        DO 1200 IEL=1,NEL
        
            IF (IFRZER(IEL)) THEN
                CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
            ENDIF
        
            DO 1190 J=1,NY1
                DO 1190 I=1,NX1
                    IF (YM1(I,J,1,IEL) /= 0.) THEN
                        TERM1 = BM1(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                        IF (IFRZER(IEL)) THEN
                            TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J) &
                            *JACM1(I,1,1,IEL)/YSM1(I)
                        ELSE
                            TERM2 = 0.
                        ENDIF
                        DPCM1(I,J,1,IEL) = DPCM1(I,J,1,IEL) &
                        + HELM1(I,J,1,IEL)*(TERM1+TERM2)
                    ENDIF
            1190 END DO
        1200 END DO
    ENDIF

    CALL DSSUM (DPCM1,NX1,NY1,NZ1)
    CALL INVCOL1 (DPCM1,NTOT)

    return
    end subroutine setprec

!=======================================================================
    subroutine chktcg1 (tol,res,h1,h2,mask,mult,imesh,isd)
!-------------------------------------------------------------------

!     Check that the tolerances are not too small for the CG-solver.
!     Important when calling the CG-solver (Gauss-Lobatto mesh) with
!     zero Neumann b.c.

!-------------------------------------------------------------------
    use size_m
    use eigen
    use input
    use mass
    COMMON  /CPRINT/ IFPRINT
    LOGICAL ::          IFPRINT
    COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT) &
    ,             W2   (LX1,LY1,LZ1,LELT)
    REAL :: RES  (LX1,LY1,LZ1,1)
    REAL :: H1   (LX1,LY1,LZ1,1)
    REAL :: H2   (LX1,LY1,LZ1,1)
    REAL :: MULT (LX1,LY1,LZ1,1)
    REAL :: MASK (LX1,LY1,LZ1,1)

    IF (EIGAA /= 0.) THEN
        ACONDNO = EIGGA/EIGAA
    ELSE
        ACONDNO = 10.
    ENDIF

!     Single or double precision???

    DELTA = 1.E-9
    X     = 1.+DELTA
    Y     = 1.
    DIFF  = ABS(X-Y)
    IF (DIFF == 0.) EPS = 1.E-6
    IF (DIFF > 0.) EPS = 1.E-13

    IF (IMESH == 1) THEN
        NL  = NELV
        VOL = VOLVM1
    ELSEIF (IMESH == 2) THEN
        NL  = NELT
        VOL = VOLTM1
    ENDIF
    NTOT1 = NX1*NY1*NZ1*NL
    CALL COPY (W1,RES,NTOT1)

    IF (IMESH == 1) THEN
        CALL COL3 (W2,BINVM1,W1,NTOT1)
        RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLVM1)
    ELSE
        CALL COL3 (W2,BINTM1,W1,NTOT1)
        RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLTM1)
    ENDIF
    RMIN   = EPS*RINIT
    IF (TOL < RMIN) THEN
        IF (NID == 0 .AND. IFPRINT) &
        WRITE (6,*) 'New CG1-tolerance (RINIT*epsm) = ',RMIN,TOL
        TOL = RMIN
    ENDIF

    CALL RONE (W1,NTOT1)
    BCNEU1 = GLSC3(W1,MASK,MULT,NTOT1)
    BCNEU2 = GLSC3(W1,W1  ,MULT,NTOT1)
    BCTEST = ABS(BCNEU1-BCNEU2)

    CALL AXHELM (W2,W1,H1,H2,IMESH,ISD)
    CALL COL2   (W2,W2,NTOT1)
    CALL COL2   (W2,BM1,NTOT1)
    BCROB  = SQRT(GLSUM(W2,NTOT1)/VOL)

    IF ((BCTEST < .1) .AND. (BCROB < (EPS*ACONDNO))) THEN
    !         OTR = GLSC3 (W1,RES,MULT,NTOT1)
        TOLMIN = RINIT*EPS*10.
        IF (TOL < TOLMIN) THEN
            TOL = TOLMIN
            IF (NID == 0 .AND. IFPRINT) &
            WRITE(6,*) 'New CG1-tolerance (Neumann) = ',TOLMIN
        ENDIF
    ENDIF

    return
    end subroutine chktcg1
!=======================================================================
    subroutine cggo(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name)
!-------------------------------------------------------------------------

!     Solve the Helmholtz equation, H*U = RHS,
!     using preconditioned conjugate gradient iteration.
!     Preconditioner: diag(H).

!------------------------------------------------------------------------
    use size_m
    use geom
    use fdmh1
    use input
    use mass
    use soln
    use tstep

    COMMON  /CPRINT/ IFPRINT, IFHZPC
    LOGICAL ::          IFPRINT, IFHZPC

    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV
    logical :: ifmcor,ifprint_hmh

    real :: x(1),f(1),h1(1),h2(1),mask(1),mult(1),binv(1)
    parameter        (lg=lx1*ly1*lz1*lelt)
    COMMON /SCRCG/ d (lg) , scalar(2)
    common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)

    parameter (maxcg=900)
    common /tdarray/ diagt(maxcg),upper(maxcg)
    common /iterhm/ niterhm
    character(4) :: name

    if (ifsplit .AND. name == 'PRES' .AND. param(42) == 0) then
        n = nx1*ny1*nz1*nelv
        call copy      (x,f,n)
        call hmh_gmres (x,h1,h2,mult,iter)
        niterhm = iter
        return
    endif
!      write(6,*) ifsplit,name,param(44),' P44 C'

! **  zero out stuff for Lanczos eigenvalue estimator
    call rzero(diagt,maxcg)
    call rzero(upper,maxcg)
    rho = 0.00

!     Initialization

    NXYZ   = NX1*NY1*NZ1
    NEL    = NELV
    VOL    = VOLVM1
    IF (IMSH == 2) NEL=NELT
    IF (IMSH == 2) VOL=VOLTM1
    n      = NEL*NXYZ

    tol=abs(tin)
    if (param(22) /= 0) tol=abs(param(22))
    if (name == 'PRES' .AND. param(21) /= 0) tol=abs(param(21))
    if (tin < 0)       tol=abs(tin)
    niter = min(maxit,maxcg)

!     Speed-up for undeformed elements and constant properties.
    if ( .NOT. ifsolv) then
        call setfast(h1,h2,imesh)
        ifsolv = .TRUE. 
    endif

!     Set up diag preconditioner.

    if (kfldfdm < 0) then
        call setprec(D,h1,h2,imsh,isd)
    elseif(param(100) /= 2) then
        call set_fdm_prec_h1b(d,h1,h2,nel)
    endif

    call copy (r,f,n)
    call rzero(x,n)
    call rzero(p,n)

!     Check for non-trivial null-space

    ifmcor = .FALSE. 
    h2max = glmax(h2  ,n)
    skmin = glmin(mask,n)
    if (skmin > 0 .AND. h2max == 0) ifmcor = .TRUE. 

    if (name == 'PRES') then
    !        call ortho (r)           ! Commented out March 15, 2011,pff
    elseif (ifmcor) then

        smean = -1./glsum(bm1,n) ! Modified 5/4/12 pff
        rmean = smean*glsc2(r,mult,n)
        call copy(x,bm1,n)
        call dssum(x,nx1,ny1,nz1)
        call add2s2(r,x,rmean,n)
        call rzero(x,n)
    endif

    krylov = 0
    rtz1=1.0
    niterhm = 0

    do iter=1,niter
    
        if (kfldfdm < 0) then  ! Jacobi Preconditioner
        !           call copy(z,r,n)
            call col3(z,r,d,n)
        else                                       ! Schwarz Preconditioner
        write (*,*) "Oops: kfldfdm"
#if 0
            if (name == 'PRES' .AND. param(100) == 2) then
                call h1_overlap_2(z,r,mask)
                call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
                call add2         (z,w,n)
            else
                call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
                if (name == 'PRES') then
                    call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
                    call add2         (z,w,n)
                endif
            endif
#endif
        endif
    
        if (name == 'PRES') then
            call ortho (z)
        elseif (ifmcor) then
            rmean = smean*glsc2(z,bm1,n)
            call cadd(z,rmean,n)
        endif
    !        write(6,*) rmean,ifmcor,' ifmcor'
    
        rtz2=rtz1
        scalar(1)=vlsc3 (z,r,mult,n)
        scalar(2)=vlsc32(r,mult,binv,n)
        call gop(scalar,w,'+  ',2)
        rtz1=scalar(1)
        rbn2=sqrt(scalar(2)/vol)
        if (iter == 1) rbn0 = rbn2
        if (param(22) < 0) tol=abs(param(22))*rbn0
        if (tin < 0)       tol=abs(tin)*rbn0

        ifprint_hmh = .FALSE. 
        if (nid == 0 .AND. ifprint .AND. param(74) /= 0) ifprint_hmh= .TRUE. 
        if (nid == 0 .AND. istep == 1)                 ifprint_hmh= .TRUE. 

        if (ifprint_hmh) &
        write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1),h2(1)


    !        Always take at least one iteration   (for projection) pff 11/23/98
#ifndef TST_WSCAL
        IF (rbn2 <= TOL .AND. (iter > 1 .OR. istep <= 5)) THEN
#else
            iter_max = param(150)
            if (name == 'PRES') iter_max = param(151)
            if (iter > iter_max) then
#endif
            !        IF (rbn2.LE.TOL) THEN
                NITER = ITER-1
            !           IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
                if (nid == 0) &
                write(6,3000) istep,name,niter,rbn2,rbn0,tol
                goto 9999
            ENDIF
        
            beta = rtz1/rtz2
            if (iter == 1) beta=0.0
            call add2s1 (p,z,beta,n)
            call axhelm (w,p,h1,h2,imsh,isd)
            call dssum  (w,nx1,ny1,nz1)
            call col2   (w,mask,n)
        
            rho0 = rho
            rho  = glsc3(w,p,mult,n)
            alpha=rtz1/rho
            alphm=-alpha
            call add2s2(x,p ,alpha,n)
            call add2s2(r,w ,alphm,n)
        
        !        Generate tridiagonal matrix for Lanczos scheme
            if (iter == 1) then
                krylov = krylov+1
                diagt(iter) = rho/rtz1
            elseif (iter <= maxcg) then
                krylov = krylov+1
                diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
                upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
            endif
            1000 enddo
            niter = iter-1
        
            if (nid == 0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
            3000 format(4x,i7,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
            3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
            3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
            9999 continue
            niterhm = niter
            ifsolv = .FALSE. 
        
        
        !     Call eigenvalue routine for Lanczos scheme:
        !          two work arrays are req'd if you want to save "diag & upper"
        
        !     if (iter.ge.3) then
        !        niter = iter-1
        !        call calc (diagt,upper,w,z,krylov,dmax,dmin)
        !        cond = dmax/dmin
        !        if (nid.eq.0) write(6,6) istep,cond,dmin,dmax,' lambda'
        !     endif
        !   6 format(i9,1p3e12.4,4x,a7)
        
        !     if (n.gt.0) write(6,*) 'quit in cggo'
        !     if (n.gt.0) call exitt
        !     call exitt
            return
            end subroutine cggo
        !=======================================================================
            function vlsc32(r,b,m,n)
            real :: r(1),b(1),m(1)
            s = 0.
            do i=1,n
                s = s + b(i)*m(i)*r(i)*r(i)
            enddo
            vlsc32 = s
            return
            end function vlsc32
        !=======================================================================
            subroutine calc (diag,upper,d,e,n,dmax,dmin)
        
            dimension diag(n),upper(n)
            dimension d(n),e(n)
        
            call copy (d,diag ,n)
            call copy (e,upper,n)
        
            do 15 l=1,n
                iter = 0
            
                1 do 12 m=l,n-1
                    dd = abs( d(m) ) + abs( d(m+1) )
                    if ( abs(e(m)) + dd == dd ) goto 2
                12 END DO
            
                m = n
                2 if ( m /= l ) then
                
                    if ( iter == 30 ) then
                        write (6,*) 'too many iterations'
                        return
                    endif
                
                    iter = iter + 1
                    g = ( d(l+1) - d(l) ) / ( 2.0 * e(l) )
                    r = sqrt( g**2 + 1.0 )
                
                ! sign is defined as a(2) * abs( a(1) )
                
                    g = d(m) - d(l) + e(l)/(g+sign(r,g))
                    s = 1.0
                    c = 1.0
                    p = 0.0
                
                    do 14 i = m-1,l,-1
                        f = s * e(i)
                        b = c * e(i)
                        if ( abs(f) >= abs(g) ) then
                            c = g/f
                            r = sqrt( c**2 + 1.0 )
                            e(i+1) = f*r
                            s = 1.0/r
                            c = c*s
                        else
                            s = f/g
                            r = sqrt( s**2 + 1.0 )
                            e(i+1) = g*r
                            c = 1.0 / r
                            s = s * c
                        endif
                    
                        g = d(i+1) - p
                        r = ( d(i) - g ) * s + 2.0 * c * b
                        p = s * r
                        d(i+1) = g + p
                        g = c*r - b
                    14 END DO
                
                    d(l) = d(l) - p
                    e(l) = g
                    e(m) = 0.0
                    goto 1
                
                endif
            
            15 END DO
        
            dmax = 0.0
            dmin = d(1)
        
            do 40 i = 1 , n
                dmax = abs( max( d(i) , dmax ) )
                dmin = abs( min( d(i) , dmin ) )
            40 END DO
        
            return
            end subroutine calc
        !-----------------------------------------------------------------------
#if 0
            subroutine fdm_h1(z,r,d,mask,mult,nel,kt,rr)
            use size_m
            use fdmh1
            include 'TOTAL'
        
            common /ctmp0/ w(lx1,ly1,lz1)
        
        
        !     Overlapping Schwarz, FDM based
        
            real :: z(nx1,ny1,nz1,1)
            real :: r(nx1,ny1,nz1,1)
            real :: d(nx1,ny1,nz1,1)
            real :: mask(nx1,ny1,nz1,1)
            real :: mult(nx1,ny1,nz1,1)
            real :: rr(nx1,ny1,nz1,1)
        
            integer :: kt(lelt,3)
        
            integer :: icalld
            save    icalld
            data    icalld /0/
        
            n1 = nx1
            n2 = nx1*nx1
            n3 = nx1*nx1*nx1
            ntot = nx1*ny1*nz1*nel
        
            if (ifbhalf) then
                call col3(rr,r,bhalf,ntot)
            else
                call copy(rr,r,ntot)
            !        call col2(rr,mult,ntot)
            endif
        !      if (nid.eq.0.and.icalld.eq.0) write(6,*) 'In fdm_h1',nel
            icalld = icalld+1
        
        
            do ie=1,nel
                if (if3d) then
                !           Transfer to wave space:
                    call mxm(fdst(1,kt(ie,1)),n1,rr(1,1,1,ie),n1,w,n2)
                    do iz=1,n1
                        call mxm(w(1,1,iz),n1,fds (1,kt(ie,2)),n1,z(1,1,iz,ie),n1)
                    enddo
                    call mxm(z(1,1,1,ie),n2,fds (1,kt(ie,3)),n1,w,n1)
                
                !           fdsolve:
                
                    call col2(w,d(1,1,1,ie),n3)
                
                !           Transfer to physical space:
                
                    call mxm(w,n2,fdst(1,kt(ie,3)),n1,z(1,1,1,ie),n1)
                    do iz=1,n1
                        call mxm(z(1,1,iz,ie),n1,fdst(1,kt(ie,2)),n1,w(1,1,iz),n1)
                    enddo
                    call mxm(fds (1,kt(ie,1)),n1,w,n1,z(1,1,1,ie),n2)
                
                else
                !           Transfer to wave space:
                    call mxm(fdst(1,kt(ie,1)),n1,rr(1,1,1,ie),n1,w,n1)
                    call mxm(w,n1,fds (1,kt(ie,2)),n1,z(1,1,1,ie),n1)
                
                !           fdsolve:
                
                    call col2(z(1,1,1,ie),d(1,1,1,ie),n2)
                
                !           Transfer to physical space:
                
                    call mxm(z(1,1,1,ie),n1,fdst(1,kt(ie,2)),n1,w,n1)
                    call mxm(fds (1,kt(ie,1)),n1,w,n1,z(1,1,1,ie),n1)
                
                endif
            enddo
        
        !     call copy(vx,rr,ntot)
        !     call copy(vy,z,ntot)
        !     call prepost(.true.)
        !     write(6,*) 'quit in fdm'
        !     call exitt
        
            if (ifbhalf) call col2(z,bhalf,ntot)
        
        !     call col2 (z,mult,ntot)
            call dssum(z,nx1,ny1,nz1)
            call col2 (z,mask,ntot)
        
            return
            end subroutine fdm_h1
#endif
        !-----------------------------------------------------------------------
            subroutine set_fdm_prec_h1A_gen
        
            use size_m
            use dxyz
            use fdmh1
            use input
            use mass
            use wz_m
        
        
            COMMON /CTMP0/ W(LX1,LX1),aa(lx1,lx1),bb(lx1,lx1)
        
            integer :: left,right
        
        !     Set up generic operators for fdm applied to H1 operator (Helmholtz)
        
        !     3 cases:   E (or P),   "D"  or "N"  for E-E bc, Dirichlet, or Neuamann.
        
        !     Since there are 2 endpoints, there are a total of 9 types.
        
        
            n  = nx1
            n2 = nx1*nx1
        
            delta = abs( zgm1(2,1) - zgm1(1,1) )
            bbh   = 0.5*delta
            aah   = 1./delta
        
            l = 0
            do right = 1,3
                do left  = 1,3
                    l = l+1
                
                    call rzero(bb,n2)
                    do i=1,nx1
                        bb(i,i) = wxm1(i)
                    enddo
                
                !        A = D^T B D
                
                    call mxm(BB,n,Dxm1 ,n,w,n)
                    call mxm(Dxtm1,n,w,n,AA,n)
                    if (left == 1) then
                    !           Internal
                        bb(1,1) = bb(1,1) + bbh
                        aa(1,1) = aa(1,1) + aah
                    elseif (left == 2) then
                    !           Dirichlet
                        bb(1,1) = 1.
                        do i=1,n
                            aa(i,1) = 0.
                            aa(1,i) = 0.
                        enddo
                        aa(1,1) = 1.
                    endif
                
                    if (right == 1) then
                    !           Internal
                        bb(n,n) = bb(n,n) + bbh
                        aa(n,n) = aa(n,n) + aah
                    elseif (right == 2) then
                    !           Dirichlet
                        bb(n,n) = 1.
                        do i=1,n
                            aa(i,n) = 0.
                            aa(n,i) = 0.
                        enddo
                        aa(n,n) = 1.
                    endif
                
                !        Scale out mass matrix, so we can precondition w/ binvhf.
                
                !        ifbhalf = .true.

                    ifbhalf = .FALSE. 
                    if (ifbhalf) call rescale_abhalf (aa,bb,w,n)
                
                !        Now, compute eigenvectors/eigenvalues
                
                    call generalev(aa,bb,dd(1,l),n,w)
                    call copy(fds(1,l),aa,n*n)
                    call transpose(fdst(1,l),n,fds(1,l),n)
                
                enddo
            enddo
            ntot = nx1*ny1*nz1*nelv
            if (ifbhalf) call copy (bhalf,binvm1,ntot)
            if (ifbhalf) call vsqrt(bhalf,ntot)
        
            return
            end subroutine set_fdm_prec_h1A_gen
        !-----------------------------------------------------------------------
            subroutine set_fdm_prec_h1A_els
        
            use size_m
            use dxyz
            use fdmh1
            use geom
            use input
            use soln
            use topol
            use wz_m
        
            COMMON /CTMP0/ W(LX1,LX1),aa(lx1,lx1),bb(lx1,lx1) &
            , mask(lx1,ly1,lz1,lelt)
            real :: mask
            character(3) :: cb
        
        
        !     Set up element specific information
        
        !     3 cases:   E (or P),   "D"  or "N"  for E-E bc, Dirichlet, or Neuamann.
        
        !     Since there are 2 endpoints, there are a total of 9 types.
        
        
            ntot = nx1*ny1*nz1*nelt
            kf0 = 1
            kf1 = 0
            if (ifheat)  kf0 = 0
            if (ifflow)  kf1 = ndim
            if (ifsplit) kf1 = ndim+1
            do kfld=kf0,kf1
                ifld = 1
                if (kfld == 0) ifld = 2
            
                if (kfld == 0)      call copy(mask, tmask,ntot)
                if (kfld == 1)      call copy(mask,v1mask,ntot)
                if (kfld == 2)      call copy(mask,v2mask,ntot)
                if (kfld == 3)      call copy(mask,v3mask,ntot)
                if (kfld == ndim+1) call copy(mask, pmask,ntot)
            
                do ie=1,nelv
                    do ifacedim = 1,ndim
                    
                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    !              Mask pointers
                    
                        ii = 2
                        jj = 2
                        kk = 2
                    
                        if (ifacedim == 1) ii = 1
                        if (ifacedim == 2) jj = 1
                        if (ifacedim == 3) kk = 1
                        k1 = ii+nx1*(jj-1)
                        if (if3d) k1 = ii+nx1*(jj-1) + nx1*nx1*(kk-1)
                    
                        if (ifacedim == 1) ii = nx1
                        if (ifacedim == 2) jj = nx1
                        if (ifacedim == 3) kk = nx1
                        k2 = ii+nx1*(jj-1)
                        if (if3d) k2 = ii+nx1*(jj-1) + nx1*nx1*(kk-1)
                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    
                        iface = 2*ifacedim-1
                        jface = iface+1
                    
                    !              Convert to preproc   :(
                        iface = eface(iface)
                        jface = eface(jface)
                    
                    !              "left" bc
                    
                        cb = cbc(iface,ie,ifld)
                        if (cb == 'E  ' .OR. cb == 'P  ' .OR. cb == 'p  ') then
                        !                 Internal
                            ic1 = 1
                        elseif (mask(k1,1,1,ie) == 0) then
                        !                 Dirichlet
                            ic1 = 2
                        else
                        !                 Neumann
                            ic1 = 3
                        endif
                    !              write(6,*) ie,iface,'cbl: ',cb,ic1,k1,mask(k1,1,1,ie)
                    
                    !              "right" bc
                    
                        cb = cbc(jface,ie,ifld)
                        if (cb == 'E  ' .OR. cb == 'P  ' .OR. cb == 'p  ') then
                        !                 Internal
                            jc1 = 1
                        elseif (mask(k2,1,1,ie) == 0) then
                        !                 Dirichlet
                            jc1 = 2
                        else
                        !                 Neumann
                            jc1 = 3
                        endif
                    !              write(6,*) ie,jface,'cbr: ',cb,jc1,k2,mask(k2,1,1,ie)
                    
                        ijc = ic1 + 3*(jc1-1)
                        ktype(ie,ifacedim,kfld) = ijc
                    
                    enddo
                enddo
            enddo
        
        !     Boundary condition issues resolved... now resolve length scales
        
        
        !     Compute mean distance between faces for each direction in ref. domain
        
            do ie = 1,nelt
                do idim=1,ndim
                    k1 = 1
                    k2 = nz1
                    if (idim == 3 .OR. ndim == 2) k2=1
                    j1 = 1
                    j2 = ny1
                    if (idim == 2) j2=1
                    i1 = 1
                    i2 = nx1
                    if (idim == 1) i2=1
                
                !        l -- face 1,  l+jump = face 2
                
                    jump = (nx1-1)*nx1**(idim-1)
                    l    = 0
                    dlm  = 0
                    wgt  = 0
                    do k=k1,k2
                        do j=j1,j2
                            do i=i1,i2
                                l = l+1
                                dl2 = (xm1(i+jump,j,k,ie)-xm1(i,j,k,ie))**2 &
                                + (ym1(i+jump,j,k,ie)-ym1(i,j,k,ie))**2 &
                                + (zm1(i+jump,j,k,ie)-zm1(i,j,k,ie))**2
                                dlm = dlm + dl2*wxm1(i)*wxm1(j)*wxm1(k)
                                wgt = wgt + wxm1(i)*wxm1(j)*wxm1(k)
                            
                            enddo
                        enddo
                    enddo
                
                    dlm             = sqrt(dlm/wgt)
                    elsize(idim,ie) = dlm/2.
                
                enddo
            !        write(6,1) ie,' elsize:',(elsize(k,ie),k=1,ndim)
            enddo
            1 format(i8,a8,1p3e15.4)
        
            return
            end subroutine set_fdm_prec_h1A_els
        !-----------------------------------------------------------------------
            subroutine set_fdm_prec_h1b(d,h1,h2,nel)
            use size_m
            use fdmh1
            use geom
            use input
            real :: d (nx1,ny1,nz1,1)
            real :: h1(nx1,ny1,nz1,1)
            real :: h2(nx1,ny1,nz1,1)
        
        !     Set up diagonal for FDM for each spectral element
        
            nxyz = nx1*ny1*nz1
            if (if3d) then
                do ie=1,nel
                    h1b = vlsum(h1(1,1,1,ie),nxyz)/nxyz
                    h2b = vlsum(h2(1,1,1,ie),nxyz)/nxyz
                    k1 = ktype(ie,1,kfldfdm)
                    k2 = ktype(ie,2,kfldfdm)
                    k3 = ktype(ie,3,kfldfdm)
                    vol = elsize(1,ie)*elsize(2,ie)*elsize(3,ie)
                    vl1 = elsize(2,ie)*elsize(3,ie)/elsize(1,ie)
                    vl2 = elsize(1,ie)*elsize(3,ie)/elsize(2,ie)
                    vl3 = elsize(1,ie)*elsize(2,ie)/elsize(3,ie)
                    do i3=1,nz1
                        do i2=1,ny1
                            do i1=1,nx1
                                den = h1b*(vl1*dd(i1,k1) + vl2*dd(i2,k2) + vl3*dd(i3,k3)) &
                                + h2b*vol
                                if (ifbhalf) den = den/vol
                                if (den /= 0) then
                                    d(i1,i2,i3,ie) = 1./den
                                else
                                    d(i1,i2,i3,ie) = 0.
                                
                                !                 write(6,3) 'd=0:'
                                !    $                 ,h1(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2),dd(i3,k3)
                                !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2,k3
                                    3 format(a4,1p4e12.4,8i8)
                                
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            else
                do ie=1,nel
                    if (ifaxis) then
                        h1b = vlsc2(h1(1,1,1,ie),ym1(1,1,1,ie),nxyz)/nxyz
                        h2b = vlsc2(h2(1,1,1,ie),ym1(1,1,1,ie),nxyz)/nxyz
                    else
                        h1b = vlsum(h1(1,1,1,ie),nxyz)/nxyz
                        h2b = vlsum(h2(1,1,1,ie),nxyz)/nxyz
                    endif
                    k1 = ktype(ie,1,kfldfdm)
                    k2 = ktype(ie,2,kfldfdm)
                    vol = elsize(1,ie)*elsize(2,ie)
                    vl1 = elsize(2,ie)/elsize(1,ie)
                    vl2 = elsize(1,ie)/elsize(2,ie)
                    i3=1
                    do i2=1,ny1
                        do i1=1,nx1
                            den = h1b*( vl1*dd(i1,k1) + vl2*dd(i2,k2) ) &
                            + h2b*vol
                            if (ifbhalf) den = den/vol
                            if (den /= 0) then
                                d(i1,i2,i3,ie) = 1./den
                            !                 write(6,3) 'dn0:'
                            !    $                 ,d(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2)
                            !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2
                            else
                                d(i1,i2,i3,ie) = 0.
                            !                 write(6,3) 'd=0:'
                            !    $                 ,h1(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2)
                            !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2
                                2 format(a4,1p3e12.4,8i8)
                            endif
                        !           write(6,1) ie,i1,i2,k1,k2,'d:',d(i1,i2,i3,ie),vol,vl1,vl2
                        !   1       format(5i3,2x,a2,1p4e12.4)
                        enddo
                    enddo
                enddo
            endif
        
            return
            end subroutine set_fdm_prec_h1b
        !-----------------------------------------------------------------------
            subroutine set_fdm_prec_h1A
            use size_m
        
            call set_fdm_prec_h1A_gen
            call set_fdm_prec_h1A_els
        
            return
            end subroutine set_fdm_prec_h1A
        !-----------------------------------------------------------------------
            subroutine generalev(a,b,lam,n,w)
        
        !     Solve the generalized eigenvalue problem  A x = lam B x
        
        !     A -- symm.
        !     B -- symm., pos. definite
        
        !     "SIZE" is included here only to deduce WDSIZE, the working
        !     precision, in bytes, so as to know whether dsygv or ssygv
        !     should be called.
        
            use size_m
            use parallel
        
            real :: a(n,n),b(n,n),lam(n),w(n,n)
            real :: aa(100),bb(100)
        
            parameter (lbw=4*lx1*ly1*lz1*lelv)
            common /bigw/ bw(lbw)
        
            lw = n*n
        !     write(6,*) 'in generalev, =',info,n,ninf
        
        !     call outmat2(a,n,n,n,'aa  ')
        !     call outmat2(b,n,n,n,'bb  ')
        
            call copy(aa,a,100)
            call copy(bb,b,100)
        
            if (ifdblas) then
                call dsygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
            else
                call ssygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
            endif
        
        !     call outmat2(a,n,n,n,'Aeig')
        !     call outmat2(lam,1,n,n,'Deig')
        
            if (info /= 0) then
            
                if (nid == 0) then
                    call outmat2(aa ,n,n,n,'aa  ')
                    call outmat2(bb ,n,n,n,'bb  ')
                    call outmat2(a  ,n,n,n,'Aeig')
                    call outmat2(lam,1,n,n,'Deig')
                endif
            
                ninf = n-info
                write(6,*) 'Error in generalev, info=',info,n,ninf
                call exitt
            endif
        
            return
            end subroutine generalev
        !-----------------------------------------------------------------------
            subroutine outmat2(a,m,n,k,name)
            use size_m
            real :: a(m,n)
            character(4) :: name
        
            n2 = min(n,8)
            write(6,2) nid,name,m,n,k
            do i=1,m
                write(6,1) nid,name,(a(i,j),j=1,n2)
            enddo
        !   1 format(i3,1x,a4,16f6.2)
            1 format(i3,1x,a4,1p8e14.5)
            2 format(/,'Matrix: ',i3,1x,a4,3i8)
            return
            end subroutine outmat2
        !-----------------------------------------------------------------------
            subroutine rescale_abhalf (a,b,w,n)
            real :: a(n,n),b(n,n),w(n)
        
        !             -1/2      -1/2
        !     Set A = B    A  B
        
        
        !     NOTE:   B is *diagonal*
        
        
            do i=1,n
                w(i) = 1./sqrt(b(i,i))
            enddo
        
            do j=1,n
                do i=1,n
                    a(i,j) = a(i,j)*w(i)*w(j)
                enddo
            enddo
        
        !     duh....  don't forget to change B ...  duh...
        
            call ident(b,n)
        
            return
            end subroutine rescale_abhalf
        !-----------------------------------------------------------------------
