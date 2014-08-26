!-----------------------------------------------------------------------
    subroutine cdscal (igeom)

!     Solve the convection-diffusion equation for passive scalar IPSCAL

    use size_m
    use geom
    use input
    use mass
    use mvgeom
    use soln
    use tstep
    COMMON  /CPRINT/ IFPRINT
    LOGICAL ::          IFPRINT
    LOGICAL ::          IFCONV

    COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT) &
    ,TB(LX1,LY1,LZ1,LELT)
    COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT) &
    ,H2(LX1,LY1,LZ1,LELT)

    parameter (ktot = lx1*ly1*lz1*lelt)
    parameter (laxt = mxprev)

    common /trthoi/ napprox(2)
    common /trthov/ approx(ktot,0:laxt)
    common /trthoc/ name4
    character(4) ::     name4

!max    include 'ORTHOT'

    napprox(1) = laxt

    nel    = nelfld(ifield)
    ntot   = nx1*ny1*nz1*nel

    if (igeom == 1) then   ! geometry at t^{n-1}
        call makeq
        call lagscal
    else                   ! geometry at t^n

        IF (IFPRINT) THEN
            IF (IFMODEL .AND. IFKEPS) THEN
                NFLDT = NFIELD - 1
                IF (IFIELD == NFLDT .AND. NID == 0) THEN
                    WRITE (6,*) ' Turbulence Model - k/epsilon solution'
                ENDIF
            ELSE
                IF (IFIELD == 2 .AND. NID == 0) THEN
                    WRITE (6,*) ' Temperature/Passive scalar solution'
                ENDIF
            ENDIF
        ENDIF
        if1=ifield-1
        write(name4,1) if1-1
        1 format('PS',i2)
        if(ifield == 2) write(name4,'(A4)') 'TEMP'

    
    !        New geometry
    
        isd = 1
        if (ifaxis .AND. ifaziv .AND. ifield == 2) isd = 2
    !        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

        do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

            INTYPE = 0
            IF (IFTRAN) INTYPE = -1
            CALL SETHLM  (H1,H2,INTYPE)
            CALL BCNEUSC (TA,-1)
            CALL ADD2    (H2,TA,NTOT)
            CALL BCDIRSC (T(1,1,1,1,IFIELD-1))
            CALL AXHELM  (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMESH,isd)
            CALL SUB3    (TB,BQ(1,1,1,1,IFIELD-1),TA,NTOT)
            CALL BCNEUSC (TA,1)
            CALL ADD2    (TB,TA,NTOT)

        !        CALL HMHOLTZ (name4,TA,TB,H1,H2
        !    $                 ,TMASK(1,1,1,1,IFIELD-1)
        !    $                 ,TMULT(1,1,1,1,IFIELD-1)
        !    $                 ,IMESH,TOLHT(IFIELD),NMXH,isd)

            if(iftmsh(ifield)) then
                call hsolve  (name4,TA,TB,H1,H2 &
                ,tmask(1,1,1,1,ifield-1) &
                ,tmult(1,1,1,1,ifield-1) &
                ,imesh,tolht(ifield),nmxh,1 &
                ,approx,napprox,bintm1)
            else
                call hsolve  (name4,TA,TB,H1,H2 &
                ,tmask(1,1,1,1,ifield-1) &
                ,tmult(1,1,1,1,ifield-1) &
                ,imesh,tolht(ifield),nmxh,1 &
                ,approx,napprox,binvm1)
            endif

            call add2    (t(1,1,1,1,ifield-1),ta,ntot)

            call cvgnlps (ifconv) ! Check convergence for nonlinear problem
            if (ifconv) goto 2000

        !        Radiation case, smooth convergence, avoid flip-flop (ER).
            CALL CMULT (TA,0.5,NTOT)
            CALL SUB2  (T(1,1,1,1,IFIELD-1),TA,NTOT)

        1000 END DO
        2000 CONTINUE
        CALL BCNEUSC (TA,1)
        CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TA,NTOT) ! no idea why... pf

    endif

    return
    end subroutine cdscal

!-----------------------------------------------------------------------
    subroutine makeuq

!     Fill up user defined forcing function and collocate will the
!     mass matrix on the Gauss-Lobatto mesh.

    use size_m
    use mass
    use soln
    use tstep

    ntot = nx1*ny1*nz1*nelfld(ifield)

    time = time-dt        ! Set time to t^n-1 for user function

    call rzero   ( bq(1,1,1,1,ifield-1) ,    ntot)
    call setqvol ( bq(1,1,1,1,ifield-1)          )
    call col2    ( bq(1,1,1,1,ifield-1) ,bm1,ntot)

    time = time+dt        ! Restore time

    return
    end subroutine makeuq
!-----------------------------------------------------------------------
    subroutine setqvol(bql)

!     Set user specified volumetric forcing function (e.g. heat source).

    use size_m
    use input
    use soln
    use tstep

    real :: bql(lx1*ly1*lz1,lelt)

#ifndef MOAB
    nel   = nelfld(ifield)
    nxyz1 = nx1*ny1*nz1
    ntot1 = nxyz1*nel

    do iel=1,nel
        igrp = igroup(iel)
        if (matype(igrp,ifield) == 1) then ! constant source within a group
            cqvol = cpgrp(igrp,ifield,3)
            call cfill (bql(1,iel),cqvol,nxyz1)
        else  !  pff 2/6/96 ............ default is to look at userq
            call nekuq (bql,iel)
        endif
    enddo

! 101 FORMAT(' Wrong material type (',I3,') for group',I3,', field',I2
!    $    ,/,' Aborting in SETQVOL.')
#else
! pulling in temperature right now, since we dont have anything else
    call userq2(bql)
#endif

    return
    end subroutine setqvol

    subroutine nekuq (bql,iel)
!------------------------------------------------------------------

!     Generate user-specified volumetric source term (temp./p.s.)

!------------------------------------------------------------------
    use size_m
    use input
    use mass
    use nekuse
    use parallel
    use soln
    use tstep

    real :: bql(lx1,ly1,lz1,lelt)

    ielg = lglel(iel)
    do 10 k=1,nz1
        do 10 j=1,ny1
            do 10 i=1,nx1
                call nekasgn (i,j,k,iel)
                qvol = 0.0
                call userq   (i,j,k,ielg)
                bql(i,j,k,iel) = qvol
    10 END DO

    return
    end subroutine nekuq
!-----------------------------------------------------------------------
    subroutine convab
!---------------------------------------------------------------

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

!---------------------------------------------------------------
    use size_m
    use mass
    use soln
    use tstep

    COMMON /SCRUZ/ TA (LX1,LY1,LZ1,LELT)

    NEL = NELFLD(IFIELD)
    NTOT1 = NX1*NY1*NZ1*NEL
    CALL CONVOP  (TA,T(1,1,1,1,IFIELD-1))
    CALL COL2    (TA,VTRANS(1,1,1,1,IFIELD),NTOT1)
    CALL SUBCOL3 (BQ(1,1,1,1,IFIELD-1),BM1,TA,NTOT1)

    return
    end subroutine convab
!-----------------------------------------------------------------------
    subroutine makeabq

!     Sum up contributions to 3rd order Adams-Bashforth scheme.

    use size_m
    use soln
    use tstep

    COMMON /SCRUZ/ TA (LX1,LY1,LZ1,LELT)

    AB0   = AB(1)
    AB1   = AB(2)
    AB2   = AB(3)
    NEL   = NELFLD(IFIELD)
    NTOT1 = NX1*NY1*NZ1*NEL

    CALL ADD3S2 (TA,VGRADT1(1,1,1,1,IFIELD-1), &
    VGRADT2(1,1,1,1,IFIELD-1),AB1,AB2,NTOT1)
    CALL COPY   (   VGRADT2(1,1,1,1,IFIELD-1), &
    VGRADT1(1,1,1,1,IFIELD-1),NTOT1)
    CALL COPY   (   VGRADT1(1,1,1,1,IFIELD-1), &
    BQ(1,1,1,1,IFIELD-1),NTOT1)
    CALL ADD2S1 (BQ(1,1,1,1,IFIELD-1),TA,AB0,NTOT1)

    return
    end subroutine makeabq
!-----------------------------------------------------------------------
    subroutine makebdq
!-----------------------------------------------------------------------

!     Add contributions to F from lagged BD terms.

!-----------------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use soln
    use tstep

    COMMON /SCRNS/ TA (LX1,LY1,LZ1,LELT) &
    ,             TB (LX1,LY1,LZ1,LELT) &
    ,             H2 (LX1,LY1,LZ1,LELT)

    NEL   = NELFLD(IFIELD)
    NTOT1 = NX1*NY1*NZ1*NEL
    CONST = 1./DT
    CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
    CALL CMULT (H2,CONST,NTOT1)

    CALL COL3  (TB,BM1,T(1,1,1,1,IFIELD-1),NTOT1)
    CALL CMULT (TB,BD(2),NTOT1)

    DO 100 ILAG=2,NBD
        IF (IFGEOM) THEN
            CALL COL3 (TA,BM1LAG(1,1,1,1,ILAG-1), &
            TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
        ELSE
            CALL COL3 (TA,BM1, &
            TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
        ENDIF
        CALL CMULT (TA,BD(ILAG+1),NTOT1)
        CALL ADD2  (TB,TA,NTOT1)
    100 END DO

    CALL COL2 (TB,H2,NTOT1)
    CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TB,NTOT1)

    return
    end subroutine makebdq
!-----------------------------------------------------------------------
    subroutine thyprk (tch,ilag)

!     Convection of a passive scalar.
!     Runge-Kutta scheme.

    use size_m
    use mass
    use soln
    use tstep

    REAL ::           TCH   (LX1,LY1,LZ1,1)
    COMMON /SCRNS/ VXN   (LX1,LY1,LZ1,LELV) &
    ,             VYN   (LX1,LY1,LZ1,LELV) &
    ,             VZN   (LX1,LY1,LZ1,LELV) &
    ,             HTMASK(LX1,LY1,LZ1,LELT) &
    ,             WORK  (LX1,LY1,LZ1,LELT)
    COMMON /CTMP1/ RK1   (LX1,LY1,LZ1,LELT) &
    ,             RK2   (LX1,LY1,LZ1,LELT) &
    ,             RK3   (LX1,LY1,LZ1,LELT) &
    ,             RK4   (LX1,LY1,LZ1,LELT)

    NTOT1 = NX1*NY1*NZ1*NELV

    CALL OPCOPY  (VXN,VYN,VZN,VX,VY,VZ)
    CALL HYPMSK1 (HTMASK)
    CALL TAUINIT (TAU,ILAG)
    CALL TCHINIT (TCH,ILAG)
    CALL VELCONV (VXN,VYN,VZN,TAU)

    DO 1000 JLAG=ILAG,1,-1
    
        DTAU   = DTLAG(JLAG)/(NTAUBD)
        DTHALF = DTAU/2.
        CRK1   = DTAU/6.
        CRK2   = DTAU/3.
    
        DO 900 ITAU=1,NTAUBD
        
            CALL FRKCONV (RK1,TCH,HTMASK)
        
            TAU = TAU + DTHALF
            CALL VELCONV (VXN,VYN,VZN,TAU)
        
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK1,-DTHALF,NTOT1)
            CALL FRKCONV (RK2,WORK,HTMASK)
        
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK2,-DTHALF,NTOT1)
            CALL FRKCONV (RK3,WORK,HTMASK)
        
            TAU = TAU + DTHALF
            CALL VELCONV (VXN,VYN,VZN,TAU)
        
            CALL COPY    (WORK,TCH,NTOT1)
            CALL ADD2S2  (WORK,RK3,-DTAU,NTOT1)
            CALL FRKCONV (RK4,WORK,HTMASK)
        
            CALL ADD2S2  (TCH,RK1,-CRK1,NTOT1)
            CALL ADD2S2  (TCH,RK2,-CRK2,NTOT1)
            CALL ADD2S2  (TCH,RK3,-CRK2,NTOT1)
            CALL ADD2S2  (TCH,RK4,-CRK1,NTOT1)
        
        900 END DO
    1000 END DO

    CALL OPCOPY (VX,VY,VZ,VXN,VYN,VZN)

    return
    end subroutine thyprk

!-----------------------------------------------------------------------
    subroutine hypmsk1 (htmask)

!     Generate mask-array for the hyperbolic system (passive scalar).

    use size_m
    use geom
    use input
    use soln
    use tstep
    REAL ::           HTMASK(LX1,LY1,LZ1,1)
    CHARACTER      CB*3
    PARAMETER (LXYZ1=LX1*LY1*LZ1)
    COMMON /CTMP1/ WORK   (LXYZ1,LELT)

    NFACES= 2*NDIM
    NEL   = NELFLD(IFIELD)
    NTOT1 = NX1*NY1*NZ1*NEL
    CALL RZERO(WORK  ,NTOT1)
    CALL RONE (HTMASK,NTOT1)

    DO 100 IE=1,NELV
        DO 100 IFACE=1,NFACES
            CB=CBC(IFACE,IE,IFIELD)
            IF (CB == 'T  ' .OR. CB == 't  ' .OR. &
            CB == 'KD ' .OR. CB == 'kd ' .OR. &
            CB == 'ED ' .OR. CB == 'ed ') THEN
                CALL FACCL3 (WORK(1,IE),VX(1,1,1,IE),UNX(1,1,IFACE,IE),IFACE)
                CALL FADDCL3(WORK(1,IE),VY(1,1,1,IE),UNY(1,1,IFACE,IE),IFACE)
                IF (IF3D) &
                CALL FADDCL3(WORK(1,IE),VZ(1,1,1,IE),UNZ(1,1,IFACE,IE),IFACE)
                CALL FCAVER (TAVER,WORK,IE,IFACE)
            
                IF (TAVER < 0.) CALL FACEV (HTMASK,IE,IFACE,0.0,NX1,NY1,NZ1)
                IF (CB == 'KWS' .OR. CB == 'EWS') &
                CALL FACEV (HTMASK,IE,IFACE,0.0,NX1,NY1,NZ1)
            ENDIF
    100 END DO

    return
    end subroutine hypmsk1

!-----------------------------------------------------------------------
    subroutine tchinit (tch,ilag)

!     Set initial conditions for subintegration

    use size_m
    use soln
    use tstep
    REAL :: TCH (LX1,LY1,LZ1,1)

    NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)
    IF (ILAG == 1) THEN
        CALL COPY (TCH,T(1,1,1,1,IFIELD-1),NTOT1)
    ELSE
        CALL COPY (TCH,TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
    ENDIF
    return
    end subroutine tchinit

!-----------------------------------------------------------------------
    subroutine lagscal
!-----------------------------------------------------------------------

!     Keep old passive scalar field(s)

!-----------------------------------------------------------------------
    use size_m
    use input
    use soln
    use tstep

    NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)

    DO 100 ILAG=NBDINP-1,2,-1
        CALL COPY (TLAG(1,1,1,1,ILAG  ,IFIELD-1), &
        TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
    100 END DO

    CALL COPY (TLAG(1,1,1,1,1,IFIELD-1),T(1,1,1,1,IFIELD-1),NTOT1)

    return
    end subroutine lagscal
!-----------------------------------------------------------------------
