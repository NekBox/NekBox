!-----------------------------------------------------------------------
    subroutine cdscal (igeom)

!     Solve the convection-diffusion equation for passive scalar IPSCAL

    use size_m
    use geom
    use input
    use mass
    use mvgeom
    include 'SOLN'
    include 'TSTEP'
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
    include 'SOLN'
    include 'TSTEP'

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
    include 'SOLN'
    include 'TSTEP'

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
    include 'PARALLEL'
    include 'SOLN'
    include 'TSTEP'

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
    include 'SOLN'
    include 'TSTEP'

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
    include 'SOLN'
    include 'TSTEP'

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
    include 'SOLN'
    include 'TSTEP'

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
    subroutine convch_old

!     Compute convective contribution using
!     operator-integrator-factor method (characteristics).

    use size_m
    use mass
    include 'SOLN'
    include 'TSTEP'

    COMMON /SCRUZ/ TCH (LX1,LY1,LZ1,LELT) &
    ,             H2  (LX1,LY1,LZ1,LELT)
    LOGICAL :: IFCHAB

    IFCHAB = .FALSE. 
    NEL    = NELFLD(IFIELD)
    NTOT1  = NX1*NY1*NZ1*NEL
    CONST  = 1./DT
    CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
    CALL CMULT (H2,CONST,NTOT1)

    DO 100 ILAG=NBD,1,-1
        IF (IFCHAB) THEN
            CALL THYPAB (TCH,ILAG)
        ELSE
            CALL THYPRK (TCH,ILAG)
        ENDIF
        CALL COL2  (TCH,BM1,NTOT1)
        CALL COL2  (TCH,H2 ,NTOT1)
        CALL CMULT (TCH,BD(ILAG+1),NTOT1)
        CALL ADD2  (BQ(1,1,1,1,IFIELD-1),TCH,NTOT1)
    100 END DO

    return
    end subroutine convch_old

!-----------------------------------------------------------------------
    subroutine thyprk (tch,ilag)

!     Convection of a passive scalar.
!     Runge-Kutta scheme.

    use size_m
    use mass
    include 'SOLN'
    include 'TSTEP'

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
    subroutine thypab (tch,ilag)

!     Convection of a passive scalar.
!     Adams-Bashforth.

    use size_m
    use mass
    include 'SOLN'
    include 'TSTEP'

    REAL ::           TCH   (LX1,LY1,LZ1,1)
    COMMON /SCRNS/ VXN   (LX1,LY1,LZ1,LELV) &
    ,             VYN   (LX1,LY1,LZ1,LELV) &
    ,             VZN   (LX1,LY1,LZ1,LELV) &
    ,             HTMASK(LX1,LY1,LZ1,LELT) &
    ,             WORK  (LX1,LY1,LZ1,LELT)
    COMMON /CTMP1/ TMP1  (LX1,LY1,LZ1,LELT) &
    ,             TMP2  (LX1,LY1,LZ1,LELT)

    NEL   = NELFLD(IFIELD)
    NTOT1 = NX1*NY1*NZ1*NEL

    CALL OPCOPY  (VXN,VYN,VZN,VX,VY,VZ)
    CALL HYPMSK1 (HTMASK)
    CALL TAUINIT (TAU,ILAG)
    CALL TCHINIT (TCH,ILAG)
    CALL VELCONV (VXN,VYN,VZN,TAU)

    DT2   = 0.
    DT1   = 0.

    DO 1000 JLAG=ILAG,1,-1
    
        DT2 = DT1
        DT1 = DTAU
        DTAU   = DTLAG(JLAG)/(NTAUBD)
    
        IF (JLAG == ILAG) THEN
            DTAU1 = 0.1*DTAU**1.5
            DTAU3 = (DTLAG(JLAG)-2.*DTAU1)/(NTAUBD-2)
        ELSE
            DTAU1 = DTAU
            DTAU3 = DTAU
        ENDIF
    
        DO 900 ITAU=1,NTAUBD
        
            IF (ITAU > 1) THEN
                DT2 = DT1
                DT1 = DTAU
            ENDIF
            IF ((JLAG == ILAG) .AND. (ITAU <= 2)) THEN
                AB1  = 1.
                AB2  = 0.
                AB3  = 0.
                DTAU = DTAU1
            ELSE
                DTAU = DTAU3
                AB3  = (DTAU/DT2)*((DTAU/3.+DT1/2.)/(DT1+DT2))
                AB2  = -(DTAU/DT1)*(0.5+(DTAU/3.+DT1/2.)/DT2)
                AB1  = 1.-AB2-AB3
            ENDIF
            TAU = TAU + DTAU
        
            CALL ADD3S2  (WORK,TMP1,TMP2,AB2,AB3,NTOT1)
            CALL COPY    (TMP2,TMP1,NTOT1)
            CALL FRKCONV (TMP1,TCH,HTMASK)
            CALL ADD2S2  (WORK,TMP1,AB1,NTOT1)
            CALL ADD2S2  (TCH,WORK,-DTAU,NTOT1)
        
            CALL VELCONV (VXN,VYN,VZN,TAU)
        
        900 END DO
    1000 END DO

    CALL OPCOPY (VX,VY,VZ,VXN,VYN,VZN)

    return
    end subroutine thypab

!-----------------------------------------------------------------------
    subroutine hypmsk1 (htmask)

!     Generate mask-array for the hyperbolic system (passive scalar).

    use size_m
    use geom
    use input
    include 'SOLN'
    include 'TSTEP'
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
    include 'SOLN'
    include 'TSTEP'
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
    include 'SOLN'
    include 'TSTEP'

    NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)

    DO 100 ILAG=NBDINP-1,2,-1
        CALL COPY (TLAG(1,1,1,1,ILAG  ,IFIELD-1), &
        TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
    100 END DO

    CALL COPY (TLAG(1,1,1,1,1,IFIELD-1),T(1,1,1,1,IFIELD-1),NTOT1)

    return
    end subroutine lagscal
!-----------------------------------------------------------------------
    subroutine outfldrq (x,txt10,ichk)
    use size_m
    INCLUDE 'TSTEP'
    real :: x(nx1,ny1,nz1,lelt)
    character(10) :: txt10

    integer :: idum,e
    save idum
    data idum /-3/

    if (idum < 0) return


    mtot = nx1*ny1*nz1*nelv
    if (nx1 > 8 .OR. nelv > 16) return
    xmin = glmin(x,mtot)
    xmax = glmax(x,mtot)

    nell = nelt
    rnel = nell
    snel = sqrt(rnel)+.1
    ne   = snel
    ne1  = nell-ne+1
    k = 1
    do ie=1,1
        ne = 0
        write(6,116) txt10,k,ie,xmin,xmax,istep,time
        do l=0,1
            write(6,117)
            do j=ny1,1,-1
                if (nx1 == 2) write(6,102) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 3) write(6,103) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 4) write(6,104) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 5) write(6,105) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 6) write(6,106) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 7) write(6,107) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
                if (nx1 == 8) write(6,118) ((x(i,j,k,e+l),i=1,nx1),e=1,1)
            enddo
        enddo
    enddo

    102 FORMAT(4(2f9.5,2x))
    103 FORMAT(4(3f9.5,2x))
    104 FORMAT(4(4f7.3,2x))
    105 FORMAT(5f9.5,10x,5f9.5)
    106 FORMAT(6f9.5,5x,6f9.5)
    107 FORMAT(7f8.4,5x,7f8.4)
    108 FORMAT(8f8.4,4x,8f8.4)
    118 FORMAT(8f12.9)

    116 FORMAT(  /,5X,'     ^              ',/, &
    &     5X,'   Y |              ',/, &
    &     5X,'     |              ',A10,/, &
    &     5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/, &
    &     5X,'       X            ','Step  =',I9,f15.5)
    117 FORMAT(' ')

    if (ichk == 1 .AND. idum > 0) call checkit(idum)
    return
    end subroutine outfldrq
!-----------------------------------------------------------------------
    subroutine cdscal_expl (igeom)

!     explicit convection-diffusion equation for passive scalar

    use size_m
    use geom
    use input
    use mass
    use mvgeom
    include 'SOLN'
    include 'TSTEP'
    common  /cprint/ ifprint
    logical ::          ifprint
    logical ::          ifconv

    common /scrns/ ta(lx1,ly1,lz1,lelt) &
    ,tb(lx1,ly1,lz1,lelt)
    common /scrvh/ h1(lx1,ly1,lz1,lelt) &
    ,h2(lx1,ly1,lz1,lelt)


!     QUESTIONABLE support for Robin BC's at this point! (5/15/08)

    nel    = nelfld(ifield)
    ntot   = nx1*ny1*nz1*nel

    if (igeom == 1) then   ! geometry at t^{n-1}

        call makeq
        call lagscal

    else                   ! geometry at t^n

        if ( .TRUE. .AND. nid == 0) &
        write (6,*) istep,ifield,' explicit step'


    !        New geometry

        isd = 1
        if (ifaxis .AND. ifmhd) isd = 2 !This is a problem if T is to be T!

        intype = 0
        if (iftran) intype = -1
        call sethlm  (h1,h2,intype)

        call bcneusc (ta,-1)       ! Modify diagonal for Robin condition
        call add2    (h2,ta ,ntot)
        call col2    (h2,BM1,ntot)

        call bcneusc (tb,1)        ! Modify rhs for flux bc
        call add2    (bq(1,1,1,1,ifield-1),tb,ntot)

        call dssum   (bq(1,1,1,1,ifield-1),nx1,ny1,nz1)
        call dssum   (h2,nx1,ny1,nz1)

        call invcol3 (t(1,1,1,1,ifield-1),bq(1,1,1,1,ifield-1),h2,ntot)

        call bcdirsc (t(1,1,1,1,ifield-1)) ! --> no mask needed

    endif                   ! geometry at t^n

    return
    end subroutine cdscal_expl
!-----------------------------------------------------------------------
    subroutine diffab  ! explicit treatment of diffusion operator

!     Eulerian scheme, add diffusion term to forcing function
!     at current time step.


    use size_m
    use input
    use mass
    include 'SOLN'
    include 'TSTEP'

    common /scruz/ ta(lx1,ly1,lz1,lelt) &
    ,h2(lx1,ly1,lz1,lelt)

    nel = nelfld(ifield)
    ntot1 = nx1*ny1*nz1*nel

    intype = 0
    if (iftran) intype = -1

    isd = 1
    if (ifaxis .AND. ifmhd) isd = 2 !This is a problem if T is to be T!

    imesh = 1
!      if (iftmsh(ifield)) imesh=2

    call rzero   (h2,ntot1)
    call axhelm  (ta,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield) &
    ,h2,imesh,isd)
    call sub2    (bq(1,1,1,1,ifield-1),ta,ntot1)

    return
    end subroutine diffab
!-----------------------------------------------------------------------
