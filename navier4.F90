!-----------------------------------------------------------------------
#if 0
    subroutine setrhs(p,h1,h2,h2inv)

!     Project rhs onto best fit in the "E" norm.

    use size_m
    use input
    use mass
    use soln
    use tstep

    REAL ::             P    (LX2,LY2,LZ2,LELV)
    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    logical :: ifdump
    save    ifdump
    data    ifdump / .FALSE. /

    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev
    REAL :: ALPHA,WORK


    integer :: icalld
    save    icalld
    data    icalld/0/

!     First call, we have no vectors to orthogonalize against.
    IF (ICALLD == 0) THEN
        icalld=icalld+1
        Nprev=0
        Mprev=param(93)
        Mprev=min(Mprev,Mxprev)
    ENDIF

!     Diag to see how much reduction in the residual is attained.

    NTOT2  = NX2*NY2*NZ2*NELV
    ALPHA1 = GLSC3(p,p,bm2inv,NTOT2)
    if (alpha1 > 0) ALPHA1 = sqrt(alpha1/volvm2)

!     Update rhs's if E-matrix has changed

    CALL UPDRHSE(P,H1,H2,H2INV,ierr)
    if (ierr == 1) Nprev=0

!     Perform Gram-Schmidt for previous rhs's.

    DO 10 I=1,Nprev
        ALPHA(i) = VLSC2(P,RHS(1,i),NTOT2)
    10 END DO

    IF (Nprev > 0) CALL gop(alpha,WORK,'+  ',Nprev)

    CALL RZERO(Pbar,NTOT2)
    DO 20 I=1,Nprev
        alphas = alpha(i)
        CALL ADD2S2(Pbar,RHS(1,i),alphas,NTOT2)
    20 END DO

    if (Nprev > 0) then
        INTETYPE = 1
        CALL CDABDTP(Pnew,Pbar,H1,H2,H2INV,INTETYPE)
        CALL SUB2   (P,Pnew,NTOT2)
    !    ................................................................
    !      Diag.
        ALPHA2 = GLSC3(p,p,bm2inv,NTOT2)
        if (alpha2 > 0) ALPHA2 = sqrt(alpha2/volvm2)
        ratio  = alpha1/alpha2
        n10=min(10,nprev)
        IF (NID == 0) WRITE(6,11)ISTEP,Nprev,(ALPHA(I),I=1,N10)
        IF (NID == 0) WRITE(6,12) ISTEP,nprev,ALPHA1,ALPHA2,ratio
        11 FORMAT(2I5,' alpha:',1p10e12.4)
        12 FORMAT(I6,i4,1p3e12.4,' alph12')
    
    !        if (alpha1.gt.0.001 .and. .not.ifdump) then
    !           IF (NID.EQ.0) WRITE(6,*) 'alph1 large ... '
    !           if (istep.gt.10) then
    !              IF (NID.EQ.0) WRITE(6,*) ' ... dumping'
    !              call prepost (.true.,'   ')
    !              ifdump = .true.
    !           else
    !              IF (NID.EQ.0) WRITE(6,*) ' ... doing nothing'
    !           endif
    !        endif
    !        if (alpha1.gt.0.1.and.istep.gt.10) then
    !           IF (NID.EQ.0) WRITE(6,*) 'alph1 too large ... aborting'
    !           call prepost (.true.,'   ')
    !           call exitt
    !        endif
    !    ................................................................
    endif


    RETURN
    end subroutine setrhs
!-----------------------------------------------------------------------
    subroutine gensoln(p,h1,h2,h2inv)

!     Reconstruct the solution to the original problem by adding back
!     the previous solutions
!     know the soln.

    use size_m
    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev
    REAL :: ALPHA,WORK

    REAL ::             P    (LX2,LY2,LZ2,LELV)
    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    NTOT2=NX2*NY2*NZ2*NELV

!     First, save current solution

    CALL COPY (Pnew,P,NTOT2)

!     Reconstruct solution

    CALL ADD2(P,Pbar,NTOT2)

!     Update the set of <p,rhs>

    CALL UPDTSET(P,H1,H2,H2INV,ierr)
    if (ierr == 1) Nprev = 0

    RETURN
    end subroutine gensoln
!-----------------------------------------------------------------------
    subroutine updtset(p,h1,h2,h2inv,IERR)

!     Update the set of rhs's and the corresponding p-set:

!        . Standard case is to add P_new, and RHS_new = E*P_new

!        . However, when Nprev=Mprev (max. allowed), we throw out
!          the old set, and set P_1 = P, RHS_1=E*P_1

!        . Other schemes are possible, e.g., let's save a bunch of
!          old vectors, perhaps chosen wisely via P.O.D.


    use size_m
    use input
    use mass
    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev

    REAL :: ALPHA,WORK

    REAL ::             P    (LX2,LY2,LZ2,LELV)
    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    NTOT2=NX2*NY2*NZ2*NELV

    IF (Nprev == Mprev) THEN
        CALL COPY(Pnew,P,NTOT2)
        Nprev=0
    ENDIF

!     Increment solution set
    Nprev = Nprev+1

    CALL COPY   (RHS(1,Nprev),Pnew,NTOT2)

!     Orthogonalize rhs against previous rhs and normalize

    CALL ECONJ (Nprev,H1,H2,H2INV,ierr)
!     CALL ECHECK(Nprev,H1,H2,H2INV,INTETYPE)

!     Save last sol'n
    CALL COPY(Pnew,P,NTOT2)

    RETURN
    end subroutine updtset
#endif
!-----------------------------------------------------------------------
    subroutine econj(kprev,h1,h2,h2inv,ierr)

!     Orthogonalize the rhs wrt previous rhs's for which we already
!     know the soln.

    use size_m
    use input
    use mass
    use soln
    use tstep

    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),Pbrr(ltot2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev
    REAL :: ALPHA,WORK
    real :: ALPHAd


    ierr  = 0
    NTOT2 = NX2*NY2*NZ2*NELV
    INTETYPE=1

!     Gram Schmidt, w re-orthogonalization

    npass=1
    if (abs(param(105)) == 2) npass=2
    do ipass=1,npass
    
        CALL CDABDTP(Pbrr,RHS(1,Kprev),H1,H2,H2INV,INTETYPE)
    
    !        Compute part of the norm
        Alphad = GLSC2(RHS(1,Kprev),Pbrr,NTOT2)
    
    !        Gram-Schmidt
        Kprev1=Kprev-1
        DO 10 I=1,Kprev1
            ALPHA(I) = VLSC2(Pbrr,RHS(1,i),NTOT2)
        10 END DO
        IF (Kprev1 > 0) CALL gop(alpha,WORK,'+  ',Kprev1)
    
        DO 20 I=1,Kprev1
            alpham = -alpha(i)
            CALL ADD2S2(RHS(1,Kprev),RHS (1,i),alpham,NTOT2)
            Alphad = Alphad - alpha(i)**2
        20 END DO
    enddo

!    .Normalize new element in P~

    if (ALPHAd <= 0.0) then
        write(6,*) 'ERROR:  alphad <= 0 in ECONJ',alphad,Kprev
        ierr = 1
        return
    endif
    ALPHAd = 1.0/SQRT(ALPHAd)
    ALPHAN = Alphad
    CALL CMULT(RHS (1,Kprev),alphan,NTOT2)

    RETURN
    end subroutine econj
!-----------------------------------------------------------------------
    subroutine chkptol
!--------------------------------------------------------------------

!     Check pressure tolerance for transient case.

!     pff 6/20/92
!     This routine has been modified for diagnostic purposes only.
!     It can be replaced with the standard nekton version.

!--------------------------------------------------------------------
    use size_m
    use input
    use mass
    use soln
    use tstep
    COMMON /CTOLPR/ DIVEX
    COMMON /CPRINT/ IFPRINT
    LOGICAL ::         IFPRINT

    COMMON /SCRUZ/  DIVV (LX2,LY2,LZ2,LELV) &
    ,              BDIVV(LX2,LY2,LZ2,LELV)

    if (ifsplit) return
    IF (param(102) == 0 .AND. (TOLPDF /= 0. .OR. ISTEP <= 5)) RETURN
    five = 5.0
    if (param(102) /= 0.0) five=param(102)

    NTOT2 = NX2*NY2*NZ2*NELV
    if (ifield == 1) then     ! avo: sub arguments?
        CALL OPDIV (BDIVV,VX,VY,VZ)
    else
        CALL OPDIV (BDIVV,BX,BY,BZ)
    endif
    CALL COL3 (DIVV,BDIVV,BM2INV,NTOT2)
    DNORM = SQRT(GLSC2(DIVV,BDIVV,NTOT2)/VOLVM2)

    if (nid == 0) WRITE (6,*) istep,' DNORM, DIVEX',DNORM,DIVEX

!     IF (istep.gt.10.and.DNORM.GT.(1.01*DIVEX).AND.DIVEX.GT.0.) then
!        if (DNORM.gt.1e-8) then
!           if (nid.eq.0) WRITE(6,*) 'DNORM-DIVEX div. ... aborting'
!           call prepost (.true.,'   ')
!           call exitt
!        else
!           if (nid.eq.0) WRITE(6,*) 'DNORM-DIVEX div. ... small'
!        endif
!     endif

!     IF (DNORM.GT.(1.2*DIVEX).AND.DIVEX.GT.0.) TOLPDF = 5.*DNORM
    IF (istep > 5 .AND. tolpdf == 0.0 .AND. &
    DNORM > (1.2*DIVEX) .AND. DIVEX > 0.) &
    TOLPDF = FIVE*DNORM

    RETURN
    end subroutine chkptol
    FUNCTION VLSC3(X,Y,B,N)
    use opctr

!     local inner product, with weight

    DIMENSION X(1),Y(1),B(1)
    REAL :: DT


    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'VLSC3 '
    endif
    isbcnt = 3*n
    dct(myrout) = dct(myrout) + dfloat(isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + dfloat(isbcnt)

    DT = 0.0
    DO 10 I=1,N
        T = X(I)*Y(I)*B(I)
        DT = DT+T
    10 END DO
    T=DT
    VLSC3 = T
    RETURN
    END FUNCTION VLSC3
!-----------------------------------------------------------------------
#if 0
    subroutine updrhse(p,h1,h2,h2inv,ierr)

!     Update rhs's if E-matrix has changed


    use size_m
    use input
    use mass
    use tstep

    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev
    COMMON /ORTHOL/ IFNEWE
    REAL :: ALPHA,WORK
    LOGICAL :: IFNEWE


    REAL ::             P    (LX2,LY2,LZ2,LELV)
    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    integer :: icalld
    save    icalld
    data    icalld/0/
    NTOT2=NX2*NY2*NZ2*NELV


!     First, we have to decide if the E matrix has changed.

    IF (icalld == 0) THEN
        icalld=1
        DTlast=DT
    ENDIF

    IFNEWE= .FALSE. 
    IF (IFMVBD) THEN
        IFNEWE= .TRUE. 
        CALL INVERS2(bm2inv,bm2,Ntot2)
    ELSEIF (DTlast /= DT) THEN
        IFNEWE= .TRUE. 
        DTlast=DT
    ENDIF
    IF (IFNEWE .AND. nid == 0) write(6,*) 'reorthogo:',nprev



!     Next, we reconstruct a new rhs set.

    IF (IFNEWE) THEN
    
    !        new idea...
        if (nprev > 0) nprev=1
        call copy(rhs,pnew,ntot2)
    
        Nprevt = Nprev
        DO 100 Iprev=1,Nprevt
        !           Orthogonalize this rhs w.r.t. previous rhs's
            CALL ECONJ (Iprev,H1,H2,H2INV,ierr)
            if (ierr == 1) then
                Nprev = 0
                return
            endif
        100 END DO
    
    ENDIF

    RETURN
    end subroutine updrhse
!-----------------------------------------------------------------------
    subroutine echeck(kprev,h1,h2,h2inv,intetype)

!     Orthogonalize the rhs wrt previous rhs's for which we already
!     know the soln.

    use size_m
    use input
    use mass
    use soln
    use tstep

    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)

    PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
    COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
    COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
    COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
    COMMON /ORTHOI/ Nprev,Mprev
    REAL :: ALPHA,WORK,GLSC2
    REAL :: Alphad


    NTOT2=NX2*NY2*NZ2*NELV

!     Compute part of the norm

    do 20 j=1,kprev
        CALL CDABDTP(Pbar,rhs(1,j),H1,H2,H2INV,INTETYPE)
        do 10 i=1,kprev
            Alphad = GLSC2(RHS(1,i),Pbar,NTOT2)
            Alphas = alphad
            if (nid == 0) then
                write(6,5) i,j,alphad,alphas,istep,kprev
                5 format(' E-check:',2i4,e16.8,g12.5,i6,i4)
            endif
        10 END DO
    20 END DO
    return
    end subroutine echeck
#endif
!-----------------------------------------------------------------------
!     THE ROUTINES BELOW ARE THE NEW Helmholtz projectors
!-----------------------------------------------------------------------
    subroutine projh(r,h1,h2,bi,vml,vmk,approx,napprox,wl,ws,name4)

!     Orthogonalize the rhs wrt previous rhs's for which we already
!     know the soln.

!     Input:   r         -- residual
!              h1,h2     -- Helmholtz arrays
!              bi        -- inverse mass matrix
!              vml,vmk   -- multiplicity and mask arrays
!              approx    -- approximation space
!              napprox   -- (1) = max vecs,  (2) = current number of vecs
!              wl        -- large work array of size lx1*ly1*lz1*nelv
!              ws        -- small work array of size 2*max vecs

    use size_m
    use input
    use mass
    use soln
    use tstep

    parameter(lt=lx1*ly1*lz1*lelt)

    real :: r(1),h1(1),h2(1),vml(1),vmk(1)
    real :: bi(1)
    real :: wl(1),ws(1)
    real :: approx(lt,0:1)
    integer :: napprox(2)
    character(4) :: name4

    n_max = napprox(1)
    n_sav = napprox(2)
    if (n_sav == 0) return
    nel =nelfld(ifield)
    ntot=nx1*ny1*nz1*nel

    vol = voltm1
    if (nel == nelv) vol = volvm1

!     Diag to see how much reduction in the residual is attained.

    alpha1 = glsc23(r,bi,vml,ntot)
    if (alpha1 > 0) alpha1 = sqrt(alpha1/vol)

!     Update approximation space if dt has changed
    call updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)


!     Perform Gram-Schmidt for previous soln's

    do i=1,n_sav
        ws(i) = vlsc3(r,approx(1,i),vml,ntot)
    enddo
    call gop    (ws,ws(n_sav+1),'+  ',n_sav)

    call cmult2   (approx(1,0),approx(1,1),ws(1),ntot)
    do i=2,n_sav
        call add2s2(approx(1,0),approx(1,i),ws(i),ntot)
    enddo

    call axhelm  (wl,approx(1,0),h1,h2,1,1)
    call col2    (wl,vmk,ntot)
    call dssum   (wl,nx1,ny1,nz1)
    call sub2    (r ,wl,ntot)
! ................................................................
!   Diag.
    alpha2 = glsc23(r,bi,vml,ntot)
    if (alpha2 > 0) alpha2 = sqrt(alpha2/vol)
    ratio  = alpha1/alpha2
    n10=min(10,n_sav)

    if (nid == 0) write(6,10) istep,name4,alpha1,alpha2,ratio,n_sav
    10 format(4X,I7,4x,a4,' alph1n',1p3e12.4,i6)

    if (nid == 0) write(6,11) istep,name4,n_sav,(ws(i),i=1,n10)
    11 format(4X,I7,4x,a4,' halpha',i6,10(1p10e12.4,/,17x))

    return
    end subroutine projh
!-----------------------------------------------------------------------
    subroutine gensh(v1,h1,h2,vml,vmk,approx,napprox,wl,ws,name4)

!     Reconstruct the solution to the original problem by adding back
!     the previous solutions

    use size_m
    use input
    use mass
    use soln
    use tstep

    common /iterhm/ niterhm

    parameter(lt=lx1*ly1*lz1*lelt)
    real :: v1(1),h1(1),h2(1),vml(1),vmk(1)
    real :: wl(1),ws(1)
    real :: approx(lt,0:1)
    integer :: napprox(2)
    character(4) :: name4

    n_max = napprox(1)
    n_sav = napprox(2)
    ntot=nx1*ny1*nz1*nelfld(ifield)

!     Reconstruct solution and save current du

    if (n_sav < n_max) then
    
        if (niterhm > 0) then      ! new vector not in space
            n_sav = n_sav+1
            call copy(approx(1,n_sav),v1,ntot)
            call add2(v1,approx(1,0),ntot)
        !           orthogonalize rhs against previous rhs and normalize
            call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)

        !           if (ierr.ne.0) n_sav = n_sav-1
            if (ierr /= 0) n_sav = 0

        else

            call add2(v1,approx(1,0),ntot)

        endif
    else
        n_sav = 1
        call add2(v1,approx(1,0),ntot)
        call copy(approx(1,n_sav),v1,ntot)
    !        normalize
        call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)
        if (ierr /= 0) n_sav = 0
    endif

    napprox(2)=n_sav

    return
    end subroutine gensh
!-----------------------------------------------------------------------
    subroutine hconj(approx,k,h1,h2,vml,vmk,ws,name4,ierr)

!     Orthonormalize the kth vector against vector set

    use size_m
    use input
    use mass
    use parallel
    use soln
    use tstep

    parameter  (lt=lx1*ly1*lz1*lelt)
    real :: approx(lt,0:1),h1(1),h2(1),vml(1),vmk(1),ws(1)
    character(4) :: name4

    ierr=0
    ntot=nx1*ny1*nz1*nelfld(ifield)

    call axhelm  (approx(1,0),approx(1,k),h1,h2,1,1)
    call col2    (approx(1,0),vmk,ntot)
    call dssum   (approx(1,0),nx1,ny1,nz1)
    call col2    (approx(1,0),vml        ,ntot)

!     Compute part of the norm   (Note:  a(0) already scaled by vml)

    alpha = glsc2(approx(1,0),approx(1,k),ntot)
    alph1 = alpha

!     Gram-Schmidt

    km1=k-1
    do i=1,km1
        ws(i) = vlsc2(approx(1,0),approx(1,i),ntot)
    enddo
    if (km1 > 0) call gop(ws,ws(k),'+  ',km1)

    do i=1,km1
        alpham = -ws(i)
        call add2s2(approx(1,k),approx(1,i),alpham,ntot)
        alpha = alpha - ws(i)**2
    enddo

!    .Normalize new element in approximation space

    eps = 1.e-7
    if (wdsize == 8) eps = 1.e-15
    ratio = alpha/alph1

    if (ratio <= 0) then
        ierr=1
        if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
        12 format(I6,1x,a4,' alpha b4 sqrt:',i4,1p2e12.4)
    elseif (ratio <= eps) then
        ierr=2
        if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
    else
        ierr=0
        alpha = 1.0/sqrt(alpha)
        call cmult(approx(1,k),alpha,ntot)
    endif

    if (ierr /= 0) then
        call axhelm  (approx(1,0),approx(1,k),h1,h2,1,1)
        call col2    (approx(1,0),vmk,ntot)
        call dssum   (approx(1,0),nx1,ny1,nz1)
        call col2    (approx(1,0),vml        ,ntot)
    
    !        Compute part of the norm   (Note:  a(0) already scaled by vml)
    
        alpha = glsc2(approx(1,0),approx(1,k),ntot)
        if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
        if (alpha <= 0) then
            ierr=3
            if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
            return
        endif
        alpha = 1.0/sqrt(alpha)
        call cmult(approx(1,k),alpha,ntot)
        ierr = 0
    endif

    return
    end subroutine hconj
!-----------------------------------------------------------------------
    subroutine updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)

!     Reorthogonalize approx if dt has changed

    use size_m
    use input
    use mass
    use tstep


    parameter  (lt=lx1*ly1*lz1*lelt)
    real :: approx(lt,0:1),h1(1),h2(1),vml(1),vmk(1),ws(1)
    integer :: napprox(2)
    character(4) :: name4

    logical :: ifupdate
    logical :: ifnewdt
    save    ifnewdt
    data    ifnewdt / .FALSE. /

    character(4) :: name_old
    save    name_old
    data    name_old/'DMIR'/

    real ::    dtold
    save    dtold
    data    dtold/0.0/

!     First, we have to decide if the dt has changed.

    ifupdate = .FALSE. 
    if (dt /= dtold) then
        dtold    = dt
        name_old = name4
        ifnewdt  = .TRUE. 
        ifupdate = .TRUE. 
    elseif (ifnewdt) then
        if (name4 == name_old) then
            ifnewdt = .FALSE. 
        else
            ifupdate = .TRUE. 
        endif
    endif
    if (ifvarp(ifield)) ifupdate = .TRUE. 
    if (iflomach)       ifupdate = .TRUE. 

    if (ifupdate) then    ! reorthogonalize
        n_sav = napprox(2)
        l     = 1
        do k=1,n_sav
        !           Orthogonalize kth vector against {v_1,...,v_k-1}
            if (k /= l) then
                ntot = nx1*ny1*nz1*nelfld(ifield)
                call copy(approx(1,l),approx(1,k),ntot)
            endif
            call hconj(approx,l,h1,h2,vml,vmk,ws,name4,ierr)
            if (ierr == 0) l=l+1
        enddo
        napprox(2)=min(l,n_sav)
    endif

    return
    end subroutine updrhsh
!-----------------------------------------------------------------------
    subroutine hmhzpf(name,u,r,h1,h2,mask,mult,imesh,tli,maxit,isd,bi)
    use ctimer
    use size_m
    use fdmh1
    use input
    use mass

    CHARACTER(4) ::    NAME
    REAL ::           U    (LX1,LY1,LZ1,1)
    REAL ::           R    (LX1,LY1,LZ1,1)
    REAL ::           H1   (LX1,LY1,LZ1,1)
    REAL ::           H2   (LX1,LY1,LZ1,1)
    REAL ::           MASK (LX1,LY1,LZ1,1)
    REAL ::           MULT (LX1,LY1,LZ1,1)
    REAL ::           bi   (LX1,LY1,LZ1,1)
    COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT) &
    ,             W2   (LX1,LY1,LZ1,LELT)

    etime1=dnekclock()

    IF (IMESH == 1) NTOT = NX1*NY1*NZ1*NELV
    IF (IMESH == 2) NTOT = NX1*NY1*NZ1*NELT

    tol = tli
    if (param(22) /= 0) tol = abs(param(22))
    CALL CHKTCG1 (TOL,R,H1,H2,MASK,MULT,IMESH,ISD)


!     Set flags for overlapping Schwarz preconditioner (pff 11/12/98)

    kfldfdm = -1
!     if (name.eq.'TEMP') kfldfdm =  0
!     if (name.eq.'VELX') kfldfdm =  1
!     if (name.eq.'VELY') kfldfdm =  2
!     if (name.eq.'VELZ') kfldfdm =  3
    if (name == 'PRES') kfldfdm =  ndim+1

    call cggo &
    (u,r,h1,h2,mask,mult,imesh,tol,maxit,isd,bi,name)
    thmhz=thmhz+(dnekclock()-etime1)


    return
    end subroutine hmhzpf
!-----------------------------------------------------------------------
    subroutine hsolve(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd &
    ,approx,napprox,bi)

!     Either std. Helmholtz solve, or a projection + Helmholtz solve

    use ctimer
    use size_m
    use fdmh1
    use input
    use mass
    use tstep

    CHARACTER(4) ::    NAME
    REAL ::           U    (LX1,LY1,LZ1,1)
    REAL ::           R    (LX1,LY1,LZ1,1)
    REAL ::           H1   (LX1,LY1,LZ1,1)
    REAL ::           H2   (LX1,LY1,LZ1,1)
    REAL ::           vmk  (LX1,LY1,LZ1,1)
    REAL ::           vml  (LX1,LY1,LZ1,1)
    REAL ::           bi   (LX1,LY1,LZ1,1)
    REAL ::           approx (1)
    integer ::        napprox(1)
    common /ctmp2/ w1   (lx1,ly1,lz1,lelt)
    common /ctmp3/ w2   (2+2*mxprev)

    logical :: ifstdh
    character(4) ::  cname

    call chcopy(cname,name,4)
    call capit (cname,4)

    ifstdh = .TRUE. 

    if ( .NOT. ifflow) ifstdh = .FALSE. 

    if (param(95) /= 0 .AND. istep > param(95)) then
        if (cname == 'PRES') ifstdh = .FALSE. 
    elseif (param(94) /= 0 .AND. istep > param(94)) then
        ifstdh = .FALSE. 
    endif

    if (param(93) == 0) ifstdh = .TRUE. 

    if (ifstdh) then
        call hmholtz(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)
    else

        n = nx1*ny1*nz1*nelfld(ifield)

        call dssum  (r,nx1,ny1,nz1)
        call col2   (r,vmk,n)
        call projh  (r,h1,h2,bi,vml,vmk,approx,napprox,w1,w2,name)
        call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)
        call gensh  (u,h1,h2,vml,vmk,approx,napprox,w1,w2,name)

    endif

    return
    end subroutine hsolve
!-----------------------------------------------------------------------
