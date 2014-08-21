#if 0
    SUBROUTINE PLAN3 (IGEOM)
!-----------------------------------------------------------------------

!     Compute pressure and velocity using consistent approximation spaces.
!     Operator splitting technique.

!-----------------------------------------------------------------------
    use size_m
    use input
    use eigen
    include 'SOLN'
    include 'TSTEP'

    COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV) &
    ,              RESV2 (LX1,LY1,LZ1,LELV) &
    ,              RESV3 (LX1,LY1,LZ1,LELV) &
    ,              DV1   (LX1,LY1,LZ1,LELV) &
    ,              DV2   (LX1,LY1,LZ1,LELV) &
    ,              DV3   (LX1,LY1,LZ1,LELV)
    COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV) &
    ,              H2    (LX1,LY1,LZ1,LELV)

    IF (IGEOM == 1) THEN
    
    !        Old geometry
    
        CALL MAKEF
    
    ELSE
    
    !        New geometry, new b.c.
    
        INTYPE = -1
        CALL SETHLM  (H1,H2,INTYPE)
        CALL CRESVIF (RESV1,RESV2,RESV3,H1,H2)

        mstep = abs(param(94))
        if (param(94) /= 0. .AND. istep >= mstep) then
            CALL OPHINVpr(DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
        !          CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
        else
            CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
        endif
        CALL OPADD2  (VX,VY,VZ,DV1,DV2,DV3)
    
    !        Default Filtering
    
    !        alpha_filt = 0.05
    !        if (param(103).ne.0.) alpha_filt=param(103)
    !        call q_filter(alpha_filt)
    
    !        CALL SSNORMD (DV1,DV2,DV3)
    
        call incomprn(vx,vy,vz,pr)
    
    ENDIF

    RETURN
    END SUBROUTINE PLAN3

    SUBROUTINE LAGPRES
!--------------------------------------------------------------------

!     Keep old pressure values

!--------------------------------------------------------------------
    use size_m
    include 'SOLN'
    include 'TSTEP'

    common /cgeom/ igeom

    IF (NBDINP == 3 .AND. igeom <= 2) THEN
        NTOT2 = NX2*NY2*NZ2*NELV
        CALL COPY (PRLAG,PR,NTOT2)
    ENDIF
    RETURN
    END SUBROUTINE LAGPRES

    subroutine cresvif (resv1,resv2,resv3,h1,h2)
!---------------------------------------------------------------------

!     Compute startresidual/right-hand-side in the velocity solver

!---------------------------------------------------------------------
    use size_m
    include 'TOTAL'
    REAL ::           RESV1 (LX1,LY1,LZ1,1)
    REAL ::           RESV2 (LX1,LY1,LZ1,1)
    REAL ::           RESV3 (LX1,LY1,LZ1,1)
    REAL ::           H1    (LX1,LY1,LZ1,1)
    REAL ::           H2    (LX1,LY1,LZ1,1)
    COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV) &
    ,             W2    (LX1,LY1,LZ1,LELV) &
    ,             W3    (LX1,LY1,LZ1,LELV)

    common /cgeom/ igeom

    NTOT1 = NX1*NY1*NZ1*NELV
    NTOT2 = NX2*NY2*NZ2*NELV
    if (igeom == 2) CALL LAGVEL
    CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
    IF (IFSTRS)  CALL BCNEUTR

    call extrapp (pr,prlag)
    call opgradt (resv1,resv2,resv3,pr)
    CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)
    CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
    CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)

    RETURN
    end subroutine cresvif

    SUBROUTINE EXTRAPP_old (PREXTR)
!--------------------------------------------------------------------

!     Pressure extrapolation

!--------------------------------------------------------------------
    use size_m
    include 'SOLN'
    include 'TSTEP'
    COMMON /CTMP0/ DPR (LX2,LY2,LZ2,LELV)
    REAL ::        PREXTR (LX2,LY2,LZ2,LELV)

    common /cgeom/ igeom

    NTOT2 = NX2*NY2*NZ2*NELV

    IF (NBD == 1 .OR. NBD == 2 .OR. igeom > 2) THEN
        CALL COPY (PREXTR,PR,NTOT2)
    ELSEIF (NBD == 3) THEN
        CONST = DTLAG(1)/DTLAG(2)
        CALL SUB3 (DPR,PR,PRLAG,NTOT2)
        CALL CMULT(DPR,CONST,NTOT2)
        CALL ADD3 (PREXTR,PR,DPR,NTOT2)
    ELSEIF (NBD > 3) THEN
        WRITE (6,*) 'Pressure extrapolation cannot be completed'
        WRITE (6,*) 'Try a lower-order temporal scheme'
        call exitt
    ENDIF
    RETURN
    end SUBROUTINE EXTRAPP_old
!-----------------------------------------------------------------------
    subroutine ophinvpr(ot1,ot2,ot3,in1,in2,in3,h1,h2,tolh,nmxi)

!     OT = (H1*A+H2*B)-1 * IN  (implicit)

    use size_m
    use input
    include 'ORTHOV'
    include 'TSTEP'
    include 'SOLN'

    REAL :: OT1 (LX1,LY1,LZ1,1)
    REAL :: OT2 (LX1,LY1,LZ1,1)
    REAL :: OT3 (LX1,LY1,LZ1,1)
    REAL :: IN1 (LX1,LY1,LZ1,1)
    REAL :: IN2 (LX1,LY1,LZ1,1)
    REAL :: IN3 (LX1,LY1,LZ1,1)
    REAL :: H1  (LX1,LY1,LZ1,1)
    REAL :: H2  (LX1,LY1,LZ1,1)


    IMESH = 1

    IF (IFSTRS) THEN
#if 0
        MATMOD = 0
        CALL HMHZSF  ('NOMG',OT1,OT2,OT3,IN1,IN2,IN3,H1,H2, &
        V1MASK,V2MASK,V3MASK,VMULT, &
        TOLH,NMXi,MATMOD)
#endif
    ELSE
        CALL hmzpf2 ('VELX',OT1,IN1,H1,H2,V1MASK,VMULT, &
        IMESH,TOLH,NMXi,1)
        CALL hmzpf2 ('VELY',OT2,IN2,H1,H2,V2MASK,VMULT, &
        IMESH,TOLH,NMXi,2)
        IF (NDIM == 3) &
        CALL hmzpf2 ('VELZ',OT3,IN3,H1,H2,V3MASK,VMULT, &
        IMESH,TOLH,NMXi,3)
    ENDIF

    return
    end subroutine ophinvpr
!-----------------------------------------------------------------------
    subroutine hmzpf2(nm,u,rhs,h1,h2,mask,mult,imsh,tol,mxit,isd)
    use size_m
    use input
    include 'MASS'
    character(4) :: nm

    REAL ::           U    (LX1,LY1,LZ1,1)
    REAL ::           RHS  (LX1,LY1,LZ1,1)
    REAL ::           H1   (LX1,LY1,LZ1,1)
    REAL ::           H2   (LX1,LY1,LZ1,1)
    REAL ::           MASK (LX1,LY1,LZ1,1)
    REAL ::           MULT (LX1,LY1,LZ1,1)

    ntot1 = nx1*ny1*nz1*nelv
    if (imsh == 2) ntot1 = nx1*ny1*nz1*nelt

    call col2   (rhs,mask,ntot1)
    call dssum  (rhs,nx1,ny1,nz1)
    call projh2 (rhs,h1,h2,mult,mask,isd,imsh)
    if (imsh == 1) then
        call hmhzpf (nm,u,rhs,h1,h2,mask,mult,imsh,tol,mxit,isd,binvm1)
    else
        call hmhzpf (nm,u,rhs,h1,h2,mask,mult,imsh,tol,mxit,isd,bintm1)
    endif
    call gensh2 (u,h1,h2,mult,mask,isd,imsh)

    return
    end subroutine hmzpf2
!-----------------------------------------------------------------------
    subroutine projh2(v1,h1,h2,vml,vmask,isd,imsh)

!     Orthogonalize the rhs wrt previous rhs's for which we already
!     know the soln.

    use size_m
    use input
    include 'MASS'
    include 'SOLN'
    include 'TSTEP'
    include 'ORTHOV'

    real :: v1(1),h1(1),h2(1),vml(1),vmask(1)
    real :: work(mxprev)

    integer :: icalld
    save    icalld
    data    icalld/0/

    ntot1=nx1*ny1*nz1*nelv
    if (imsh == 2) ntot1=nx1*ny1*nz1*nelt

    if (icalld == 0) then ! First call, no vectors to orthogonalize against.
        call izero(nprev,ndim)
        mprev=param(93)
        mprev=min(mprev,mxprev)
        if (mprev == 0) mprev = mxprev
        if (nid == 0) write(6,*) 'this is mprev:',mprev,mxprev
    endif

!     Diag to see how much reduction in the residual is attained.
    if (imsh == 1) then
        alpha1 = glsc23(v1,binvm1,vml,ntot1)
        alpha1 = sqrt(alpha1/volvm1)
    else
        alpha1 = glsc23(v1,bintm1,vml,ntot1)
        alpha1 = sqrt(alpha1/voltm1)
    endif

!     if (icalld.eq.0.and.nid.eq.0)
!    $     write(6,*) 'alpha1:',alpha1,volvm1,ntot1
!     if (icalld.eq.0.and.nid.eq.0)
!    $     write(6,*) 'binvm1:',binvm1(1,1,1,1),vml(1),v1(1)


    call updrhsh2(h1,h2,vml,vmask,isd,imsh) ! Update rhs's if matrix has changed


    call rzero(alpha,mxprev) !  Gram-Schmidt for previous soln's
    ioff = 1
    do i=1,nprev(isd)
        alpha(i) = vlsc3(v1,sln(ioff,isd),vml,ntot1)
        ioff = ioff + ntot1
    enddo

    if (nprev(isd) > 0) then
        call gop(alpha,work,'+  ',nprev(isd))
        call cmult2(vbar(1,isd),sln(1,isd),alpha(1),ntot1)

        do i=2,nprev(isd)
            ioff = ntot1*(i-1)+1
            call add2s2 (vbar(1,isd),sln(ioff,isd),alpha(i),ntot1)
        enddo
    !        alphmn = glmin (vbar(1,isd),ntot1)
    !        alphmx = glmax (vbar(1,isd),ntot1)
        call axhelm    (vnew(1,isd),vbar(1,isd),H1,H2,1,1)
    !        alp1mn = glmin (vnew(1,isd),ntot1)
    !        alp1mx = glmax (vnew(1,isd),ntot1)
        call col2      (vnew(1,isd),vmask,ntot1)
    !        alp2mn = glmin (vnew(1,isd),ntot1)
    !        alp2mx = glmax (vnew(1,isd),ntot1)
        call dssum     (vnew(1,isd),nx1,ny1,nz1)
    !        alp3mn = glmin (vnew(1,isd),ntot1)
    !        alp3mx = glmax (vnew(1,isd),ntot1)
        call sub2      (v1,vnew(1,isd),ntot1)
    else
        call rzero     (vnew(1,isd),ntot1)
        call rzero     (vbar(1,isd),ntot1)
    endif

!     if (nid.eq.0) write(6,90) istep,alphmn,alphmx
!    $              ,alp1mn,alp1mx,alp2mn,alp2mx,alp3mn,alp3mx
!  90 format(i4,1p8e11.3,' xx')

!     Diag. ............................................................
    if (imsh == 1) then
        alpha2 = glsc23(v1,binvm1,vml,ntot1)
        alpha2 = sqrt(alpha2/volvm1)
    else
        alpha2 = glsc23(v1,bintm1,vml,ntot1)
        alpha2 = sqrt(alpha2/voltm1)
    endif
    ratio  = alpha1/alpha2
    n10=min(10,nprev(isd))
    if (nid == 0) write(6,10) istep,isd,alpha1,alpha2,ratio,nprev(isd)
    10 format(i8,i3,1p3e12.4,i4,' alph1x')
    if (nid == 0) write(6,11) istep,nprev(isd),(alpha(I),I=1,n10)
    11 format(i6,' halpha',i4,10(1p10e12.4,/,17x))

!     alphmn = glmax(v1,ntot1)
!     alphmx = glmin(v1,ntot1)
!     if (nid.eq.0) write(6,10) istep,alphmn,alphmx,ratio,nprev(isd)

!     Diag. .............................................................

    icalld=icalld+1

    return
    end subroutine projh2
!-----------------------------------------------------------------------
    subroutine gensh2(v1,h1,h2,vml,vmask,isd,imsh)

!     Reconstruct the solution to the original problem by adding back
!     the previous solutions

    use size_m
    use input
    include 'MASS'
    include 'SOLN'
    include 'TSTEP'
    include 'ORTHOV'
    real :: v1(1),h1(1),h2(1),vml(1),vmask(1)

    ntot1=nx1*ny1*nz1*nelv
    if (imsh == 2) ntot1=nx1*ny1*nz1*nelt

    call copy (vnew(1,isd),v1,ntot1)            !  Save current solution
    call add2(v1,vbar(1,isd),ntot1)             !  Reconstruct solution
    call updtseth2(v1,h1,h2,vml,vmask,isd,imsh) !  Update {SLN}

    return
    end subroutine gensh2
!-----------------------------------------------------------------------
    subroutine updtseth2(v1,h1,h2,vml,vmask,isd,imsh)

!     Update the set of rhs's and the corresponding p-set:

!        . Standard case is to add P_new, and RHS_new = E*P_new

!        . However, when nprev=mprev (max. allowed), we throw out
!          the old set, and set P_1 = P, RHS_1=E*P_1

!        . Other schemes are possible, e.g., let's save a bunch of
!          old vectors, perhaps chosen wisely via P.O.D.

    use size_m
    use input
    include 'MASS'
    include 'SOLN'
    include 'TSTEP'
    include 'ORTHOV'
    real :: v1(1),h1(1),h2(1),vml(1),vmask(1)

    ntot1=nx1*ny1*nz1*nelv
    if (imsh == 2) ntot1=nx1*ny1*nz1*nelt
         
    if (nprev(isd) == mprev) then
        call copy(vnew(1,isd),v1,ntot1)
        nprev(isd)=0
    endif

!     Increment solution set
    nprev(isd) = nprev(isd)+1
    ioff = ntot1*(nprev(isd)-1)+1
    call copy(sln(ioff,isd),vnew(1,isd),ntot1)

!     Orthogonalize rhs against previous rhs and normalize
    call hconj2(nprev(isd),h1,h2,vml,vmask,isd,imsh)

!     Save last sol'n
    call copy(vnew(1,isd),v1,ntot1)

    return
    end subroutine updtseth2
!-----------------------------------------------------------------------
    subroutine hconj2(kprev,h1,h2,vml,vmask,isd,imsh)

!     Orthonormalize the last saved vector against vector set

    use size_m
    use input
    include 'MASS'
    include 'SOLN'
    include 'TSTEP'
    include 'ORTHOV'
    real :: h1(1),h2(1),vml(1),vmask(1)
    real :: work(mxprev)

    ntot1 = nx1*ny1*nz1*nelv
    if (imsh == 2) ntot1 = nx1*ny1*nz1*nelt

    kprev1 = kprev-1
    i1     = kprev1*ntot1 + 1

    call axhelm  (vbar(1,isd),sln(i1,isd),h1,h2,1,1)
    call col2    (vbar(1,isd),vmask,ntot1)
    call dssum   (vbar(1,isd),nx1,ny1,nz1)
    call col2    (vbar(1,isd),vml   ,ntot1) ! Compute part of the norm
    alphad=glsc2 (vbar(1,isd),sln(i1,isd),ntot1)

    do i=1,kprev1                    ! Gram-Schmidt
        ioff = (i-1)*ntot1 + 1
        alpha(i) = vlsc2(vbar(1,isd),sln(ioff,isd),ntot1)
    enddo
    if (kprev1 > 0) call gop(alpha,work,'+  ',kprev1)

    do i=1,kprev1
        alpham = -alpha(i)
        ioff = (i-1)*ntot1 + 1
        call add2s2(sln(i1,isd),sln(ioff,isd),alpham,ntot1)
        alphad = alphad - alpha(i)**2
    enddo

!    .Normalize new element in P~
    alphad = 1.0/sqrt(alphad)
    call cmult(sln(i1,isd),alphad,ntot1)

    return
    end subroutine hconj2
!-----------------------------------------------------------------------
    subroutine updrhsh2(h1,h2,vml,vmask,isd,imsh)

!     Update rhs's if A-matrix has changed

    use size_m
    use input
    include 'MASS'
    include 'ORTHOV'
    include 'TSTEP'


    real :: vml(1),h1(1),h2(1),vmask(1)

    real ::    dtold
    save    dtold
    data    dtold/0.0/

!     First, we have to decide if the E matrix has changed.

    if (dt == dtold) return
    dtold = dt
    call izero(nprev,ldim)
    return

    do iprev=1,nprev(isd) ! Orthogonalize this rhs w.r.t. previous rhs's
        call hconj2(iprev,h1,h2,vml,vmask,isd,imsh)
    enddo

    return
    end subroutine updrhsh2
#endif
!-----------------------------------------------------------------------
