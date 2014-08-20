!-----------------------------------------------------------------------

!     To do:

!        Differing BC's imposed for ophinv, incomprn, etc.

!        1-shot Fast solver for Helmholtz and pressure


!-----------------------------------------------------------------------
    subroutine induct (igeom)

!     Solve the convection-diffusion equation for the B-field, with
!     projection onto a div-free space.


    use size_m
    include 'INPUT'
    include 'EIGEN'
    include 'SOLN'
    include 'TSTEP'
    include 'MASS'

    common /scrns/  resv1 (lx1,ly1,lz1,lelv) &
    ,              resv2 (lx1,ly1,lz1,lelv) &
    ,              resv3 (lx1,ly1,lz1,lelv) &
    ,              dv1   (lx1,ly1,lz1,lelv) &
    ,              dv2   (lx1,ly1,lz1,lelv) &
    ,              dv3   (lx1,ly1,lz1,lelv)
    common /scrvh/  h1    (lx1,ly1,lz1,lelv) &
    ,              h2    (lx1,ly1,lz1,lelv)

    ifield = ifldmhd

    if (igeom == 1) then  ! old geometry, old velocity

        call makebsource_mhd

    else

        call lagbfield
        call lagvel

        call elsasserh(igeom)

        call vol_flow        ! check for fixed flow rate


    endif

    return
    end subroutine induct
!--------------------------------------------------------------------
    subroutine lagbfield

!     Keep old B-field(s)

    use size_m
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'

    do ilag=nbdinp-1,2,-1
        call opcopy &
        ( bxlag(1,ilag  ),bylag(1,ilag  ),bzlag(1,ilag  ) &
        , bxlag(1,ilag-1),bylag(1,ilag-1),bzlag(1,ilag-1) )
    enddo
    call opcopy (bxlag,bylag,bzlag,bx,by,bz)

    return
    end subroutine lagbfield
!--------------------------------------------------------------------
    subroutine makebsource_mhd

!     Make rhs for induction equation

    use ctimer
    use size_m
    include 'SOLN'
    include 'MASS'
    include 'INPUT'
    include 'TSTEP'

    if (icalld == 0) tbmhd=0.0
    icalld = icalld+1
    nbmhd  = icalld
    etime1 = dnekclock()

    ifield = 1
    call makeuf

    ifield = ifldmhd
    call makeufb
    if (ifaxis) then
    !        do ifield = 2,3
        do ifield = 2,npscal+1
            call makeuq  ! nonlinear terms
        enddo
        ifield = ifldmhd
    endif

    if (ifnav .AND. ( .NOT. ifchar)) then
        call advab_elsasser_fast
    endif
    if (ifchar) then
        write(6,*) 'No IFCHAR for MHD, yet.'
        call exitt
    endif

    ifield = 1
    if (iftran)                     call makeabf
    call makebdf

    ifield = ifldmhd
    if (iftran)                     call makextb
    call makebdfb

    tbmhd=tbmhd+(dnekclock()-etime1)
    return
    end subroutine makebsource_mhd
!--------------------------------------------------------------------
    subroutine makeufb

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

    use size_m
    include 'SOLN'
    include 'MASS'
    include 'TSTEP'

    time = time-dt
    call nekuf   (bmx,bmy,bmz)
    call opcolv2 (bmx,bmy,bmz,vtrans(1,1,1,1,ifield),bm1)
    time = time+dt

    return
    end subroutine makeufb
!--------------------------------------------------------------------
    subroutine makextb

!     Add extrapolation terms to magnetic source terms

!     (nek5 equivalent for velocity is "makeabf")

    use size_m
    include 'INPUT'
    include 'SOLN'
    include 'MASS'
    include 'TSTEP'

    common /scrns/ ta1 (lx1,ly1,lz1,lelv) &
    ,             ta2 (lx1,ly1,lz1,lelv) &
    ,             ta3 (lx1,ly1,lz1,lelv)

    ntot1 = nx1*ny1*nz1*nelv

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)
    call add3s2 (ta1,bbx1,bbx2,ab1,ab2,ntot1)
    call add3s2 (ta2,bby1,bby2,ab1,ab2,ntot1)
    call copy   (bbx2,bbx1,ntot1)
    call copy   (bby2,bby1,ntot1)
    call copy   (bbx1,bmx,ntot1)
    call copy   (bby1,bmy,ntot1)
    call add2s1 (bmx,ta1,ab0,ntot1)
    call add2s1 (bmy,ta2,ab0,ntot1)
    if (ndim == 3) then
        call add3s2 (ta3,bbz1,bbz2,ab1,ab2,ntot1)
        call copy   (bbz2,bbz1,ntot1)
        call copy   (bbz1,bmz,ntot1)
        call add2s1 (bmz,ta3,ab0,ntot1)
    endif

    return
    end subroutine makextb
!--------------------------------------------------------------------
    subroutine makebdfb

!     Add contributions to magnetic source from lagged BD terms.

    use size_m
    include 'SOLN'
    include 'MASS'
    include 'GEOM'
    include 'INPUT'
    include 'TSTEP'

    COMMON /SCRNS/ TA1(LX1,LY1,LZ1,LELV) &
    ,             TA2(LX1,LY1,LZ1,LELV) &
    ,             TA3(LX1,LY1,LZ1,LELV) &
    ,             TB1(LX1,LY1,LZ1,LELV) &
    ,             TB2(LX1,LY1,LZ1,LELV) &
    ,             TB3(LX1,LY1,LZ1,LELV) &
    ,             H2 (LX1,LY1,LZ1,LELV)

    NTOT1 = NX1*NY1*NZ1*NELV
    CONST = 1./DT
    CALL CMULT2(H2,vtrans(1,1,1,1,ifield),CONST,NTOT1)
    CALL OPCOLV3c (TB1,TB2,TB3,BX,BY,BZ,BM1,bd(2))

    DO 100 ILAG=2,NBD
        IF (IFGEOM) THEN
            CALL OPCOLV3c(TA1,TA2,TA3,BXLAG (1,ILAG-1), &
            BYLAG (1,ILAG-1), &
            BZLAG (1,ILAG-1), &
            BM1LAG(1,1,1,1,ILAG-1),bd(ilag+1))
        ELSE
            CALL OPCOLV3c(TA1,TA2,TA3,BXLAG (1,ILAG-1), &
            BYLAG (1,ILAG-1), &
            BZLAG (1,ILAG-1), &
            BM1                   ,bd(ilag+1))
        ENDIF
        CALL OPADD2  (TB1,TB2,TB3,TA1,TA2,TA3)
    100 END DO
    CALL OPADD2col (BMX,BMY,BMZ,TB1,TB2,TB3,h2)

    return
    end subroutine makebdfb
!--------------------------------------------------------------------
    subroutine cresvib(resv1,resv2,resv3,h1,h2)

!     Account for inhomogeneous Dirichlet boundary contributions
!     in rhs of induction eqn.
!                                               n
!     Also, subtract off best estimate of grad p

    use size_m
    include 'TOTAL'
    real ::           resv1 (lx1,ly1,lz1,1)
    real ::           resv2 (lx1,ly1,lz1,1)
    real ::           resv3 (lx1,ly1,lz1,1)
    real ::           h1    (lx1,ly1,lz1,1)
    real ::           h2    (lx1,ly1,lz1,1)
    common /scruz/ w1    (lx1,ly1,lz1,lelv) &
    ,             w2    (lx1,ly1,lz1,lelv) &
    ,             w3    (lx1,ly1,lz1,lelv)


    call bcdirvc (bx,by,bz,b1mask,b2mask,b3mask)
    call extrapp (pm,pmlag)
    call opgradt (resv1,resv2,resv3,pm)
    call opadd2  (resv1,resv2,resv3,bmx,bmy,bmz)
!     call opcopy  (resv1,resv2,resv3,bmx,bmy,bmz)
    call ophx    (w1,w2,w3,bx,by,bz,h1,h2)
    call opsub2  (resv1,resv2,resv3,w1,w2,w3)

    return
    end subroutine cresvib
!--------------------------------------------------------------------
    subroutine cresvibp(resv1,resv2,resv3,h1,h2)

!     Account for inhomogeneous Dirichlet boundary contributions
!     in rhs of momentum eqn.
!                                               n
!     Also, subtract off best estimate of grad p

    use size_m
    include 'TOTAL'
    real ::           resv1 (lx1,ly1,lz1,1)
    real ::           resv2 (lx1,ly1,lz1,1)
    real ::           resv3 (lx1,ly1,lz1,1)
    real ::           h1    (lx1,ly1,lz1,1)
    real ::           h2    (lx1,ly1,lz1,1)
    common /scruz/ w1    (lx1,ly1,lz1,lelv) &
    ,             w2    (lx1,ly1,lz1,lelv) &
    ,             w3    (lx1,ly1,lz1,lelv)

    call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)
    if (ifstrs)  call bcneutr
    call extrapp (pr,prlag)
    call opgradt (resv1,resv2,resv3,pr)
    call opadd2  (resv1,resv2,resv3,bfx,bfy,bfz)
!     call opcopy  (resv1,resv2,resv3,bfx,bfy,bfz)
    call ophx    (w1,w2,w3,vx,vy,vz,h1,h2)
    call opsub2  (resv1,resv2,resv3,w1,w2,w3)

    return
    end subroutine cresvibp
!--------------------------------------------------------------------
    subroutine incomprn (ux,uy,uz,up)

!     Project U onto the closest incompressible field

!     Input:  U     := (ux,uy,uz)

!     Output: updated values of U, iproj, proj; and
!             up    := pressure currection req'd to impose div U = 0


!     Dependencies: ifield ==> which "density" (vtrans) is used.

!     Notes  1.  up is _not_ scaled by bd(1)/dt.  This should be done
!                external to incompr().

!            2.  up accounts _only_ for the perturbation pressure,
!                not the current pressure derived from extrapolation.


    use ctimer
    use size_m
    include 'TOTAL'

    common /scrns/ w1    (lx1,ly1,lz1,lelv) &
    ,             w2    (lx1,ly1,lz1,lelv) &
    ,             w3    (lx1,ly1,lz1,lelv) &
    ,             dv1   (lx1,ly1,lz1,lelv) &
    ,             dv2   (lx1,ly1,lz1,lelv) &
    ,             dv3   (lx1,ly1,lz1,lelv) &
    ,             dpp    (lx2,ly2,lz2,lelv)
    common /scrvh/ h1    (lx1,ly1,lz1,lelv) &
    ,             h2    (lx1,ly1,lz1,lelv)
    common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

    parameter(nset = 1 + lbelv/lelv)
    common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
    common /orthbi/ nprv(2)
    logical :: ifprjp

    ifprjp= .FALSE.    ! Project out previous pressure solutions?
    istart=param(95)
    if (istep >= istart .AND. istart /= 0) ifprjp= .TRUE. 

    if (icalld == 0) tpres=0.0
    icalld = icalld+1
    npres  = icalld
    etime1 = dnekclock()

    ntot1  = nx1*ny1*nz1*nelv
    ntot2  = nx2*ny2*nz2*nelv
    intype = 1

    call rzero   (h1,ntot1)
    call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
    call invers2 (h2inv,h2,ntot1)

    call opdiv   (dpp,ux,uy,uz)

    bdti = -bd(1)/dt
    call cmult   (dpp,bdti,ntot2)

    call add2col2(dpp,bm2,usrdiv,ntot2) ! User-defined divergence.

    call ortho   (dpp)

    i = 1 + ifield/ifldmhd
    if (ifprjp)   call setrhsp  (dpp,h1,h2,h2inv,pset(1,i),nprv(i))
    scaledt = dt/bd(1)
    scaledi = 1./scaledt
    call cmult(dpp,scaledt,ntot2)        ! scale for tol
    call esolver  (dpp,h1,h2,h2inv,intype)
    call cmult(dpp,scaledi,ntot2)
    if (ifprjp)   call gensolnp (dpp,h1,h2,h2inv,pset(1,i),nprv(i))

    call add2(up,dpp,ntot2)

    call opgradt  (w1 ,w2 ,w3 ,dpp)
    call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
    dtb  = dt/bd(1)
    call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )

    if (ifmhd)  call chkptol	! to avoid repetition

    tpres=tpres+(dnekclock()-etime1)

    return
    end subroutine incomprn
!-----------------------------------------------------------------------
    subroutine extrapp(p,plag)

!     Pressure extrapolation

    use size_m
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'

    real ::  p    (lx2,ly2,lz2,1) &
    ,plag (lx2,ly2,lz2,1)

    common /cgeom/ igeom

    ntot2 = nx2*ny2*nz2*nelv

    if (nbd == 3 .AND. igeom <= 2) then

        const = dtlag(1)/dtlag(2)

        do i=1,ntot2
            pnm1          = p   (i,1,1,1)
            pnm2          = plag(i,1,1,1)
            p   (i,1,1,1) = pnm1 + const*(pnm1-pnm2)
            plag(i,1,1,1) = pnm1
        enddo

    elseif (nbd > 3) then
        WRITE (6,*) 'Pressure extrapolation cannot be completed'
        WRITE (6,*) 'Try a lower-order temporal scheme'
        call exitt
    endif
    return
    end subroutine extrapp
!-----------------------------------------------------------------------
    subroutine opzero(ux,uy,uz)
    use size_m
    include 'TOTAL'
    real :: ux(1),uy(1),uz(1)

    n = nx1*ny1*nz1*nelfld(ifield)
    call rzero(ux,n)
    call rzero(uy,n)
    if (if3d) call rzero(uz,n)

    return
    end subroutine opzero
!-----------------------------------------------------------------------
    subroutine opnorm(unorm,ux,uy,uz,type3)
    use size_m
    include 'TOTAL'
    character(3) :: type3
    real :: ux(1),uy(1),uz(1)
    real :: un(3),wn(3)

    n = nx1*ny1*nz1*nelfld(ifield)
    if (type3 == 'L2 ') then
        if (if3d) then
            un(1) = vlsc3(ux,ux,bm1,n)
            un(2) = vlsc3(uy,uy,bm1,n)
            un(3) = vlsc3(uz,uz,bm1,n)
            un(1) = un(1) + un(2) + un(3)
            unorm = glsum(un(1),1)
            if (unorm > 0) unorm = sqrt(unorm/volvm1)
        else
            un(1) = vlsc3(ux,ux,bm1,n)
            un(2) = vlsc3(uy,uy,bm1,n)
            un(1) = un(1) + un(2)
            unorm = glsum(un(1),1)
            if (unorm > 0) unorm = sqrt(unorm/volvm1)
        endif
    endif

    return
    end subroutine opnorm
!-----------------------------------------------------------------------
    subroutine lorentz_force (lf,b1,b2,b3,w1,w2)

!     Compute Lorentz force

!     Input:  B     := (b1,b2,b3)

!     Output: lf(1,ldim)

!     Work arrays: w1(ltot) and w2(ltot)

!     The output will not be continuous.  However, it will be in
!     the form appropriate for incorporation as a body force term
!     in the variational formulation of the Navier-Stokes equations.

!     (i.e.,   rhs(NS) = rhs(NS) + B*lf,  where B is the mass matrix)

    use size_m

    real :: lf(lx1*ly1*lz1*lelv,ldim)
    real :: b1(lx1*ly1*lz1*lelv)
    real :: b2(lx1*ly1*lz1*lelv)
    real :: b3(lx1*ly1*lz1*lelv)

    call curl(lf,b1,b2,b3, .FALSE. ,w1,w2)

    ntot = nx1*ny1*nz1*nelv

    do i=1,ntot
        c1 = lf(i,2)*b3(i) - lf(i,3)*b2(i)
        c2 = lf(i,3)*b1(i) - lf(i,1)*b3(i)
        c3 = lf(i,1)*b2(i) - lf(i,2)*b1(i)
        lf(i,1) = c1
        lf(i,2) = c2
        lf(i,3) = c3
    enddo

    return
    end subroutine lorentz_force
!-----------------------------------------------------------------------
    subroutine curl(vort,u,v,w,ifavg,work1,work2)

    use size_m
    include 'TOTAL'

    logical :: ifavg

    parameter(lt=lx1*ly1*lz1*lelv)
    real :: vort(lt,3),work1(1),work2(1),u(1),v(1),w(1)

    ntot  = nx1*ny1*nz1*nelv
    if (if3d) then
    !        work1=dw/dy ; work2=dv/dz
        call dudxyz(work1,w,rym1,sym1,tym1,jacm1,1,2)
        call dudxyz(work2,v,rzm1,szm1,tzm1,jacm1,1,3)
        call sub3(vort(1,1),work1,work2,ntot)
    !        work1=du/dz ; work2=dw/dx
        call dudxyz(work1,u,rzm1,szm1,tzm1,jacm1,1,3)
        call dudxyz(work2,w,rxm1,sxm1,txm1,jacm1,1,1)
        call sub3(vort(1,2),work1,work2,ntot)
    !        work1=dv/dx ; work2=du/dy
        call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
        call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
        call sub3(vort(1,3),work1,work2,ntot)
    else
    !        work1=dv/dx ; work2=du/dy
        call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
        call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
        call sub3(vort(1,3),work1,work2,ntot)
    endif

!    Avg at bndry

    if (ifavg) then
        ifielt = ifield
        ifield = 1
        if (if3d) then
            do idim=1,ndim
                call col2  (vort(1,idim),bm1,ntot)
                call dssum (vort(1,idim),nx1,ny1,nz1)
                call col2  (vort(1,idim),binvm1,ntot)
            enddo
        else
            call col2  (vort(1,3),bm1,ntot)    ! NOTE:  This differs from
            call dssum (vort(1,3),nx1,ny1,nz1) ! "comp_vort", which returns
            call col2  (vort(1,3),binvm1,ntot) ! vorticity as 1st entry in vort
        endif
        ifield = ifielt
    endif

    return
    end subroutine curl
!-----------------------------------------------------------------------
    subroutine lorentz_force2(lf,b1,b2,b3)

!     Compute Lorentz force

!     Input:  B     := (b1,b2,b3)

!     Output: lf(1,ldim)

!     Work arrays: w1(ltot) and w2(ltot)

!     The output will not be continuous.  However, it will be in
!     the form appropriate for incorporation as a body force term
!     in the variational formulation of the Navier-Stokes equations.

!     (i.e.,   rhs(NS) = rhs(NS) + B*lf,  where B is the mass matrix)

    use size_m

    real :: lf(lx1*ly1*lz1*ldim,lelt)
    real :: b1(lx1*ly1*lz1,lelt)
    real :: b2(lx1*ly1*lz1,lelt)
    real :: b3(lx1*ly1*lz1,lelt)

    integer :: e

    do e = 1,nelt   ! NOTE:  the order is different from v. 1
        call lorentz_force_e(lf(1,e),b1(1,e),b2(1,e),b3(1,e),e)
    enddo

    return
    end subroutine lorentz_force2
!-----------------------------------------------------------------------
    subroutine lorentz_force_e(lf,b1,b2,b3,e)

!     Compute Lorentz force field for a single element

!     Input:  B     := (b1,b2,b3)

!     Output: lf(1,ldim)

!     Work arrays: cb(lxyzd) and w2(lxyzd)

!     The output will not be continuous.  However, it will be in
!     the form appropriate for incorporation as a body force term
!     in the variational formulation of the Navier-Stokes equations.
!     (i.e.,   rhs(NS) = rhs(NS) + B*lf,  where B is the mass matrix)

!     Dealiasing strategy:

!       e      e  -1   ~ T  ~e   ~ ~
!     lf  = ( B  )     I    B  ( I j x I b )
!       i              ~                    i

!            ~
!     Here,  I is the interpolant from N to M, where M = 3/2 N.

!            ~                                  -1
!            j is a special curl  (sans Jacobian   )

!            b is the B-field on the N-points

!            ~e
!            B  is the local mass matrix on the M-points (sans Jacabian)

!             e
!            B  is the local mass matrix on the N-points (with Jacabian)

!     The last B is req'd to compensate for the subsequent multiplication
!     by B that takes place once the rhs is formed.



    use size_m
    use dealias
    include 'GEOM'
    include 'WZ'

    real :: lf(lx1*ly1*lz1,3)
    real :: b1(lx1*ly1*lz1)
    real :: b2(lx1*ly1*lz1)
    real :: b3(lx1*ly1*lz1)
    integer :: d,e

    integer :: icalld
    save    icalld
    data    icalld /0/

    common /ctmp1x/ lfd(lxd*lyd*lzd,3) &
    ,  bd (lxd*lyd*lzd,3) &
    ,  cb (lx1*ly1*lz1,3) &
    ,  cbd(lxd*lyd*lzd,3)
    real :: lfd

    if (icalld == 0) then
        write(6,*) 'CALL SET PROJ',nx1,nxd
        call setmap (nx1,nxd)   ! Set up interpolation operators
        call setproj(nx1,nxd)   ! Set up interpolation operators
        icalld = icalld + 1
    endif

    call spec_curl_e (cb,b1,b2,b3                   & !local curl, w/o Jacobian
    , rxm1(1,1,1,e),rym1(1,1,1,e),rzm1(1,1,1,e) &
    , sxm1(1,1,1,e),sym1(1,1,1,e),szm1(1,1,1,e) &
    , txm1(1,1,1,e),tym1(1,1,1,e),tzm1(1,1,1,e) )

    do d=1,ndim                                   ! interpolate to M points
        call specmp(cbd(1,d),nxd,cb(1,d),nx1,im1d,im1dt,bd)
    enddo
    call specmp(bd(1,1),nxd,b1,nx1,im1d,im1dt,lfd)
    call specmp(bd(1,2),nxd,b2,nx1,im1d,im1dt,lfd)
    call specmp(bd(1,3),nxd,b3,nx1,im1d,im1dt,lfd)

    nxyzd = nxd*nyd*nzd
    do i=1,nxyzd
        lfd(i,1) = cbd(i,2)*bd(i,3) - cbd(i,3)*bd(i,2) ! Curl B x B
        lfd(i,2) = cbd(i,3)*bd(i,1) - cbd(i,1)*bd(i,3)
        lfd(i,3) = cbd(i,1)*bd(i,2) - cbd(i,2)*bd(i,1)
    enddo

!     Project back and simultaneous collocate with local quadrature weights

!                                 ~        ~        ~
!        P := Pz x Py x Px  =  Iz*B  x  Iy*B  x  Ix*B


!           ~            M
!     where B = diag (rho  )
!                        i

    call specmp(lf(1,1),nx1,lfd(1,1),nxd,pmd1,pmd1t,cbd)
    call specmp(lf(1,2),nx1,lfd(1,2),nxd,pmd1,pmd1t,cbd)
    call specmp(lf(1,3),nx1,lfd(1,3),nxd,pmd1,pmd1t,cbd)

!     Finally, divide by local mass matrix in anticipation of subsequent
!     multiply by BM1.

    nxyz = nx1*ny1*nz1
    do i=1,nxyz
    !        scale = 1./(w3m1(i,1,1)*jacm1(i,1,1,e))
        scale = 1./jacm1(i,1,1,e)
        lf(i,1) = scale*lf(i,1)
        lf(i,2) = scale*lf(i,2)
        lf(i,3) = scale*lf(i,3)
    enddo

    return
    end subroutine lorentz_force_e
!-----------------------------------------------------------------------
    subroutine spec_curl_e (cb,b1,b2,b3,rx,ry,rz,sx,sy,sz,tx,ty,tz)

!     local curl, multiplied by Jacobian

    use size_m
    include 'DXYZ'

    real :: cb(lx1*ly1*lz1,3)  ! Output J*curl B  (J:=Jacobian)
    real :: b1(1),b2(1),b3(1)  ! Input B-field

    real :: rx(1),ry(1),rz(1)  ! Metrics
    real :: sx(1),sy(1),sz(1)
    real :: tx(1),ty(1),tz(1)

    common /ctmp0x/ br(lx1*ly1*lz1),bs(lx1*ly1*lz1),bt(lx1*ly1*lz1)

!              / db3     db2 \
!     cb1 = J  | ---  -  --- |    ! Keep J ( J:=Jacobian)
!              \ dy      dz  /


!              / db1     db3 \
!     cb2 = J  | ---  -  --- |    ! Keep J ( J:=Jacobian)
!              \ dz      dx  /


!              / db2     db1 \
!     cb3 = J  | ---  -  --- |    ! Keep J ( J:=Jacobian)
!              \ dx      dy  /


!     Note:

!      db2      db2   dr     db2   ds     db2   dt
!    J --- =    --- J --  +  --- J --  +  --- J --
!      dz       dr    dz     ds    dz     dt    dz

!              etc.



    nxyz = nx1*ny1*nz1

    N=nx1-1
    call local_grad3(br,bs,bt,b1,N,1,dxm1,dxtm1)
    do i=1,nxyz
        cb(i,2) =  (br(i)*rz(i)+bs(i)*sz(i)+bt(i)*tz(i))
        cb(i,3) = -(br(i)*ry(i)+bs(i)*sy(i)+bt(i)*ty(i))
    enddo

    call local_grad3(br,bs,bt,b2,N,1,dxm1,dxtm1)
    do i=1,nxyz
        cb(i,1) = -(br(i)*rz(i)+bs(i)*sz(i)+bt(i)*tz(i))
        cb(i,3) =  (br(i)*rx(i)+bs(i)*sx(i)+bt(i)*tx(i)) &
        +  cb(i,3)
    enddo

    call local_grad3(br,bs,bt,b3,N,1,dxm1,dxtm1)
    do i=1,nxyz
        cb(i,1) =  (br(i)*ry(i)+bs(i)*sy(i)+bt(i)*ty(i)) &
        +  cb(i,1)
        cb(i,2) = -(br(i)*rx(i)+bs(i)*sx(i)+bt(i)*tx(i)) &
        +  cb(i,2)
    enddo

    return
    end subroutine spec_curl_e
!-----------------------------------------------------------------------
    subroutine specx(b,nb,a,na,ba,ab,w)

    use size_m
    include 'INPUT'
    real :: b(1),a(1)
    real :: w(1)

    n=na*na*na
    do i=1,n
        b(i) = a(i)
    enddo

    return
    end subroutine specx
!-----------------------------------------------------------------------
    subroutine phys_to_elsasser(u1,u2,u3,b1,b2,b3,n)

    real :: u1(1),u2(1),u3(1),b1(1),b2(1),b3(1)

    do i=1,n
    
        zpx = u1(i) + b1(i)
        zpy = u2(i) + b2(i)
        zpz = u3(i) + b3(i)
    
        zmx = u1(i) - b1(i)
        zmy = u2(i) - b2(i)
        zmz = u3(i) - b3(i)
    
        u1(i) = zpx
        u2(i) = zpy
        u3(i) = zpz
    
        b1(i) = zmx
        b2(i) = zmy
        b3(i) = zmz
    
    enddo

    return
    end subroutine phys_to_elsasser
!-----------------------------------------------------------------------
    subroutine elsasser_to_phys(u1,u2,u3,b1,b2,b3,n)

    real :: u1(1),u2(1),u3(1),b1(1),b2(1),b3(1)

    do i=1,n
    
        zpx = 0.5*( u1(i) + b1(i) )
        zpy = 0.5*( u2(i) + b2(i) )
        zpz = 0.5*( u3(i) + b3(i) )
    
        zmx = 0.5*( u1(i) - b1(i) )
        zmy = 0.5*( u2(i) - b2(i) )
        zmz = 0.5*( u3(i) - b3(i) )
    
        u1(i) = zpx
        u2(i) = zpy
        u3(i) = zpz
    
        b1(i) = zmx
        b2(i) = zmy
        b3(i) = zmz
    
    enddo

    return
    end subroutine elsasser_to_phys
!-----------------------------------------------------------------------
    subroutine phys_to_elsasser2(u1,b1,n)

    real :: u1(1),b1(1)

    do i=1,n
        zpx = u1(i) + b1(i)
        zmx = u1(i) - b1(i)
        u1(i) = zpx
        b1(i) = zmx
    enddo

    return
    end subroutine phys_to_elsasser2
!-----------------------------------------------------------------------
    subroutine elsasser_to_phys2(u1,b1,n)

    real :: u1(1),b1(1)

    do i=1,n
        zpx = 0.5*( u1(i) + b1(i) )
        zmx = 0.5*( u1(i) - b1(i) )
        u1(i) = zpx
        b1(i) = zmx
    enddo

    return
    end subroutine elsasser_to_phys2
!-----------------------------------------------------------------------
    subroutine elsasserh(igeom)


!     Solve MHD in Elsasser variables


    use size_m
    include 'INPUT'
    include 'EIGEN'
    include 'SOLN'
    include 'TSTEP'
    include 'MASS'
    include 'GEOM'

    common /scrnt/  besv1 (lbx1,lby1,lbz1,lbelv) &
    ,              besv2 (lbx1,lby1,lbz1,lbelv) &
    ,              besv3 (lbx1,lby1,lbz1,lbelv)
    COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV) &
    ,              RESV2 (LX1,LY1,LZ1,LELV) &
    ,              RESV3 (LX1,LY1,LZ1,LELV) &
    ,              DV1   (LX1,LY1,LZ1,LELV) &
    ,              DV2   (LX1,LY1,LZ1,LELV) &
    ,              DV3   (LX1,LY1,LZ1,LELV)
    COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV) &
    ,              H2    (LX1,LY1,LZ1,LELV)

    n  = nx1*ny1*nz1*nelv

!     New geometry, new velocity

    intype = -1

    ifield = 1
    call sethlm   (h1,h2,intype)
    call cresvibp (resv1,resv2,resv3,h1,h2)

    ifield = ifldmhd
    call sethlm   (h1,h2,intype)
    call cresvib  (besv1,besv2,besv3,h1,h2)


    ifield = 1
    call sethlm   (h1,h2,intype)

    call ophinv_pr(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxh)

    call opadd2   (vx,vy,vz,dv1,dv2,dv3)

    if (param(103) > 0) alpha_filt=param(103)      ! Optional Filtering
    if (param(103) > 0) call q_filter(alpha_filt)

    call incomprn (vx,vy,vz,pr) ! project U onto div-free space

    ifield = ifldmhd
    call sethlm   (h1,h2,intype)

    call ophinv_pr(dv1,dv2,dv3,besv1,besv2,besv3,h1,h2,tolhv,nmxh)
    call opadd2   (bx,by,bz,dv1,dv2,dv3)


!     if (param(103).gt.0) call q_filter(alpha_filt)

    call incomprn (bx,by,bz,pm) ! project B onto div-free space

    return
    end subroutine elsasserh
!--------------------------------------------------------------------
    subroutine compute_cfl(cfl,u,v,w,dt)

!     Given velocity field (u,v,w) and dt, compute current CFL number.

    use size_m
    include 'GEOM'
    include 'INPUT'
    include 'WZ'
    include 'SOLN'

    real :: u(nx1,ny1,nz1,nelv),v(nx1,ny1,nz1,nelv),w(nx1,ny1,nz1,nelv)

!     Store the inverse jacobian to speed up this operation

    common /cfldx/ dri(lx1),dsi(ly1),dti(lz1)

    integer :: e

    integer :: icalld
    save    icalld
    data    icalld /0/

    if (icalld == 0) then
        icalld=1
        call getdr(dri,zgm1(1,1),nx1)
        call getdr(dsi,zgm1(1,2),ny1)
        if (if3d) call getdr(dti,zgm1(1,3),nz1)
    endif

    cfl = 0.
    l   = 0

    if (if3d) then
        nxyz = nx1*ny1*nz1
        do e=1,nelv
            do k=1,nz1
                do j=1,ny1
                    do i=1,nx1
                        l = l+1
                        ur = ( u(i,j,k,e)*rxm1(i,j,k,e) &
                        +   v(i,j,k,e)*rym1(i,j,k,e) &
                        +   w(i,j,k,e)*rzm1(i,j,k,e) ) * jacmi(l,1)
                        us = ( u(i,j,k,e)*sxm1(i,j,k,e) &
                        +   v(i,j,k,e)*sym1(i,j,k,e) &
                        +   w(i,j,k,e)*szm1(i,j,k,e) ) * jacmi(l,1)
                        ut = ( u(i,j,k,e)*txm1(i,j,k,e) &
                        +   v(i,j,k,e)*tym1(i,j,k,e) &
                        +   w(i,j,k,e)*tzm1(i,j,k,e) ) * jacmi(l,1)
                         
                        cflr = abs(dt*ur*dri(i))
                        cfls = abs(dt*us*dsi(j))
                        cflt = abs(dt*ut*dti(k))
                         
                        cflm = cflr + cfls + cflt
                        cfl  = max(cfl,cflm)

                        cflf(i,j,k,e) = cflm
                         
                    enddo
                enddo
            enddo
        enddo
    else
        nxyz = nx1*ny1
        do e=1,nelv
            do j=1,ny1
                do i=1,nx1
                    l = l+1
                    ur = ( u(i,j,1,e)*rxm1(i,j,1,e) &
                    +   v(i,j,1,e)*rym1(i,j,1,e) ) * jacmi(l,1)
                    us = ( u(i,j,1,e)*sxm1(i,j,1,e) &
                    +   v(i,j,1,e)*sym1(i,j,1,e) ) * jacmi(l,1)

                    cflr = abs(dt*ur*dri(i))
                    cfls = abs(dt*us*dsi(j))

                    cflm = cflr + cfls
                    cfl  = max(cfl,cflm)

                    cflf(i,j,1,e) = cflm

                enddo
            enddo
        enddo
    endif

    cfl = glmax(cfl,1)

    return
    end subroutine compute_cfl
!-----------------------------------------------------------------------
    subroutine getdr(dri,zgm1,nx1)
    real :: dri(nx1),zgm1(nx1)

    dri(1) = zgm1(2) - zgm1(1)   !  Compute 1/Dx
    do i=2,nx1-1
        dri(i) = 0.5*( zgm1(i+1) - zgm1(i-1) )
    enddo
    dri(nx1) = zgm1(nx1) - zgm1(nx1-1)

    call invcol1(dri,nx1)

    return
    end subroutine getdr
!-----------------------------------------------------------------------
    subroutine ophinv_pr(o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)

!     Ok = (H1*A+H2*B)-1 * Ik  (implicit)

    use size_m
    include 'INPUT'
    include 'MASS'
    include 'SOLN'
!max    include 'ORTHOV'
    include 'VPROJ'
    include 'TSTEP'

    parameter (ktot = lx1*ly1*lz1*lelt)

    common /vrthoi/ mprev,nprev(ldim)
    common /vrthov/ vbar(ktot,ldim),vnew(ktot,ldim) &
                 , sln (ktot*mxprev,ldim)
    common /vrthos/ alpha(mxprev)
    common /vrthol/ ifproj
    logical         ifproj

    real :: o1 (lx1,ly1,lz1,1) , o2 (lx1,ly1,lz1,1) , o3 (lx1,ly1,lz1,1)
    real :: i1 (lx1,ly1,lz1,1) , i2 (lx1,ly1,lz1,1) , i3 (lx1,ly1,lz1,1)
    real :: h1 (lx1,ly1,lz1,1) , h2 (lx1,ly1,lz1,1)

    ifproj = .FALSE. 
    if (param(94) > 0) ifproj = .TRUE. 

    if ( .NOT. ifproj .OR. .NOT. if3d) then
        if (ifield == 1) call ophinvm &
        (o1,o2,o3,i1,i2,i3,v1mask,v2mask,v3mask,h1,h2,tolh,nmxhi)
        if (ifield == ifldmhd) call ophinvm &
        (o1,o2,o3,i1,i2,i3,b1mask,b2mask,b3mask,h1,h2,tolh,nmxhi)
        return
    endif

    do i=1,2*ndim
        mtmp        = param(93)
        ivproj(1,i) = min(mxprev,mtmp) - 1
    enddo

    imesh = 1

    if (ifstrs) then
#if 0
        matmod = 0
        call hmhzsf  ('nomg',o1,o2,o3,i1,i2,i3,h1,h2, &
        v1mask,v2mask,v3mask,vmult, &
        tolh,nmxhi,matmod)
#endif
    else
        if (ifield == 1) then
            call hsolve ('velx',o1,i1,h1,h2,v1mask,vmult &
            ,imesh,tolh,nmxhi,1 &
            ,vproj(1,1),ivproj(1,1),binvm1)
            call hsolve ('vely',o2,i2,h1,h2,v2mask,vmult &
            ,imesh,tolh,nmxhi,2 &
            ,vproj(1,2),ivproj(1,2),binvm1)
            if (if3d) &
            call hsolve ('velz',o3,i3,h1,h2,v3mask,vmult &
            ,imesh,tolh,nmxhi,3 &
            ,vproj(1,3),ivproj(1,3),binvm1)
        else  ! B-field
            call hsolve (' Bx ',o1,i1,h1,h2,b1mask,vmult &
            ,imesh,tolh,nmxhi,1 &
            ,vproj(1,4),ivproj(1,4),binvm1)
            call hsolve (' By ',o2,i2,h1,h2,b2mask,vmult &
            ,imesh,tolh,nmxhi,2 &
            ,vproj(1,5),ivproj(1,5),binvm1)
            if (if3d) &
            call hsolve (' Bz ',o3,i3,h1,h2,b3mask,vmult &
            ,imesh,tolh,nmxhi,3 &
            ,vproj(1,6),ivproj(1,6),binvm1)
        endif
    endif

    return
    end subroutine ophinv_pr
!--------------------------------------------------------------------
    subroutine ophinvm(o1,o2,o3,i1,i2,i3,m1,m2,m3,h1,h2,tolh,nmxhi)

!     Ok  = (H1*A+H2*B)-1 * Ik   (implicit)

    use size_m
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'
    real :: o1(lx1,ly1,lz1,1), o2(lx1,ly1,lz1,1), o3(lx1,ly1,lz1,1)
    real :: i1(lx1,ly1,lz1,1), i2(lx1,ly1,lz1,1), i3(lx1,ly1,lz1,1)
    real :: m1(lx1,ly1,lz1,1), m2(lx1,ly1,lz1,1), m3(lx1,ly1,lz1,1)
    real :: h1(lx1,ly1,lz1,1), h2(lx1,ly1,lz1,1)

    imesh = 1

    if (ifstrs) then
#if 0
        matmod = 0
        call hmhzsf  ('NOMG',o1,o2,o3,i1,i2,i3,h1,h2, &
        m1,m2,m3,vmult,tolh,nmxhi,matmod)
#endif
    elseif (ifield == 1) then
        call hmholtz ('VELX',o1,i1,h1,h2,m1,vmult,imesh,tolh,nmxhi,1)
        call hmholtz ('VELY',o2,i2,h1,h2,m2,vmult,imesh,tolh,nmxhi,2)
        if (ndim == 3) &
        call hmholtz ('VELZ',o3,i3,h1,h2,m3,vmult,imesh,tolh,nmxhi,3)
    elseif (ifield == ifldmhd) then
        call hmholtz (' BX ',o1,i1,h1,h2,m1,vmult,imesh,tolh,nmxhi,1)
        call hmholtz (' BY ',o2,i2,h1,h2,m2,vmult,imesh,tolh,nmxhi,2)
        if (ndim == 3) &
        call hmholtz (' BZ ',o3,i3,h1,h2,m3,vmult,imesh,tolh,nmxhi,3)
    endif

    return
    end subroutine ophinvm
!--------------------------------------------------------------------
#if 0
    subroutine set_ifbcor
    use size_m
    include 'GEOM'
    include 'INPUT'
!     include 'TSTEP'   ! ifield?

    common  /nekcb/ cb
    character cb*3

    ifbcor = .TRUE. 
    ifield = ifldmhd

    nface  = 2*ndim
    do iel=1,nelv
        do ifc=1,nface
            cb = cbc(ifc,iel,ifield)
            if  (cb == 'ndd' .OR. cb == 'dnd' .OR. cb == 'ddn') &
            ifbcor = .FALSE. 
        enddo
    enddo

    call gllog(ifbcor , .FALSE. )

    if (nid == 0)  write (6,*) 'IFBCOR   =',ifbcor

    return
    end subroutine set_ifbcor
#endif
!--------------------------------------------------------------------
!      subroutine set_ifbcor_old(ifbcor)
!c
!c     This is a quick hack for the rings problem - it is not general,
!c     but will also work fine for the periodic in Z problem
!c
!      use size_m
!      include 'TOTAL'

!      logical ifbcor

!      integer e,f
!      character*1 cb1(3)

!      itest = 0

!      do e=1,nelv
!      do f=1,2*ndim
!         call chcopy(cb1,cbc(f,e,ifldmhd),3)
!         if (cb1(3).eq.'n'.or.cb1(3).eq.'N') itest = 1
!      enddo
!      enddo

!      itest = iglmax(itest,1)

!      ifbcor = .true.  ! adjust mean pressure to rmv hydrostatic mode


!      if (itest.eq.1) ifbcor = .false.  ! do not adjust mean pressure

!      return
!      end
!--------------------------------------------------------------------
    subroutine setrhsp(p,h1,h2,h2inv,pset,nprev)

!     Project soln onto best fit in the "E" norm.

    use size_m
    include 'INPUT'
    include 'MASS'
    include 'SOLN'
    include 'TSTEP'

    real :: p    (lx2,ly2,lz2,lelv)
    real :: h1   (lx1,ly1,lz1,lelv)
    real :: h2   (lx1,ly1,lz1,lelv)
    real :: h2inv(lx1,ly1,lz1,lelv)
    real :: pset (lx2*ly2*lz2*lelv,mxprev)

    parameter (ltot2=lx2*ly2*lz2*lelv)
    common /orthox/ pbar(ltot2),pnew(ltot2)
    common /orthos/ alpha(mxprev),work(mxprev)

    if (nprev == 0) return

!     Diag to see how much reduction in the residual is attained.
    ntot2  = nx2*ny2*nz2*nelv
    alpha1 = glsc3(p,p,bm2inv,ntot2)
    if (alpha1 > 0) then
        alpha1 = sqrt(alpha1/volvm2)
    else
        return
    endif

!     CALL UPDRHSE(P,H1,H2,H2INV,ierr) ! update rhs's if E-matrix has changed
!     if (ierr.eq.1) Nprev=0           ! Doesn't happen w/ new formulation

    do i=1,nprev  ! Perform Gram-Schmidt for previous soln's.
        alpha(i) = vlsc2(p,pset(1,i),ntot2)
    enddo
    call gop(alpha,work,'+  ',nprev)

    call rzero(pbar,ntot2)
    do i=1,nprev
        call add2s2(pbar,pset(1,i),alpha(i),ntot2)
    enddo

    intetype = 1
    call cdabdtp(pnew,pbar,h1,h2,h2inv,intetype)
    call sub2   (p,pnew,ntot2)

!    ................................................................
    alpha2 = glsc3(p,p,bm2inv,ntot2) ! Diagnostics
    if (alpha2 > 0) then
        alpha2 = sqrt(alpha2/volvm2)
        ratio  = alpha1/alpha2
        n10=min(10,nprev)
        if (nid == 0) write(6,11) istep,nprev,(alpha(i),i=1,n10)
        if (nid == 0) write(6,12) istep,nprev,alpha1,alpha2,ratio
        11 format(2i5,' alpha:',1p10e12.4)
        12 format(i6,i4,1p3e12.4,' alph12')
    endif
!    ................................................................

    return
    end subroutine setrhsp
!-----------------------------------------------------------------------
    subroutine gensolnp(p,h1,h2,h2inv,pset,nprev)

!     Reconstruct the solution to the original problem by adding back
!     the previous solutions

    use size_m
    include 'INPUT'

    real :: p    (lx2,ly2,lz2,lelv)
    real :: h1   (lx1,ly1,lz1,lelv)
    real :: h2   (lx1,ly1,lz1,lelv)
    real :: h2inv(lx1,ly1,lz1,lelv)
    real :: pset (lx2*ly2*lz2*lelv,mxprev)

    parameter (ltot2=lx2*ly2*lz2*lelv)
    common /orthox/ pbar(ltot2),pnew(ltot2)
    common /orthos/ alpha(mxprev),work(mxprev)

    mprev=param(93)
    mprev=min(mprev,mxprev)

    ntot2=nx2*ny2*nz2*nelv

    if (nprev < mprev) then
        nprev = nprev+1
        call copy  (pset(1,nprev),p,ntot2)        ! Save current solution
        call add2  (p,pbar,ntot2)                 ! Reconstruct solution.
        call econjp(pset,nprev,h1,h2,h2inv,ierr)  ! Orthonormalize set
    else                                         !          (uses pnew).
        nprev = 1
        call add2  (p,pbar,ntot2)                 ! Reconstruct solution.
        call copy  (pset(1,nprev),p,ntot2)        ! Save current solution
        call econjp(pset,nprev,h1,h2,h2inv,ierr)  !   and orthonormalize.
    endif

    return
    end subroutine gensolnp
!-----------------------------------------------------------------------
    subroutine econjp(pset,nprev,h1,h2,h2inv,ierr)

!     Orthogonalize the soln wrt previous soln's for which we already
!     know the soln.

    use size_m
    include 'INPUT'

    real :: p    (lx2,ly2,lz2,lelv)
    real :: h1   (lx1,ly1,lz1,lelv)
    real :: h2   (lx1,ly1,lz1,lelv)
    real :: h2inv(lx1,ly1,lz1,lelv)
    real :: pset (lx2*ly2*lz2*lelv,mxprev)

    parameter (ltot2=lx2*ly2*lz2*lelv)
    common /orthox/ pbar(ltot2),pnew(ltot2)
    common /orthos/ alpha(mxprev),work(mxprev)

    ierr  = 0

    ntot2=nx2*ny2*nz2*nelv


!     Gram Schmidt, w re-orthogonalization

    npass=1
!     if (abs(param(105)).eq.2) npass=2
    do ipass=1,npass

        intetype=1
        call cdabdtp(pnew,pset(1,nprev),h1,h2,h2inv,intetype)
        alphad = glsc2(pnew,pset(1,nprev),ntot2) ! compute part of the norm

        nprev1 = nprev-1
        do i=1,nprev1   !   Gram-Schmidt
            alpha(i) = vlsc2(pnew,pset(1,i),ntot2)
        enddo
        if (nprev1 > 0) call gop(alpha,work,'+  ',nprev1)

        do i=1,nprev1
            alpham = -alpha(i)
            call add2s2(pset(1,nprev),pset(1,i),alpham,ntot2)
            alphad = alphad - alpha(i)**2
        enddo
    enddo

!    .Normalize new element in P~

    if (alphad <= 0) then
        write(6,*) 'ERROR:  alphad <= 0 in econjp',alphad,nprev
        ierr = 1
        return
    endif
    alphad = 1./sqrt(alphad)
    call cmult(pset(1,nprev),alphad,ntot2)

    return
    end subroutine econjp
!-----------------------------------------------------------------------
    subroutine advab_elsasser_fast

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

    use size_m
    include 'INPUT'
    include 'GEOM'
    include 'SOLN'
    include 'MASS'
    include 'TSTEP'

    parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
    common /scrns/ wk(2*ltd) &
    , fx(lxy),fy(lxy),fz(lxy) &
    , gx(lxy),gy(lxy),gz(lxy) &
    , zr(ltd),zs(ltd),zt(ltd) &
    , tr(ltd,3),zp(ltd,3),zm(ltd,3)

    integer :: e
    integer :: icalld
    save    icalld
    data    icalld /0/

    if (icalld == 0) call set_dealias_rx
    icalld = icalld + 1

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    if (icalld == 1 .AND. nid == 0) write(6,*) 'inside fast',nxyz1,nxyzd

    if (if3d) then
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

        !           write(6,*) nid,' inside fast',e,nxyz1

            call add3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,1),wk,nx1,nxd,if3d,0) ! 0 --> forward

            call add3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,2),wk,nx1,nxd,if3d,0)

            call add3(wk,vz(1,1,1,e),bz(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,3),wk,nx1,nxd,if3d,0)

            call sub3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,1),wk,nx1,nxd,if3d,0)

            call sub3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,2),wk,nx1,nxd,if3d,0)

            call sub3(wk,vz(1,1,1,e),bz(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,3),wk,nx1,nxd,if3d,0)

            do i=1,nxyzd  ! Convert convector (zm) to r-s-t coordinates
                tr(i,1)= &
                rx(i,1,e)*zm(i,1)+rx(i,2,e)*zm(i,2)+rx(i,3,e)*zm(i,3)
                tr(i,2)= &
                rx(i,4,e)*zm(i,1)+rx(i,5,e)*zm(i,2)+rx(i,6,e)*zm(i,3)
                tr(i,3)= &
                rx(i,7,e)*zm(i,1)+rx(i,8,e)*zm(i,2)+rx(i,9,e)*zm(i,3)
            enddo


            call grad_rst(zr,zs,zt,zp(1,1),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fx,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zp(1,2),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fy,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zp(1,3),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fz,wk,nx1,nxd,if3d,1) ! Project back to coarse


            do i=1,nxyzd  ! Convert convector (zp) to r-s-t coordinates
                tr(i,1)= &
                rx(i,1,e)*zp(i,1)+rx(i,2,e)*zp(i,2)+rx(i,3,e)*zp(i,3)
                tr(i,2)= &
                rx(i,4,e)*zp(i,1)+rx(i,5,e)*zp(i,2)+rx(i,6,e)*zp(i,3)
                tr(i,3)= &
                rx(i,7,e)*zp(i,1)+rx(i,8,e)*zp(i,2)+rx(i,9,e)*zp(i,3)
            enddo


            call grad_rst(zr,zs,zt,zm(1,1),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gx,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zm(1,2),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gy,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zm(1,3),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gz,wk,nx1,nxd,if3d,1) ! Project back to coarse

            tmp = -0.5 ! vtrans() assumed to be 1.0 !
            do i=1,nxyz1
                bfx(i,1,1,e) = bfx(i,1,1,e) + tmp*( fx(i) + gx(i) )
                bmx(i,1,1,e) = bmx(i,1,1,e) + tmp*( fx(i) - gx(i) )

                bfy(i,1,1,e) = bfy(i,1,1,e) + tmp*( fy(i) + gy(i) )
                bmy(i,1,1,e) = bmy(i,1,1,e) + tmp*( fy(i) - gy(i) )

                bfz(i,1,1,e) = bfz(i,1,1,e) + tmp*( fz(i) + gz(i) )
                bmz(i,1,1,e) = bmz(i,1,1,e) + tmp*( fz(i) - gz(i) )
            enddo

        enddo

    else ! 2D
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call add3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,1),wk,nx1,nxd,if3d,0) ! 0 --> forward

            call add3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,2),wk,nx1,nxd,if3d,0)

            call sub3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,1),wk,nx1,nxd,if3d,0)

            call sub3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,2),wk,nx1,nxd,if3d,0)

            do i=1,nxyzd  ! Convert convector (zm) to r-s-t coordinates
                tr(i,1)= &
                rx(i,1,e)*zm(i,1)+rx(i,2,e)*zm(i,2)
                tr(i,2)= &
                rx(i,3,e)*zm(i,1)+rx(i,4,e)*zm(i,2)
            enddo

            call grad_rst(zr,zs,zt,zp(1,1),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(fx,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zp(1,2),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(fy,wk,nx1,nxd,if3d,1) ! Project back to coarse

            do i=1,nxyzd  ! Convert convector (zp) to r-s-t coordinates
                tr(i,1)= &
                rx(i,1,e)*zp(i,1)+rx(i,2,e)*zp(i,2)
                tr(i,2)= &
                rx(i,3,e)*zp(i,1)+rx(i,4,e)*zp(i,2)
            enddo

            call grad_rst(zr,zs,zt,zm(1,1),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(gx,wk,nx1,nxd,if3d,1) ! Project back to coarse

            call grad_rst(zr,zs,zt,zm(1,2),nxd,if3d)
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(gy,wk,nx1,nxd,if3d,1) ! Project back to coarse

            tmp = -0.5 ! vtrans() assumed to be 1.0 !
            do i=1,nxyz1
                bfx(i,1,1,e) = bfx(i,1,1,e) + tmp*( fx(i) + gx(i) )
                bmx(i,1,1,e) = bmx(i,1,1,e) + tmp*( fx(i) - gx(i) )

                bfy(i,1,1,e) = bfy(i,1,1,e) + tmp*( fy(i) + gy(i) )
                bmy(i,1,1,e) = bmy(i,1,1,e) + tmp*( fy(i) - gy(i) )
            enddo
        enddo
    endif

    return
    end subroutine advab_elsasser_fast
!-----------------------------------------------------------------------
    subroutine set_dealias_rx

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

    use size_m
    include 'INPUT'
    include 'GEOM'
    include 'TSTEP' ! for istep

    common /dealias1/ zd(lxd),wd(lxd)
    integer :: e

    integer :: ilstep
    save    ilstep
    data    ilstep /-1/

    if ( .NOT. ifgeom .AND. ilstep > 1) return  ! already computed
    if (ifgeom .AND. ilstep == istep)  return  ! already computed
    ilstep = istep

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    call zwgl (zd,wd,nxd)  ! zwgl -- NOT zwgll!

    if (if3d) then
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),rzm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,5,e),sym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,6,e),szm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,7,e),txm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,8,e),tym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,9,e),tzm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd

            l = 0
            do k=1,nzd
                do j=1,nyd
                    do i=1,nxd
                        l = l+1
                        w = wd(i)*wd(j)*wd(k)
                        do ii=1,9
                            rx(l,ii,e) = w*rx(l,ii,e)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    else ! 2D
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),sxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd

            l = 0
            do j=1,nyd
                do i=1,nxd
                    l = l+1
                    w = wd(i)*wd(j)
                    do ii=1,4
                        rx(l,ii,e) = w*rx(l,ii,e)
                    enddo
                enddo
            enddo
        enddo

    endif

    return
    end subroutine set_dealias_rx
!-----------------------------------------------------------------------
    subroutine cfl_check

    use size_m
    include 'INPUT'
    include 'SOLN'
    include 'MASS'
    include 'TSTEP'

    common /scrns/ ta1 (lx1,ly1,lz1,lelv) &
    ,             ta2 (lx1,ly1,lz1,lelv) &
    ,             ta3 (lx1,ly1,lz1,lelv) &
    ,             tb1 (lx1,ly1,lz1,lelv) &
    ,             tb2 (lx1,ly1,lz1,lelv) &
    ,             tb3 (lx1,ly1,lz1,lelv)


    call opcopy(ta1,ta2,ta3,vx,vy,vz)
    call opcopy(tb1,tb2,tb3,bx,by,bz)

    ntot1 = nx1*ny1*nz1*nelv
    call phys_to_elsasser(ta1,ta2,ta3,tb1,tb2,tb3,ntot1) ! crude, but effective

    call compute_cfl(cflp,ta1,ta2,ta3,dt)   !  vx = U+B
    call compute_cfl(cflm,tb1,tb2,tb3,dt)   !  bx = U-B

    courno = max(cflp,cflm)

    if (nid == 0) write(6,1) istep,time,dt,cflp,cflm
    1 format(i9,1p4e15.7,' CFL')

    return
    end subroutine cfl_check
!--------------------------------------------------------------------
