    subroutine load_fld(string)

    use size_m
    use input
    use restart

    character(1) :: string(1),fout(132),BLNK
    character(6) :: ext
    DATA BLNK/' '/

    call blank  (initc(1),132)

    L1=0
    DO 100 I=1,132
        IF (STRING(I) == BLNK) GOTO 200
        L1=I
    100 END DO
    200 CONTINUE
    LEN=L1

    call chcopy (initc(1),string,len)
    call setics

    return
    end subroutine load_fld
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine find_lam3(lam,aa,w,ndim,ierr)
    real :: aa(ndim,ndim),lam(ndim),w(ndim,ndim),lam2

!     Use cubic eqn. to compute roots

    common /ecmnr/ a,b,c,d,e,f,f2,ef,df,r0,r1
    common /ecmni/ nr
    common /ecmnl/ iffout,ifdefl
    logical ::        iffout,ifdefl


    iffout = .FALSE. 
    ierr = 0

!     2D case....


    if (ndim == 2) then
        a = aa(1,1)
        b = aa(1,2)
        c = aa(2,1)
        d = aa(2,2)
        aq = 1.
        bq = -(a+d)
        cq = a*d-c*b
    
        call quadratic(x1,x2,aq,bq,cq,ierr)
    
        lam(1) = min(x1,x2)
        lam(2) = max(x1,x2)
    
        return
    endif

!     Else ...  3D case....
!                                    a d e
!     Get symmetric 3x3 matrix       d b f
!                                    e f c

    a = aa(1,1)
    b = aa(2,2)
    c = aa(3,3)
    d = 0.5*(aa(1,2)+aa(2,1))
    e = 0.5*(aa(1,3)+aa(3,1))
    f = 0.5*(aa(2,3)+aa(3,2))
    ef = e*f
    df = d*f
    f2 = f*f


!     Use cubic eqn. to compute roots

!     ax = a-x
!     bx = b-x
!     cx = c-x
!     y = ax*(bx*cx-f2) - d*(d*cx-ef) + e*(df-e*bx)

    a1 = -(a+b+c)
    a2 =  (a*b+b*c+a*c) - (d*d+e*e+f*f)
    a3 =  a*f*f + b*e*e + c*d*d - a*b*c - 2*d*e*f

    call cubic  (lam,a1,a2,a3,ierr)
    call sort   (lam,w,3)

    return
    end subroutine find_lam3
!-----------------------------------------------------------------------
    subroutine quadratic(x1,x2,a,b,c,ierr)

!     Stable routine for computation of real roots of quadratic

    ierr = 0
    x1 = 0.
    x2 = 0.

    if (a == 0.) then
        if (b == 0) then
            if (c /= 0) then
            !              write(6,10) x1,x2,a,b,c
                ierr = 1
            endif
            return
        endif
        ierr = 2
        x1 = -c/b
    !        write(6,11) x1,a,b,c
        return
    endif

    d = b*b - 4.*a*c
    if (d < 0) then
        ierr = 1
    !        write(6,12) a,b,c,d
        return
    endif
    if (d > 0) d = sqrt(d)

    if (b > 0) then
        x1 = -2.*c / ( d+b )
        x2 = -( d+b ) / (2.*a)
    else
        x1 =  ( d-b ) / (2.*a)
        x2 = -2.*c / ( d-b )
    endif

    10 format('ERROR: Both a & b zero in routine quadratic NO ROOTS.' &
    ,1p5e12.4)
    11 format('ERROR: a = 0 in routine quadratic, only one root.' &
    ,1p5e12.4)
    12 format('ERROR: negative discriminate in routine quadratic.' &
    ,1p5e12.4)

    return
    end subroutine quadratic
!-----------------------------------------------------------------------
    subroutine cubic(xo,ai1,ai2,ai3,ierr)
    real :: xo(3),ai1,ai2,ai3
    complex*16 x(3),a1,a2,a3,q,r,d,arg,t1,t2,t3,theta,sq,a13

!     Compute real solutions to cubic root eqn. (Num. Rec. v. 1, p. 146)
!     pff/Sang-Wook Lee  Jan 19 , 2004

!     Assumption is that all x's are *real*

    real*8 :: twopi
    save   twopi
    data   twopi /6.283185307179586476925286766/

    ierr = 0

    zero = 0.
    a1   = cmplx(ai1,zero)
    a2   = cmplx(ai2,zero)
    a3   = cmplx(ai3,zero)

    q = (a1*a1 - 3*a2)/9.
    if (q == 0) goto 999

    r = (2*a1*a1*a1 - 9*a1*a2 + 27*a3)/54.

    d = q*q*q - r*r

!     if (d.lt.0) goto 999

    arg   = q*q*q
    arg   = sqrt(arg)
    arg   = r/arg

    if (abs(arg) > 1) goto 999
    theta = acos(abs(arg))

    t1    = theta / 3.
    t2    = (theta + twopi) / 3.
    t3    = (theta + 2.*twopi) / 3.

    sq  = -2.*sqrt(q)
    a13 = a1/3.
    x(1) = sq*cos(t1) - a13
    x(2) = sq*cos(t2) - a13
    x(3) = sq*cos(t3) - a13

    xo(1) = real(x(1))
    xo(2) = real(x(2))
    xo(3) = real(x(3))

    return

    999 continue   ! failed
    ierr = 1
    call rzero(x,3)

    return
    end subroutine cubic
!-----------------------------------------------------------------------
    subroutine comp_gije(gije,u,v,w,e)

!                                         du_i
!     Compute the gradient tensor G_ij := ----  ,  for element e
!                                         du_j

    use size_m
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

    real :: gije(lx1*ly1*lz1,ldim,ldim)
    real :: u   (lx1*ly1*lz1)
    real :: v   (lx1*ly1*lz1)
    real :: w   (lx1*ly1*lz1)

    real :: ur  (lx1*ly1*lz1)
    real :: us  (lx1*ly1*lz1)
    real :: ut  (lx1*ly1*lz1)

    integer :: e

    n    = nx1-1      ! Polynomial degree
    nxyz = nx1*ny1*nz1

    if (if3d) then     ! 3D CASE

        do k=1,3
            if (k == 1) call local_grad3(ur,us,ut,u,n,1,dxm1,dxtm1)
            if (k == 2) call local_grad3(ur,us,ut,v,n,1,dxm1,dxtm1)
            if (k == 3) call local_grad3(ur,us,ut,w,n,1,dxm1,dxtm1)

            do i=1,nxyz
                dj = jacmi(i,e)

            ! d/dx
                gije(i,k,1) = dj*( &
                ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e))
            ! d/dy
                gije(i,k,2) = dj*( &
                ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e))
            ! d/dz
                gije(i,k,3) = dj*( &
                ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e))

            enddo
        enddo

    elseif (ifaxis) then   ! AXISYMMETRIC CASE
        if(nid == 0) write(6,*) &
        'ABORT: comp_gije no axialsymmetric support for now'
        call exitt
    else              ! 2D CASE
#if 0
        do k=1,2
            if (k == 1) call local_grad2(ur,us,u,n,1,dxm1,dxtm1)
            if (k == 2) call local_grad2(ur,us,v,n,1,dxm1,dxtm1)
            do i=1,nxyz
                dj = jacmi(i,e)
            ! d/dx
                gije(i,k,1)=dj*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))
            ! d/dy
                gije(i,k,2)=dj*(ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e))
            enddo
        enddo
#endif
    endif

    return
    end subroutine comp_gije
!-----------------------------------------------------------------------
    subroutine filter_s0(scalar,wght,ncut,name5) ! filter scalar field

    use size_m
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

    real :: scalar(1)
    character(5) :: name5

    parameter (l1=lx1*lx1)
    real :: intdv(l1),intuv(l1),intdp(l1),intup(l1),intv(l1),intp(l1)
    save intdv    ,intuv    ,intdp    ,intup    ,intv    ,intp

    common /ctmp0/ intt
    common /screv/ wk1,wk2
    common /scrvh/ zgmv,wgtv,zgmp,wgtp,tmax(100)

    real :: intt (lx1,lx1)
    real :: wk1  (lx1,lx1,lx1,lelt)
    real :: wk2  (lx1,lx1,lx1)
    real :: zgmv (lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)


    integer :: icall
    save    icall
    data    icall /0/

    logical :: ifdmpflt

    imax = nid
    imax = iglmax(imax,1)
    jmax = iglmax(imax,1)

!      if (icall.eq.0) call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
    call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)

    icall = 1

    call filterq(scalar,intv,nx1,nz1,wk1,wk2,intt,if3d,fmax)
    fmax = glmax(fmax,1)

    if (nid == 0) write(6,1) istep,fmax,name5
    1 format(i8,' sfilt:',1pe12.4,a10)

    return
    end subroutine filter_s0
!-----------------------------------------------------------------------
    subroutine intpts_setup(tolin,ih)

! setup routine for interpolation tool
! tolin ... stop point seach interation if 1-norm of the step in (r,s,t)
!           is smaller than tolin

    use size_m
    use geom

    common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

    tol = tolin
    if (tolin < 0) tol = 1e-13 ! default tolerance

    n       = lx1*ly1*lz1*lelt
    npt_max = 256
    nxf     = 2*nx1 ! fine mesh for bb-test
    nyf     = 2*ny1
    nzf     = 2*nz1
    bb_t    = 0.1 ! relative size to expand bounding boxes by

    if(nidd == 0) write(6,*) 'initializing intpts(), tol=', tol
    call findpts_setup(ih,nekcomm,npp,ndim, &
    xm1,ym1,zm1,nx1,ny1,nz1, &
    nelt,nxf,nyf,nzf,bb_t,n,n, &
    npt_max,tol)

    return
    end subroutine intpts_setup
!-----------------------------------------------------------------------
    subroutine intpts(fieldin,nfld,pts,n,fieldout,ifot,ifpts,ih)

! simple wrapper to interpolate input field at given points

! in:
! fieldin ... input field(s) to interpolate (lelt*lxyz,nfld)
! nfld    ... number of fields in fieldin
! pts     ... packed list of interpolation points
!             pts:=[x(1)...x(n),y(1)...y(n),z(1)...z(n)]
! n       ... local number of interpolation points
! ifto    ... transpose output (n,nfld)^T = (nfld,n)
! itpts   ... find interpolation points
! ih      ... interpolation handle

! out:
! fieldout ... packed list of interpolated values (n,nfld)

    use size_m

    real ::    fieldin(1),fieldout(1)
    real ::    pts(1)

    real ::    dist(lpart) ! squared distance
    real ::    rst(lpart*ldim)
    integer :: rcode(lpart),elid(lpart),proc(lpart)

    common /intp_r/ rst,dist
    common /intp_i/ rcode,elid,proc

    integer :: nn(2)
    logical :: ifot,ifpts

    if(nid == 0) write(6,*) 'call intpts'

    if(n > lpart) then
        write(6,*) &
        'ABORT: intpts() n>lpart, increase lpart in SIZE', n, lpart
        call exitt
    endif

! locate points (iel,iproc,r,s,t)
    nfail = 0
    if(ifpts) then
        if(nid == 0) write(6,*) 'call findpts'
        call findpts(ih,rcode,1, &
        proc,1, &
        elid,1, &
        rst,ndim, &
        dist,1, &
        pts(    1),1, &
        pts(  n+1),1, &
        pts(2*n+1),1,n)
        do in=1,n
        ! check return code
            if(rcode(in) == 1) then
                if(dist(in) > 1e-12) then
                    nfail = nfail + 1
                    if (nfail <= 5) write(6,'(a,1p4e15.7)') &
                    ' WARNING: point on boundary or outside the mesh xy[z]d^2: ', &
                    (pts(in+k*n),k=0,ndim-1),dist(in)
                endif
            elseif(rcode(in) == 2) then
                nfail = nfail + 1
                if (nfail <= 5) write(6,'(a,1p3e15.7)') &
                ' WARNING: point not within mesh xy[z]: !', &
                (pts(in+k*n),k=0,ndim-1)
            endif
        enddo
    endif

! evaluate inut field at given points
    ltot = lelt*lx1*ly1*lz1
    do ifld = 1,nfld
        iin    = (ifld-1)*ltot + 1
        iout   = (ifld-1)*n + 1
        is_out = 1
        if(ifot) then ! transpose output
            iout   = ifld
            is_out = nfld
        endif
        call findpts_eval(ih,fieldout(iout),is_out, &
        rcode,1, &
        proc,1, &
        elid,1, &
        rst,ndim,n, &
        fieldin(iin))
    enddo

    nn(1) = iglsum(n,1)
    nn(2) = iglsum(nfail,1)
    if(nid == 0) then
        write(6,1) nn(1),nn(2)
        1 format('   total number of points = ',i12,/,'   failed = ' &
        ,i12,/,' done :: intpts')
    endif

    return
    end subroutine intpts
!-----------------------------------------------------------------------
    subroutine intpts_done(ih)

    call findpts_free(ih)

    return
    end subroutine intpts_done
!-----------------------------------------------------------------------
    subroutine tens3d1(v,u,f,ft,nv,nu)  ! v = F x F x F x u

!     Note: this routine assumes that nx1=ny1=nz1

    use size_m
    use input

    parameter (lw=4*lx1*lx1*lz1)
    common /ctensor/ w1(lw),w2(lw)

    real :: v(nv,nv,nv),u(nu,nu,nu)
    real :: f(1),ft(1)

    if (nu*nu*nv > lw) then
        write(6,*) nid,nu,nv,lw,' ERROR in tens3d1. Increase lw.'
        call exitt
    endif

    if (if3d) then
        nuv = nu*nv
        nvv = nv*nv
        call mxm(f,nv,u,nu,w1,nu*nu)
        k=1
        l=1
        do iz=1,nu
            call mxm(w1(k),nv,ft,nu,w2(l),nv)
            k=k+nuv
            l=l+nvv
        enddo
        call mxm(w2,nvv,ft,nu,v,nv)
    else
        call mxm(f ,nv,u,nu,w1,nu)
        call mxm(w1,nv,ft,nu,v,nv)
    endif
    return
    end subroutine tens3d1
!-----------------------------------------------------------------------
    subroutine build_1d_filt(fh,fht,trnsfr,nx,nid)

!     This routing builds a 1D filter with transfer function diag()

!     Here, nx = number of points

    real :: fh(nx,nx),fht(nx,nx),trnsfr(nx)

    parameter (lm=40)
    parameter (lm2=lm*lm)
    common /cfiltr/ phi(lm2),pht(lm2),diag(lm2),rmult(lm),Lj(lm) &
    , zpts(lm)
    real :: Lj

    common /cfilti/ indr(lm),indc(lm),ipiv(lm)

    if (nx > lm) then
        write(6,*) 'ABORT in set_filt:',nx,lm
        call exitt
    endif

    call zwgll(zpts,rmult,nx)

    kj = 0
    n  = nx-1
    do j=1,nx
        z = zpts(j)
        call legendre_poly(Lj,z,n)
        kj = kj+1
        pht(kj) = Lj(1)
        kj = kj+1
        pht(kj) = Lj(2)
        do k=3,nx
            kj = kj+1
            pht(kj) = Lj(k)-Lj(k-2)
        enddo
    enddo
    call transpose (phi,nx,pht,nx)
    call copy      (pht,phi,nx*nx)
    call gaujordf  (pht,nx,nx,indr,indc,ipiv,ierr,rmult)

    call rzero(diag,nx*nx)
    k=1
    do i=1,nx
        diag(k) = trnsfr(i)
        k = k+(nx+1)
    enddo

    call mxm  (diag,nx,pht,nx,fh,nx)      !          -1
    call mxm  (phi ,nx,fh,nx,pht,nx)      !     V D V

    call copy      (fh,pht,nx*nx)
    call transpose (fht,nx,fh,nx)

    do k=1,nx*nx
        pht(k) = 1.-diag(k)
    enddo
    np1 = nx+1
    if (nid == 0) then
        write(6,6) 'flt amp',(pht (k),k=1,nx*nx,np1)
        write(6,6) 'flt trn',(diag(k),k=1,nx*nx,np1)
        6 format(a8,16f7.4,6(/,8x,16f7.4))
    endif

    return
    end subroutine build_1d_filt
!-----------------------------------------------------------------------
    subroutine comp_sije(gije)

!     Compute symmetric part of a tensor G_ij for element e

    use size_m
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

    real :: gije(lx1*ly1*lz1,ldim,ldim)

    nxyz = nx1*ny1*nz1

    k = 1

    do j=1,ndim
        do i=k,ndim
            do l=1,nxyz
                gije(l,i,j) = 0.5*(gije(l,i,j)+gije(l,j,i))
                gije(l,j,i) = gije(l,i,j)
            enddo
        enddo
        k = k + 1
    enddo

    return
    end subroutine comp_sije
!-----------------------------------------------------------------------
    subroutine map2reg(ur,n,u,nel)

!     Map scalar field u() to regular n x n x n array ur

    use size_m
    real :: ur(1),u(lx1*ly1*lz1,1)

    integer :: e

    ldr = n**ndim

    k=1
    do e=1,nel
        if (ndim == 2) call map2reg_2di_e(ur(k),n,u(1,e),nx1)
        if (ndim == 3) call map2reg_3di_e(ur(k),n,u(1,e),nx1)
        k = k + ldr
    enddo

    return
    end subroutine map2reg
!-----------------------------------------------------------------------
    subroutine map2reg_2di_e(uf,n,uc,m) ! Fine, uniform pt

    real :: uf(n,n),uc(m,m)

    parameter (l=50)
    common /cmap2d/ j(l*l),jt(l*l),w(l*l),z(l)

    integer :: mo,no
    save    mo,no
    data    mo,no / 0,0 /

    if (m > l) call exitti('map2reg_2di_e memory 1$',m)
    if (n > l) call exitti('map2reg_2di_e memory 2$',n)

    if (m /= mo .OR. n /= no ) then

        call zwgll (z,w,m)
        call zuni  (w,n)

        call gen_int_gz(j,jt,w,n,z,m)

    endif

    call mxm(j,n,uc,m,w ,m)
    call mxm(w,n,jt,m,uf,n)

    return
    end subroutine map2reg_2di_e
!-----------------------------------------------------------------------
    subroutine map2reg_3di_e(uf,n,uc,m) ! Fine, uniform pt

    real :: uf(n,n,n),uc(m,m,m)

    parameter (l=50)
    common /cmap3d/ j(l*l),jt(l*l),v(l*l*l),w(l*l*l),z(l)

    integer :: mo,no
    save    mo,no
    data    mo,no / 0,0 /

    if (m > l) call exitti('map2reg_3di_e memory 1$',m)
    if (n > l) call exitti('map2reg_3di_e memory 2$',n)

    if (m /= mo .OR. n /= no ) then

        call zwgll (z,w,m)
        call zuni  (w,n)

        call gen_int_gz(j,jt,w,n,z,m)

    endif

    mm = m*m
    mn = m*n
    nn = n*n

    call mxm(j,n,uc,m,v ,mm)
    iv=1
    iw=1
    do k=1,m
        call mxm(v(iv),n,jt,m,w(iw),n)
        iv = iv+mn
        iw = iw+nn
    enddo
    call mxm(w,nn,jt,m,uf,n)

    return
    end subroutine map2reg_3di_e
!-----------------------------------------------------------------------
    subroutine gen_int_gz(j,jt,g,n,z,m)

!     Generate interpolater from m z points to n g points

!        j   = interpolation matrix, mapping from z to g
!        jt  = transpose of interpolation matrix
!        m   = number of points on z grid
!        n   = number of points on g grid

    real :: j(n,m),jt(m,n),g(n),z(m)

    mpoly  = m-1
    do i=1,n
        call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
    enddo

    call transpose(j,n,jt,m)

    return
    end subroutine gen_int_gz
!-----------------------------------------------------------------------
    subroutine zuni(z,np)

!     Generate equaly spaced np points on the interval [-1:1]

    real :: z(1)

    dz = 2./(np-1)
    z(1) = -1.
    do i = 2,np-1
        z(i) = z(i-1) + dz
    enddo
    z(np) = 1.

    return
    end subroutine zuni
!-----------------------------------------------------------------------
    subroutine gen_re2_xyz
    use size_m
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

    parameter (lv=2**ldim,lblock=1000)
    common /scrns/ xyz(lv,ldim,lblock),wk(lv*ldim*lblock)
    common /scruz/ igr(lblock)

    integer :: e,eb,eg,ierr,wdsiz2

    integer :: isym2pre(8)   ! Symmetric-to-prenek vertex ordering
    save    isym2pre
    data    isym2pre / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

    real*4 :: buf (50)  ! nwds * 2 for double precision
    real ::   buf2(25)  ! double precsn
    equivalence (buf,buf2)


    nxs = nx1-1
    nys = ny1-1
    nzs = nz1-1
     
    wdsiz2=4
    if(wdsize == 8) wdsiz2=8
    nblock = lv*ldim*lblock   !memory size of data blocks

    ierr=0

    do eb=1,nelgt,lblock
        nemax = min(eb+lblock-1,nelgt)
        call rzero(xyz,nblock)
        call izero(igr,lblock)
        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            l   = 0
            if (mid == nid .AND. if3d) then ! fill owning processor
                igr(kb) = igroup(e)
                do k=0,1
                    do j=0,1
                        do i=0,1
                            l=l+1
                            li=isym2pre(l)
                            xyz(li,1,kb) = xm1(1+i*nxs,1+j*nys,1+k*nzs,e)
                            xyz(li,2,kb) = ym1(1+i*nxs,1+j*nys,1+k*nzs,e)
                            xyz(li,3,kb) = zm1(1+i*nxs,1+j*nys,1+k*nzs,e)
                        enddo
                    enddo
                enddo
            elseif (mid == nid) then    ! 2D
                igr(kb) = igroup(e)
                do j=0,1
                    do i=0,1
                        l =l+1
                        li=isym2pre(l)
                        xyz(li,1,kb) = xm1(1+i*nxs,1+j*nys,1,e)
                        xyz(li,2,kb) = ym1(1+i*nxs,1+j*nys,1,e)
                    enddo
                enddo
            endif
        enddo
        call  gop(xyz,wk,'+  ',nblock)  ! Sum across all processors
        call igop(igr,wk,'+  ',lblock)  ! Sum across all processors

        if (nid == 0 .AND. ierr == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                if(wdsiz2 == 8) then
                    rgrp=igr(kb)
                    call byte_write(rgrp,2,ierr)
                                     
                    if(if3d) then
                        buf2(1) = xyz(1,1,kb)
                        buf2(2) = xyz(2,1,kb)
                        buf2(3) = xyz(3,1,kb)
                        buf2(4) = xyz(4,1,kb)

                        buf2(5) = xyz(5,1,kb)
                        buf2(6) = xyz(6,1,kb)
                        buf2(7) = xyz(7,1,kb)
                        buf2(8) = xyz(8,1,kb)

                        buf2(9) = xyz(1,2,kb)
                        buf2(10)= xyz(2,2,kb)
                        buf2(11)= xyz(3,2,kb)
                        buf2(12)= xyz(4,2,kb)
                                          
                        buf2(13)= xyz(5,2,kb)
                        buf2(14)= xyz(6,2,kb)
                        buf2(15)= xyz(7,2,kb)
                        buf2(16)= xyz(8,2,kb)

                        buf2(17)= xyz(1,3,kb)
                        buf2(18)= xyz(2,3,kb)
                        buf2(19)= xyz(3,3,kb)
                        buf2(20)= xyz(4,3,kb)

                        buf2(21)= xyz(5,3,kb)
                        buf2(22)= xyz(6,3,kb)
                        buf2(23)= xyz(7,3,kb)
                        buf2(24)= xyz(8,3,kb)

                        if(ierr == 0) call byte_write(buf,48,ierr)
                    else
                        buf2(1) = xyz(1,1,kb)
                        buf2(2) = xyz(2,1,kb)
                        buf2(3) = xyz(3,1,kb)
                        buf2(4) = xyz(4,1,kb)

                        buf2(5) = xyz(1,2,kb)
                        buf2(6) = xyz(2,2,kb)
                        buf2(7) = xyz(3,2,kb)
                        buf2(8) = xyz(4,2,kb)
                            
                        if(ierr == 0) call byte_write(buf,16,ierr)
                    endif
                else  !!!! 4byte precision !!!!
                    call byte_write(igr(kb),1,ierr)
                    if (if3d) then

                        buf(1)  = xyz(1,1,kb)
                        buf(2)  = xyz(2,1,kb)
                        buf(3)  = xyz(3,1,kb)
                        buf(4)  = xyz(4,1,kb)

                        buf(5)  = xyz(5,1,kb)
                        buf(6)  = xyz(6,1,kb)
                        buf(7)  = xyz(7,1,kb)
                        buf(8)  = xyz(8,1,kb)

                        buf(9)  = xyz(1,2,kb)
                        buf(10) = xyz(2,2,kb)
                        buf(11) = xyz(3,2,kb)
                        buf(12) = xyz(4,2,kb)
                                          
                        buf(13) = xyz(5,2,kb)
                        buf(14) = xyz(6,2,kb)
                        buf(15) = xyz(7,2,kb)
                        buf(16) = xyz(8,2,kb)

                        buf(17) = xyz(1,3,kb)
                        buf(18) = xyz(2,3,kb)
                        buf(19) = xyz(3,3,kb)
                        buf(20) = xyz(4,3,kb)

                        buf(21) = xyz(5,3,kb)
                        buf(22) = xyz(6,3,kb)
                        buf(23) = xyz(7,3,kb)
                        buf(24) = xyz(8,3,kb)

                        if(ierr == 0) call byte_write(buf,24,ierr)

                    else ! 2D
                        buf(1)  = xyz(1,1,kb)
                        buf(2)  = xyz(2,1,kb)
                        buf(3)  = xyz(3,1,kb)
                        buf(4)  = xyz(4,1,kb)

                        buf(5)  = xyz(1,2,kb)
                        buf(6)  = xyz(2,2,kb)
                        buf(7)  = xyz(3,2,kb)
                        buf(8)  = xyz(4,2,kb)

                        if(ierr == 0) call byte_write(buf,8,ierr)
                    endif
                endif
                              
            enddo
        endif
    enddo
    call err_chk(ierr,'Error writing to newre2.re2 in gen_re2_xyz$')

    return
    end subroutine gen_re2_xyz
!-----------------------------------------------------------------------
    subroutine gen_re2_curve(imid)

!     This routine is complex because we must first count number of
!     nontrivial curved sides.

!     A two pass strategy is used:  first count, then write

    use size_m
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

    integer :: e,eb,eg,wdsiz2
    character(1) :: cc(4)

    real*4 :: buf (16)  ! nwds * 2 for double precision
    real ::   buf2( 8)  ! double precsn
    equivalence (buf,buf2)

    parameter (lblock=500)
    common /scrns/ vcurve(5,12,lblock),wk(5*12*lblock)
    common /scruz/ icurve(12,lblock)

    wdsiz2=4
    if(wdsize == 8) wdsiz2=8

    nblock = lv*ldim*lblock   !memory size of data blocks

    if (imid > 0) then

    !        imid = 0  ! No midside node defs
    !        imid = 1  ! Midside defs where current curve sides don't exist
    !        imid = 2  ! All nontrivial midside node defs

        if (imid == 2) call blank(ccurve,12*lelt)

        do e=1,nelt
            call gen_rea_midside_e(e)
        enddo

    endif
    nedge = 4 + 8*(ndim-2)

    ncurvn = 0
    do e=1,nelt
        do i=1,nedge
            if (ccurve(i,e) /= ' ') ncurvn = ncurvn+1
        enddo
    enddo
    ncurvn = iglsum(ncurvn,1)
     
    ierr=0

    if(nid == 0 .AND. wdsiz2 == 8) then
        rcurvn = ncurvn
        call byte_write(rcurvn,2,ierr)
    elseif(nid == 0) then
        call byte_write(ncurvn,1,ierr)
    endif

    do eb=1,nelgt,lblock

        nemax = min(eb+lblock-1,nelgt)
        call izero(icurve,12*lblock)
        call rzero(vcurve,60*lblock)

        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            if (mid == nid) then ! fill owning processor
                do i=1,nedge
                    icurve(i,kb) = 0
                    if (ccurve(i,e) == 'C') icurve(i,kb) = 1
                    if (ccurve(i,e) == 's') icurve(i,kb) = 2
                    if (ccurve(i,e) == 'm') icurve(i,kb) = 3
                    call copy(vcurve(1,i,kb),curve(1,i,e),5)
                enddo
            endif
        enddo
        call igop(icurve,wk,'+  ',12*lblock)  ! Sum across all processors
        call  gop(vcurve,wk,'+  ',60*lblock)  ! Sum across all processors

        if (nid == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                do i=1,nedge
                    ii = icurve(i,kb)   ! equivalenced to s4
                    if (ii /= 0) then
                        if (ii == 1) cc(1)='C'
                        if (ii == 2) cc(1)='s'
                        if (ii == 3) cc(1)='m'

                        if(wdsiz2 == 8) then

                            buf2(1) = eg
                            buf2(2) = i
                            call copy  (buf2(3),vcurve(1,i,kb),5)!real*8 write
                            call blank(buf2(8),8)
                            call chcopy(buf2(8),cc,4)
                            iz = 16
                        else
                            call icopy(buf(1),eg,1)
                            call icopy(buf(2), i,1)
                            call copyX4(buf(3),vcurve(1,i,kb),5) !real*4 write
                            call blank(buf(8),4)
                            call chcopy(buf(8),cc,4)
                            iz = 8
                        endif

                        if(ierr == 0) call byte_write(buf,iz,ierr)
                    endif
                enddo
            enddo
        endif

    enddo
    call err_chk(ierr,'Error writing to newre2.re2 in gen_re2_curve$')

    return
    end subroutine gen_re2_curve
!-----------------------------------------------------------------------
    subroutine gen_re2_bc (ifld)

    use size_m
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

    integer :: e,eb,eg,wdsiz2

    parameter (lblock=500)
    common /scrns/ vbc(5,6,lblock),wk(5*6*lblock)
    common /scruz/ ibc(6,lblock)

    character(1) :: s4(4)
    character(3) :: s3
    integer ::     i4
    equivalence(i4,s4)
    equivalence(s3,s4)

    character(1) :: chtemp
    save        chtemp
    data        chtemp /' '/   ! For mesh bcs

    real*4 :: buf (16)  ! nwds * 2 for double precision
    real ::   buf2( 8)  ! double precsn
    equivalence (buf,buf2)


    nface = 2*ndim
    ierr = 0
    nbc  = 0
    rbc  = 0

    wdsiz2=4
    if(wdsize == 8) wdsiz2=8

    nlg = nelg(ifld)

    if (ifld == 1 .AND. .NOT. ifflow) then ! NO B.C.'s for this field
        if(nid == 0 .AND. wdsiz2 == 4) call byte_write(nbc,1,ierr)
        if(nid == 0 .AND. wdsiz2 == 8) call byte_write(rbc,2,ierr)
        call err_chk(ierr,'Error writing to newre2.re2 in gen_re2_bc$')
        return
    endif

    do ii = 1,nelt
        do jj = 1,nface
            if(cbc(jj,ii,ifld) /= 'E  ')nbc=nbc+1
        enddo
    enddo
    call igop(nbc,wk,'+  ', 1       )  ! Sum across all processors
    if(nid == 0 .AND. wdsiz2 == 8) then
        rbc = nbc
        call byte_write(rbc,2,ierr)
    elseif(nid == 0) then
        call byte_write(nbc,1,ierr)
    endif

    do eb=1,nlg,lblock
        nemax = min(eb+lblock-1,nlg)
        call izero(ibc, 6*lblock)
        call rzero(vbc,30*lblock)
        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            if (mid == nid) then ! fill owning processor
                do i=1,nface
                    i4 = 0
                    call chcopy(s4,cbc(i,e,ifld),3)
                    ibc(i,kb) = i4
                    call copy(vbc(1,i,kb),bc(1,i,e,ifld),5)
                enddo
            endif
        enddo
        call igop(ibc,wk,'+  ', 6*lblock)  ! Sum across all processors
        call  gop(vbc,wk,'+  ',30*lblock)  ! Sum across all processors

        if (nid == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                do i=1,nface
                    i4 = ibc(i,kb)   ! equivalenced to s4
                    if (s3 /= 'E  ') then

                        if(wdsiz2 == 8) then
                            buf2(1)=eg
                            buf2(2)=i
                            call copy    (buf2(3),vbc(1,i,eg),5)
                            call blank   (buf2(8),8)
                            call chcopy  (buf2(8),s3,3)
                            if(nlg >= 1000000) then
                                call icopy(i_vbc,vbc(1,i,eg),1)
                                buf2(3)=i_vbc
                            endif
                            iz=16
                        else
                            call icopy   (buf(1),eg,1)
                            call icopy   (buf(2),i,1)
                            call copyX4  (buf(3),vbc(1,i,eg),5)
                            call blank   (buf(8),4)
                            if(nlg >= 1000000)call icopy(buf(3),vbc(1,i,eg),1)
                            call chcopy  (buf(8),s3,3)
                            iz=8
                        endif

                        if(ierr == 0) call byte_write (buf,iz,ierr)
                                             
                    endif
                enddo
            enddo
        endif
    enddo
    call err_chk(ierr,'Error writing to newre2.re2 in gen_re2_bc$')

    return
    end subroutine gen_re2_bc
!-----------------------------------------------------------------------
    subroutine gen_rea_xyz
    use size_m
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

    parameter (lv=2**ldim,lblock=1000)
    common /scrns/ xyz(lv,ldim,lblock),wk(lv*ldim*lblock)
    common /scruz/ igr(lblock)

    integer :: e,eb,eg
    character(1) :: letapt

    integer :: isym2pre(8)   ! Symmetric-to-prenek vertex ordering
    save    isym2pre
    data    isym2pre / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

    letapt = 'a'
    numapt = 1

    nxs = nx1-1
    nys = ny1-1
    nzs = nz1-1
    nblock = lv*ldim*lblock

    if (nid == 0) &
    write(10,'(i12,i3,i12,'' NEL,NDIM,NELV'')') nelgt,ndim,nelgv

    do eb=1,nelgt,lblock
        nemax = min(eb+lblock-1,nelgt)
        call rzero(xyz,nblock)
        call izero(igr,lblock)
        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            l   = 0
            if (mid == nid .AND. if3d) then ! fill owning processor
                igr(kb) = igroup(e)
                do k=0,1
                    do j=0,1
                        do i=0,1
                            l=l+1
                            li=isym2pre(l)
                            xyz(li,1,kb) = xm1(1+i*nxs,1+j*nys,1+k*nzs,e)
                            xyz(li,2,kb) = ym1(1+i*nxs,1+j*nys,1+k*nzs,e)
                            xyz(li,3,kb) = zm1(1+i*nxs,1+j*nys,1+k*nzs,e)
                        enddo
                    enddo
                enddo
            elseif (mid == nid) then    ! 2D
                igr(kb) = igroup(e)
                do j=0,1
                    do i=0,1
                        l =l+1
                        li=isym2pre(l)
                        xyz(li,1,kb) = xm1(1+i*nxs,1+j*nys,1,e)
                        xyz(li,2,kb) = ym1(1+i*nxs,1+j*nys,1,e)
                    enddo
                enddo
            endif
        enddo
        call  gop(xyz,wk,'+  ',nblock)  ! Sum across all processors
        call igop(igr,wk,'+  ',lblock)  ! Sum across all processors

        if (nid == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                write(10,'(a15,i12,a2,i5,a1,a10,i6)') &
                '      ELEMENT  ',eg,' [',numapt,letapt,']    GROUP',igr(kb)

                if (if3d) then

                    write(10,'(4g15.7)')(xyz(ic,1,kb),ic=1,4)
                    write(10,'(4g15.7)')(xyz(ic,2,kb),ic=1,4)
                    write(10,'(4g15.7)')(xyz(ic,3,kb),ic=1,4)

                    write(10,'(4g15.7)')(xyz(ic,1,kb),ic=5,8)
                    write(10,'(4g15.7)')(xyz(ic,2,kb),ic=5,8)
                    write(10,'(4g15.7)')(xyz(ic,3,kb),ic=5,8)

                else ! 2D

                    write(10,'(4g15.7)')(xyz(ic,1,kb),ic=1,4)
                    write(10,'(4g15.7)')(xyz(ic,2,kb),ic=1,4)

                endif

            enddo
        endif
    enddo

    return
    end subroutine gen_rea_xyz
!-----------------------------------------------------------------------
    subroutine gen_rea_curve(imid)

!     This routine is complex because we must first count number of
!     nontrivial curved sides.

!     A two pass strategy is used:  first count, then write

    use size_m
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

    integer :: e,eb,eg
    character(1) :: cc

    parameter (lblock=500)
    common /scrns/ vcurve(5,12,lblock),wk(5*12*lblock)
    common /scruz/ icurve(12,lblock)

    if (imid > 0) then

    !        imid = 0  ! No midside node defs
    !        imid = 1  ! Midside defs where current curve sides don't exist
    !        imid = 2  ! All nontrivial midside node defs

        if (imid == 2) call blank(ccurve,12*lelt)

        do e=1,nelt
            call gen_rea_midside_e(e)
        enddo

    endif

    nedge = 4 + 8*(ndim-2)

    ncurvn = 0
    do e=1,nelt
        do i=1,nedge
            if (ccurve(i,e) /= ' ') ncurvn = ncurvn+1
        enddo
    enddo
    ncurvn = iglsum(ncurvn,1)

    if (nid == 0) then
        WRITE(10,*)' ***** CURVED SIDE DATA *****'
        WRITE(10,'(I10,A20,A33)') ncurvn,' Curved sides follow', &
        ' IEDGE,IEL,CURVE(I),I=1,5, CCURVE'
    endif

    do eb=1,nelgt,lblock

        nemax = min(eb+lblock-1,nelgt)
        call izero(icurve,12*lblock)
        call rzero(vcurve,60*lblock)

        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            if (mid == nid) then ! fill owning processor
                do i=1,nedge
                    icurve(i,kb) = 0
                    if (ccurve(i,e) == 'C') icurve(i,kb) = 1
                    if (ccurve(i,e) == 's') icurve(i,kb) = 2
                    if (ccurve(i,e) == 'm') icurve(i,kb) = 3
                    call copy(vcurve(1,i,kb),curve(1,i,e),5)
                enddo
            endif
        enddo
        call igop(icurve,wk,'+  ',12*lblock)  ! Sum across all processors
        call  gop(vcurve,wk,'+  ',60*lblock)  ! Sum across all processors

        if (nid == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                do i=1,nedge
                    ii = icurve(i,kb)   ! equivalenced to s4
                    if (ii /= 0) then
                        if (ii == 1) cc='C'
                        if (ii == 2) cc='s'
                        if (ii == 3) cc='m'
                        if (nelgt < 1000) then
                            write(10,'(i3,i3,5g14.6,1x,a1)') i,eg, &
                            (vcurve(k,i,kb),k=1,5),cc
                        elseif (nelgt < 1000000) then
                            write(10,'(i2,i6,5g14.6,1x,a1)') i,eg, &
                            (vcurve(k,i,kb),k=1,5),cc
                        else
                            write(10,'(i2,i12,5g14.6,1x,a1)') i,eg, &
                            (vcurve(k,i,kb),k=1,5),cc
                        endif
                    endif
                enddo
            enddo
        endif

    enddo

    return
    end subroutine gen_rea_curve
!-----------------------------------------------------------------------
    subroutine gen_rea_bc (ifld)

    use size_m
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

    integer :: e,eb,eg

    parameter (lblock=500)
    common /scrns/ vbc(5,6,lblock),wk(5*6*lblock)
    common /scruz/ ibc(6,lblock)

    character(1) :: s4(4)
    character(3) :: s3
    integer ::     i4
    equivalence(i4,s4)
    equivalence(s3,s4)

    character(1) :: chtemp
    save        chtemp
    data        chtemp /' '/   ! For mesh bcs

    nface = 2*ndim

    nlg = nelg(ifld)

    if (ifld == 1 .AND. .NOT. ifflow) then ! NO B.C.'s for this field
        if (nid == 0) write(10,*) &
        ' ***** NO FLUID   BOUNDARY CONDITIONS *****'
        return
    elseif (ifld == 1 .AND. nid == 0) then ! NO B.C.'s for this field
        write(10,*) ' *****    FLUID   BOUNDARY CONDITIONS *****'
    elseif (ifld >= 2 .AND. nid == 0) then ! NO B.C.'s for this field
        write(10,*) ' *****    THERMAL BOUNDARY CONDITIONS *****'
    endif

    do eb=1,nlg,lblock
        nemax = min(eb+lblock-1,nlg)
        call izero(ibc, 6*lblock)
        call rzero(vbc,30*lblock)
        kb = 0
        do eg=eb,nemax
            mid = gllnid(eg)
            e   = gllel (eg)
            kb  = kb+1
            if (mid == nid) then ! fill owning processor
                do i=1,nface
                    i4 = 0
                    call chcopy(s4,cbc(i,e,ifld),3)
                    ibc(i,kb) = i4
                    call copy(vbc(1,i,kb),bc(1,i,e,ifld),5)
                enddo
            endif
        enddo
        call igop(ibc,wk,'+  ', 6*lblock)  ! Sum across all processors
        call  gop(vbc,wk,'+  ',30*lblock)  ! Sum across all processors

        if (nid == 0) then
            kb = 0
            do eg=eb,nemax
                kb  = kb+1

                do i=1,nface
                    i4 = ibc(i,kb)   ! equivalenced to s4

                !                 chtemp='   '
                !                 if (ifld.eq.1 .or. (ifld.eq.2 .and. .not. ifflow))
                !    $               chtemp = cbc(i,kb,0)

                    if (nlg < 1000) then
                        write(10,'(a1,a3,2i3,5g14.6)') &
                        chtemp,s3,eg,i,(vbc(ii,i,kb),ii=1,5)
                    elseif (nlg < 100000) then
                        write(10,'(a1,a3,i5,i1,5g14.6)') &
                        chtemp,s3,eg,i,(vbc(ii,i,kb),ii=1,5)
                    elseif (nlg < 1000000) then
                        write(10,'(a1,a3,i6,5g14.6)') &
                        chtemp,s3,eg,(vbc(ii,i,kb),ii=1,5)
                    else
                        write(10,'(a1,a3,i12,5g18.11)') &
                        chtemp,s3,eg,(vbc(ii,i,kb),ii=1,5)
                    endif
                enddo
            enddo
        endif
    enddo

    return
    end subroutine gen_rea_bc
!-----------------------------------------------------------------------
    subroutine gen_rea_midside_e(e)

    use size_m
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

    common /scrns/ x3(27),y3(27),z3(27),xyz(3,3)
    character(1) :: ccrve(12)
    integer :: e,edge

    integer :: e3(3,12)
    save    e3
    data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1 &
    , 19,20,21,   21,24,27,   27,26,25,   25,22,19 &
    ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /

    real :: len

    call chcopy(ccrve,ccurve(1,e),12)

    call map2reg(x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
    call map2reg(y3,3,ym1(1,1,1,e),1)
    if (if3d) call map2reg(z3,3,zm1(1,1,1,e),1)

!     Take care of spherical curved face defn
    if (ccurve(5,e) == 's') then
        call chcopy(ccrve(1),'ssss',4) ! face 5
        call chcopy(ccrve(5),' ',1)    ! face 5
    endif
    if (ccurve(6,e) == 's') then
        call chcopy(ccrve(5),'ssss',4) ! face 6
    endif

    tol   = 1.e-4
    tol2  = tol**2
    nedge = 4 + 8*(ndim-2)

    do i=1,nedge
        if (ccrve(i) == ' ') then
            do j=1,3
                xyz(1,j)=x3(e3(j,i))
                xyz(2,j)=y3(e3(j,i))
                xyz(3,j)=z3(e3(j,i))
            enddo
            len = 0.
            h   = 0.
            do j=1,ndim
                xmid = .5*(xyz(j,1)+xyz(j,3))
                h    = h   + (xyz(j,2)-xmid)**2
                len  = len + (xyz(j,3)-xyz(j,1))**2
            enddo
            if (h > tol2*len) ccurve(i,e) = 'm'
            if (h > tol2*len) call copy(curve(1,i,e),xyz(1,2),ndim)
        endif
    enddo

    return
    end subroutine gen_rea_midside_e
!-----------------------------------------------------------------------
    subroutine g2gi_buf2v(v,buf,ndim,nel,nxyz)

    real*4 :: buf(1)
    real :: v(nel*nxyz,1)

    do i = 1,nel
        k = (i-1) * ndim*nxyz
        jj = (i-1)*nxyz
        do j = 1,ndim
            call copy4r(v(jj+1,j),buf(k+1),nxyz)
            k = k + nxyz
        enddo
    enddo

    return
    end subroutine g2gi_buf2v
!-----------------------------------------------------------------------
    subroutine g2gi_v2buf(buf,v,ndim,nel,nxyz)

    real*4 :: buf(1)
    real :: v(nel*nxyz,1)

    do i = 1,nel
        k = (i-1) * ndim*nxyz
        jj = (i-1)*nxyz
        do j = 1,ndim
            call copyx4(buf(k+1),v(jj+1,j),nxyz)
            k = k + nxyz
        enddo
    enddo

    return
    end subroutine g2gi_v2buf
!-----------------------------------------------------------------------
    subroutine buffer_in(buffer,npp,npoints,nbuf)
            
    use size_m
    use parallel

    real ::    buffer(ldim,nbuf)

    ierr = 0
    if(nid == 0) then
        write(6,*) 'reading hpts.in'
        open(50,file='hpts.in',status='old',err=100)
        read(50,*,err=100) npoints
        goto 101
        100 ierr = 1
        101 continue
    endif
    ierr=iglsum(ierr,1)
    if(ierr > 0) then
        write(6,*) 'Cannot open hpts.in in subroutine hpts()'
        call exitt
    endif
          
    call bcast(npoints,isize)
    if(npoints > lhis*np) then
        if(nid == 0) write(6,*) 'ABORT: Too many pts to read in hpts()!'
        call exitt
    endif
    if(nid == 0) write(6,*) 'found ', npoints, ' points'


    npass =  npoints/nbuf +1  !number of passes to cover all pts
    n0    =  mod(npoints,nbuf)!remainder
    if(n0 == 0) then
        npass = npass-1
        n0    = nbuf
    endif

    len = wdsize*ndim*nbuf
    if (nid > 0 .AND. nid < npass) msg_id=irecv(nid,buffer,len)
    call nekgsync
          
    npp=0
    if(nid == 0) then
        i1 = nbuf
        do ipass = 1,npass
            if(ipass == npass) i1 = n0
            do i = 1,i1
                read(50,*) (buffer(j,i),j=1,ndim)
            enddo
            if(ipass < npass)call csend(ipass,buffer,len,ipass,0)
        enddo
        close(50)
        npp = n0
        open(50,file='hpts.out',status='new')
        write(50,'(A)') &
        '# time  vx  vy  [vz]  pr  T  PS1  PS2  ...'
    elseif (nid < npass)  then !processors receiving data
        call msgwait(msg_id)
        npp=nbuf
    endif

    return
    end subroutine buffer_in
!-----------------------------------------------------------------------
    subroutine hpts_in(pts,npts,npoints)
!                        npts=local count; npoints=total count

    use size_m
    use parallel

    parameter (lt2=2*lx1*ly1*lz1*lelt)
    common /scrns/ xyz(ldim,lt2)
    common /scruz/ mid(lt2)  ! Target proc id
    real ::    pts(ldim,npts)

    if (lt2 > npts) then

        call buffer_in(xyz,npp,npoints,lt2)
        if(npoints > np*npts) then
            if(nid == 0)write(6,*)'ABORT in hpts(): npoints > NP*lhis!!'
            if(nid == 0)write(6,*)'Change SIZE: ',np,npts,npoints
            call exitt
        endif

        npmax = (npoints/npts)
        if(mod(npoints,npts) == 0) npmax=npmax+1

        if(nid > 0 .AND. npp > 0) then
            npts_b = lt2*(nid-1)               ! # pts  offset(w/o 0)
            nprc_b = npts_b/npts               ! # proc offset(w/o 0)

            istart = mod(npts_b,npts)          ! istart-->npts pts left
            ip     = nprc_b + 1                ! PID offset
            icount = istart                    ! point offset
        elseif(nid == 0) then
            npts0   = mod1(npoints,lt2)        ! Node 0 pts
            npts_b  = npoints - npts0          ! # pts before Node 0
            nprc_b  = npts_b/npts

            istart  = mod(npts_b,npts)
            ip      = nprc_b + 1
            icount  = istart
        endif

        do i =1,npp
            icount = icount + 1
            if(ip > npmax) ip = 0
            mid(i) = ip
            if (icount == npts) then
                ip     = ip+1
                icount = 0
            endif
        enddo

        call crystal_tuple_transfer &
        (cr_h,npp,lt2,mid,1,pts,0,xyz,ldim,1)

        call copy(pts,xyz,ldim*npp)
    else
        call buffer_in(pts,npp,npoints,npts)
    endif
    npts = npp


    return
    end subroutine hpts_in
!-----------------------------------------------------------------------
    subroutine hpts_out(fieldout,nflds,nfldm,npoints,nbuff)

    use size_m
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

    real :: buf(nfldm,nbuff),fieldout(nfldm,nbuff)

    len = wdsize*nfldm*nbuff


    npass = npoints/nbuff + 1
    il = mod(npoints,nbuff)
    if(il == 0) then
        il = nbuff
        npass = npass-1
    endif

    do ipass = 1,npass

        call nekgsync

        if(ipass < npass) then
            if(nid == 0) then
                call crecv(ipass,buf,len)
                do ip = 1,nbuff
                    write(50,'(1p20E15.7)') time, &
                    (buf(i,ip), i=1,nflds)
                enddo
            elseif(nid == ipass) then
                call csend(ipass,fieldout,len,0,nid)
            endif

        else  !ipass == npass

            if(nid == 0) then
                do ip = 1,il
                    write(50,'(1p20E15.7)') time, &
                    (fieldout(i,ip), i=1,nflds)
                enddo
            endif

        endif
    enddo

    return
    end subroutine hpts_out
!-----------------------------------------------------------------------
