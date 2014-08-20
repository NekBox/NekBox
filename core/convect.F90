!-----------------------------------------------------------------------

!    Stability limits:

!    AB3:    .7236                     w/safety (1.2):   .603

!    RK3:    1.73   (sqrt 3)           w/safety (1.2):   1.44

!    RK4:    2.828                     w/safety (1.2):   2.36

!    SEM Safety factor:  1.52 for N=3
!                     <  1.20 for N=16
!                     ~  1.16 for N=256

!-----------------------------------------------------------------------
    subroutine setup_convect(igeom)
    use size_m
    include 'TOTAL'
    logical :: ifnew

    common /cchar/ ct_vx(0:lorder+1) ! time for each slice in c_vx()

    if (igeom == 1) return
    if (param(99) < 0) return ! no dealiasing

    if (ifchar) then

        nelc = nelv
        if (ifmhd) nelc = max(nelv,nelfld(ifldmhd))
        if (ifmhd) call exitti('no characteristics for mhd yet$',istep)

        ifnew = .TRUE. 
        if (igeom > 2) ifnew = .FALSE. 

!max        call set_conv_char(ct_vx,c_vx,vx,vy,vz,nelc,time,ifnew)

    else

        if ( .NOT. ifpert) then
            if (ifcons) then
!                call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
            else
                call set_convect_new  (vxd,vyd,vzd,vx,vy,vz)
            endif
        endif

    endif

!     write(6,*) istep,' conv',ifnew,igeom,' continu? ',time
!     read(5,*) dum

    return
    end subroutine setup_convect
!-----------------------------------------------------------------------
    subroutine char_conv(p0,u,ulag,msk,c,cs,gsl)


!     Convect over last NBD steps using characteristics scheme

!     NOTE:  Here, we assume that ulag is stored by time-slice first,
!            then by field number (this is opposite to prior Nek5000)


    use size_m
    include 'TOTAL'
    real ::    p0(1),u(1),ulag(1),msk(1),c(1),cs(0:1)
    integer :: gsl

    common /scrns/ ct  (lxd*lyd*lzd*lelv*ldim)

    common /scrvh/ bmsk(lx1*ly1*lz1*lelv) &
    , u1  (lx1*ly1*lz1*lelv)

    common /scrmg/ r1  (lx1*ly1*lz1*lelv) &
    , r2  (lx1*ly1*lz1*lelv) &
    , r3  (lx1*ly1*lz1*lelv) &
    , r4  (lx1*ly1*lz1*lelv)
          

    nelc = nelv            ! number of elements in convecting field
    if (ifield == ifldmhd) nelc = nelfld(ifield)

    nc  = cs(0)            ! number of stored convecting fields

    ln  = lx1*ly1*lz1*lelt
    n   = nx1*ny1*nz1*nelfld(ifield)
    m   = nxd*nyd*nzd*nelc*ndim

    if (ifield == ifldmhd) then
        call col3(bmsk,bintm1,msk,n)
    elseif (ifield == 1) then
        call col3(bmsk,binvm1,msk,n)
    else ! if (ifield == 2) then
        call col3(bmsk,bintm1,msk,n)
    endif

    call char_conv1 &
    (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs(1),nc,ct,u1,r1,r2,r3,r4)

    return
    end subroutine char_conv
!-----------------------------------------------------------------------
    subroutine char_conv1 &
    (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs,nc,ct,u1,r1,r2,r3,r4)

    use size_m
    include 'INPUT'
    include 'TSTEP'

    real ::    p0(n),u(n),ulag(ln,1),bmsk(n),c(m,0:nc),cs(0:nc)

    real ::    ct(m),u1(n),r1(n),r2(n),r3(n),r4(n) ! work arrays

    integer :: gsl


!     Convect over last NBD steps using characteristics scheme

!              n-q                                      n-1
!     Given u(t    , X ),  q = 1,2,...,nbd, compute  phi     ( := p0 )

!        n-1       nbd   ~n-q
!     phi     :=  sum    u
!                  q=1

!          ~n-q             ~n-q
!     each u     satisfies  u    := v  such that


!     dv
!     -- + C.grad v = 0  t \in [t^n-q,t^n],   v(t^n-q,X) = u(t^n-q,X)
!     dt


    tau = time-vlsum(dtlag,nbd)         ! initialize time for u^n-k
    call int_vel (ct,tau,c,m,nc,cs,nid) ! ct(t) = sum w_k c(.,k)

    call rzero(p0,n)

    do ilag = nbd,1,-1

        um = 0
        if (ilag == 1) then
            do i=1,n
                p0(i) = p0(i)+bd(ilag+1)*u(i)
                um=max(um,u(i))
            enddo
        else
            do i=1,n
                p0(i) = p0(i)+bd(ilag+1)*ulag(i,ilag-1)
                um=max(um,ulag(i,ilag-1))
            enddo
        endif

    !        write(6,1) istep,ilag,bd(ilag),bd(ilag+1),um
    ! 1      format(i5,i4,1p3e14.5,' bdf')

        dtau = dtlag(ilag)/ntaubd
        do itau = 1,ntaubd ! ntaubd=number of RK4 substeps (typ. 1 or 2)

            tau1 = tau + dtau

            c1 = 1.
            c2 = -dtau/2.
            c3 = -dtau
            th = tau+dtau/2.

            call conv_rhs(r1,p0,ct,bmsk,gsl)         !  STAGE 1

            call int_vel (ct,th,c,m,nc,cs,nid)       !  STAGE 2
            call add3s12 (u1,p0,r1,c1,c2,n)
            call conv_rhs(r2,u1,ct,bmsk,gsl)

            call add3s12 (u1,p0,r2,c1,c2,n)          !  STAGE 3
            call conv_rhs(r3,u1,ct,bmsk,gsl)

            call int_vel (ct,tau1,c,m,nc,cs,nid)     !  STAGE 4
            call add3s12 (u1,p0,r3,c1,c3,n)
            call conv_rhs(r4,u1,ct,bmsk,gsl)

            c1 = -dtau/6.
            c2 = -dtau/3.
            do i=1,n
                p0(i) = p0(i)+c1*(r1(i)+r4(i))+c2*(r2(i)+r3(i))
            enddo
            tau = tau1
        enddo
    enddo

    return
    end subroutine char_conv1
!-----------------------------------------------------------------------
    subroutine int_vel(c_t,t0,c,n,nc,ct,nid)

!     Interpolate convecting velocity field c(t_1,...,t_nconv) to
!     time t0 and return result in c_t.

!     Ouput:   c_t = sum wt_k * ct_i(k)

!     Here, t0 is the time of interest

    real :: c_t(n),c(n,0:nc),ct(0:nc)

    parameter (lwtmax=10)
    real :: wt(0:lwtmax)

    if (nc > lwtmax) then
        write(6,*) nid,'ERROR int_vel: lwtmax too small',lwtmax,m0
        call exitt
    endif

    no = nc-1
    call fd_weights_full(t0,ct(0),no,0,wt)  ! interpolation weights

    call rzero(c_t,n)
    do j=1,n
        do i=0,no
            c_t(j) = c_t(j) + wt(i)*c(j,i)
        enddo
    enddo

    return
    end subroutine int_vel
!-----------------------------------------------------------------------
    subroutine conv_rhs (du,u,c,bmsk,gsl)

    use size_m
    include 'TOTAL'

!     apply convecting field c(1,ndim) to scalar field u(1)

    real :: du(1),u(1),c(1),bmsk(1)
    integer :: gsl

    logical :: ifconv

!     ifconv = .false.
    ifconv = .TRUE. 

    n = nx1*ny1*nz1*nelv

    if (ifconv) then

        if (ifcons) then
            if (if3d     ) call convop_cons_3d (du,u,c,nx1,nxd,nelv)
            if ( .NOT. if3d) call convop_cons_2d (du,u,c,nx1,nxd,nelv)
        else
            if (if3d     ) call convop_fst_3d  (du,u,c,nx1,nxd,nelv)
            if ( .NOT. if3d) call convop_fst_2d  (du,u,c,nx1,nxd,nelv)
        endif

        call gs_op(gsl,du,1,1,0)  !  +

    else
        call rzero   (du,n)
        return
    endif

    do i=1,n
        du(i) = bmsk(i)*du(i)  ! Binv * msk
    enddo

    return
    end subroutine conv_rhs
!-----------------------------------------------------------------------
    subroutine convop_fst_3d(du,u,c,mx,md,nel)

    use size_m

!     apply convecting field c to scalar field u

    real :: du(mx*mx*mx,nel)
    real ::  u(mx*mx*mx,nel)
    real ::  c(md*md*md,nel,3)
    parameter (ldd=lxd*lyd*lzd)
    common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)

    logical :: if3d,ifd
    integer :: e

    if3d = .TRUE. 
    ifd  = .FALSE. 
    if (md /= mx) ifd= .TRUE. 
    ifd= .TRUE. 

    nrstd = md**3
    call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')

    do e=1,nel
        call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
        if (ifd) then    ! dealiased
            do i=1,nrstd
            !              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
                ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
        else
            do i=1,nrstd
            !              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
                du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
        endif
    enddo

    return
    end subroutine convop_fst_3d
!-----------------------------------------------------------------------
    subroutine convop_fst_2d(du,u,c,mx,md,nel)

    use size_m

!     apply convecting field c to scalar field u

    real :: du(mx*mx,nel)
    real ::  u(mx*mx,nel)
    real ::  c(md*md,nel,2)
    parameter (ldd=lxd*lyd*lzd)
    common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)

    logical :: if3d,ifd
    integer :: e

    if3d = .FALSE. 
    ifd  = .FALSE. 
    if (md /= mx) ifd= .TRUE. 
    ifd= .TRUE. 

    nrstd = md**2
    call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')

    do e=1,nel
        call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
        if (ifd) then    ! dealiased
            do i=1,nrstd
            !              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
                ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
        else
            do i=1,nrstd
            !              C has the mass matrix factored in per (4.8.5), p. 227, DFM.
                du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
        endif
    enddo

    return
    end subroutine convop_fst_2d
!-----------------------------------------------------------------------
    subroutine grad_rstd(ur,us,ut,u,mx,md,if3d,ju) ! GLL->GL grad

    use size_m
    use dxyz

    real ::    ur(1),us(1),ut(1),u(1),ju(1)
    logical :: if3d

    parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
    common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg) &
    , wkd(lwkd)
    real :: jgl,jgt

    call intp_rstd(ju,u,mx,md,if3d,0) ! 0 = forward

    m0 = md-1
    call get_dgl_ptr (ip,md,md)
    if (if3d) then
        call local_grad3(ur,us,ut,ju,m0,1,dg(ip),dgt(ip))
    else
        call local_grad2(ur,us   ,ju,m0,1,dg(ip),dgt(ip))
    endif

    return
    end subroutine grad_rstd
!-----------------------------------------------------------------------
    subroutine intp_rstd(ju,u,mx,md,if3d,idir) ! GLL->GL interpolation

!     GLL interpolation from mx to md.

!     If idir ^= 0, then apply transpose operator  (md to mx)

    use size_m

    real ::    ju(1),u(1)
    logical :: if3d

    parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
    common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg) &
    , wkd(lwkd)
    real :: jgl,jgt

    parameter (ld=2*lxd)
    common /ctmp0/ w(ld**ldim,2)

    call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
    call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')

    ldw = 2*(ld**ldim)

    call get_int_ptr (i,mx,md)

    if (idir == 0) then
        call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
    else
        call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
    endif

    return
    end subroutine intp_rstd
!-----------------------------------------------------------------------
    subroutine gen_int(jgl,jgt,mp,np,w)

!     Generate interpolation from np GLL points to mp GL points

!        jgl  = interpolation matrix, mapping from velocity nodes to pressure
!        jgt  = transpose of interpolation matrix
!        w    = work array of size (np+mp)

!        np   = number of points on GLL grid
!        mp   = number of points on GL  grid


    real :: jgl(mp,np),jgt(np*mp),w(1)

    iz = 1
    id = iz + np

    call zwgll (w(iz),jgt,np)
    call zwgl  (w(id),jgt,mp)

    n  = np-1
    do i=1,mp
        call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
        do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
        enddo
    enddo

    call transpose(jgt,np,jgl,mp)

    return
    end subroutine gen_int
!-----------------------------------------------------------------------
    subroutine gen_dgl(dgl,dgt,mp,np,w)

!     Generate derivative from np GL points onto mp GL points

!        dgl  = interpolation matrix, mapping from velocity nodes to pressure
!        dgt  = transpose of interpolation matrix
!        w    = work array of size (3*np+mp)

!        np   = number of points on GLL grid
!        mp   = number of points on GL  grid



    real :: dgl(mp,np),dgt(np*mp),w(1)


    iz = 1
    id = iz + np

    call zwgl  (w(iz),dgt,np)  ! GL points
    call zwgl  (w(id),dgt,mp)  ! GL points

    ndgt = 2*np
    ldgt = mp*np
    call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')

    n  = np-1
    do i=1,mp
        call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
        do j=1,np
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
        enddo
    enddo

    call transpose(dgt,np,dgl,mp)

    return
    end subroutine gen_dgl
!-----------------------------------------------------------------------
    subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
    use size_m            ! need nid
    character(5) ::  avar5,lvar5
    character(10) :: sub_name10

    if (n > m) then
        write(6,1) nid,n,m,avar5,lvar5,sub_name10
        1 format(i8,' ERROR: :',2i12,2(1x,a5),1x,a10)
        call exitti('lim_chk problem. $',n)
    endif

    return
    end subroutine lim_chk
!-----------------------------------------------------------------------
    subroutine get_int_ptr (ip,mx,md) ! GLL-->GL pointer

!     Get pointer to jgl() for interpolation pair (mx,md)

    use size_m

    parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
    common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg) &
    , wkd(lwkd)
    real :: jgl,jgt

    parameter (ld=2*lxd)
    common /igrad/ pd    (0:ld*ld) &
    , pdg   (0:ld*ld) &
    , pjgl  (0:ld*ld)
    integer :: pd , pdg , pjgl

    ij = md + ld*(mx-1)
    ip = pjgl(ij)

    if (ip == 0) then
    
        nstore   = pjgl(0)
        pjgl(ij) = nstore+1
        nstore   = nstore + md*mx
        pjgl(0)  = nstore
        ip       = pjgl(ij)
    
        nwrkd = mx + md
        call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
        call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
    
        call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
    endif

    return
    end subroutine get_int_ptr
!-----------------------------------------------------------------------
    subroutine get_dgl_ptr (ip,mx,md)

!     Get pointer to GL-GL interpolation dgl() for pair (mx,md)

    use size_m

    parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
    common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg) &
    , wkd(lwkd)
    real :: jgl,jgt

    parameter (ld=2*lxd)
    common /jgrad/ pd    (0:ld*ld) &
    , pdg   (0:ld*ld) &
    , pjgl  (0:ld*ld)
    integer :: pd , pdg , pjgl

    ij = md + ld*(mx-1)
    ip = pdg (ij)

    if (ip == 0) then

        nstore   = pdg (0)
        pdg (ij) = nstore+1
        nstore   = nstore + md*mx
        pdg (0)  = nstore
        ip       = pdg (ij)
    
        nwrkd = mx + md
        call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
        call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
    
        call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)
    endif

    return
    end subroutine get_dgl_ptr
!-----------------------------------------------------------------------
#if 0
    subroutine set_conv_char(ct,c,ux,uy,uz,nelc,tau,ifnew)
    use size_m
    include 'TSTEP'

    real :: ct(0:1)               ! time stamps for saved field (0=#flds)
    real :: c(1)                  ! saved vel. fields, dealiased etc.
    real :: ux(1),uy(1),uz(1)     ! input vel. field
    integer :: nelc               ! number of elements in conv. field
    logical :: ifnew              ! =true if shifting stack of fields

    numr      = lxd*lyd*lzd*lelv*ldim*(lorder+1)
    denr      = nxd*nyd*nzd*nelv*ndim
    nconv_max = numr/denr
    if (nconv_max < nbdinp+1) &
    call exitti( &
    'ABORT: not enough memory for characteristics scheme!$', &
    nconv_max)

    nc = ct(0)

    m  = nxd*nyd*nzd*nelc*ndim

!     write(6,*) nelc,ifnew,' set conv_char',istep,nc,nconv_max
    call set_ct_cvx &
    (ct,c,m,ux,uy,uz,tau,nc,nconv_max,nelc,ifnew)

    nc = min (nc,nbdinp)
    ct(0) = nc  ! store current count

    return
    end subroutine set_conv_char
#endif
!-----------------------------------------------------------------------
    subroutine set_ct_cvx(ct,c,m,u,v,w,tau,nc,mc,nelc,ifnew)
    use size_m
    include 'INPUT'  ! ifcons

    real :: ct(0:1),c(m,1)
    real :: u(1),v(1),w(1)
    logical :: ifnew

    if (ifnew) then

    !        Shift existing convecting fields
    !        Note:  "1" entry is most recent

        nc = nc+1
        nc = min(nc,mc)
        ct(0) = nc

        do i=nc,2,-1
            call copy(c(1,i),c(1,i-1),m)
            ct(i) = ct(i-1)
        enddo
    endif

!     Save time and map the current velocity to rst coordinates.

    ix = 1
    iy = ix + nxd*nyd*nzd*nelc
    iz = iy + nxd*nyd*nzd*nelc

    if (ifcons) then
!        call set_convect_cons(c(ix,1),c(iy,1),c(iz,1),u,v,w)
    else
        call set_convect_new (c(ix,1),c(iy,1),c(iz,1),u,v,w)
    endif

    ct(1) = tau

    return
    end subroutine set_ct_cvx
!-----------------------------------------------------------------------
    subroutine grad_rst(ur,us,ut,u,md,if3d) ! Gauss-->Gauss grad

    use size_m
    use dxyz

    real ::    ur(1),us(1),ut(1),u(1)
    logical :: if3d

    parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
    common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg) &
    , wkd(lwkd)
    real :: jgl,jgt

    m0 = md-1
    call get_dgl_ptr (ip,md,md)
    if (if3d) then
        call local_grad3(ur,us,ut,u,m0,1,dg(ip),dgt(ip))
    else
        call local_grad2(ur,us   ,u,m0,1,dg(ip),dgt(ip))
    endif

    return
    end subroutine grad_rst
!-----------------------------------------------------------------------
    subroutine convect_new(bdu,u,ifuf,cx,cy,cz,ifcf)

!     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians

    use size_m
    include 'TOTAL'

    real :: bdu(1),u(1),cx(1),cy(1),cz(1)
    logical :: ifuf,ifcf            ! u and/or c already on fine mesh?

    parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
    common /scrcv/ fx(ltd),fy(ltd),fz(ltd) &
    , ur(ltd),us(ltd),ut(ltd) &
    , tr(ltd,3),uf(ltd)

    integer :: e

    call set_dealias_rx

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    nxyzu = nxyz1
    if (ifuf) nxyzu = nxyzd

    nxyzc = nxyz1
    if (ifcf) nxyzc = nxyzd

    iu = 1    ! pointer to scalar field u
    ic = 1    ! pointer to vector field C
    ib = 1    ! pointer to scalar field Bdu


    do e=1,nelv

        if (ifcf) then

            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)

        else  ! map coarse velocity to fine mesh (C-->F)

            call intp_rstd(fx,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
            call intp_rstd(fy,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
            if (if3d) call intp_rstd(fz,cz(ic),nx1,nxd,if3d,0) ! 0 --> forward

            if (if3d) then  ! Convert convector F to r-s-t coordinates

                do i=1,nxyzd
                    tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
                    tr(i,2)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
                    tr(i,3)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
                enddo

            else

                do i=1,nxyzd
                    tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
                    tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
                enddo

            endif

        endif

        if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),nxd,if3d)
        else
            call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward
            call grad_rst(ur,us,ut,uf,nxd,if3d)
        endif

        if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
            enddo
        else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
                uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
            enddo
        endif
        call intp_rstd(bdu(ib),uf,nx1,nxd,if3d,1) ! Project back to coarse

        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1

    enddo

    return
    end subroutine convect_new
!-----------------------------------------------------------------------
#if 0
    subroutine convect_cons(bdu,u,ifuf,cx,cy,cz,ifcf)

!     Compute dealiased form:  J^T Bf *div. JC Ju w/ correct Jacobians

!     conservative form


    use size_m
    include 'TOTAL'

    real :: bdu(1),u(1),cx(1),cy(1),cz(1)

    logical :: ifuf,ifcf            ! u and/or c already on fine mesh?

    parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
    common /scrcv/ uf(ltd),cf(ltd),cu(ltd) &
    , cr(ltd),cs(ltd),ct(ltd)


    integer :: e

    call set_dealias_rx

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    nxyzu = nxyz1
    if (ifuf) nxyzu = nxyzd

    nxyzc = nxyz1
    if (ifcf) nxyzc = nxyzd

    iu = 1    ! pointer to scalar field u
    ic = 1    ! pointer to vector field C
    ib = 1    ! pointer to scalar field Bdu

    do e=1,nelv

        call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward

        call rzero(cu,nxyzd)
        do i=1,ndim

            if (ifcf) then  ! C is already on fine mesh

                call exitt  ! exit for now

            else  ! map coarse velocity to fine mesh (C-->F)

                if (i == 1) call intp_rstd(cf,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
                if (i == 2) call intp_rstd(cf,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
                if (i == 3) call intp_rstd(cf,cz(ic),nx1,nxd,if3d,0) ! 0 --> forward

                call col2(cf,uf,nxyzd)   !  collocate C and u on fine mesh

                call grad_rst(cr,cs,ct,cf,nxd,if3d)  ! d/dr (C_i*u)

                if (if3d) then

                    do j=1,nxyzd
                        cu(j)=cu(j) &
                        +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+3,e)+ct(j)*rx(j,i+6,e)
                    enddo

                else  ! 2D

                    do j=1,nxyzd
                        cu(j)=cu(j) &
                        +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+2,e)
                    enddo

                endif
            endif
        enddo

        call intp_rstd(bdu(ib),cu,nx1,nxd,if3d,1) ! Project back to coarse

        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1

    enddo

    return
    end subroutine convect_cons
!-----------------------------------------------------------------------
    subroutine set_convect_cons(cx,cy,cz,ux,uy,uz)

!     Put vx,vy,vz on fine mesh (for conservation form)


    use size_m
    include 'TOTAL'

    parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

    real :: cx(ltd,1),cy(ltd,1),cz(ltd,1)
    real :: ux(lxy,1),uy(lxy,1),uz(lxy,1)

    integer :: e

    call set_dealias_rx

    do e=1,nelv    ! Map coarse velocity to fine mesh (C-->F)

        call intp_rstd(cx(1,e),ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
        call intp_rstd(cy(1,e),uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
        if (if3d) call intp_rstd(cz(1,e),uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward

    enddo

    return
    end subroutine set_convect_cons
#endif
!-----------------------------------------------------------------------
    subroutine set_convect_new(cr,cs,ct,ux,uy,uz)

!     Put vxd,vyd,vzd into rst form on fine mesh

!     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).

    use size_m
    include 'TOTAL'

    parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

    real :: cr(ltd,1),cs(ltd,1),ct(ltd,1)
    real :: ux(lxy,1),uy(lxy,1),uz(lxy,1)

    common /scrcv/ fx(ltd),fy(ltd),fz(ltd) &
    , ur(ltd),us(ltd),ut(ltd) &
    , tr(ltd,3),uf(ltd)

    integer :: e

    call set_dealias_rx

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    ic = 1    ! pointer to vector field C

    do e=1,nelv

    !        Map coarse velocity to fine mesh (C-->F)

        call intp_rstd(fx,ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
        call intp_rstd(fy,uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
        if (if3d) call intp_rstd(fz,uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward

    !        Convert convector F to r-s-t coordinates

        if (if3d) then

            do i=1,nxyzd
                cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
                cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
                ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
            enddo

        else

            do i=1,nxyzd
                cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
                cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
            enddo

        endif
    enddo

    return
    end subroutine set_convect_new
!-----------------------------------------------------------------------
    subroutine set_char_mask(mask,u,v,w) ! mask for hyperbolic system

    use size_m
    use geom
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'
    integer :: msk(0:1)
    character      cb*3
    parameter (lxyz1=lx1*ly1*lz1)
    common /ctmp1/ work(lxyz1,lelt)

    real :: mask(lxyz1,1),u(lxyz1,1),v(lxyz1,1),w(lxyz1,1)

    integer :: e,f

    nfaces= 2*ndim
    ntot1 = nx1*ny1*nz1*nelv
    call rzero (work,ntot1)
    call rone  (mask,NTOT1)

    ifldv = 1
    do 100 e=1,nelv
        do 100 f=1,nfaces
            cb=cbc(f,e,ifldv)
            if (cb(1:1) == 'v' .OR. cb(1:1) == 'V') then

                call faccl3 (work(1,e),u(1,e),unx(1,1,f,e),f)
                call faddcl3(work(1,e),v(1,e),uny(1,1,f,e),f)
                if (if3d) &
                call faddcl3(work(1,e),w(1,e),unz(1,1,f,e),f)

                call fcaver (vaver,work,e,f)

                if (vaver < 0) call facev (mask,e,f,0.0,nx1,ny1,nz1)
            endif
            if (cb(1:2) == 'ws' .OR. cb(1:2) == 'WS') &
            call facev (mask,e,f,0.0,nx1,ny1,nz1)
    100 END DO

    return
    end subroutine set_char_mask
!-----------------------------------------------------------------------
    subroutine advchar

!     Compute convective contribution using
!     operator-integrator-factor method (characteristics).

    use ctimer
    use size_m
    include 'MASS'
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'
    include 'PARALLEL'


    common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()

    common /scruz/ phx  (lx1*ly1*lz1*lelt) &
    ,             phy  (lx1*ly1*lz1*lelt) &
    ,             phz  (lx1*ly1*lz1*lelt) &
    ,             hmsk (lx1*ly1*lz1*lelt)

#ifndef NOTIMER
    if (icalld == 0) tadvc=0.0
    icalld=icalld+1
    nadvc=icalld
    etime1=dnekclock()
#endif

    dti = 1./dt
    n   = nx1*ny1*nz1*nelv

    call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system

    call char_conv(phx,vx,vxlag,hmsk,c_vx,ct_vx,gsh_fld(1))
    call char_conv(phy,vy,vylag,hmsk,c_vx,ct_vx,gsh_fld(1))

    if (if3d) then

        call char_conv(phz,vz,vzlag,hmsk,c_vx,ct_vx,gsh_fld(1))

        do i=1,n
            h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
            bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
            bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
            bfz(i,1,1,1) = bfz(i,1,1,1)+phz(i)*h2i
        enddo

    else

        do i=1,n
            h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
            bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
            bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
        enddo

    endif

#ifndef NOTIMER
    tadvc=tadvc+(dnekclock()-etime1)
#endif
    return
    end subroutine advchar
!-----------------------------------------------------------------------
    subroutine convch

!     Compute convective contribution using
!     operator-integrator-factor method (characteristics).

    use ctimer
    use size_m
    include 'MASS'
    include 'INPUT'
    include 'SOLN'
    include 'TSTEP'
    include 'PARALLEL'

    common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()

    common /scruz/ phi  (lx1*ly1*lz1*lelt) &
    ,             hmsk (lx1*ly1*lz1*lelt)

    if (icalld == 0) tadvc=0.0
    icalld=icalld+1
    nadvc=icalld
    etime1=dnekclock()

    n   = nx1*ny1*nz1*nelv
    dti = 1./dt

    if (ifield == 2) then  ! set convecting velocity and mask
    !        call setup_convect(1)
        call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system
    endif


    call char_conv(phi,t(1,1,1,1,ifield-1),tlag(1,1,1,1,1,ifield-1) &
    ,hmsk,c_vx,ct_vx,gsh_fld(1))

!     pmax = glamax(phi,n)
!     qmax = glamax(vtrans(1,1,1,1,2),n)
!     write(6,*) istep,dti,pmax,' pmax'

    do i=1,n
        bq(i,1,1,1,ifield-1) = bq(i,1,1,1,ifield-1) &
        + phi(i)*bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)*dti
    enddo

    tadvc=tadvc+(dnekclock()-etime1)

    return
    end subroutine convch
!-----------------------------------------------------------------------
    subroutine convop_cons_3d(du,u,c,mx,md,nel) ! Conservation form

!     Apply convecting field c to scalar field u, conservation form d/dxj cj phi

!     Assumes that current convecting field is on dealias mesh, in c()

    use size_m
    use dealias
    use geom

    real :: du(mx*mx*mx,nel)
    real ::  u(mx*mx*mx,nel)
    real ::  c(md*md*md,nel,3)
    parameter (ldd=lxd*lyd*lzd)
    common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
    real :: ju

    logical :: if3d,ifd
    integer :: e

    if3d = .TRUE. 
    ifd  = .FALSE. 
    if (md /= mx) ifd= .TRUE. 
    ifd= .TRUE. 

    nrstd = md**3
    call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')

    do e=1,nel

        call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss points!
        call rzero     (ud,nrstd)

        do j=1,ndim
            do i=1,nrstd
                tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss) mesh

            j0 = j+0
            j3 = j+3
            j6 = j+6
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
                ud(i)=ud(i) &
                +rx(i,j0,e)*ur(i)+rx(i,j3,e)*us(i)+rx(i,j6,e)*ut(i)
            enddo
        enddo

        call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose

    enddo

    return
    end subroutine convop_cons_3d
!-----------------------------------------------------------------------
    subroutine convop_cons_2d(du,u,c,mx,md,nel) ! Conservation form

!     Apply convecting field c to scalar field u, conservation form d/dxj cj phi

!     Assumes that current convecting field is on dealias mesh, in c()

    use size_m
    use geom
    include 'TSTEP'


    real :: du(mx*mx,nel)
    real ::  u(mx*mx,nel)
    real ::  c(md*md,nel,2)
    parameter (ldd=lxd*lyd*lzd)
    common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
    real :: ju

    logical :: if3d,ifd
    integer :: e

    if3d = .FALSE. 
    ifd  = .FALSE. 
    if (md /= mx) ifd= .TRUE. 

    nrstd = md**2
    call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')

    if (nid == 0 .AND. istep < 3) write(6,*) 'convp_cons',istep

    do e=1,nel

        call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss points!
        call rzero     (ud,nrstd)

    !        call outmat(c(1,e,1),md,md,'fine u',e)
    !        call outmat(c(1,e,2),md,md,'fine v',e)
    !        call outmat(ju      ,md,md,'fine T',e)

        do j=1,ndim
            do i=1,nrstd
                tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss) mesh

            j0 = j+0
            j2 = j+2
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mesh
                ud(i)=ud(i)+rx(i,j0,e)*ur(i)+rx(i,j2,e)*us(i)
            enddo
        enddo

        call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose

    enddo

!     call exitti('convop_cons_2d$',istep)

    return
    end subroutine convop_cons_2d
!-----------------------------------------------------------------------
