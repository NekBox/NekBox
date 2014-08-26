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
#if 0
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
#endif
!-----------------------------------------------------------------------
