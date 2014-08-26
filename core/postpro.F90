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
