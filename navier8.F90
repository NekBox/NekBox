!-----------------------------------------------------------------------

    subroutine set_vert(glo_num,ngv,nx,nel,vertex,ifcenter)

!     Given global array, vertex, pointing to hex vertices, set up
!     a new array of global pointers for an nx^ndim set of elements.

    use size_m
    use input

    integer*8 :: glo_num(1),ngv
    integer :: vertex(1),nx
    logical :: ifcenter

    if (if3d) then
        call setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
    else
!max        call setvert2d(glo_num,ngv,nx,nel,vertex,ifcenter)
    endif

!     Check for single-element periodicity 'p' bc
    nz = 1
    if (if3d) nz = nx
    call check_p_bc(glo_num,nx,nx,nz,nel)

    if(nid == 0) write(6,*) 'call usrsetvert'
    call usrsetvert(glo_num,nel,nx,nx,nx)
    if(nid == 0) write(6,'(A,/)') ' done :: usrsetvert'

    return
    end subroutine set_vert

!-----------------------------------------------------------------------
!      subroutine test_h1_crs
!      use size_m
!      use domain
!      common /scrxxt/ x(lcr*lelv),b(lcr*lelv)
!      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
!      real x,b
!      ntot=nelv*nxyz_c
!      do i=1,12
!         call rzero(b,ntot)
!         if(mp.eq.1) then
!            b(i)=1
!         else if(mid.eq.0) then
!            if(i.gt.8) b(i-8)=1
!         else
!            if(i.le.8) b(i)=1
!         endif
!         call hsmg_coarse_solve(x,b)
!         print *, 'Column ',i,':',(x(j),j=1,ntot)
!      enddo
!      return
!      end
!-----------------------------------------------------------------------

subroutine set_up_h1_crs
  use kinds, only : DP
  use size_m
  use domain
  use geom
  use input
  use parallel
  use tstep
  implicit none

  integer :: mid, mp, nekcomm, nekgroup, nekreal
  common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

  common /ivrtx/ vertex ((2**ldim)*lelt)
  integer :: vertex

  integer :: gs_handle
  integer :: null_space,e

  character(3) :: cb
  real(DP) :: cmlt, mask
  common /scrxxt/ cmlt(lcr,lelv),mask(lcr,lelv)
  integer, allocatable :: ia(:,:,:), ja(:,:,:)
  real :: z

  integer :: key(2),aa(2)
  real(DP), allocatable :: a(:)
  integer :: iwork
  common /scrch/ iwork(2,lx1*ly1*lz1*lelv)
  common /scrns/ w(7*lx1*ly1*lz1*lelv)
  integer :: w
  real :: wr(1)
  equivalence (wr,w)
 
  real(DP), allocatable :: h1(:), h2(:), w1(:), w2(:)

  integer*8 :: ngv

  real(DP) :: t0
  integer :: nxc, ncr, ntot, nzc, nfaces, ie, iface, nz, n
  real(DP), external :: dnekclock

  allocate(ia(lcr,lcr,lelv), ja(lcr,lcr,lelv))
  allocate(a(27*lx1*ly1*lz1*lelv))
  allocate(h1(lx1*ly1*lz1*lelv),h2(lx1*ly1*lz1*lelv))
  allocate(w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelv))
  t0 = dnekclock()

!     nxc is order of coarse grid space + 1, nxc=2, linear, 3=quad,etc.
!     nxc=param(82)
!     if (nxc.gt.lxc) then
!        nxc=lxc
!        write(6,*) 'WARNING :: coarse grid space too large',nxc,lxc
!     endif
!     if (nxc.lt.2) nxc=2

  nxc     = 2
  nx_crs  = nxc

  if(nid == 0) write(6,*) 'setup h1 coarse grid, nx_crs=', nx_crs

  ncr     = nxc**ndim
  nxyz_c  = ncr

!   Set SEM_to_GLOB

  call get_vertex
  call set_vert(se_to_gcrs,ngv,nxc,nelv,vertex, .TRUE. )

!   Set mask
  z=0
  ntot=nelv*nxyz_c
  nzc=1
  if (if3d) nzc=nxc
  call rone(mask,ntot)
  call rone(cmlt,ntot)
  nfaces=2*ndim
!   ifield=1			!c? avo: set in set_overlap through 'TSTEP'?

  if (ifield == 1) then
      do ie=1,nelv
          do iface=1,nfaces
              cb=cbc(iface,ie,ifield)
              if (cb == 'O  '  .OR.  cb == 'ON '  .OR.  cb == 'MM '  .OR. &
              cb == 'mm '  .OR.  cb == 'ms '  .OR.  cb == 'MS ') &
              call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
          enddo
      enddo
  elseif (ifield == ifldmhd) then   ! no ifmhd ?avo?
      do ie=1,nelv
          do iface=1,nfaces
              cb=cbc(iface,ie,ifield)
              if (cb == 'ndd'  .OR.  cb == 'dnd'  .OR.  cb == 'ddn') &
              call facev(mask,ie,iface,z,nxc,nxc,nzc)
          enddo
      enddo
  endif

!   Set global index of dirichlet nodes to zero; xxt will ignore them

  call gs_setup(gs_handle,se_to_gcrs,ntot,nekcomm,mp)
  call gs_op   (gs_handle,mask,1,2,0)  !  "*"
  call gs_op   (gs_handle,cmlt,1,1,0)  !  "+"
  call gs_free (gs_handle)
  call set_jl_crs_mask(ntot,mask,se_to_gcrs)

  call invcol1(cmlt,ntot)

!   Setup local SEM-based Neumann operators (for now, just full...)

!    if (param(51).eq.1) then     ! old coarse grid
!       nxyz1=nx1*ny1*nz1
!       lda = 27*nxyz1*lelt
!       ldw =  7*nxyz1*lelt
!       call get_local_crs(a,lda,nxc,h1,h2,w,ldw)
!    else
!      NOTE: a(),h1,...,w2() must all be large enough
  n = nx1*ny1*nz1*nelv
  call rone (h1,n)
  call rzero(h2,n)
  call get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)
!    endif

  call set_mat_ij(ia,ja,ncr,nelv)
  null_space=0
  if (ifield == 1) then
      if (ifvcor)  null_space=1
  elseif (ifield == ifldmhd) then
      if (ifbcor)  null_space=1
  endif

  nz=ncr*ncr*nelv
  call crs_setup(xxth(ifield),nekcomm,mp, ntot,se_to_gcrs, &
  nz,ia,ja,a, null_space)
!   call crs_stats(xxth(ifield))

  t0 = dnekclock()-t0
  if (nid == 0) then
      write(6,*) 'done :: setup h1 coarse grid ',t0, ' sec'
      write(6,*) ' '
  endif

  return
end subroutine set_up_h1_crs

!-----------------------------------------------------------------------
    subroutine set_jl_crs_mask(n, mask, se_to_gcrs)
    real :: mask(1)
    integer*8 :: se_to_gcrs(1)
    do i=1,n
        if(mask(i) < 0.1) se_to_gcrs(i)=0
    enddo
    return
    end subroutine set_jl_crs_mask
!-----------------------------------------------------------------------
    subroutine set_mat_ij(ia,ja,n,ne)
    integer :: n,ne
    integer :: ia(n,n,ne), ja(n,n,ne)

    integer :: i,j,ie
    do ie=1,ne
        do j=1,n
            do i=1,n
                ia(i,j,ie)=(ie-1)*n+i-1
                ja(i,j,ie)=(ie-1)*n+j-1
            enddo
        enddo
    enddo
    return
    end subroutine set_mat_ij
!-----------------------------------------------------------------------

    subroutine ituple_sort(a,lda,n,key,nkey,ind,aa)

!     Use Heap Sort (p 231 Num. Rec., 1st Ed.)

    integer :: a(lda,n),aa(lda)
    integer :: ind(1),key(nkey)
    logical :: iftuple_ialtb

    dO 10 j=1,n
        ind(j)=j
    10 END DO

    if (n <= 1) return
    L=n/2+1
    ir=n
    100 continue
    if (l > 1) then
        l=l-1
    !           aa  = a  (l)
        call icopy(aa,a(1,l),lda)
        ii  = ind(l)
    else
    !           aa =   a(ir)
        call icopy(aa,a(1,ir),lda)
        ii = ind(ir)
    !           a(ir) =   a( 1)
        call icopy(a(1,ir),a(1,1),lda)
        ind(ir) = ind( 1)
        ir=ir-1
        if (ir == 1) then
        !              a(1) = aa
            call icopy(a(1,1),aa,lda)
            ind(1) = ii
            return
        endif
    endif
    i=l
    j=l+l
    200 continue
    if (j <= ir) then
        if (j < ir) then
            if (iftuple_ialtb(a(1,j),a(1,j+1),key,nkey)) j=j+1
        endif
        if (iftuple_ialtb(aa,a(1,j),key,nkey)) then
        !              a(i) = a(j)
            call icopy(a(1,i),a(1,j),lda)
            ind(i) = ind(j)
            i=j
            j=j+j
        else
            j=ir+1
        endif
        GOTO 200
    endif
!        a(i) = aa
    call icopy(a(1,i),aa,lda)
    ind(i) = ii
    GOTO 100
    end subroutine ituple_sort

!-----------------------------------------------------------------------

    subroutine tuple_sort(a,lda,n,key,nkey,ind,aa)

!     Use Heap Sort (p 231 Num. Rec., 1st Ed.)

    real :: a(lda,n),aa(lda)
    integer :: ind(1),key(nkey)
    logical :: iftuple_altb

    dO 10 j=1,n
        ind(j)=j
    10 END DO

    if (n <= 1) return
    L=n/2+1
    ir=n
    100 continue
    if (l > 1) then
        l=l-1
    !           aa  = a  (l)
        call copy(aa,a(1,l),lda)
        ii  = ind(l)
    else
    !           aa =   a(ir)
        call copy(aa,a(1,ir),lda)
        ii = ind(ir)
    !           a(ir) =   a( 1)
        call copy(a(1,ir),a(1,1),lda)
        ind(ir) = ind( 1)
        ir=ir-1
        if (ir == 1) then
        !              a(1) = aa
            call copy(a(1,1),aa,lda)
            ind(1) = ii
            return
        endif
    endif
    i=l
    j=l+l
    200 continue
    if (j <= ir) then
        if (j < ir) then
        !              if ( a(j).lt.a(j+1) ) j=j+1
            if (iftuple_altb(a(1,j),a(1,j+1),key,nkey)) j=j+1
        endif
    !           if (aa.lt.a(j)) then
        if (iftuple_altb(aa,a(1,j),key,nkey)) then
        !              a(i) = a(j)
            call copy(a(1,i),a(1,j),lda)
            ind(i) = ind(j)
            i=j
            j=j+j
        else
            j=ir+1
        endif
        GOTO 200
    endif
!        a(i) = aa
    call copy(a(1,i),aa,lda)
    ind(i) = ii
    GOTO 100
    end subroutine tuple_sort

!-----------------------------------------------------------------------

    logical function iftuple_ialtb(a,b,key,nkey)
    integer :: a(1),b(1)
    integer :: key(nkey)

    do i=1,nkey
        k=key(i)
        if (a(k) < b(k)) then
            iftuple_ialtb = .TRUE. 
            return
        elseif (a(k) > b(k)) then
            iftuple_ialtb = .FALSE. 
            return
        endif
    enddo
    iftuple_ialtb = .FALSE. 
    return
    end function iftuple_ialtb

!-----------------------------------------------------------------------

    logical function iftuple_altb(a,b,key,nkey)
    real :: a(1),b(1)
    integer :: key(nkey)

    do i=1,nkey
        k=key(i)
        if (a(k) < b(k)) then
            iftuple_altb = .TRUE. 
            return
        elseif (a(k) > b(k)) then
            iftuple_altb = .FALSE. 
            return
        endif
    enddo
    iftuple_altb = .FALSE. 
    return
    end function iftuple_altb

!-----------------------------------------------------------------------

    logical function iftuple_ianeb(a,b,key,nkey)
    integer :: a(1),b(1)
    integer :: key(nkey)

    do i=1,nkey
        k=key(i)
        if (a(k) /= b(k)) then
            iftuple_ianeb = .TRUE. 
            return
        endif
    enddo
    iftuple_ianeb = .FALSE. 
    return
    end function iftuple_ianeb

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    subroutine bindec(bin_in)
    integer :: bin_in,d,b,b2

    keep  = bin_in
    d  = bin_in
    b2 = 1
    b  = 0
    do l=1,12
        b  = b + b2*mod(d,10)
        d  = d/10
        b2 = b2*2
        if (d == 0) goto 1
    enddo
    1 continue
    bin_in = b
    return
    end subroutine bindec

!-----------------------------------------------------------------------
    subroutine get_local_A_tet(a,x,y,z,kt,ie)

!     Generate local tetrahedral matrix


    real :: a(4,4), g(4,4)
    real :: x(4),y(4),z(4)

    11 continue
    x23 = x(2) - x(3)
    y23 = y(2) - y(3)
    z23 = z(2) - z(3)
    x34 = x(3) - x(4)
    y34 = y(3) - y(4)
    z34 = z(3) - z(4)
    x41 = x(4) - x(1)
    y41 = y(4) - y(1)
    z41 = z(4) - z(1)
    x12 = x(1) - x(2)
    y12 = y(1) - y(2)
    z12 = z(1) - z(2)

    xy234 = x34*y23 - x23*y34
    xy341 = x34*y41 - x41*y34
    xy412 = x12*y41 - x41*y12
    xy123 = x12*y23 - x23*y12
    xz234 = x23*z34 - x34*z23
    xz341 = x41*z34 - x34*z41
    xz412 = x41*z12 - x12*z41
    xz123 = x23*z12 - x12*z23
    yz234 = y34*z23 - y23*z34
    yz341 = y34*z41 - y41*z34
    yz412 = y12*z41 - y41*z12
    yz123 = y12*z23 - y23*z12

    g(1,1) = -(x(2)*yz234 + y(2)*xz234 + z(2)*xy234)
    g(2,1) = -(x(3)*yz341 + y(3)*xz341 + z(3)*xy341)
    g(3,1) = -(x(4)*yz412 + y(4)*xz412 + z(4)*xy412)
    g(4,1) = -(x(1)*yz123 + y(1)*xz123 + z(1)*xy123)
    g(1,2) = yz234
    g(2,2) = yz341
    g(3,2) = yz412
    g(4,2) = yz123
    g(1,3) = xz234
    g(2,3) = xz341
    g(3,3) = xz412
    g(4,3) = xz123
    g(1,4) = xy234
    g(2,4) = xy341
    g(3,4) = xy412
    g(4,4) = xy123

!        vol36 = 1/(36*volume) = 1/(6*determinant)

    det = x(1)*yz234 + x(2)*yz341 + x(3)*yz412 + x(4)*yz123
    vol36 = 1.0/(6.0*det)
    if (vol36 < 0) then
        write(6,*) 'Error: tetrahedron not right-handed',ie
        write(6,1) 'x',(x(k),k=1,4)
        write(6,1) 'y',(y(k),k=1,4)
        write(6,1) 'z',(z(k),k=1,4)
        1 format(a1,1p4e15.5)

    !        call exitt                 ! Option 1

        xx = x(1)                  ! Option 2
        x(1) = x(2)                !  -- this is the option that
        x(2) = xx                  !     actually works. 11/25/07

        xx = y(1)
        y(1) = y(2)
        y(2) = xx

        xx = z(1)
        z(1) = z(2)
        z(2) = xx

        goto 11

    !        call rzero(a,16)           ! Option 3
    !        return

    !        vol36 = abs(vol36)         ! Option 4

    endif

    do j=1,4
        do i=1,4
            a(i,j)=vol36*(g(i,2)*g(j,2)+g(i,3)*g(j,3)+g(i,4)*g(j,4))
        enddo
    enddo

    return
    end subroutine get_local_A_tet

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)

!     -  Spectral interpolation from A to B via tensor products
!     -  scratch arrays: w(na*na*nb + nb*nb*na)

!     5/3/00  -- this routine replaces specmp in navier1.f, which
!                has a potential memory problem


    logical :: if3d

    real :: b(nb,nb,nb),a(na,na,na)
    real :: w(ldw)

    ltest = na*nb
    if (if3d) ltest = na*na*nb + nb*na*na
    if (ldw < ltest) then
        write(6,*) 'ERROR specmp:',ldw,ltest,if3d
        call exitt
    endif

    if (if3d) then
        nab = na*nb
        nbb = nb*nb
        call mxm(ba,nb,a,na,w,na*na)
        k=1
        l=na*na*nb + 1
        do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
        enddo
        l=na*na*nb + 1
        call mxm(w(l),nbb,ab,na,b,nb)
    else
        call mxm(ba,nb,a,na,w,na)
        call mxm(w,nb,ab,na,b,nb)
    endif
    return
    end subroutine specmpn

!-----------------------------------------------------------------------

    subroutine irank(A,IND,N)

!     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.

    integer :: A(1),IND(1)

    if (n <= 1) return
    DO 10 J=1,N
        IND(j)=j
    10 END DO

    if (n == 1) return
    L=n/2+1
    ir=n
    100 continue
    IF (l > 1) THEN
        l=l-1
        indx=ind(l)
        q=a(indx)
    ELSE
        indx=ind(ir)
        q=a(indx)
        ind(ir)=ind(1)
        ir=ir-1
        if (ir == 1) then
            ind(1)=indx
            return
        endif
    ENDIF
    i=l
    j=l+l
    200 continue
    IF (J <= IR) THEN
        IF (J < IR) THEN
            IF ( A(IND(j)) < A(IND(j+1)) ) j=j+1
        ENDIF
        IF (q < A(IND(j))) THEN
            IND(I)=IND(J)
            I=J
            J=J+J
        ELSE
            J=IR+1
        ENDIF
        GOTO 200
    ENDIF
    IND(I)=INDX
    GOTO 100
    end subroutine irank
!-----------------------------------------------------------------------

    subroutine ifacev_redef(a,ie,iface,val,nx,ny,nz)

!     Assign the value VAL to face(IFACE,IE) of array A.
!     IFACE is the input in the pre-processor ordering scheme.

    use size_m
    integer :: a(nx,ny,nz,lelt),val
    call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
    do 100 iz=kz1,kz2
        do 100 iy=ky1,ky2
            do 100 ix=kx1,kx2
                a(ix,iy,iz,ie)=val
    100 END DO
    return
    end subroutine ifacev_redef

!-----------------------------------------------------------------------

    subroutine out_se1(se2crs,nx,name)

    use size_m
    integer :: se2crs(nx,nx,1)
    character(4) :: name

    write(6,*)
    write(6,*) 'out_se',nx,name
    do ie=nelv-1,1,-2
        write(6,*)
        do j=nx,1,-1
            if(nx == 4) then
                write(6,4) name,((se2crs(i,j,k+ie),i=1,nx),k=0,1)
            elseif(nx == 3) then
                write(6,3) name,((se2crs(i,j,k+ie),i=1,nx),k=0,1)
            else
                write(6,2) name,((se2crs(i,j,k+ie),i=1,nx),k=0,1)
            endif
        enddo
    enddo

    4 format(a4,5x,2(4i5,3x))
    3 format(a4,5x,2(3i5,3x))
    2 format(a4,5x,2(2i5,3x))

    return
    end subroutine out_se1

!-----------------------------------------------------------------------

    subroutine out_se0(se2crs,nx,nel,name)

    use size_m
    integer :: se2crs(nx,nx,1)
    character(4) :: name

    write(6,*)
    write(6,*) 'out_se',nx,name,nel
    do ie=nel-3,1,-4
        write(6,*)
        do j=nx,1,-1
            if(nx == 4) then
                write(6,4) name,((se2crs(i,j,k+ie),i=1,nx),k=0,3)
            elseif(nx == 3) then
                write(6,3) name,((se2crs(i,j,k+ie),i=1,nx),k=0,3)
            else
                write(6,2) name,((se2crs(i,j,k+ie),i=1,nx),k=0,3)
            endif
        enddo
    enddo

    4 format(a4,5x,4(4i5,3x))
    3 format(a4,5x,4(3i5,3x))
    2 format(a4,5x,4(2i5,3x))

    return
    end subroutine out_se0

!-----------------------------------------------------------------------
    subroutine set_h1_basis_bilin

    use size_m
    use domain
    use wz_m

    do ix=1,nx1
        h1_basis(ix) = 0.5*(1.0-zgm1(ix,1))
        h1_basis(ix+nx1) = 0.5*(1.0+zgm1(ix,1))
    enddo
    call transpose(h1_basist,2,h1_basis,lx1)

    return
    end subroutine set_h1_basis_bilin

!-----------------------------------------------------------------------
    subroutine get_local_crs_galerkin(a,ncl,nxc,h1,h2,w1,w2)

!     This routine generates Nelv submatrices of order ncl using
!     Galerkin projection

    use size_m

    real ::    a(ncl,ncl,1),h1(1),h2(1)
    real ::    w1(nx1*ny1*nz1,nelv),w2(nx1*ny1*nz1,nelv)

    parameter (lcrd=lx1**ldim)
    common /ctmp1z/ b(lcrd,8)

    integer :: e

    do j=1,ncl
        call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
    enddo

    isd  = 1
    imsh = 1

    nxyz = nx1*ny1*nz1
    do j = 1,ncl
        do e = 1,nelv
            call copy(w1(1,e),b(1,j),nxyz)
        enddo

        call axhelm (w2,w1,h1,h2,imsh,isd)        ! A^e * bj

        do e = 1,nelv
            do i = 1,ncl
                a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * A^e * bj
            enddo
        enddo

    enddo

    return
    end subroutine get_local_crs_galerkin
!-----------------------------------------------------------------------
    subroutine gen_crs_basis(b,j) ! bi- tri-linear

    use size_m
    real :: b(nx1,ny1,nz1)

    real :: z0(lx1),z1(lx1)
    real :: zr(lx1),zs(lx1),zt(lx1)

    integer :: p,q,r

    call zwgll(zr,zs,nx1)

    do i=1,nx1
        z0(i) = .5*(1-zr(i))  ! 1-->0
        z1(i) = .5*(1+zr(i))  ! 0-->1
    enddo

    call copy(zr,z0,nx1)
    call copy(zs,z0,nx1)
    call copy(zt,z0,nx1)

    if (mod(j,2) == 0)                        call copy(zr,z1,nx1)
    if (j == 3 .OR. j == 4 .OR. j == 7 .OR. j == 8) call copy(zs,z1,nx1)
    if (j > 4)                               call copy(zt,z1,nx1)

    if (ndim == 3) then
        do r=1,nx1
            do q=1,nx1
                do p=1,nx1
                    b(p,q,r) = zr(p)*zs(q)*zt(r)
                enddo
            enddo
        enddo
    else
        do q=1,nx1
            do p=1,nx1
                b(p,q,1) = zr(p)*zs(q)
            enddo
        enddo
    endif

    return
    end subroutine gen_crs_basis
!-----------------------------------------------------------------------
    subroutine get_vertex
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
    use zper

    common /ivrtx/ vertex ((2**ldim)*lelt)
    integer :: vertex

    integer :: icalld
    save    icalld
    data    icalld  /0/

    if (icalld > 0) return
    icalld = 1

    if (ifgtp) then
        call gen_gtp_vertex    (vertex, ncrnr)
    else
        call get_vert
    endif

    return
    end subroutine get_vertex
!-----------------------------------------------------------------------
    subroutine assign_gllnid(gllnid,iunsort,nelgt,nelgv,np)

    integer :: gllnid(1),iunsort(1),nelgt,np
    integer :: e,eg


    log2p = log2(np)
    np2   = 2**log2p
    if (np2 == np .AND. nelgv == nelgt) then   ! std power of 2 case

        npstar = ivlmax(gllnid,nelgt)+1
        nnpstr = npstar/np
        do eg=1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
        enddo

        return

    elseif (np2 == np) then   ! std power of 2 case, conjugate heat xfer

    !        Assign fluid elements
        npstar = max(np,ivlmax(gllnid,nelgv)+1)
        nnpstr = npstar/np
        do eg=1,nelgv
            gllnid(eg) = gllnid(eg)/nnpstr
        enddo

    !        Assign solid elements
        nelgs  = nelgt-nelgv  ! number of solid elements
        npstar = max(np,ivlmax(gllnid(nelgv+1),nelgs)+1)
        nnpstr = npstar/np
        do eg=nelgv+1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
        enddo

        return

    elseif (nelgv /= nelgt) then
        call exitti &
        ('Conjugate heat transfer requires P=power of 2.$',np)
    endif


!  Below is the code for P a non-power of two:

!  Split the sorted gllnid array (read from .map file)
!  into np contiguous partitions.

!  To load balance the partitions in case of mod(nelgt,np)>0
!  add 1 contiguous entry out of the sorted list to NODE_i
!  where i = np-mod(nelgt,np) ... np


    nel   = nelgt/np       ! number of elements per processor
    nmod  = mod(nelgt,np)  ! bounded between 1 ... np-1
    npp   = np - nmod      ! how many paritions of size nel
     
! sort gllnid
    call isort(gllnid,iunsort,nelgt)

! setup partitions of size nel
    k   = 0
    do ip = 0,npp-1
        do e = 1,nel
            k = k + 1
            gllnid(k) = ip
        enddo
    enddo
! setup partitions of size nel+1
    if(nmod > 0) then
        do ip = npp,np-1
            do e = 1,nel+1
                k = k + 1
                gllnid(k) = ip
            enddo
        enddo
    endif

! unddo sorting to restore initial ordering by
! global element number
    call iswapt_ip(gllnid,iunsort,nelgt)

    return
    end subroutine assign_gllnid
!-----------------------------------------------------------------------
    subroutine get_vert
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
    use zper

    common /ivrtx/ vertex ((2**ldim),lelt)
    integer :: vertex

    integer :: e,eg

    integer :: icalld
    save    icalld
    data    icalld  /0/
    if (icalld > 0) return
    icalld = 1

    ncrnr = 2**ndim

    if (ifmoab) then
#ifdef MOAB
        call nekMOAB_loadConn (vertex, nelgt, ncrnr)
#endif
    else
        call get_vert_map(vertex, ncrnr, nelgt, '.map', ifgfdm)
    endif

    return
    end subroutine get_vert
!-----------------------------------------------------------------------
    subroutine get_vert_map(vertex, nlv, nel, suffix, ifgfdm)
    use size_m
    use input
    use parallel
    logical :: ifgfdm
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
    integer :: vertex(nlv,1)
    character(4) :: suffix

    parameter(mdw=2+2**ldim)
    parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
    common /scrns/ wk(mdw,ndw)   ! room for long ints, if desired
    integer :: wk,e,eg,eg0,eg1

    character(132) :: mapfle
    character(1) ::   mapfle1(132)
    equivalence  (mapfle,mapfle1)

    iok = 0
    if (nid == 0) then
        lfname = ltrunc(reafle,132) - 4
        call blank (mapfle,132)
        call chcopy(mapfle,reafle,lfname)
        call chcopy(mapfle1(lfname+1),suffix,4)
        open(unit=80,file=mapfle,status='old',err=99)
        read(80,*,err=99) neli,nnzi
        iok = 1
    endif
    99 continue
    iok = iglmax(iok,1)
    if (iok == 0) goto 999     ! Mapfile not found

    if (nid == 0) then
        neli = iglmax(neli,1)   ! communicate to all procs
    else
        neli = 0
        neli = iglmax(neli,1)   ! communicate neli to all procs
    endif

    npass = 1 + (neli/ndw)
    if (npass > np) then
        if (nid == 0) write(6,*) npass,np,neli,ndw,'Error get_vert_map'
        call exitt
    endif

    len = 4*mdw*ndw
    if (nid > 0 .AND. nid < npass) msg_id=irecv(nid,wk,len)
    call nekgsync

    if (nid == 0) then
        eg0 = 0
        do ipass=1,npass
            eg1 = min(eg0+ndw,neli)
            m   = 0
            do eg=eg0+1,eg1
                m = m+1
                read(80,*,end=998) (wk(k,m),k=2,mdw)
                if( .NOT. ifgfdm)  gllnid(eg) = wk(2,m)  !proc map,  must still be divided
                wk(1,m)    = eg
            enddo
            if (ipass < npass) call csend(ipass,wk,len,ipass,0) !send to ipass
            eg0 = eg1
        enddo
        close(80)
        ntuple = m
    elseif (nid < npass) then
        call msgwait(msg_id)
        ntuple = ndw
    else
        ntuple = 0
    endif

!     Distribute and assign partitions
    if ( .NOT. ifgfdm) then             ! gllnid is already assigned for gfdm
        lng = isize*neli
        call bcast(gllnid,lng)
        call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch

    !       if(nid.eq.0) then
    !         write(99,*) (gllnid(i),i=1,nelgt)
    !       endif
    !       call exitt
    endif

    nelt=0 !     Count number of elements on this processor
    nelv=0
    do eg=1,neli
        if (gllnid(eg) == nid) then
            if (eg <= nelgv) nelv=nelv+1
            if (eg <= nelgt) nelt=nelt+1
        endif
    enddo
    if (np <= 64) write(6,*) nid,nelv,nelt,nelgv,nelgt,' NELV'

!     NOW: crystal route vertex by processor id

    do i=1,ntuple
        eg=wk(1,i)
        wk(2,i)=gllnid(eg)        ! processor id for element eg
    enddo

    key = 2  ! processor id is in wk(2,:)
    call crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)

    if ( .NOT. ifgfdm) then            ! no sorting for gfdm?
        key = 1  ! Sort tuple list by eg := wk(1,:)
        nkey = 1
        call crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
    endif

    iflag = 0
    if (ntuple /= nelt) then
        write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
        write(6,*) 'Check that .map file and .rea file agree'
        iflag=1
    else
        nv = 2**ndim
        do e=1,nelt
            call icopy(vertex(1,e),wk(3,e),nv)
        enddo
    endif

    iflag = iglmax(iflag,1)
    if (iflag > 0) then
        do mid=0,np-1
            call nekgsync
            if (mid == nid) &
            write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FB'
            call nekgsync
        enddo
        call nekgsync
        call exitt
    endif

    return

    999 continue
    if (nid == 0) write(6,*) 'ABORT: Could not find map file ',mapfle
    call exitt

    998 continue
    if (nid == 0) write(6,*)ipass,npass,eg0,eg1,mdw,m,eg,'get v fail'
    call exitt0  ! Emergency exit

    return
    end subroutine get_vert_map
!-----------------------------------------------------------------------
    subroutine irank_vecn(ind,nn,a,m,n,key,nkey,aa)

!     Compute rank of each unique entry a(1,i)

!     Output:   ind(i)  i=1,...,n    (global) rank of entry a(*,i)
!               nn  = max(rank)
!               a(j,i) is permuted

!     Input:    a(j,i) j=1,...,m;  i=1,...,n
!               m      :   leading dim. of v  (ldv must be .ge. m)
!               key    :   sort key
!               nkey   :

!     Although not mandatory, this ranking procedure is probably
!     most effectively employed when the keys are pre-sorted. Thus,
!     the option is provided to sort vi() prior to the ranking.


    integer :: ind(n),a(m,n)
    integer :: key(nkey),aa(m)
    logical :: iftuple_ianeb,a_ne_b

    nk = min(nkey,m)
    call ituple_sort(a,m,n,key,nk,ind,aa)

!     Find unique a's
    call icopy(aa,a,m)
    nn     = 1
    ind(1) = nn

    do i=2,n
        a_ne_b = iftuple_ianeb(aa,a(1,i),key,nk)
        if (a_ne_b) then
            call icopy(aa,a(1,i),m)
            nn = nn+1
        endif
        ind(i) = nn ! set ind() to rank
    enddo

    return
    end subroutine irank_vecn
!-----------------------------------------------------------------------
    subroutine gbtuple_rank(tuple,m,n,nmax,cr_h,nid,np,ind)

!     Return a unique rank for each matched tuple set. Global.  Balanced.

!     tuple is destroyed.

!     By "balanced" we mean that none of the tuple entries is likely to
!     be much more uniquely populated than any other, so that any of
!     the tuples can serve as an initial (parallel) sort key

!     First two slots in tuple(:,i) assumed empty

    integer :: ind(nmax),tuple(m,nmax),cr_h

    parameter (mmax=40)
    integer :: key(mmax),wtuple(mmax)

    if (m > mmax) then
        write(6,*) nid,m,mmax,' gbtuple_rank fail'
        call exitt
    endif

    do i=1,n
        tuple(1,i) = mod(tuple(3,i),np) ! destination processor
        tuple(2,i) = i                  ! return location
    enddo

    ni= n
    ky=1  ! Assumes crystal_new already called
    call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

    nimx = iglmax(ni,1)
    if (ni > nmax)   write(6,*) ni,nmax,n,'cr_xfer problem, A'
    if (nimx > nmax) call exitt

    nkey = m-2
    do k=1,nkey
        key(k) = k+2
    enddo

    call irank_vecn(ind,nu,tuple,m,ni,key,nkey,wtuple)! tuple re-ordered,
! but contents same

    nu_tot   = igl_running_sum(nu) ! running sum over P processors
    nu_prior = nu_tot - nu

    do i=1,ni
        tuple(3,i) = ind(i) + nu_prior  ! global ranking
    enddo

    call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

    nk = 1  ! restore to original order, local rank: 2; global: 3
    ky = 2
    call ituple_sort(tuple,m,n,ky,nk,ind,wtuple)


    return
    end subroutine gbtuple_rank
!-----------------------------------------------------------------------
    subroutine setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)

!     setup unique ids for dssum
!     note:
!     total number of unique vertices, edges and faces has to be smaller
!     than 2**31 (integer-4 limit).
!     if nelgt < 2**31/12 we're ok for sure (independent of N)!

    use ctimer
    use size_m
    use geom
    use parallel
    use topol

    integer*8 :: glo_num(1),ngv
    integer :: vertex(0:1,0:1,0:1,1),nx
    logical :: ifcenter

    integer ::  edge(0:1,0:1,0:1,3,lelt),enum(12,lelt),fnum(6,lelt)
    common  /scrmg/ edge,enum,fnum

    parameter (nsafe=8)  ! OFTEN, nsafe=2 suffices
    integer :: etuple(4,12*lelt*nsafe),ftuple(5,6,lelt*nsafe)
    integer :: ind(12*lelt*nsafe)
    common  /scrns/ ind,etuple
    equivalence  (etuple,ftuple)

    integer :: gvf(4),facet(4),aa(3),key(3),e
    logical :: ifij
          
    integer*8 :: igv,ig0
    integer*8 :: ngvv,ngve,ngvs,ngvi,ngvm
    integer*8 :: n_on_edge,n_on_face,n_in_interior
    integer*8 :: i8glmax

    ny   = nx
    nz   = nx
    nxyz = nx*ny*nz

    key(1)=1
    key(2)=2
    key(3)=3

!     Assign hypercube ordering of vertices
!     -------------------------------------

!     Count number of unique vertices
    nlv  = 2**ndim
    ngvv = iglmax(vertex,nlv*nel)

    do e=1,nel
        do k=0,1
            do j=0,1
                do i=0,1
                !           Local to global node number (vertex)
                    il  = 1 + (nx-1)*i + nx*(nx-1)*j + nx*nx*(nx-1)*k
                    ile = il + nx*ny*nz*(e-1)
                    glo_num(ile)   = vertex(i,j,k,e)
                enddo
            enddo
        enddo
    enddo
    ngv  = ngvv

    if (nx == 2) return

!     Assign global vertex numbers to SEM nodes on each edge
!     ------------------------------------------------------

!     Assign edge labels by bounding vertices.
    do e=1,nel
        do k=0,1
            do j=0,1
                do i=0,1
                    edge(i,j,k,1,e) = vertex(i,j,k,e)  ! r-edge
                    edge(j,i,k,2,e) = vertex(i,j,k,e)  ! s-edge
                    edge(k,i,j,3,e) = vertex(i,j,k,e)  ! t-edge
                enddo
            enddo
        enddo
    enddo

!     Sort edges by bounding vertices.
    do i=0,12*nel-1
        if (edge(0,i,0,1,1) > edge(1,i,0,1,1)) then
            kswap = edge(0,i,0,1,1)
            edge(0,i,0,1,1) = edge(1,i,0,1,1)
            edge(1,i,0,1,1) = kswap
        endif
        etuple(3,i+1) = edge(0,i,0,1,1)
        etuple(4,i+1) = edge(1,i,0,1,1)
    enddo

!     Assign a number (rank) to each unique edge
    m    = 4
    n    = 12*nel
    nmax = 12*lelt*nsafe  ! nsafe for crystal router factor of safety
    call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
    do i=1,12*nel
        enum(i,1) = etuple(3,i)
    enddo
    n_unique_edges = iglmax(enum,12*nel)

    n_on_edge = nx-2
    ngve      = n_unique_edges*n_on_edge
    do e=1,nel
        iedg_loc = 0
    
    !        Edges 1-4
        do k=0,1
            do j=0,1
                igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
                i0  = nx*(nx-1)*j + nx*nx*(nx-1)*k
                i0e = i0 + nxyz*(e-1)
                if (glo_num(i0e+1) < glo_num(i0e+nx)) then
                    do i=2,nx-1                                   ! std forward case
                        glo_num(i0e+i) = igv + i-1
                    enddo
                else
                    do i=2,nx-1                                   ! backward case
                        glo_num(i0e+i) = igv + 1 + n_on_edge-(i-1)
                    enddo
                endif
                iedg_loc = iedg_loc + 1
            enddo
        enddo
    
    !        Edges 5-8
        do k=0,1
            do i=0,1
                igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
                i0  = 1+(nx-1)*i + nx*nx*(nx-1)*k
                i0e = i0 + nxyz*(e-1)
                if (glo_num(i0e) < glo_num(i0e+nx*(nx-1))) then
                    do j=2,nx-1                                   ! std forward case
                        glo_num(i0e+(j-1)*nx) = igv + j-1
                    enddo
                else
                    do j=2,nx-1                                   ! backward case
                        glo_num(i0e+(j-1)*nx) = igv + 1 + n_on_edge-(j-1)
                    enddo
                endif
                iedg_loc = iedg_loc + 1
            enddo
        enddo
    
    !        Edges 9-12
        do j=0,1
            do i=0,1
                igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
                i0  = 1 + (nx-1)*i + nx*(nx-1)*j
                i0e = i0 + nxyz*(e-1)
                if (glo_num(i0e) < glo_num(i0e+nx*nx*(nx-1))) then
                    do k=2,nx-1                                   ! std forward case
                        glo_num(i0e+(k-1)*nx*nx) = igv + k-1
                    enddo
                else
                    do k=2,nx-1                                   ! backward case
                        glo_num(i0e+(k-1)*nx*nx) = igv + 1 + n_on_edge-(k-1)
                    enddo
                endif
                iedg_loc = iedg_loc + 1
            enddo
        enddo
    enddo
    ngv   = ngv + ngve

!     Asign global node numbers on the interior of each face
!     ------------------------------------------------------

!     Assign faces by 3-tuples

!     (The following variables all take the symmetric
!     notation of IFACE as arguments:)

!     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
!                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.

!                        3+-----+4    ^ Y
!                        /  2  /|     |
!     Edge 1 extends    /     / |     |
!       from vertex   7+-----+8 +2    +----> X
!       1 to 2.        |  4  | /     /
!                      |     |/     /
!                     5+-----+6    Z
!                         3

    nfaces=ndim*2
    ncrnr =2**(ndim-1)
    do e=1,nel
        do ifac=1,nfaces
            do icrn=1,ncrnr
                i                  = icface(icrn,ifac)-1
                facet(icrn)        = vertex(i,0,0,e)
            enddo
            call isort(facet,ind,ncrnr)
            call icopy(ftuple(3,ifac,e),facet,ncrnr-1)
        enddo
    enddo

!     Assign a number (rank) to each unique face
    m    = 5
    n    = 6*nel
    nmax = 6*lelt*nsafe  ! nsafe for crystal router factor of safety
    call gbtuple_rank(ftuple,m,n,nmax,cr_h,nid,np,ind)
    do i=1,6*nel
        fnum(i,1) = ftuple(3,i,1)
    enddo
    n_unique_faces = iglmax(fnum,6*nel)

    call dsset (nx,ny,nz)
    do e=1,nel
        do iface=1,nfaces
            i0 = skpdat(1,iface)
            i1 = skpdat(2,iface)
            is = skpdat(3,iface)
            j0 = skpdat(4,iface)
            j1 = skpdat(5,iface)
            js = skpdat(6,iface)
        
        !        On each face, count from minimum global vertex number,
        !        towards smallest adjacent vertex number.  e.g., suppose
        !        the face is defined by the following global vertex numbers:
        
        
        !                    11+--------+81
        !                      |c      d|
        !                      |        |
        !                      |        |
        !                      |a      b|
        !                    15+--------+62
        
        !        We would count from c-->a, then towards d.
        
            gvf(1) = glo_num(i0+nx*(j0-1)+nxyz*(e-1))
            gvf(2) = glo_num(i1+nx*(j0-1)+nxyz*(e-1))
            gvf(3) = glo_num(i0+nx*(j1-1)+nxyz*(e-1))
            gvf(4) = glo_num(i1+nx*(j1-1)+nxyz*(e-1))
        
            call irank(gvf,ind,4)
        
        !        ind(1) tells which element of gvf() is smallest.
        
            ifij = .FALSE. 
            if (ind(1) == 1) then
                idir =  1
                jdir =  1
                if (gvf(2) < gvf(3)) ifij = .TRUE. 
            elseif (ind(1) == 2) then
                idir = -1
                jdir =  1
                if (gvf(1) < gvf(4)) ifij = .TRUE. 
            elseif (ind(1) == 3) then
                idir =  1
                jdir = -1
                if (gvf(4) < gvf(1)) ifij = .TRUE. 
            elseif (ind(1) == 4) then
                idir = -1
                jdir = -1
                if (gvf(3) < gvf(2)) ifij = .TRUE. 
            endif
        
            if (idir < 0) then
                it=i0
                i0=i1
                i1=it
                is=-is
            endif
        
            if (jdir < 0) then
                jt=j0
                j0=j1
                j1=jt
                js=-js
            endif
        
            nxx = nx*nx
            n_on_face = (nx-2)*(ny-2)
            ngvs  = n_unique_faces*n_on_face
            ig0 = ngv + n_on_face*(fnum(iface,e)-1)
            if (ifij) then
                k=0
                l=0
                do j=j0,j1,js
                    do i=i0,i1,is
                        k=k+1
                    !              this is a serious kludge to stay on the face interior
                        if (k > nx .AND. k < nxx-nx .AND. &
                        mod(k,nx) /= 1 .AND. mod(k,nx) /= 0) then
                        !                 interior
                            l = l+1
                            glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
                        endif
                    enddo
                enddo
            else
                k=0
                l=0
                do i=i0,i1,is
                    do j=j0,j1,js
                        k=k+1
                    !              this is a serious kludge to stay on the face interior
                        if (k > nx .AND. k < nxx-nx .AND. &
                        mod(k,nx) /= 1 .AND. mod(k,nx) /= 0) then
                        !                 interior
                            l = l+1
                            glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
    ngv   = ngv + ngvs

!     Finally,  number interiors (only ifcenter=.true.)
!     -------------------------------------------------

    n_in_interior = (nx-2)*(ny-2)*(nz-2)
    ngvi = n_in_interior*nelgt
    if (ifcenter) then
        do e=1,nel
            ig0 = ngv + n_in_interior*(lglel(e)-1)
            l = 0
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        l = l+1
                        glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = ig0+l
                    enddo
                enddo
            enddo
        enddo
        ngv = ngv + ngvi
    else
        do e=1,nel
            l = 0
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        l = l+1
                        glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = 0
                    enddo
                enddo
            enddo
        enddo
    endif

!     Quick check on maximum #dofs:
    m    = nxyz*nelt
    ngvm = i8glmax(glo_num,m)
    ngvv = ngvv + ngve + ngvs  ! number of unique ids w/o interior
    ngvi = ngvi + ngvv         ! total number of unique ids
    if (nid == 0) write(6,1) nx,ngvv,ngvi,ngv,ngvm
    1 format('   setvert3d:',i4,4i12)

    return
    end subroutine setvert3d
!-----------------------------------------------------------------------
    subroutine check_p_bc(glo_num,nx,ny,nz,nel)

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

    integer*8 :: glo_num(nx,ny,nz,nel)
    integer*8 :: gmn

    integer :: e,f,fo,ef,efo
    integer :: eface0(6)
    save    eface0
    data    eface0 / 4,2,1,3,5,6 /

    ifld = 2
    if (ifflow) ifld = 1

    nface=2*ndim
    do e=1,nelt
        do f=1,nface,2
            fo  = f+1
            ef  = eface0(f)
            efo = eface0(fo)
            if (cbc(ef,e,ifld) == 'p  ' .AND. cbc(efo,e,ifld) == 'p  ') then
                if (f == 1) then  ! r=-/+1
                    do k=1,nz
                        do j=1,ny
                            gmn = min(glo_num(1,j,k,e),glo_num(nx,j,k,e))
                            glo_num(1 ,j,k,e) = gmn
                            glo_num(nx,j,k,e) = gmn
                        enddo
                    enddo
                elseif (f == 3) then  ! s=-/+1
                    do k=1,nz
                        do i=1,nx
                            gmn = min(glo_num(i,1,k,e),glo_num(i,ny,k,e))
                            glo_num(i,1 ,k,e) = gmn
                            glo_num(i,ny,k,e) = gmn
                        enddo
                    enddo
                else
                    do j=1,ny
                        do i=1,nx
                            gmn = min(glo_num(i,j,1,e),glo_num(i,j,nz,e))
                            glo_num(i,j,1 ,e) = gmn
                            glo_num(i,j,nz,e) = gmn
                        enddo
                    enddo
                endif
            endif
        enddo
    enddo

    return
    end subroutine check_p_bc
!-----------------------------------------------------------------------
