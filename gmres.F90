!-----------------------------------------------------------------------
    subroutine uzawa_gmres(res,h1,h2,h2inv,intype,iter)

!     Solve the pressure equation by right-preconditioned
!     GMRES iteration.
!     intype =  0  (steady)
!     intype =  1  (explicit)
!     intype = -1  (implicit)

    use size_m
    use gmres
    include 'TOTAL'
    common  /ctolpr/ divex
    common  /cprint/ ifprint
    logical ::          ifprint
    real ::             res  (lx2*ly2*lz2*lelv)
    real ::             h1   (lx1,ly1,lz1,lelv)
    real ::             h2   (lx1,ly1,lz1,lelv)
    real ::             h2inv(lx1,ly1,lz1,lelv)

    common /scrmg/    wp (lx2,ly2,lz2,lelv)

    common /ctmp0/   wk1(lgmres),wk2(lgmres)
    common /cgmres1/ y(lgmres)

    real :: alpha, l, temp
    integer :: j,m

    logical :: iflag
    save    iflag
    data    iflag / .FALSE. /
    real ::    norm_fac
    save    norm_fac

    real*8 :: etime1,dnekclock

    if( .NOT. iflag) then
        iflag= .TRUE. 
        call uzawa_gmres_split0(ml,mu,bm2,bm2inv,nx2*ny2*nz2*nelv)
        norm_fac = 1./sqrt(volvm2)
    endif

    etime1 = dnekclock()
    etime_p = 0.
    divex = 0.
    iter  = 0
    m = lgmres

    call chktcg2(tolps,res,iconv)
    if (param(21) > 0 .AND. tolps > abs(param(21))) &
    tolps = abs(param(21))
!     if (param(21).lt.0) tolps = abs(param(21))
    if (istep == 0) tolps = 1.e-4
    tolpss = tolps

    ntot2  = nx2*ny2*nz2*nelv

    iconv = 0
    call rzero(x,ntot2)

    do while(iconv == 0 .AND. iter < 100)

        if(iter == 0) then
        !      -1
            call col3(r,ml,res,ntot2)             ! r = L  res
        !           call copy(r,res,ntot2)
        else
        ! pdate residual
            call copy(r,res,ntot2)                ! r = res
            call cdabdtp(w,x,h1,h2,h2inv,intype)  ! w = A x
            call add2s2(r,w,-1.,ntot2)            ! r = r - w
        !      -1
            call col2(r,ml,ntot2)                 ! r = L   r
        endif
    !            ______
        gamma(1) = sqrt(glsc2(r,r,ntot2))        ! gamma  = \/ (r,r)
    !      1
        if(iter == 0) then
            div0 = gamma(1)*norm_fac
            if (param(21) < 0) tolpss=abs(param(21))*div0
        endif

    ! heck for lucky convergence
        rnorm = 0.
        if(gamma(1) == 0.) goto 9000
        temp = 1./gamma(1)
        call cmult2(v(1,1),r,temp,ntot2)         ! v  = r / gamma
    !  1            1
        do j=1,m
            iter = iter+1
        !       -1
            call col3(w,mu,v(1,j),ntot2)          ! w  = U   v
        !           j
                        
            etime2 = dnekclock()
            if(param(43) == 1) then
                call uzprec(z(1,j),w,h1,h2,intype,wp)
            else                                  !       -1
                call hsmg_solve(z(1,j),w)          ! z  = M   w
            !              call copy(z(1,j),w,ntot2)          ! z  = M   w
            endif
            etime_p = etime_p + dnekclock()-etime2
                 
            call cdabdtp(w,z(1,j),                 & ! w = A z
            h1,h2,h2inv,intype)      !        j
                 
        !      -1
            call col2(w,ml,ntot2)                 ! w = L   w

        !           !modified Gram-Schmidt
        !           do i=1,j
        !              h(i,j)=glsc2(w,v(1,i),ntot2)       ! h    = (w,v )
        !                                                 !  i,j       i
        !              call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
        !           enddo                                 !          i,j  i


        !           2-PASS GS, 1st pass:

            do i=1,j
                h(i,j)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
            enddo                                 !  i,j       i

            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
                call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
            enddo                                 !          i,j  i


        !           2-PASS GS, 2nd pass:
        
        !           do i=1,j
        !              wk1(i)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
        !           enddo                                 !  i,j       i
        !                                                 !
        !           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
        
        !           do i=1,j
        !              call add2s2(w,v(1,i),-wk1(i),ntot2)! w = w - h    v
        !              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
        !           enddo


        ! pply Givens rotations to new column
            do i=1,j-1
                temp = h(i,j)
                h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
                h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
        !            ______
            alpha = sqrt(glsc2(w,w,ntot2))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha == 0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

        !            call outmat(h,m,j,' h    ',j)
                        
            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint .AND. nid == 0) &
            write (6,66) iter,tolpss,rnorm,div0,ratio,istep
            66 format(i5,1p4e12.5,i8,' Divergence')

#ifndef TST_WSCAL
            if (rnorm < tolpss) goto 900  !converged
#else
            if (iter > param(151)-1) goto 900
#endif
            if (j == m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,ntot2)   ! v    = w / alpha
        !  j+1
        enddo
        900 iconv = 1
        1000 continue
    ! ack substitution
    !     -1
    !c = H   gamma
        do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
                temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
        enddo
    ! um up Arnoldi vectors
        do i=1,j
            call add2s2(x,z(1,i),c(i),ntot2)     ! x = x + c  z
        !          i  i
        enddo
    !        if(iconv.eq.1) call dbg_write(x,nx2,ny2,nz2,nelv,'esol',3)
    enddo
    9000 continue

    divex = rnorm
!     iter = iter - 1

!     DIAGNOSTICS
!      call copy   (w,x,ntot2)
    call ortho  (w) ! Orthogonalize wrt null space, if present
!      call copy(r,res,ntot2) !r = res
!      call cdabdtp(r,w,h1,h2,h2inv,intype)  ! r = A w
!      do i=1,ntot2
!         r(i) = res(i) - r(i)               ! r = res - r
!      enddo
!      call uzawa_gmres_temp(r,bm2inv,ntot2)
!                                               !            ______
!      gamma(1) = sqrt(glsc2(r,r,ntot2)/volvm2) ! gamma  = \/ (r,r)
!                                               !      1
!      print *, 'GMRES end resid:',gamma(1)
!     END DIAGNOSTICS
    call copy(res,x,ntot2)

    call ortho (res)  ! Orthogonalize wrt null space, if present

    etime1 = dnekclock()-etime1
    if (nid == 0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p, &
    etime1
!     call flush_hack
    9999 format(i10,' U-PRES gmres:',I7,1p5e12.4)

    return
    end subroutine uzawa_gmres

!-----------------------------------------------------------------------

    subroutine uzawa_gmres_split0(l,u,b,binv,n)
    integer :: n
    real :: l(n),u(n),b(n),binv(n)
    integer :: i
    do i=1,n
        l(i)=sqrt(binv(i))
        u(i)=sqrt(b(i))
        if(abs(u(i)*l(i)-1.0) > 1e-13) print *, i, u(i)*l(i)
    enddo
    return
    end subroutine uzawa_gmres_split0

!-----------------------------------------------------------------------
    subroutine uzawa_gmres_split(l,u,b,binv,n)
    integer :: n
    real :: l(n),u(n),b(n),binv(n)
    integer :: i
    do i=1,n
    !        l(i)=sqrt(binv(i))
    !        u(i)=sqrt(b(i))

    !        u(i)=sqrt(b(i))
    !        l(i)=1./u(i)

    !        l(i)=sqrt(binv(i))
        l(i)=1.
        u(i)=1./l(i)


    !        if(abs(u(i)*l(i)-1.0).gt.1e-13)write(6,*) i,u(i)*l(i),' gmr_sp'
    enddo
    return
    end subroutine uzawa_gmres_split

!-----------------------------------------------------------------------
    subroutine uzawa_gmres_temp(a,b,n)
    integer :: n
    real :: a(n),b(n)
    integer :: i
    do i=1,n
        a(i)=sqrt(b(i))*a(i)
    enddo
    return
    end subroutine uzawa_gmres_temp
!-----------------------------------------------------------------------
    subroutine ax(w,x,h1,h2,n)
    use size_m
    include 'TOTAL'


!     w = A*x for pressure iteration


    integer :: n
    real :: w(n),x(n),h1(n),h2(n)

    imsh = 1
    isd  = 1
    call axhelm (w,x,h1,h2,imsh,isd)
    call dssum  (w,nx1,ny1,nz1)
    call col2   (w,pmask,n)

    return
    end subroutine ax
!-----------------------------------------------------------------------
    subroutine hmh_gmres(res,h1,h2,wt,iter)

!     Solve the Helmholtz equation by right-preconditioned
!     GMRES iteration.

         
    use size_m
    use fdmh1
    use gmres
    include 'TOTAL'
    common  /ctolpr/ divex
    common  /cprint/ ifprint
    logical ::          ifprint
    real ::             res  (lx1*ly1*lz1*lelv)
    real ::             h1   (lx1,ly1,lz1,lelv)
    real ::             h2   (lx1,ly1,lz1,lelv)
    real ::             wt   (lx1,ly1,lz1,lelv)

    common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

    common /cgmres1/ y(lgmres)
    common /ctmp0/   wk1(lgmres),wk2(lgmres)
    real :: alpha, l, temp
    integer :: outer

    logical :: iflag,if_hyb
    save    iflag,if_hyb
!     data    iflag,if_hyb  /.false. , .true. /
    data    iflag,if_hyb  / .FALSE. , .FALSE. /
    real ::    norm_fac
    save    norm_fac

    real*8 :: etime1,dnekclock


    n = nx1*ny1*nz1*nelv

    etime1 = dnekclock()
    etime_p = 0.
    divex = 0.
    iter  = 0
    m     = lgmres

    if( .NOT. iflag) then
        iflag= .TRUE. 
        call uzawa_gmres_split(ml,mu,bm1,binvm1,nx1*ny1*nz1*nelv)
        norm_fac = 1./sqrt(volvm1)
    endif

    if (param(100) /= 2) call set_fdm_prec_h1b(d,h1,h2,nelv)

    call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
    if (param(21) > 0 .AND. tolps > abs(param(21))) &
    tolps = abs(param(21))
    if (istep == 0) tolps = 1.e-4
    tolpss = tolps

    iconv = 0
    call rzero(x,n)

    outer = 0
    do while (iconv == 0 .AND. iter < 500)
        outer = outer+1

        if(iter == 0) then               !      -1
            call col3(r,ml,res,n)         ! r = L  res
        !           call copy(r,res,n)
        else
        ! pdate residual
            call copy  (r,res,n)                  ! r = res
            call ax    (w,x,h1,h2,n)              ! w = A x
            call add2s2(r,w,-1.,n)                ! r = r - w
        !      -1
            call col2(r,ml,n)                     ! r = L   r
        endif
    !            ______
        gamma(1) = sqrt(glsc3(r,r,wt,n))         ! gamma  = \/ (r,r)
    !      1
        if(iter == 0) then
            div0 = gamma(1)*norm_fac
            if (param(21) < 0) tolpss=abs(param(21))*div0
        endif

    ! heck for lucky convergence
        rnorm = 0.
        if(gamma(1) == 0.) goto 9000
        temp = 1./gamma(1)
        call cmult2(v(1,1),r,temp,n)             ! v  = r / gamma
    !  1            1
        do j=1,m
            iter = iter+1
        !       -1
            call col3(w,mu,v(1,j),n)              ! w  = U   v
        !           j

        ! . . . . . Overlapping Schwarz + coarse-grid . . . . . . .

            etime2 = dnekclock()

        !           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
            if (ifmgrid) then
                call h1mg_solve(z(1,j),w,if_hyb)   ! z  = M   w
            else                                  !  j
                write(*,*) "Oops: ifmgrid"
#if 0
                kfldfdm = ndim+1
                if (param(100) == 2) then
                    call h1_overlap_2 (z(1,j),w,pmask)
                else
                    call fdm_h1 &
                    (z(1,j),w,d,pmask,vmult,nelv,ktype(1,1,kfldfdm),wk)
                endif
                call crs_solve_h1 (wk,w)           ! z  = M   w
                call add2         (z(1,j),wk,n)    !  j
#endif
            endif


            call ortho        (z(1,j)) ! Orthogonalize wrt null space, if present
            etime_p = etime_p + dnekclock()-etime2
        ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

                 
            call ax  (w,z(1,j),h1,h2,n)           ! w = A z
        !        j
                 
        !      -1
            call col2(w,ml,n)                     ! w = L   w

        !           !modified Gram-Schmidt

        !           do i=1,j
        !              h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
        !                                                 !  i,j       i

        !              call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
        !           enddo                                 !          i,j  i

        !           2-PASS GS, 1st pass:

            do i=1,j
                h(i,j)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
            enddo                                 !  i,j       i

            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
                call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !          i,j  i


        !           2-PASS GS, 2nd pass:
        
        !           do i=1,j
        !              wk1(i)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
        !           enddo                                 !  i,j       i
        !                                                 !
        !           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
        
        !           do i=1,j
        !              call add2s2(w,v(1,i),-wk1(i),n)    ! w = w - h    v
        !              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
        !           enddo

        ! pply Givens rotations to new column
            do i=1,j-1
                temp = h(i,j)
                h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
                h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
        !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha == 0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint .AND. nid == 0) &
            write (6,66) iter,tolpss,rnorm,div0,ratio,istep
            66 format(i5,1p4e12.5,i8,' Divergence')

#ifndef TST_WSCAL
            if (rnorm < tolpss) goto 900  !converged
#else
            if (iter > param(151)-1) goto 900
#endif
            if (j == m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
        !  j+1
        enddo
        900 iconv = 1
        1000 continue
    ! ack substitution
    !     -1
    !c = H   gamma
        do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
                temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
        enddo
    ! um up Arnoldi vectors
        do i=1,j
            call add2s2(x,z(1,i),c(i),n)     ! x = x + c  z
        enddo                               !          i  i
    !        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
    enddo
    9000 continue

    divex = rnorm
    call copy(res,x,n)

    call ortho   (res) ! Orthogonalize wrt null space, if present

    etime1 = dnekclock()-etime1
    if (nid == 0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p, &
    etime1,if_hyb
!     call flush_hack
    9999 format(i9,' PRES gmres:',i5,1p5e12.4,1x,l4)

    if (outer <= 2) if_hyb = .FALSE. 

    return
    end subroutine hmh_gmres
!-----------------------------------------------------------------------
    subroutine set_overlap2

!     Sets up the gather scatter and the SEM operators

    use size_m
    include 'TOTAL'
     
    common /c_is1/ glo_num(lxs*lys*lzs*lelv)
    integer*8 :: glo_num
    common /ivrtx/ vertex ((2**ldim)*lelt)
    common /handle/ gsh_dd
    integer :: vertex,gsh_dd

    mz = ndim-2
    nx = nx1+2
    ny = ny1+2
    nz = nz1+2*mz
    call setupds_no_crn(gsh_dd,nx,ny,nz,nelv,nelgv,vertex,glo_num)
    call swap_lengths ! Set up Matrices for FDM
    call gen_fast_g

    return
    end subroutine set_overlap2
!-----------------------------------------------------------------------
    subroutine h1_overlap_2(u,v,mask)

!     Local overlapping Schwarz solves with overlap of 2

    use size_m
    include 'TOTAL'


    common /cwork1/ v1(lxs,lys,lzs,lelt)
    common /handle/ gsh_dd
    integer :: gsh_dd

    real :: u(lx1,lx1,lz1,1),v(1),mask(1)
    integer :: e

    n = nx1*ny1*nz1*nelfld(ifield)

    call dd_swap_vals(v1,v,gsh_dd)

    iz1 = 0
    if (if3d) iz1=1
    do ie=1,nelfld(ifield)
        do iz=1,nz1
            do iy=1,ny1
                do ix=1,nx1
                    u(ix,iy,iz,ie) = v1(ix+1,iy+1,iz+iz1,ie)
                enddo
            enddo
        enddo
    enddo

    call dssum (u,nx1,ny1,nz1)
    call col2  (u,mask,n)


    return
    end subroutine h1_overlap_2
!-----------------------------------------------------------------------
    subroutine dd_swap_vals(v1,v0,gsh_dd)

    use size_m
    include 'TOTAL'

    common /work1/ w1(lxs,lys,lzs)            & ! work arrarys for locals
    ,w2(lxs,lys,lzs)
    integer :: gsh_dd

    real :: v1(lxs,lys,lzs,lelt)
    real :: v0(lx1,ly1,lz1,lelt)

    integer :: e

    mz = ndim-2

    nx = nx1+2
    ny = ny1+2
    nz = nz1+2*mz

    n  = nx1*ny1*nz1*nelv
    m  = (nx1+2)*(ny1+2)*(nz1+2*mz)*nelv

    do e=1,nelv
        call rzero_g        (v1,e,nx,ny,nz)
        call fill_interior_g(v1,v0,e,nx1,nz1,mz,nelv)    ! v0      --> v1(int)
        call dface_ext_g    (v1,2,e,nx,nz)               ! v1(int) --> v1(face)
    enddo

!                 ~ ~T
!     This is the Q Q  part

    call gs_op(gsh_dd,v1,1,1,0)  ! 1 ==> +         ! swap v1 & add vals

    do e =1,nelv
        call dface_add1si_g  (v1,-1.,2,e,nx,nz)
        call fastdm1_g       (v1(1,1,1,e),e,w1,w2)
        call s_face_to_int2_g(v1,-1.,2,e,nx,nz)
    enddo

!     Exchange/add elemental solutions

    call gs_op              (gsh_dd,v1,1,1,0)
    do e =1,nelv
        call s_face_to_int2_g(v1,1.,2,e,nx,nz)
    enddo

    return
    end subroutine dd_swap_vals
!-----------------------------------------------------------------------
    subroutine gen_fast_g

!     Generate fast diagonalization matrices for each element

    use size_m
    use input
    include 'PARALLEL'
    include 'SOLN'
    include 'WZ'

    parameter (lxss=lxs*lxs)
    common /fastg/  sr(lxss,2,lelv),ss(lxss,2,lelv),st(lxss,2,lelv) &
    ,  df(lxs*lys*lzs,lelv)

    common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4) &
    , llr(lelt),lls(lelt),llt(lelt) &
    , lmr(lelt),lms(lelt),lmt(lelt) &
    , lrr(lelt),lrs(lelt),lrt(lelt)
    real :: lr ,ls ,lt
    real :: llr,lls,llt
    real :: lmr,lms,lmt
    real :: lrr,lrs,lrt

    integer :: lbr,rbr,lbs,rbs,lbt,rbt
         

    call load_semhat_weighted   !   Fills the SEMHAT arrays
         

    ierr = 0
    do ie=1,nelv

        call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,3,ierr)
    
    !        Set up matrices for each element.
    
        call set_up_fast_1D_sem_g( sr(1,1,ie),lr,nr ,lbr,rbr &
        ,llr(ie),lmr(ie),lrr(ie),ie)
        call set_up_fast_1D_sem_g( ss(1,1,ie),ls,ns ,lbs,rbs &
        ,lls(ie),lms(ie),lrs(ie),ie)
        if (if3d) then
            call set_up_fast_1D_sem_g( st(1,1,ie),lt,nt ,lbt,rbt &
            ,llt(ie),lmt(ie),lrt(ie),ie)
        endif
    
    !        Set up diagonal inverse
    
        if (if3d) then
            eps = 1.e-5 * (vlmax(lr(2),nr-2) &
            +  vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            l   = 1
            do k=1,nt
                do j=1,ns
                    do i=1,nr
                        diag = lr(i) + ls(j) + lt(k)
                        if (diag > eps) then
                            df(l,ie) = 1.0/diag
                        else
                        !                 write(6,3) ie,'Reset Eig in gen fast:',i,j,k,l
                        !    $                         ,eps,diag,lr(i),ls(j),lt(k)
                        !   3             format(i6,1x,a21,4i5,1p5e12.4)
                            df(l,ie) = 0.0
                        endif
                        l = l+1
                    enddo
                enddo
            enddo
        else
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            l   = 1
            do j=1,ns
                do i=1,nr
                    diag = lr(i) + ls(j)
                                   
                    if (diag > eps) then
                        df(l,ie) = 1.0/diag
                                          
                    else
                    !                 write(6,2) ie,'Reset Eig in gen fast:',i,j,l
                    !    $                         ,eps,diag,lr(i),ls(j)
                    !   2             format(i6,1x,a21,3i5,1p4e12.4)
                        df(l,ie) = 0.0
                    endif
                    l = l+1
                enddo
            enddo
        endif
    
    !        Next element ....
    
    enddo
     
    return
    end subroutine gen_fast_g
!-----------------------------------------------------------------------
    subroutine set_up_fast_1D_sem_g(s,lam,n,lbc,rbc,ll,lm,lr,ie)
    use size_m
    include 'SEMHAT'

    parameter (lr3=2*lxs*lxs)
    common /fast1dsem/ g(lr3),w(lr3)

    real :: g,w
    real :: s(1),lam(1),ll,lm,lr
    integer :: lbc,rbc

    integer :: bb0,bb1,eb0,eb1,n,n1
    logical :: l,r

    n=nx1-1
! cs on E are from normal vel component
    if(lbc == 2 .OR. lbc == 3) then !wall,sym  - Neumann
        eb0=0
        bb0=0
    else !outflow,element                      - Dirichlet
        eb0=1
        bb0=1
    endif
    if(rbc == 2 .OR. rbc == 3) then !wall,sym  - Neumann
        eb1=n
        bb1=n
    else !outflow,element                      - Dirichlet
        eb1=n-1
        bb1=n-1
    endif
    l = (lbc == 0)
    r = (rbc == 0)

!     calculate A tilde operator
    call set_up_fast_1D_sem_op_a(s,eb0,eb1,l,r,ll,lm,lr,ah)
!     calculate B tilde operator
    call set_up_fast_1D_sem_op_b(g,bb0,bb1,l,r,ll,lm,lr,bh)

    n=n+3
!     call outmat   (s,n,n,'  A   ',ie)
!     call outmat   (g,n,n,'  B   ',ie)
    call generalev(s,g,lam,n,w)
    if( .NOT. l) call row_zero(s,n,n,1)
    if( .NOT. r) call row_zero(s,n,n,n)
    call transpose(s(n*n+1),n,s,n) ! compute the transpose of s
!     call outmat   (s,n,n,'  S   ',ie)
!     call outmat   (s(n*n+1),n,n,'  St  ',1)
!     call exitt


    return
    end subroutine set_up_fast_1D_sem_g
!-----------------------------------------------------------------------
    subroutine set_up_fast_1D_sem_op_a(g,b0,b1,l &
    ,r,ll,lm,lr,ah)
!            -1 T
!     G = J B  J

!     gives the inexact restriction of this matrix to
!     an element plus two node on either side

!     g - the output matrix
!     b0, b1 - the range for Bhat indices for the element
!              (enforces boundary conditions)
!     l, r - whether there is a left or right neighbor
!     ll,lm,lr - lengths of left, middle, and right elements
!     ah - hat matrix for A or B

!     result is inexact because:
!        neighbor's boundary condition at far end unknown
!        length of neighbor's neighbor unknown
!        (these contribs should be small for large N and
!         elements of nearly equal size)

    use size_m
    real :: g(0:lx1+1,0:lx1+1)
    real :: ah(0:lx1-1,0:lx1-1)
    real :: ll,lm,lr
    integer :: b0,b1
    logical :: l,r

    real :: bl,bm,br
    integer :: n

    n =nx1-1


!     compute the weight of A hat

    if (l) bl = 2. /ll
    bm = 2. /lm
    if (r) br = 2. /lr

    call rzero(g,(n+3)*(n+3))
    do j=1,n+1
        do i=1,n+1
            g(i,j) = ah(i-1,j-1)*bm
        enddo
    enddo

    if (l) then
        g(0,0) = g(0,0) + bl*ah(n-1,n-1)
        g(0,1) = g(0,1) + bl*ah(n-1,n  )
        g(1,0) = g(1,0) + bl*ah(n  ,n-1)
        g(1,1) = g(1,1) + bl*ah(n  ,n  )
    elseif (b0 == 0) then          !Neumann BC
        g(0,0) = 1.
    else                           !Dirichlet BC
        g(0,0)     = 1.
        g(1,1)     = 1.
        do i=2,n+1
            g(i,1)  = 0.
            g(1,i)  = 0.
        enddo
    endif

    if (r) then
        g(n+1,n+1) = g(n+1,n+1) + br*ah(0,0)
        g(n+1,n+2) = g(n+1,n+2) + br*ah(0,1)
        g(n+2,n+1) = g(n+2,n+1) + br*ah(0,1)
        g(n+2,n+2) = g(n+2,n+2) + br*ah(1,1)
    elseif (b1 == n) then          !Neumann BC
        g(n+2,n+2)=1.
    else                           !Dirichlet BC
        g(n+2,n+2)   = 1.
        g(n+1,n+1)   = 1.
        do i=1,n
            g(i,n+1)  = 0.
            g(n+1,i)  = 0.
        enddo
    endif

    return
    end subroutine set_up_fast_1D_sem_op_a
!-----------------------------------------------------------------------
    subroutine set_up_fast_1D_sem_op_b(g,b0,b1,l &
    ,r,ll,lm,lr,bh)
!            -1 T
!     G = J B  J

!     gives the inexact restriction of this matrix to
!     an element plus two node on either side

!     g - the output matrix
!     b0, b1 - the range for Bhat indices for the element
!              (enforces boundary conditions)
!     l, r - whether there is a left or right neighbor
!     ll,lm,lr - lengths of left, middle, and right elements
!     bh - hat matrix for  B

!     result is inexact because:
!        neighbor's boundary condition at far end unknown
!        length of neighbor's neighbor unknown
!        (these contribs should be small for large N and
!         elements of nearly equal size)

    use size_m
    real :: g(0:lx1+1,0:lx1+1)
    real :: bh(0:lx1-1)
    real :: ll,lm,lr
    integer :: b0,b1
    logical :: l,r

    real :: bl,bm,br
    integer :: n

    n =nx1-1

!     compute the weight of B hat

    if (l) bl = ll / 2.
    bm = lm / 2.
    if (r) br = lr / 2.

    call rzero(g,(n+3)*(n+3))
    do i=1,n+1
        g(i,i) = bh(i-1)*bm
    enddo

    if (l) then
        g(0,0) = g(0,0) + bl*bh(n-1)
        g(1,1) = g(1,1) + bl*bh(n  )
    elseif (b0 == 0) then          !Neumann BC
        g(0,0) = 1.
    else                           !Dirichlet BC
        g(0,0)     = 1.
        g(1,1)     = 1.
        do i=2,n+1
            g(i,1)  = 0.
            g(1,i)  = 0.
        enddo
    endif

    if (r) then
        g(n+1,n+1) = g(n+1,n+1) + br*bh(0)
        g(n+2,n+2) = g(n+2,n+2) + br*bh(1)
    elseif (b1 == n) then          !Neumann BC
        g(n+2,n+2)=1.
    else                           !Dirichlet BC
        g(n+2,n+2)   = 1.
        g(n+1,n+1)   = 1.
        do i=1,n
            g(i,n+1)  = 0.
            g(n+1,i)  = 0.
        enddo
    endif

    return
    end subroutine set_up_fast_1D_sem_op_b
!-----------------------------------------------------------------------
    subroutine fill_interior_g(v1,v,e,nx,nz,iz1,nel)

    real :: v1(nx+2,nx+2,nz+2*iz1,nel) ! iz1=ndim-2
    real :: v (nx  ,nx  ,nz      ,nel)
    integer :: e

    ny = nx

    do iz=1,nz
        do iy=1,ny
            do ix=1,nx
                v1(ix+1,iy+1,iz+iz1,e) = v(ix,iy,iz,e)
            enddo
        enddo
    enddo

    return
    end subroutine fill_interior_g
!-----------------------------------------------------------------------
    subroutine dface_ext_g(x,t,e,nx,nz)
!     Extend interior to face of element

    use size_m
    use input
    real :: x(nx,nx,nz,1)
    integer :: e,t

    ny = nx

    if (if3d) then
        do iz=2,nz-1
            do ix=2,nx-1
                x(ix,1 ,iz,e) = x(ix, 1+t,iz,e)
                x(ix,ny,iz,e) = x(ix,ny-t,iz,e)
            enddo
        enddo
          
        do iz=2,nz-1
            do iy=2,ny-1
                x(1 ,iy,iz,e) = x( 1+t,iy,iz,e)
                x(nx,iy,iz,e) = x(nx-t,iy,iz,e)
            enddo
        enddo
          
        do iy=2,ny-1
            do ix=2,nx-1
                x(ix,iy,1 ,e) = x(ix,iy, 1+t,e)
                x(ix,iy,nz,e) = x(ix,iy,nz-t,e)
            enddo
        enddo
    else
        do ix=2,nx-1
            x(ix,1 ,1,e) = x(ix, 1+t,1,e)
            x(ix,ny,1,e) = x(ix,ny-t,1,e)
        enddo
        do iy=2,ny-1
            x(1 ,iy,1,e) = x( 1+t,iy,1,e)
            x(nx,iy,1,e) = x(nx-t,iy,1,e)
        enddo
    endif
          
    return
    end subroutine dface_ext_g
!-----------------------------------------------------------------------
    subroutine dface_add1si_g(x,c,t,e,nx,nz)
!     Scale interior and add to face of element

    use size_m
    use input
    real :: x(nx,nx,nz,1)

    integer :: e,t

    ny = nx

    if (if3d) then

        do iz=2,nz-1
            do ix=2,nx-1
                x(ix,1 ,iz,e) = x(ix,1 ,iz,e) + c*x(ix, 1+t,iz,e)
                x(ix,ny,iz,e) = x(ix,ny,iz,e) + c*x(ix,ny-t,iz,e)
            enddo
        enddo

        do iz=2,nz-1
            do iy=2,ny-1
                x(1 ,iy,iz,e) = x(1 ,iy,iz,e) + c*x( 1+t,iy,iz,e)
                x(nx,iy,iz,e) = x(nx,iy,iz,e) + c*x(nx-t,iy,iz,e)
            enddo
        enddo

        do iy=2,ny-1
            do ix=2,nx-1
                x(ix,iy,1 ,e) = x(ix,iy,1 ,e) + c*x(ix,iy, 1+t,e)
                x(ix,iy,nz,e) = x(ix,iy,nz,e) + c*x(ix,iy,nz-t,e)
            enddo
        enddo

    else  ! 2D

        do ix=2,nx-1
            x(ix,1 ,1,e) = x(ix,1 ,1,e) + c*x(ix, 1+t,1,e)
            x(ix,ny,1,e) = x(ix,ny,1,e) + c*x(ix,ny-t,1,e)
        enddo
        do iy=2,ny-1
            x(1  ,iy,1,e) = x(1 ,iy,1,e) + c*x( 1+t,iy,1,e)
            x(nx ,iy,1,e) = x(nx,iy,1,e) + c*x(nx-t,iy,1,e)
        enddo

    endif
          
    return
    end subroutine dface_add1si_g
!-----------------------------------------------------------------------
    subroutine fastdm1_g(R,ie,w1,w2)

!     Fast diagonalization solver for FEM on mesh 1

    use size_m

    parameter (lxss=lxs*lxs)
    common /fastg/  sr(lxss,2,lelv),ss(lxss,2,lelv),st(lxss,2,lelv) &
    ,  df(lxs*lys*lzs,lelv)

    parameter (lxyz = lxs*lys*lzs)

    real :: r(1),w1(1),w2(1)

    nx = nx1+2

!      T
!     S  r
    call tensr3 (w1,nx,r ,nx,sr(1,2,ie),ss(1,1,ie),st(1,1,ie),w2)


!      -1 T
!     D  S  r

    call col2   (w1,df(1,ie),lxyz)


!        -1 T
!     S D  S  r

    call tensr3 (r ,nx,w1,nx,sr(1,1,ie),ss(1,2,ie),st(1,2,ie),w2)

    return
    end subroutine fastdm1_g
!-----------------------------------------------------------------------
    subroutine s_face_to_int2_g(x,c,t,e,nx,nz)

!     Scale face and add to interior of element

    use size_m
    use input
    real :: x(nx,nx,nz,1)
    integer :: t,e

    ny=nx


    if (if3d) then

        do iz=2,nz-1
            do ix=2,nx-1
                x(ix, 1+t,iz,e) = c*x(ix,1 ,iz,e) + x(ix, 1+t,iz,e)
                x(ix,ny-t,iz,e) = c*x(ix,ny,iz,e) + x(ix,ny-t,iz,e)
            enddo
        enddo

        do iz=2,nz-1
            do iy=2,ny-1
                x( 1+t,iy,iz,e) = c*x(1 ,iy,iz,e) + x( 1+t,iy,iz,e)
                x(nx-t,iy,iz,e) = c*x(nx,iy,iz,e) + x(nx-t,iy,iz,e)
            enddo
        enddo

        do iy=2,ny-1
            do ix=2,nx-1
                x(ix,iy, 1+t,e) = c*x(ix,iy,1 ,e) + x(ix,iy, 1+t,e)
                x(ix,iy,nz-t,e) = c*x(ix,iy,nz,e) + x(ix,iy,nz-t,e)
            enddo
        enddo

    else
    !     2D
        do ix=2,nx-1
            x(ix, 1+t,1,e) = c*x(ix,1 ,1,e) + x(ix, 1+t,1,e)
            x(ix,ny-t,1,e) = c*x(ix,ny,1,e) + x(ix,ny-t,1,e)
        enddo
        do iy=2,ny-1
            x( 1+t,iy,1,e) = c*x(1 ,iy,1,e) + x( 1+t,iy,1,e)
            x(nx-t,iy,1,e) = c*x(nx,iy,1,e) + x(nx-t,iy,1,e)
        enddo
    endif
          
    return
    end subroutine s_face_to_int2_g
!-----------------------------------------------------------------------
    subroutine outfldr_g(x,txt10,nx,nz,ichk)
    use size_m
    INCLUDE 'TSTEP'
    real :: x(nx,nx,nz,lelt)
    character(10) :: txt10

    integer :: idum,e
    save idum
    data idum /3/

    if (idum < 0) return

    ny = nx


    mtot = nx*ny*nz*nelv
    if (nx > 8 .OR. nelv > 16) return
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
        write(6,116) txt10,k,ie,xmin,xmax,istep,time,nx
        do l=12,0,-4
            write(6,117)
            do j=ny,1,-1
                if (nx == 2) write(6,102) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 3) write(6,103) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 4) write(6,104) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 5) write(6,105) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 6) write(6,106) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 7) write(6,107) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 8) write(6,118) ((x(i,j,k,e+l),i=1,nx),e=1,4)
            enddo
        enddo
    enddo

    102 FORMAT(4(2f9.5,2x))
    103 FORMAT(4(3f9.5,2x))
    104 FORMAT(4(4f7.3,2x))
    105 FORMAT(5f9.5,10x,5f9.5)
    106 FORMAT(6f7.1,5x,6f7.1,5x,6f7.1,5x,6f7.1)
    107 FORMAT(7f8.4,5x,7f8.4)
    108 FORMAT(8f8.4,4x,8f8.4)
    118 FORMAT(8f12.9)

    116 FORMAT(  /,5X,'     ^              ',/, &
    &     5X,'   Y |              ',/, &
    &     5X,'     |              ',A10,/, &
    &     5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/, &
    &     5X,'       X            ','Step  =',I9,f15.5,i2)
    117 FORMAT(' ')

    if (ichk == 1 .AND. idum > 0) call checkit(idum)
    return
    end subroutine outfldr_g
!-----------------------------------------------------------------------
    subroutine outfldi_g(x,txt10,nx,nz,ichk)
    use size_m
    INCLUDE 'TSTEP'
    integer :: x(nx,nx,nz,lelt)
    character(10) :: txt10

    integer :: idum,e,xmin,xmax
    save idum
    data idum /3/

    if (idum < 0) return

    ny = nx


    mtot = nx*ny*nz*nelv
    if (nx > 8 .OR. nelv > 16) return
    xmin = iglmin(x,mtot)
    xmax = iglmax(x,mtot)

    nell = nelt
    rnel = nell
    snel = sqrt(rnel)+.1
    ne   = snel
    ne1  = nell-ne+1
    k = 1
    do ie=1,1
        ne = 0
        write(6,116) txt10,k,ie,xmin,xmax,istep,time,nx
        do l=12,0,-4
            write(6,117)
            do j=ny,1,-1
                if (nx == 2) write(6,102) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 3) write(6,103) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 4) write(6,104) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 5) write(6,105) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 6) write(6,106) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 7) write(6,107) ((x(i,j,k,e+l),i=1,nx),e=1,4)
                if (nx == 8) write(6,118) ((x(i,j,k,e+l),i=1,nx),e=1,4)
            enddo
        enddo
    enddo

    102 FORMAT(4(2i9,2x))
    103 FORMAT(4(3i9,2x))
    104 FORMAT(4(4i7,2x))
    105 FORMAT(5i9,10x,5i9)
    106 FORMAT(6i7,5x,6i7,5x,6i7,5x,6i7)
    107 FORMAT(7i8,5x,7i8)
    108 FORMAT(8i8,4x,8i8)
    118 FORMAT(8i12)

    116 FORMAT(  /,5X,'     ^              ',/, &
    &     5X,'   Y |              ',/, &
    &     5X,'     |              ',A10,/, &
    &     5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2i12,/, &
    &     5X,'       X            ','Step  =',I9,f15.5,i2)
    117 FORMAT(' ')

    if (ichk == 1 .AND. idum > 0) call checkit(idum)
    return
    end subroutine outfldi_g
!-----------------------------------------------------------------------
    subroutine setupds_no_crn(gs_h,nx,ny,nz,nel,melg,vertex,glo_num)
    use size_m
    use input
    include 'PARALLEL'
    include 'NONCON'
    integer ::   gs_h,vertex(1),e
    integer*8 :: ngv,glo_num(nx,ny,nz,nel)
          

    common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

!     set up the global numbering
    call set_vert(glo_num,ngv,nx,nel,vertex, .FALSE. )

!     zero out corners
    mz1 = max(1,nz-1)
    do e=1,nel
        do k=1,nz,mz1
            do j=1,ny,ny-1
                do i=1,nx,nx-1
                    glo_num(i,j,k,e) = 0
                enddo
            enddo
        enddo
    enddo


    ntot = nx*ny*nz*nel

    t0 = dnekclock()
    call gs_setup(gs_h,glo_num,ntot,nekcomm,mp) ! initialize gs code
    t1 = dnekclock()

    et = t1-t0

!     call gs_chkr(glo_num)

    if (nid == 0) then
        write(6,1) et,nx,nel,ntot,ngv,gs_h
        1 format('   gs_init time',1pe11.4,' seconds ',i3,4i10)
    endif

    return
    end subroutine setupds_no_crn
!-----------------------------------------------------------------------
    subroutine rzero_g(a,e,nx,ny,nz)

    real :: a(nx,ny,nz,1)
    integer :: e

    do  i=1,nx
        do  j=1,ny
            do  k=1,nz
                a(i,j,k,e) = 0.0
            enddo
        enddo
    enddo

    return
    end subroutine rzero_g
!-----------------------------------------------------------------------

