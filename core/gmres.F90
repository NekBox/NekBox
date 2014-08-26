!-----------------------------------------------------------------------
    subroutine uzawa_gmres(res,h1,h2,h2inv,intype,iter)

!     Solve the pressure equation by right-preconditioned
!     GMRES iteration.
!     intype =  0  (steady)
!     intype =  1  (explicit)
!     intype = -1  (implicit)

    use size_m
    use gmres
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
#if 0
    subroutine set_up_fast_1D_sem_g(s,lam,n,lbc,rbc,ll,lm,lr,ie)
    use size_m
    use semhat

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
#endif
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

