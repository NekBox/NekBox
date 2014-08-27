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
!> \brief Solve the Helmholtz equation by right-preconditioned
!! GMRES iteration.
subroutine hmh_gmres(res,h1,h2,wt,iter)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, lgmres
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use gmres, only : c, s, h, gamma, ml, mu, x, r, w, v, z
  use input, only : param, ifmgrid
  use mass, only : bm1, binvm1, volvm1
  use soln, only : pmask, vmult
  use tstep, only : tolps, istep
  implicit none

  real(DP) ::             res  (lx1*ly1*lz1*lelv)
  real(DP) ::             h1   (lx1,ly1,lz1,lelv)
  real(DP) ::             h2   (lx1,ly1,lz1,lelv)
  real(DP) ::             wt   (lx1,ly1,lz1,lelv)
  integer :: iter

  real(DP) :: divex
  logical ::          ifprint
  common  /cprint/ ifprint
  real(DP) :: d, wk
  common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

  real(DP) :: wk1, wk2
  common /ctmp0/   wk1(lgmres),wk2(lgmres)

  real(DP) :: alpha, l, temp
  integer :: outer

  logical, save :: iflag = .false., if_hyb = .false.
!   data    iflag,if_hyb  /.false. , .true. /
  real, save ::  norm_fac

  real(DP) :: etime1, etime2, etime_p, dnekclock
  real(DP) :: rnorm, tolpss, div0, ratio
  real(DP), external :: glsc3, vlsc3

  integer :: m, n
  integer :: i, j, k, iconv

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
!   call flush_hack
  9999 format(i9,' PRES gmres:',i5,1p5e12.4,1x,l4)

  if (outer <= 2) if_hyb = .FALSE. 

  return
end subroutine hmh_gmres
!-----------------------------------------------------------------------

