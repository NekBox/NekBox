!> \file gmres.F90 \copybrief hmh_gmres()


!-----------------------------------------------------------------------
!> \brief w = A*x for pressure iteration
subroutine ax(w,x,h1,h2,n)
  use kinds, only : DP
  use soln, only : pmask
  implicit none

  integer :: n
  real(DP) :: w(n),x(n),h1(n),h2(n)
  integer :: imsh, isd

  imsh = 1
  isd  = 1
  call axhelm (w,x,h1,h2,imsh,isd)
  call dssum  (w)
  w = w * reshape(pmask, (/ n /))

  return
end subroutine ax

!-----------------------------------------------------------------------
!> \brief Solve the Helmholtz equation by right-preconditioned
!! GMRES iteration.
subroutine hmh_gmres(res,h1,h2,wt,iter)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, lelv, lgmres
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use ctimer, only : tgmres, ngmres, gmres_flop, gmres_mop, dnekclock
  use input, only : param, ifmgrid, ifprint
  use geom, only : volvm1
  use soln, only : pmask, vmult
  use tstep, only : tolps, istep
  use hsmg_routines, only : h1mg_solve
  implicit none

  real(DP) ::             res  (lx1*ly1*lz1*lelv)
  real(DP) ::             h1   (lx1,ly1,lz1,lelv)
  real(DP) ::             h2   (lx1,ly1,lz1,lelv)
  real(DP) ::             wt   (lx1,ly1,lz1,lelv)
  integer :: iter

  real(DP), allocatable :: x(:),r(:),w(:)
  real(DP), allocatable :: h(:,:),gamma(:)
  real(DP), allocatable :: c(:),s(:) ! store the Givens rotations
  real(DP), allocatable :: v(:,:) ! stores the orthogonal Krylov subspace basis
  real(DP), allocatable :: z(:,:) ! Z = M**(-1) V

  real(DP), allocatable, save :: ml(:), mu(:)
  real(DP) :: etime

  real(DP) :: divex
  real(DP), allocatable :: d(:)

  real(DP) :: wk1(lgmres)

  real(DP) :: alpha, l, temp
  integer :: outer

  logical, save :: iflag = .false., if_hyb = .false.
  real(DP), save ::  norm_fac

  real(DP) :: etime1, etime2, etime_p
  real(DP) :: rnorm, tolpss, div0, ratio
  real(DP), external :: glsc3, vlsc3

  integer :: m, n
  integer :: i, j, k, iconv

  ngmres = ngmres + 1
  etime = dnekclock()

  !> \todo move these allocations to where they are needed
  allocate(x(lx2*ly2*lz2*lelv))
  allocate(r(lx2*ly2*lz2*lelv))
  allocate(w(lx2*ly2*lz2*lelv)) 
  allocate(h(lgmres,lgmres))
  allocate(gamma(lgmres+1))
  allocate(c(lgmres), s(lgmres))

  allocate(v(lx2*ly2*lz2*lelv,lgmres)) ! verified
  allocate(z(lx2*ly2*lz2*lelv,lgmres)) ! verified


!  if (.not. allocated(ml)) allocate(ml(lx2*ly2*lz2*lelv))
!  if (.not. allocated(mu)) allocate(mu(lx2*ly2*lz2*lelv)) 

  !> \todo check if these inits are nessesary
  x = 0._dp; w = 0._dp; h = 0._dp; gamma=0._dp; c = 0._dp; s = 0._dp
  v = 0._dp; z = 0._dp;


  n = nx1*ny1*nz1*nelv

  etime1 = dnekclock()
  etime_p = 0.
  divex = 0.
  iter  = 0
  m     = lgmres

  if( .NOT. iflag) then
      iflag= .TRUE. 
!      call uzawa_gmres_split(ml,mu,n)
      norm_fac = 1./sqrt(volvm1)
  endif

  !> \todo Do we need this?
  allocate(d(lx1*ly1*lz1*lelv))
  if (param(100) /= 2) call set_fdm_prec_h1b(d,h1,h2,nelv)
  deallocate(d)

  call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
  if (param(21) > 0 .AND. tolps > abs(param(21))) &
  tolps = abs(param(21))
  if (istep == 0) tolps = 1.e-4
  tolpss = tolps

  iconv = 0

  outer = 0
  do while (iconv == 0 .AND. iter < 500)
      outer = outer+1

      if(iter == 0) then               !      -1
          r = res          ! r = L  res
      !           call copy(r,res,n)
      else
      ! pdate residual
          gmres_mop = gmres_mop + n
          call copy  (r,res,n)                  ! r = res
          etime = etime - dnekclock()
          call ax    (w,x,h1,h2,n)              ! w = A x
          etime = etime + dnekclock()
          r = r - w  
      !      -1
      !    r = r * ml ! r = L   r
      endif
  !            ______
      gmres_mop  = gmres_mop  + 2*n
      gmres_flop = gmres_flop + 3*n
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
      v(:,1) = r * temp  ! v  = r / gamma
  !  1            1
      do j=1,m
          iter = iter+1
      !       -1
          w =  v(:,j)               ! w  = U   v
      !           j

      ! . . . . . Overlapping Schwarz + coarse-grid . . . . . . .

          etime2 = dnekclock()

      !           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
          if (ifmgrid) then
              etime = etime - dnekclock()
              call h1mg_solve(z(1,j),w,if_hyb)   ! z  = M   w
              etime = etime + dnekclock()
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

               
          etime = etime - dnekclock()
          call ax  (w,z(1,j),h1,h2,n)           ! w = A z
          etime = etime + dnekclock()
      !        j
               
      !      -1
      !    w = w * ml  ! w = L   w

      !           !modified Gram-Schmidt

      !           do i=1,j
      !              h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
      !                                                 !  i,j       i

      !              call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
      !           enddo                                 !          i,j  i

      !           2-PASS GS, 1st pass:
          gmres_mop  = gmres_mop  + 2*n + j*n
          gmres_flop = gmres_flop + 3*j*n
          do i=1,j
              h(i,j)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
          enddo                                 !  i,j       i

          call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs

          gmres_mop  = gmres_mop  + (j+1)*n 
          gmres_flop = gmres_flop + j*2*n
          do i=1,j
              w = w - v(:,i) * h(i,j) ! w = w - h    v
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
          gmres_mop  = gmres_mop  + 2*n
          gmres_flop = gmres_flop + 3*n
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
          gmres_mop  = gmres_mop  + 2*n
          gmres_flop = gmres_flop + n
          v(:,j+1) = w * temp  ! v    = w / alpha
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
          x = x + z(:,i) * c(i)  ! x = x + c  z
      enddo                               !          i  i
  !        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
  enddo
  9000 continue

  divex = rnorm
  gmres_mop = gmres_mop + n
  call copy(res,x,n)

  call ortho   (res) ! Orthogonalize wrt null space, if present

  etime1 = dnekclock()-etime1
  if (nid == 0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p, &
  etime1,if_hyb
!   call flush_hack
  9999 format(i9,' PRES gmres:',i5,1p5e12.4,1x,l4)

  if (outer <= 2) if_hyb = .FALSE. 

  tgmres = tgmres + (dnekclock() - etime)

  return
end subroutine hmh_gmres
!-----------------------------------------------------------------------

