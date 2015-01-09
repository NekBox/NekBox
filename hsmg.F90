!-----------------------------------------------------------------------
!> \file hsmg.F90
!! \brief Module containing hybrid Schwarz multi-grid preconditioners
!!
!!  To do:
!!  1)  Why does hsmg_schwarz_toext2d not zero out a, whereas 3d does??  DONE
!!  2)  Convert all nelv refs to nelfld(ifield) or (nelmg?)  DONE
!!  3)  Define mg_schwarz_wt for up to and including mg_h1_lmax   DONE
!!  4)  MAKE CERTAIN common /hsmgw/ is LARGE enough in hsmg_tnsr and  DONE
!!      elsewhere!
!!  5)  Devise and implement UNIT tests, now, so that you can test
!!      pieces of the setup code in stages.
!!  6)  Start developing and testing, in a linear fashion, the SETUP driver.
!!  7)  Make certain dssum flags declared for all levels  DONE
!!  8)  Need TWO masks for each level:  one for A*x, and one for Schwarz!
!!      NO -- one is fine.
!!  9)  Modify axml so addition comes after dssum.  DONE
!-----------------------------------------------------------------------

! Some relevant parameters

! param(41):
!     0 - use additive SEMG
!     1 - use hybrid SEMG (not yet working... but coming soon!)

! param(42):   navier0.f, fasts.f
!     0 - use GMRES for iterative solver, use non-symmetric weighting
!     1 - use PCG for iterative solver, do not use weighting

! param(43):   uzawa_gmres.f, navier6.f
!     0 - use additive multilevel scheme (requires param(42).eq.0)
!     1 - use original 2 level scheme

! param(44):   fast3d.f, navier6.f
!     0 - base top-level additive Schwarz on restrictions of E
!     1 - base top-level additive Schwarz on restrictions of A

!----------------------------------------------------------------------

module hsmg_routines

  private
  public :: h1mg_setup, h1mg_solve
  public :: hsmg_setup

contains

!----------------------------------------------------------------------
!> \brief Sets up hybrid Schwarz multi-grid preconditioner
subroutine hsmg_setup()
  use hsmg, only : mg_fld
  use tstep, only : ifield
  implicit none

  mg_fld = 1
  if (ifield > 1) mg_fld = 2
  if (ifield == 1) call hsmg_index_0 ! initialize index sets

  call hsmg_setup_mg_nx  ! set nx values for each level of multigrid
  call hsmg_setup_semhat ! set spectral element hat matrices
  call hsmg_setup_intp
  call hsmg_setup_dssum  ! set direct stiffness summation handles
  call hsmg_setup_wtmask ! set restriction weight matrices and bc masks
  call hsmg_setup_fdm    ! set up fast diagonalization method
  call hsmg_setup_schwarz_wt( .FALSE. )
  call hsmg_setup_solve  ! set up the solver
!   call hsmg_setup_dbg

  return
end subroutine hsmg_setup

!----------------------------------------------------------------------
!> \brief Sets up Poisson preconditioner
subroutine h1mg_setup()
  use size_m, only : nx1, ny1, nz1, nelt
  use hsmg, only : mg_h1_lmax
  use input, only : param
  implicit none

  integer :: p_msk, n, l

  param(59) = 1
  call geom_reset(1)  ! Recompute g1m1 etc. with deformed only

  n = nx1*ny1*nz1*nelt

  call h1mg_setup_mg_nx
  call h1mg_setup_semhat ! SEM hat matrices for each level
  call hsmg_setup_intp   ! Interpolation operators
  call h1mg_setup_dssum  ! set direct stiffness summation handles
  call h1mg_setup_wtmask ! set restriction weight matrices and bc masks
  call h1mg_setup_fdm    ! set up fast diagonalization method
  call h1mg_setup_schwarz_wt( .FALSE. )
  call hsmg_setup_solve  ! set up the solver

  l=mg_h1_lmax
  !> \todo Is it nessesary to set h1 and h2 here?  They aren't inited until
  !! later
!  call mg_set_h1  (p_h1 ,l)
!  call mg_set_h2  (p_h2 ,l)
!  call mg_set_gb  (p_g,p_b,l)
  call mg_set_msk (p_msk,l)

  return
end subroutine h1mg_setup

!----------------------------------------------------------------------
!> \brief Solve preconditioner: z = M rhs, where \f$ M \approx A^{-1} \f$
!!
!! Assumes that preprocessing has been completed via h1mg_setup()
subroutine h1mg_solve(z,rhs,if_hybrid)  
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use hsmg, only : mg_h1_lmax, mg_h1_n, p_mg_msk, mg_imask, mg_fld ! Same array space as HSMG
  use tstep, only : nelfld, ifield
  implicit none

  real(DP), intent(out) :: z(*)      !>!< approximate solution to A z = rhs
  real(DP), intent(in)  :: rhs(*)    !>!< right hand side to Poisson equation
  logical,  intent(in)  :: if_hybrid !>!< Use hybrid or normal?
       
  integer, parameter :: lt=lx1*ly1*lz1*lelt
  real(DP), allocatable :: e(:),w(:),r(:)
  integer :: p_msk

  real(DP) :: op, om, sigma
  integer :: nel, l, n, is, im, i1, i

  nel   = nelfld(ifield)

  op    =  1.                                     ! Coefficients for h1mg_ax
  om    = -1.
  sigma =  1.
  if (if_hybrid) sigma = 2./3.

  l     = mg_h1_lmax
  n     = mg_h1_n(l,mg_fld)
  is    = 1                                       ! solve index
  call h1mg_schwarz(z,rhs,sigma,l)                ! z := sigma W M  rhs
!               Schwarz
  allocate(r(lt))
  call copy(r,rhs,n)                              ! r  := rhs
!max    if (if_hybrid) call h1mg_axm(r,z,op,om,l,w)     ! r  := rhs - A z
!  l

  allocate(e(2*lt)); e = 0_dp
  do l = mg_h1_lmax-1,2,-1                        ! DOWNWARD Leg of V-cycle
      is = is + n
      n  = mg_h1_n(l,mg_fld)
  !          T
      call h1mg_rstr(r,l, .TRUE. )                   ! r   :=  J r
  !  l         l+1
  !        OVERLAPPING Schwarz exchange and solve:
      call h1mg_schwarz(e(is),r,sigma,l)           ! e := sigma W M       r
  !  l            Schwarz l

!max        if(if_hybrid)call h1mg_axm(r,e(is),op,om,l,w)! r  := r - A e
  !  l           l
  enddo
  is = is+n
!         T
  call h1mg_rstr(r,1, .FALSE. )                     ! r  :=  J  r
!  l         l+1
  p_msk = p_mg_msk(l,mg_fld)
  call h1mg_mask(r,mg_imask(p_msk),nel)           !        -1
  call hsmg_coarse_solve ( e(is) , r )            ! e  := A   r
  call h1mg_mask(e(is),mg_imask(p_msk),nel)       !  1     1   1
  deallocate(r)

  allocate(w(lt)); w = 0_dp
  do l = 2,mg_h1_lmax-1                           ! UNWIND.  No smoothing.
      im = is
      is = is - n
      n  = mg_h1_n(l,mg_fld)
      call hsmg_intp (w,e(im),l-1)                 ! w   :=  J e
      i1=is-1                                      !            l-1
      do i=1,n
          e(i1+i) = e(i1+i) + w(i)                  ! e   :=  e  + w
      enddo                                        !  l       l
  enddo

  l  = mg_h1_lmax
  n  = mg_h1_n(l,mg_fld)
  im = is  ! solve index
  call hsmg_intp(w,e(im),l-1)                     ! w   :=  J e
  do i = 1,n                                      !            l-1
      z(i) = z(i) + w(i)                           ! z := z + w
  enddo
  deallocate(w,e)

  call dsavg(z) ! Emergency hack --- to ensure continuous z!

  return
end subroutine h1mg_solve

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PRIVATE IMPLEMENTATIONS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine hsmg_setup_semhat
  use input, only : if3d
  use hsmg, only : mg_zh, mg_lmax, mg_ah, mg_bh, mg_dh, mg_dht, mg_zh
  use hsmg, only : mg_nx, mg_nh, mg_nhz, mg_nz
  use semhat, only : ah, bh, ch, dh, zh, dph, jph, bgl, zgl, dgl, jgl, wh
  implicit none

  integer :: n,l
!   generate the SEM hat matrices for each level
!   top level
  n = mg_nx(mg_lmax)
  call generate_semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
  call copy(mg_zh(1,mg_lmax),zgl,n-1) !top level based on gl points
  mg_nh(mg_lmax)=n-1
  mg_nhz(mg_lmax)=n-1
  if( .NOT. if3d) mg_nhz(mg_lmax)=1
!   lower levels
  do l=1,mg_lmax-1
      n = mg_nx(l)
      if(n > 1) then
          call generate_semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
          call copy(mg_ah(1,l),ah,(n+1)*(n+1))
          call copy(mg_bh(1,l),bh,n+1)
          call copy(mg_dh(1,l),dh,(n+1)*(n+1))
          call transpose(mg_dht(1,l),n+1,dh,n+1)
          call copy(mg_zh(1,l),zh,n+1)
      else
          mg_zh(1,l) = -1.
          mg_zh(2,l) =  1.
      endif
      mg_nh(l)=n+1
      mg_nhz(l)=mg_nz(l)+1
  enddo
end subroutine hsmg_setup_semhat

!----------------------------------------------------------------------
subroutine hsmg_setup_intp
  use hsmg, only : mg_lmax, mg_nh, mg_jh, mg_zh, mg_jht!, mg_jhfc, mg_jhfct
  implicit none

  integer :: l,nf,nc

  do l=1,mg_lmax-1

      nf=mg_nh(l+1)
      nc=mg_nh(l)

  !        Standard multigrid coarse-to-fine interpolation
      call hsmg_setup_intpm( &
      mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
      call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)

  !        Fine-to-coarse interpolation for variable-coefficient operators
!      call hsmg_setup_intpm( &
!      mg_jhfc(1,l),mg_zh(1,l),mg_zh(1,l+1),nc,nf)
!      call transpose(mg_jhfct(1,l),nf,mg_jhfc(1,l),nc)
  !        call outmat(mg_jhfc(1,l),nc,nf,'MG_JHFC',l)

  enddo
end subroutine hsmg_setup_intp

!----------------------------------------------------------------------
subroutine hsmg_setup_intpm(jh,zf,zc,nf,nc)
  use kinds, only : DP
  use size_m, only : lx1
  implicit none
  integer :: nf,nc
  real(DP) :: jh(nf,nc),zf(*),zc(*)
  real(DP) :: w(2*lx1+2)

  integer :: i, j
  do i=1,nf
      call fd_weights_full(zf(i),zc,nc-1,1,w)
      do j=1,nc
          jh(i,j)=w(j)
      enddo
  enddo
  return
end subroutine hsmg_setup_intpm

!----------------------------------------------------------------------
subroutine hsmg_setup_dssum
  use kinds, only : i8
  use size_m, only : lx1, ly1, lz1, lelv, nelv, ndim
  use input, only : if3d
  use hsmg, only : mg_lmax, mg_nh, mg_nhz, mg_fld
  use hsmg, only : mg_gsh_handle, mg_gsh_schwarz_handle
  use mesh, only : vertex
  use parallel, only : nelgv
  implicit none

  integer, parameter :: lxyz=(lx1+2)*(ly1+2)*(lz1+2)

  integer(i8), allocatable :: glo_num(:)
  integer :: nx,ny,nz, ncrnr
  integer :: l
       

!     set up direct stiffness summation for each level
  ncrnr = 2**ndim
  call get_vert()

  allocate(glo_num(lxyz*lelv))

  do l=1,mg_lmax-1
      nx=mg_nh(l)
      ny=mg_nh(l)
      nz=mg_nhz(l)
      call setupds(mg_gsh_handle(l,mg_fld),nx,ny,nz &
      ,nelv,nelgv,vertex,glo_num)
      nx=nx+2
      ny=ny+2
      nz=nz+2
      if( .NOT. if3d) nz=1
      call setupds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz &
      ,nelv,nelgv,vertex,glo_num)
  enddo
end subroutine hsmg_setup_dssum

!----------------------------------------------------------------------
subroutine h1mg_setup_wtmask
  use kinds, only : DP
  use size_m, only : ndim, ldim, nelv, lelv
  use hsmg, only : mg_mask_index, mg_lmax, mg_rstr_wt_index, mg_nh, mg_nhz
  use hsmg, only : mg_fld, lmgs, lmg_rwt, mg_rstr_wt
  implicit none

  real(DP), allocatable :: work(:)
  integer :: i,l, itmp

  allocate(work(maxval(mg_nh(1:mg_lmax))*maxval(mg_nh(1:mg_lmax))*maxval(mg_nhz(1:mg_lmax))* nelv))

  i = mg_mask_index(mg_lmax,mg_fld-1)
  do l=1,mg_lmax
      mg_rstr_wt_index(l,mg_fld)=i
      mg_mask_index   (l,mg_fld)=i
      i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelv
      if(i > lmgs*lmg_rwt*2*ldim*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
          call exitt
      endif
      call hsmg_setup_rstr_wt( mg_rstr_wt(mg_rstr_wt_index(l,mg_fld)) &
                             , mg_nh(l),mg_nh(l),mg_nhz(l),l, work)
!      call hsmg_setup_mask( mg_mask(mg_mask_index(l,mg_fld)) &
!                          , mg_nh(l),mg_nh(l),mg_nhz(l),l, work)
  enddo
  mg_mask_index(l,mg_fld)=i

end subroutine h1mg_setup_wtmask

!----------------------------------------------------------------------
subroutine hsmg_setup_wtmask
  use kinds, only : DP
  use size_m, only : ndim, ldim, nelv, lelv
  use hsmg, only : mg_mask_index, mg_lmax, mg_rstr_wt_index, mg_nh, mg_nhz
  use hsmg, only : lmgs, lmg_rwt, mg_rstr_wt, mg_fld
  implicit none

  real(DP), allocatable :: work(:)
  integer :: i,l, itmp

  allocate(work(maxval(mg_nh(1:mg_lmax-1)) &
               *maxval(mg_nh(1:mg_lmax-1)) &
               *maxval(mg_nhz(1:mg_lmax-1))* nelv))

  i = mg_mask_index(mg_lmax,mg_fld-1)
  do l=1,mg_lmax-1
      mg_rstr_wt_index(l,mg_fld)=i
      mg_mask_index   (l,mg_fld)=i
      i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelv
      if(i > lmgs*lmg_rwt*2*ldim*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
          call exitt
      endif
      call hsmg_setup_rstr_wt( mg_rstr_wt(mg_rstr_wt_index(l,mg_fld)) &
                             , mg_nh(l),mg_nh(l),mg_nhz(l),l,work)
!      call hsmg_setup_mask( mg_mask(mg_mask_index(l,mg_fld)) &
!                          , mg_nh(l),mg_nh(l),mg_nhz(l),l,work)
  enddo
  mg_mask_index(l,mg_fld)=i
end subroutine hsmg_setup_wtmask

!----------------------------------------------------------------------
subroutine hsmg_intp(uf,uc,l) ! l is coarse level
  use kinds, only : DP
  use hsmg, only : mg_nh, mg_jh, mg_jht
  implicit none
  real(DP) :: uf(*),uc(*)
  integer :: l
  call hsmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
  return
end subroutine hsmg_intp

!----------------------------------------------------------------------
!> \brief     computes v = [A (x) A] u or v = [A (x) A (x) A] u
subroutine hsmg_tnsr(v,nv,u,nu,A,At)
  use kinds, only : DP
  use input, only : if3d
  implicit none

  integer :: nv,nu
  real(DP) :: v(*),u(*),A(*),At(*)

  if ( .NOT. if3d) then
!max        call hsmg_tnsr2d(v,nv,u,nu,A,At)
  else
      call hsmg_tnsr3d(v,nv,u,nu,A,At,At)
  endif
  return
end subroutine hsmg_tnsr

!----------------------------------------------------------------------
!> \brief computes:  v = [C (x) B (x) A] u .
subroutine hsmg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
  use kinds, only : DP
  use size_m
  implicit none

  integer :: nv,nu
  real(DP) :: v(nv*nv*nv,nelv),u(nu*nu*nu,nelv),A(*),Bt(*),Ct(*)

  integer, parameter :: lwk=(lx1+2)*(ly1+2)*(lz1+2)
  real(DP) :: work(0:lwk-1),work2(0:lwk-1)
  integer :: ie, i

  do ie=1,nelv
      call mxm(A,nv,u(1,ie),nu,work,nu*nu)
      do i=0,nu-1
          call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
  enddo
  return
end subroutine hsmg_tnsr3d

!----------------------------------------------------------------------
!> \brief computes  v = [C (x) B (x) A] u
subroutine hsmg_tnsr3d_el(v,nv,u,nu,A,Bt,Ct)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1
  implicit none

  integer :: nv,nu
  real(DP) :: v(nv*nv*nv),u(nu*nu*nu),A(*),Bt(*),Ct(*)

  integer, parameter :: lwk=(lx1+2)*(ly1+2)*(lz1+2)
  real(DP) :: work(0:lwk-1),work2(0:lwk-1)
  integer :: i

  call mxm(A,nv,u,nu,work,nu*nu)
  do i=0,nu-1
      call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
  enddo
  call mxm(work2,nv*nv,Ct,nu,v,nv)

  return
end subroutine hsmg_tnsr3d_el

!----------------------------------------------------------------------
subroutine hsmg_dssum(u,l)
  use kinds, only : DP
  use ctimer, only : dnekclock, etime1, tdadd, ifsync
  use hsmg, only : mg_gsh_handle, mg_fld
  implicit none
  real(DP) :: u(*)
  integer :: l

  if (ifsync) call nekgsync()
#ifndef NOTIMER
  etime1=dnekclock()
#endif
  call gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)
#ifndef NOTIMER
  tdadd =tdadd + dnekclock()-etime1
#endif

  return
end subroutine hsmg_dssum

!----------------------------------------------------------------------
subroutine hsmg_dsprod(u,l)
  use kinds, only : DP
  use ctimer, only : ifsync
  use hsmg, only : mg_gsh_handle, mg_fld
  implicit none
  real(DP) :: u(1)
  integer :: l

  if (ifsync) call nekgsync()

  call gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)
  return
end subroutine hsmg_dsprod

!----------------------------------------------------------------------
subroutine hsmg_schwarz_dssum(u,l)
  use kinds, only : DP
  use ctimer, only : ifsync, etime1, dnekclock, tdadd
  use hsmg, only : mg_gsh_schwarz_handle, mg_fld
  implicit none
  real(DP) :: u(1)
  integer :: l

  if (ifsync) call nekgsync()
#ifndef NOTIMER
  etime1=dnekclock()
#endif
  call gs_op(mg_gsh_schwarz_handle(l,mg_fld),u,1,1,0)
#ifndef NOTIMER
  tdadd =tdadd + dnekclock()-etime1
#endif
  return
end subroutine hsmg_schwarz_dssum

!----------------------------------------------------------------------
subroutine hsmg_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
  use kinds, only : DP
  use size_m, only : nelv
  use input, only : if3d
  implicit none

  integer :: l1,l2,nx,ny,nz
  real(DP) :: arr1(nx,ny,nz,nelv),arr2(nx,ny,nz,nelv)
  real(DP) :: f1,f2
        
  integer :: i,j,k,ie,i0,i1
  i0=2
  i1=nx-1
        
  if( .NOT. if3d) then
      do ie=1,nelv
          do j=i0,i1
              arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie) &
              +f2*arr2(l2+1 ,j,1,ie)
              arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie) &
              +f2*arr2(nx-l2,j,1,ie)
          enddo
          do i=i0,i1
              arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie) &
              +f2*arr2(i,l2+1 ,1,ie)
              arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie) &
              +f2*arr2(i,nx-l2,1,ie)
          enddo
      enddo
  else
      do ie=1,nelv
          do k=i0,i1
              do j=i0,i1
                  arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie) &
                  +f2*arr2(l2+1 ,j,k,ie)
                  arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie) &
                  +f2*arr2(nx-l2,j,k,ie)
              enddo
          enddo
          do k=i0,i1
              do i=i0,i1
                  arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie) &
                  +f2*arr2(i,l2+1 ,k,ie)
                  arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie) &
                  +f2*arr2(i,nx-l2,k,ie)
              enddo
          enddo
          do j=i0,i1
              do i=i0,i1
                  arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie) &
                  +f2*arr2(i,j,l2+1 ,ie)
                  arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie) &
                  +f2*arr2(i,j,nx-l2,ie)
              enddo
          enddo
      enddo
  endif
  return
end subroutine hsmg_extrude

!----------------------------------------------------------------------
subroutine h1mg_schwarz(e,r,sigma,l)
  use kinds, only : DP
  use hsmg, only : mg_h1_n, mg_fld
  implicit none

  real(DP) :: e(*),r(*)
  real(DP), intent(in) :: sigma
  integer :: l, n

  n = mg_h1_n(l,mg_fld)

  call h1mg_schwarz_part1 (e,r,l)
  call hsmg_schwarz_wt    (e,l)          ! e  := W e
  e(1:n) = e(1:n) * sigma

  return
end subroutine h1mg_schwarz

!----------------------------------------------------------------------
subroutine h1mg_schwarz_part1 (e,r,l)
  use kinds, only : DP
  use size_m, only : nelv
  use input, only : if3d
  use hsmg, only : mg_h1_n, p_mg_msk, mg_imask, mg_nh, mg_fld
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: e(*),r(*)

  integer :: enx,eny,enz,pm, n, i, l
  real(DP) :: zero, one, onem
  real(DP), allocatable :: work(:)

  zero =  0
  one  =  1
  onem = -1

  n  = mg_h1_n (l,mg_fld)
  pm = p_mg_msk(l,mg_fld)

  enx=mg_nh(l)+2
  eny=mg_nh(l)+2
  enz=mg_nh(l)+2

  call h1mg_mask  (r,mg_imask(pm),nelfld(ifield))  ! Zero Dirichlet nodes

  allocate(work(2*enx*eny*enz*nelv))

  if (if3d) then ! extended array
      call hsmg_schwarz_toext3d(work,r,mg_nh(l))
  else
!max        call hsmg_schwarz_toext2d(mg_work,r,mg_nh(l))
  endif

  if( .NOT. if3d) enz=1
  i = enx*eny*enz*nelv+1
     
!     exchange interior nodes
  call hsmg_extrude(work,0,zero,work,2,one,enx,eny,enz)
  call hsmg_schwarz_dssum(work,l)
  call hsmg_extrude(work,0,one ,work,2,onem,enx,eny,enz)

  call hsmg_fdm(work(i),work,l) ! Do the local solves

!     Sum overlap region (border excluded)
  call hsmg_extrude(work,0,zero,work(i),0,one ,enx,eny,enz)
  call hsmg_schwarz_dssum(work(i),l)
  call hsmg_extrude(work(i),0,one ,work,0,onem,enx,eny,enz)
  call hsmg_extrude(work(i),2,one,work(i),0,one,enx,eny,enz)

  if( .NOT. if3d) then ! Go back to regular size array
!max        call hsmg_schwarz_toreg2d(e,mg_work(i),mg_nh(l))
  else
      call hsmg_schwarz_toreg3d(e,work(i),mg_nh(l))
  endif

  call hsmg_dssum(e,l)                           ! sum border nodes
  call h1mg_mask (e,mg_imask(pm),nelfld(ifield)) ! apply mask

  return
end subroutine h1mg_schwarz_part1

!----------------------------------------------------------------------
subroutine hsmg_schwarz_toext3d(a,b,n)
  use kinds, only : DP
  use size_m, only : nelv
  implicit none
  integer :: n
  real(DP) :: a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
        
  integer :: i,j,k,ie
  a = 0._dp
  do ie=1,nelv
      do k=1,n
          do j=1,n
              do i=1,n
                  a(i,j,k,ie)=b(i,j,k,ie)
              enddo
          enddo
      enddo
  enddo
  return
end subroutine hsmg_schwarz_toext3d

!----------------------------------------------------------------------
subroutine hsmg_schwarz_toreg3d(b,a,n)
  use kinds, only : DP
  use size_m, only : nelv
  implicit none

  integer :: n
  real(DP) :: a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
        
  integer :: i,j,k,ie
  do ie=1,nelv
      do k=1,n
          do j=1,n
              do i=1,n
                  b(i,j,k,ie)=a(i,j,k,ie)
              enddo
          enddo
      enddo
  enddo
  return
end subroutine hsmg_schwarz_toreg3d

!----------------------------------------------------------------------
subroutine h1mg_setup_fdm()
  use size_m, only : ndim, ldim, nelv, lelv
  use hsmg, only : mg_fast_d_index, mg_fast_s_index, mg_nx, mg_bh, mg_ah
  use hsmg, only : mg_fast_d, mg_fast_s, lmg_fastd, lmg_fasts, mg_nh, mg_fld
  use hsmg, only : mg_lmax
  implicit none
        
  integer :: l,i,j,nl, itmp
  i = mg_fast_s_index(mg_lmax,mg_fld-1)
  j = mg_fast_d_index(mg_lmax,mg_fld-1)
  do l=2,mg_lmax
      mg_fast_s_index(l,mg_fld)=i
      nl = mg_nh(l)+2
      i=i+nl*nl*2*ndim*nelv
      if(i > lmg_fasts*2*ldim*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
          call exitt
      endif
      mg_fast_d_index(l,mg_fld)=j
      j=j+(nl**ndim)*nelv
      if(j > lmg_fastd*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
          call exitt
      endif
      call hsmg_setup_fast( &
      mg_fast_s(mg_fast_s_index(l,mg_fld)) &
      ,mg_fast_d(mg_fast_d_index(l,mg_fld)) &
      ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
  enddo
  mg_fast_s_index(l,mg_fld)=i
  mg_fast_d_index(l,mg_fld)=j
  return
end subroutine h1mg_setup_fdm

!----------------------------------------------------------------------
subroutine hsmg_setup_fdm()
  use size_m, only : ndim, ldim, nelv, lelv
  use hsmg, only : mg_fast_d_index, mg_fast_s_index, mg_nx, mg_bh, mg_ah
  use hsmg, only : mg_fast_d, mg_fast_s, lmg_fastd, lmg_fasts, mg_nh
  use hsmg, only : mg_lmax, mg_fld
  implicit none
        
  integer :: l,i,j,nl, itmp
  i = mg_fast_s_index(mg_lmax,mg_fld-1)
  j = mg_fast_d_index(mg_lmax,mg_fld-1)
  do l=2,mg_lmax-1
      mg_fast_s_index(l,mg_fld)=i
      nl = mg_nh(l)+2
      i=i+nl*nl*2*ndim*nelv
      if(i > lmg_fasts*2*ldim*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
          call exitt
      endif
      mg_fast_d_index(l,mg_fld)=j
      j=j+(nl**ndim)*nelv
      if(j > lmg_fastd*lelv) then
          itmp = i/(2*ldim*lelv)
          write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
          call exitt
      endif
      call hsmg_setup_fast( &
      mg_fast_s(mg_fast_s_index(l,mg_fld)) &
      ,mg_fast_d(mg_fast_d_index(l,mg_fld)) &
      ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
  enddo
  mg_fast_s_index(l,mg_fld)=i
  mg_fast_d_index(l,mg_fld)=j
  return
end subroutine hsmg_setup_fdm

!----------------------------------------------------------------------
!> \brief not sure    
subroutine hsmg_setup_fast(s,d,nl,ah,bh,n)
  use kinds, only : DP
  use hsmg, only : lr, llr, lrr, lmr, ls, lls, lms, lrs, lt, llt, lmt, lrt
  use size_m, only : nid, ndim, nelv
  use input, only : if3d
  implicit none

  integer :: nl
  real(DP) :: s(nl*nl,2,ndim,nelv)
  real(DP) :: d(nl**ndim,nelv)
  real(DP) :: ah(1),bh(1)
  integer :: n
      
  integer :: i,j,k
  integer :: ie,il,nr,ns,nt
  integer :: lbr,rbr,lbs,rbs,lbt,rbt,two
  integer :: ierr, ierrmx
  integer, external :: iglmax
  real(DP) :: eps,diag

  two  = 2
  ierr = 0
  do ie=1,nelv
      call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
      nr=nl
      ns=nl
      nt=nl
      call hsmg_setup_fast1d(s(1,1,1,ie),lr,nr,lbr,rbr &
      ,llr(ie),lmr(ie),lrr(ie),ah,bh,n,ie)
      call hsmg_setup_fast1d(s(1,1,2,ie),ls,ns,lbs,rbs &
      ,lls(ie),lms(ie),lrs(ie),ah,bh,n,ie)
      if(if3d) call hsmg_setup_fast1d(s(1,1,3,ie),lt,nt,lbt,rbt &
      ,llt(ie),lmt(ie),lrt(ie),ah,bh,n,ie)
      il=1
      if( .NOT. if3d) then
          eps = 1.e-5*(maxval(lr(2:nr-1)) + maxval(ls(2:ns-1)))
          do j=1,ns
              do i=1,nr
                  diag = lr(i)+ls(j)
                  if (diag > eps) then
                      d(il,ie) = 1.0/diag
                  else
                  !                 write(6,2) ie,'Reset Eig in hsmg setup fast:',i,j,l
                  !    $                         ,eps,diag,lr(i),ls(j)
                  !    2 format(i6,1x,a21,3i5,1p4e12.4)
                      d(il,ie) = 0.0
                  endif
                  il=il+1
              enddo
          enddo
      else
          eps = 1.e-5 * (maxval(lr(2:nr-1)) &
          + maxval(ls(2:ns-1)) + maxval(lt(2:nt-1)))
          do k=1,nt
              do j=1,ns
                  do i=1,nr
                      diag = lr(i)+ls(j)+lt(k)
                      if (diag > eps) then
                          d(il,ie) = 1.0/diag
                      else
                      !                 write(6,3) ie,'Reset Eig in hsmg setup fast:',i,j,k,l
                      !    $                         ,eps,diag,lr(i),ls(j),lt(k)
                      !    3 format(i6,1x,a21,4i5,1p5e12.4)
                          d(il,ie) = 0.0
                      endif
                      il=il+1
                  enddo
              enddo
          enddo
      endif
  enddo

  ierrmx = iglmax(ierr,1)
  if (ierrmx > 0) then
      if (ierr > 0) write(6,*) nid,ierr,' BC FAIL'
      call exitti('A INVALID BC FOUND in genfast$',ierrmx)
  endif

  return
end subroutine hsmg_setup_fast

!----------------------------------------------------------------------
subroutine hsmg_setup_fast1d(s,lam,nl,lbc,rbc,ll,lm,lr,ah,bh,n,ie)
  use kinds, only : DP          
  use size_m, only : lx1
  implicit none

  integer :: nl,lbc,rbc,n, ie
  real(DP) :: s(nl,nl,2),lam(nl),ll,lm,lr
  real(DP) :: ah(0:n,0:n),bh(0:n)

  integer, parameter :: lxm=lx1+2
  real(DP) :: b(2*lxm*lxm),w(2*lxm*lxm)
        
  call hsmg_setup_fast1d_a(s,lbc,rbc,ll,lm,lr,ah,n)
  call hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
          
!   if (nid.eq.0) write(6,*) 'THIS is generalev call',nl,lbc
  call generalev(s,b,lam,nl,w)
  if(lbc > 0) call row_zero(s,nl,nl,1)
  if(lbc == 1) call row_zero(s,nl,nl,2)
  if(rbc > 0) call row_zero(s,nl,nl,nl)
  if(rbc == 1) call row_zero(s,nl,nl,nl-1)
        
  call transpose(s(1,1,2),nl,s,nl)
  return
end subroutine hsmg_setup_fast1d

!----------------------------------------------------------------------
subroutine hsmg_setup_fast1d_a(a,lbc,rbc,ll,lm,lr,ah,n)
  use kinds, only : DP
  implicit none

  integer :: lbc,rbc,n
  real(DP) :: a(0:n+2,0:n+2),ll,lm,lr
  real(DP) :: ah(0:n,0:n)
        
  real(DP) :: fac
  integer :: i,j,i0,i1
  i0=0
  if(lbc == 1) i0=1
  i1=n
  if(rbc == 1) i1=n-1
       
  a = 0._dp 
  fac = 2.0/lm
  a(1,1)=1.0
  a(n+1,n+1)=1.0
  do j=i0,i1
      do i=i0,i1
          a(i+1,j+1)=fac*ah(i,j)
      enddo
  enddo
  if(lbc == 0) then
      fac = 2.0/ll
      a(0,0)=fac*ah(n-1,n-1)
      a(1,0)=fac*ah(n  ,n-1)
      a(0,1)=fac*ah(n-1,n  )
      a(1,1)=a(1,1)+fac*ah(n  ,n  )
  else
      a(0,0)=1.0
  endif
  if(rbc == 0) then
      fac = 2.0/lr
      a(n+1,n+1)=a(n+1,n+1)+fac*ah(0,0)
      a(n+2,n+1)=fac*ah(1,0)
      a(n+1,n+2)=fac*ah(0,1)
      a(n+2,n+2)=fac*ah(1,1)
  else
      a(n+2,n+2)=1.0
  endif
  return
end subroutine hsmg_setup_fast1d_a

!----------------------------------------------------------------------
subroutine hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
  use kinds, only : DP
  implicit none

  integer :: lbc,rbc,n
  real(DP) :: b(0:n+2,0:n+2),ll,lm,lr
  real(DP) :: bh(0:n)
        
  real(DP) :: fac
  integer :: i,i0,i1
  i0=0
  if(lbc == 1) i0=1
  i1=n
  if(rbc == 1) i1=n-1
        
  b = 0._dp
  fac = 0.5*lm
  b(1,1)=1.0
  b(n+1,n+1)=1.0
  do i=i0,i1
      b(i+1,i+1)=fac*bh(i)
  enddo
  if(lbc == 0) then
      fac = 0.5*ll
      b(0,0)=fac*bh(n-1)
      b(1,1)=b(1,1)+fac*bh(n  )
  else
      b(0,0)=1.0
  endif
  if(rbc == 0) then
      fac = 0.5*lr
      b(n+1,n+1)=b(n+1,n+1)+fac*bh(0)
      b(n+2,n+2)=fac*bh(1)
  else
      b(n+2,n+2)=1.0
  endif
  return
end subroutine hsmg_setup_fast1d_b

!----------------------------------------------------------------------
!> \brief clobbers r
subroutine hsmg_fdm(e,r,l)
  use kinds, only : DP
  use hsmg, only : mg_fast_s, mg_fast_d, mg_fast_s_index, mg_fast_d_index
  use hsmg, only : mg_nh, mg_fld
  implicit none
  real(DP) :: e(*), r(*)
  integer :: l
  call hsmg_do_fast(e,r, &
  mg_fast_s(mg_fast_s_index(l,mg_fld)), &
  mg_fast_d(mg_fast_d_index(l,mg_fld)), &
  mg_nh(l)+2)
  return
end subroutine hsmg_fdm

!----------------------------------------------------------------------
!> \brief clobbers r
subroutine hsmg_do_fast(e,r,s,d,nl)
  use kinds, only : DP
  use size_m, only : ndim, nelv
  use input, only : if3d
  implicit none

  integer :: nl
  real(DP) :: e(nl**ndim,nelv)
  real(DP) :: r(nl**ndim,nelv)
  real(DP) :: s(nl*nl,2,ndim,nelv)
  real(DP) :: d(nl**ndim,nelv)
        
  integer :: ie,nn,i
  nn=nl**ndim
  if( .NOT. if3d) then
#if 0
      do ie=1,nelv
          call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl &
          ,s(1,2,1,ie),s(1,1,2,ie))
          do i=1,nn
              r(i,ie)=d(i,ie)*e(i,ie)
          enddo
          call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl &
          ,s(1,1,1,ie),s(1,2,2,ie))
      enddo
#endif
  else
      do ie=1,nelv
          call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl &
          ,s(1,2,1,ie),s(1,1,2,ie),s(1,1,3,ie))
          do i=1,nn
              r(i,ie)=d(i,ie)*e(i,ie)
          enddo
          call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl &
          ,s(1,1,1,ie),s(1,2,2,ie),s(1,2,3,ie))
      enddo
  endif
  return
end subroutine hsmg_do_fast

!----------------------------------------------------------------------
!> \brief u = wt .* u
subroutine hsmg_do_wt(u,wt,nx,ny,nz)
  use kinds, only : DP
  use size_m, only : nelv, ndim
  use input, only : if3d
  implicit none

  integer :: nx,ny,nz
  real(DP) :: u(nx,ny,nz,nelv)
  real(DP) :: wt(nx,nz,2,ndim,nelv)
        
  integer :: i, j, k, ie

!   if (nx.eq.2) then
!      do e=1,nelv
!         call outmat(wt(1,1,1,1,e),nx,nz,'wt 1-1',e)
!         call outmat(wt(1,1,2,1,e),nx,nz,'wt 2-1',e)
!         call outmat(wt(1,1,1,2,e),nx,nz,'wt 1-2',e)
!         call outmat(wt(1,1,2,2,e),nx,nz,'wt 2-2',e)
!      enddo
!      call exitti('hsmg_do_wt quit$',nelv)
!   endif

  if ( .NOT. if3d) then
      do ie=1,nelv
          do j=1,ny
              u( 1,j,1,ie)=u( 1,j,1,ie)*wt(j,1,1,1,ie)
              u(nx,j,1,ie)=u(nx,j,1,ie)*wt(j,1,2,1,ie)
          enddo
          do i=2,nx-1
              u(i, 1,1,ie)=u(i, 1,1,ie)*wt(i,1,1,2,ie)
              u(i,ny,1,ie)=u(i,ny,1,ie)*wt(i,1,2,2,ie)
          enddo
      enddo
  else
      do ie=1,nelv
          do k=1,nz
              do j=1,ny
                  u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
                  u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
              enddo
          enddo
          do k=1,nz
              do i=2,nx-1
                  u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
                  u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
              enddo
          enddo
          do j=2,ny-1
              do i=2,nx-1
                  u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
                  u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
              enddo
          enddo
      enddo
  endif
  return
end subroutine hsmg_do_wt

!----------------------------------------------------------------------
subroutine hsmg_setup_rstr_wt(wt,nx,ny,nz,l,w)
  use kinds, only : DP
  use size_m, only : nelv, ndim
  use input, only : if3d
  implicit none

  integer, intent(in) :: nx,ny,nz,l
  real(DP), intent(out) :: w(nx,ny,nz,nelv)
  real(DP), intent(out) :: wt(nx,nz,2,ndim,nelv)
        
  integer :: ie, i, j, k
! nit border nodes to 1
  w = 0._dp
!     print *, 'Setup rstr wt: ',nx,ny,nz,nelv
  if ( .NOT. if3d) then
      do ie=1,nelv
          do i=1,nx
              w(i,1,1,ie)=1.0
              w(i,ny,1,ie)=1.0
          enddo
          do j=1,ny
              w(1,j,1,ie)=1.0
              w(nx,j,1,ie)=1.0
          enddo
      enddo
  else
      do ie=1,nelv
          do j=1,ny
              do i=1,nx
                  w(i,j,1,ie)=1.0
                  w(i,j,nz,ie)=1.0
              enddo
          enddo
          do k=1,nz
              do i=1,nx
                  w(i,1,k,ie)=1.0
                  w(i,ny,k,ie)=1.0
              enddo
          enddo
          do k=1,nz
              do j=1,ny
                  w(1,j,k,ie)=1.0
                  w(nx,j,k,ie)=1.0
              enddo
          enddo
      enddo
  endif
  call hsmg_dssum(w,l)
! nvert the count w to get the weight wt
  if ( .NOT. if3d) then
      do ie=1,nelv
          do j=1,ny
              wt(j,1,1,1,ie)=1.0/w(1,j,1,ie)
              wt(j,1,2,1,ie)=1.0/w(nx,j,1,ie)
          enddo
          do i=1,nx
              wt(i,1,1,2,ie)=1.0/w(i,1,1,ie)
              wt(i,1,2,2,ie)=1.0/w(i,ny,1,ie)
          enddo
      enddo
  else
      do ie=1,nelv
          do k=1,nz
              do j=1,ny
                  wt(j,k,1,1,ie)=1.0/w(1,j,k,ie)
                  wt(j,k,2,1,ie)=1.0/w(nx,j,k,ie)
              enddo
          enddo
          do k=1,nz
              do i=1,nx
                  wt(i,k,1,2,ie)=1.0/w(i,1,k,ie)
                  wt(i,k,2,2,ie)=1.0/w(i,ny,k,ie)
              enddo
          enddo
          do j=1,ny
              do i=1,nx
                  wt(i,j,1,3,ie)=1.0/w(i,j,1,ie)
                  wt(i,j,2,3,ie)=1.0/w(i,j,nz,ie)
              enddo
          enddo
      enddo
  endif
  return
end subroutine hsmg_setup_rstr_wt

!----------------------------------------------------------------------
subroutine hsmg_setup_schwarz_wt(ifsqrt)
  use size_m, only : ndim, nelv, ldim, lelv
  use input, only : if3d
  use hsmg, only : mg_schwarz_wt_index, mg_lmax, mg_fld
  use hsmg, only : mg_nh, lmg_swt, mg_schwarz_wt
  implicit none

  logical :: ifsqrt
        
  integer :: l,i,nl,nlz, itmp

  i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)
  do l=2,mg_lmax-1
      mg_schwarz_wt_index(l,mg_fld)=i
      nl = mg_nh(l)
      nlz = mg_nh(l)
      if( .NOT. if3d) nlz=1
      i=i+nl*nlz*4*ndim*nelv
      if(i > lmg_swt*4*ldim*lelv) then
          itmp = i/(4*ldim*lelv)
          write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
          call exitt
      endif

      call h1mg_setup_schwarz_wt_1( &
      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)

  enddo
  mg_schwarz_wt_index(l,mg_fld)=i

  return
end subroutine hsmg_setup_schwarz_wt

!----------------------------------------------------------------------
subroutine h1mg_setup_schwarz_wt(ifsqrt)
  use size_m, only : ndim, ldim, nelv, lelv
  use hsmg, only : mg_schwarz_wt_index, mg_schwarz_wt, mg_lmax, mg_fld
  use hsmg, only : mg_nh, mg_nhz, lmg_swt
  implicit none

  logical :: ifsqrt
        
  integer :: l,i,nl,nlz, itmp

  i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)
  do l=2,mg_lmax

      mg_schwarz_wt_index(l,mg_fld)=i
      nl  = mg_nh(l)
      nlz = mg_nhz(l)
      i   = i+nl*nlz*4*ndim*nelv

      if (i > lmg_swt*4*ldim*lelv) then
          itmp = i/(4*ldim*lelv)
          write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
          call exitt
      endif

      call h1mg_setup_schwarz_wt_1( &
      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)

  enddo

  mg_schwarz_wt_index(l,mg_fld)=i

  return
end subroutine h1mg_setup_schwarz_wt

!----------------------------------------------------------------------
subroutine hsmg_schwarz_wt(e,l)
  use kinds, only : DP
  use input, only : if3d
  use hsmg, only : mg_schwarz_wt, mg_schwarz_wt_index, mg_fld, mg_nh
  implicit none
  real(DP) :: e(*)
  integer :: l
          
#if 0
  if( .NOT. if3d) call hsmg_schwarz_wt2d( &
  e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
#endif
  if(if3d) call hsmg_schwarz_wt3d( &
  e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
  return
end subroutine hsmg_schwarz_wt

!----------------------------------------------------------------------
subroutine hsmg_schwarz_wt3d(e,wt,n)
  use kinds, only : DP
  use size_m, only : nelv
  implicit none

  integer :: n
  real(DP) :: e(n,n,n,nelv)
  real(DP) :: wt(n,n,4,3,nelv)
  integer :: ie,i,j,k
  do ie=1,nelv
      do k=1,n
          do j=1,n
              e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
              e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
              e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
              e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
          enddo
      enddo
      do k=1,n
          do i=3,n-2
              e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
              e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
              e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
              e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
          enddo
      enddo
      do j=3,n-2
          do i=3,n-2
              e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
              e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
              e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
              e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
          enddo
      enddo
  enddo
  return
end subroutine hsmg_schwarz_wt3d

!----------------------------------------------------------------------
subroutine hsmg_coarse_solve(e,r)
  use kinds, only : DP
  use ctimer, only : icalld, ncrsl, tcrsl, ifsync, etime1, dnekclock
  use parallel, only : xxth
  use tstep, only : ifield
  implicit none
  real(DP) :: e(1),r(1)

  if (icalld == 0) then ! timer info
      ncrsl=0
      tcrsl=0.0
  endif
  icalld = 1

  if (ifsync) call nekgsync()

  ncrsl  = ncrsl  + 1
#ifndef NOTIMER
  etime1=dnekclock()
#endif

  call crs_solve(xxth(ifield),e,r)

#ifndef NOTIMER
  tcrsl=tcrsl+dnekclock()-etime1
#endif

  return
end subroutine hsmg_coarse_solve

!----------------------------------------------------------------------
subroutine hsmg_setup_solve
  use size_m, only : nelv, lelv
  use hsmg, only : mg_lmax, mg_nh, mg_nhz, lmg_solve, mg_fld, mg_solve_index
  implicit none
        
  integer :: l,i, itmp

  i = mg_solve_index(mg_lmax+1,mg_fld-1)
  do l=1,mg_lmax
      mg_solve_index(l,mg_fld)=i
      i=i+mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
      if(i > lmg_solve*lelv) then
          itmp = i/lelv
          write(6,*) 'lmg_solve too small',i,itmp,lmg_solve,l
          call exitt
      endif
  enddo
  mg_solve_index(l,mg_fld)=i

  return
end subroutine hsmg_setup_solve

!----------------------------------------------------------------------
subroutine hsmg_setup_mg_nx()
  use size_m, only : lx2, nx1, ly1, lz1, lx1, nid
  use input, only : if3d
  use hsmg, only : mg_lmax, mg_nz, mg_ny, mg_nx
  implicit none

  integer :: p82, mgnx1, mgnx2, i

  integer, save :: mgn2(10) = (/ 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/)

!   if (param(82).eq.0) param(82)=2  ! nek default
!   if (np.eq.1)        param(82)=2  ! single proc. too slow
  p82 = 2                          ! potentially variable nxc
  mg_lmax = 3
!   mg_lmax = 4
  if (lx1 == 4) mg_lmax = 2
!   if (param(79).ne.0) mg_lmax = param(79)
  mgnx1    = p82-1 !1
  mg_nx(1) = mgnx1
  mg_ny(1) = mgnx1
  mg_nz(1) = mgnx1
  if ( .NOT. if3d) mg_nz(1) = 0

  mgnx2 = 2*(lx2/4) + 1
  if (lx1 == 5)  mgnx2 = 3
!   if (lx1.eq.6)  mgnx2 = 3
  if (lx1 <= 10) mgnx2 = mgn2(nx1)
  if (lx1 == 8)  mgnx2 = 4
  if (lx1 == 8)  mgnx2 = 3

!   mgnx2 = min(3,mgnx2)

  mg_nx(2) = mgnx2
  mg_ny(2) = mgnx2
  mg_nz(2) = mgnx2
  if ( .NOT. if3d) mg_nz(2) = 0

  mg_nx(3) = mgnx2+1
  mg_ny(3) = mgnx2+1
  mg_nz(3) = mgnx2+1
  if ( .NOT. if3d) mg_nz(3) = 0

  mg_nx(mg_lmax) = lx1-1
  mg_ny(mg_lmax) = ly1-1
  mg_nz(mg_lmax) = lz1-1

  if (nid == 0) write(*,*) 'mg_nx:',(mg_nx(i),i=1,mg_lmax)
  if (nid == 0) write(*,*) 'mg_ny:',(mg_ny(i),i=1,mg_lmax)
  if (nid == 0) write(*,*) 'mg_nz:',(mg_nz(i),i=1,mg_lmax)

  return
end subroutine hsmg_setup_mg_nx

!----------------------------------------------------------------------
!> \brief initialize index sets
subroutine hsmg_index_0 
  use hsmg, only : mg_rstr_wt_index, mg_solve_index, lmgn, lmgs
  use hsmg, only : mg_fast_s_index, mg_fast_d_index, mg_schwarz_wt_index
  use hsmg, only : mg_mask_index
  implicit none

  integer :: n
  n = lmgn*(lmgs+1)

  mg_mask_index       = 0
  mg_rstr_wt_index    = 0
  mg_solve_index      = 0
  mg_fast_s_index     = 0
  mg_fast_d_index     = 0
  mg_schwarz_wt_index = 0
        
  return
end subroutine hsmg_index_0

!-----------------------------------------------------------------------
subroutine h1mg_mask(w,mask,nel)
  use kinds, only : DP
  implicit none

  real(DP) ::  w(*)
  integer :: mask(*)        ! Pointer to Dirichlet BCs
  integer :: nel

  integer :: e, im
        
  do e=1,nel
      im = mask(e)
      call mg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
  enddo

  return
end subroutine h1mg_mask

!----------------------------------------------------------------------
subroutine mg_mask_e(w,mask) ! Zero out Dirichlet conditions
  use kinds, only : DP
  implicit none

  real(DP) :: w(1)
  integer :: mask(0:1)
  integer :: n, i
  n=mask(0)
  do i=1,n
  !        write(6,*) i,mask(i),n,' MG_MASK'
      w(mask(i)) = 0.
  enddo

  return
end subroutine mg_mask_e

!-----------------------------------------------------------------------
!> \brief compute  v = [A (x) A] u  or  v = [A (x) A (x) A] u
subroutine hsmg_tnsr1(v,nv,nu,A,At)
  use kinds, only : DP
  use input, only : if3d
  implicit none

  integer :: nv,nu
  real(DP) :: v(1),A(1),At(1)

  if ( .NOT. if3d) then
!max        call hsmg_tnsr1_2d(v,nv,nu,A,At)
  else
      call hsmg_tnsr1_3d(v,nv,nu,A,At,At)
  endif
  return
end subroutine hsmg_tnsr1

!-------------------------------------------------------T--------------
!> \brief compute v = [C (x) B (x) A] u
subroutine hsmg_tnsr1_3d(v,nv,nu,A,Bt,Ct)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nelv
  implicit none

  integer :: nv,nu
  real(DP) :: v(*),A(*),Bt(*),Ct(*)

  integer, parameter :: lwk=(lx1+2)*(ly1+2)*(lz1+2)
  real(DP) :: work(0:lwk-1),work2(0:lwk-1)
  integer :: e,e0,ee,es
  integer :: nu3, nv3, iu, iv, i

  e0=1
  es=1
  ee=nelv

  if (nv > nu) then
      e0=nelv
      es=-1
      ee=1
  endif

  nu3 = nu**3
  nv3 = nv**3

  do e=e0,ee,es
      iu = 1 + (e-1)*nu3
      iv = 1 + (e-1)*nv3
      call mxm(A,nv,v(iu),nu,work,nu*nu)
      do i=0,nu-1
          call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
  enddo

  return
end subroutine hsmg_tnsr1_3d

!------------------------------------------   T  -----------------------
!> \brief r =J r,   l is coarse level
subroutine h1mg_rstr(r,l,ifdssum) 
  use kinds, only : DP
  use hsmg, only : mg_rstr_wt, mg_rstr_wt_index, mg_nh, mg_nhz, mg_jht
  use hsmg, only : mg_jh, mg_fld
  implicit none
  logical :: ifdssum
  real(DP) :: r(1)
  integer :: l

  call hsmg_do_wt(r,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld)) &
  ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

  call hsmg_tnsr1(r,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))

  if (ifdssum) call hsmg_dssum(r,l)

  return
end subroutine h1mg_rstr

!-----------------------------------------------------------------------
subroutine h1mg_setup_mg_nx()
  use size_m, only : lx2, nx1, lx1, ly1, lz1, ldimt1, nid
  use input, only : if3d
  use hsmg, only : mg_h1_n, mg_lmax, mg_h1_lmax, mg_nz, mg_ny, mg_nx
  use tstep, only : nelfld
  implicit none

  integer :: mgn2(10) = (/ 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/)
  integer :: p82, mgnx1, mgnx2, i, ifld, l

!   if (param(82).eq.0) param(82)=2  ! nek default
!   if (np.eq.1)        param(82)=2  ! single proc. too slow
  p82 = 2                          ! potentially variable nxc
  mg_h1_lmax = 3
!   mg_h1_lmax = 4
  if (lx1 == 4) mg_h1_lmax = 2
!   if (param(79).ne.0) mg_h1_lmax = param(79)
  mgnx1    = p82-1 !1
  mg_nx(1) = mgnx1
  mg_ny(1) = mgnx1
  mg_nz(1) = mgnx1
  if ( .NOT. if3d) mg_nz(1) = 0

  mgnx2 = 2*(lx2/4) + 1
  if (lx1 == 5)  mgnx2 = 3
!   if (lx1.eq.6)  mgnx2 = 3
  if (lx1 <= 10) mgnx2 = mgn2(nx1)
  if (lx1 == 8)  mgnx2 = 4
  if (lx1 == 8)  mgnx2 = 3

  mgnx2 = min(3,mgnx2)  ! This choice seems best (9/24/12)

  mg_nx(2) = mgnx2
  mg_ny(2) = mgnx2
  mg_nz(2) = mgnx2
  if ( .NOT. if3d) mg_nz(2) = 0

  mg_nx(3) = mgnx2+1
  mg_ny(3) = mgnx2+1
  mg_nz(3) = mgnx2+1
  if ( .NOT. if3d) mg_nz(3) = 0

  mg_nx(mg_h1_lmax) = lx1-1
  mg_ny(mg_h1_lmax) = ly1-1
  mg_nz(mg_h1_lmax) = lz1-1

  if (nid == 0) write(*,*) 'h1_mg_nx:',(mg_nx(i),i=1,mg_h1_lmax)
  if (nid == 0) write(*,*) 'h1_mg_ny:',(mg_ny(i),i=1,mg_h1_lmax)
  if (nid == 0) write(*,*) 'h1_mg_nz:',(mg_nz(i),i=1,mg_h1_lmax)

  do ifld=1,ldimt1
      do l=1,mg_lmax
          mg_h1_n(l,ifld)=(mg_nx(l)+1) &
          *(mg_ny(l)+1) &
          *(mg_nz(l)+1)*nelfld(ifld)
      enddo
  enddo
        
  return
end subroutine h1mg_setup_mg_nx

!----------------------------------------------------------------------
!> \brief SEM hat matrices for each level
subroutine h1mg_setup_semhat 
  use hsmg, only : mg_h1_lmax, mg_nx, mg_ah, mg_bh, mg_dh, mg_dht, mg_zh
  use hsmg, only : mg_nh, mg_nhz, mg_nz
  use semhat, only : ah, bh, ch, dh, zh, dph, jph, bgl, zgl, dgl, jgl, wh
  implicit none

  integer :: l, n

  do l=1,mg_h1_lmax
      n = mg_nx(l)     ! polynomial order
      call generate_semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
      call copy(mg_ah(1,l),ah,(n+1)*(n+1))
      call copy(mg_bh(1,l),bh,n+1)
      call copy(mg_dh(1,l),dh,(n+1)*(n+1))
      call transpose(mg_dht(1,l),n+1,dh,n+1)
      call copy(mg_zh(1,l),zh,n+1)

      mg_nh(l)=n+1
      mg_nhz(l)=mg_nz(l)+1

  enddo
end subroutine h1mg_setup_semhat

!----------------------------------------------------------------------
subroutine h1mg_setup_dssum
  use kinds, only : i8
  use size_m, only : lx1, ly1, lz1, lelt, nelv, ndim
  use input, only : if3d
  use hsmg, only : mg_lmax, mg_nh, mg_nhz, mg_fld
  use hsmg, only : mg_gsh_handle, mg_gsh_schwarz_handle
  use mesh, only : vertex
  use parallel, only : nelgv
  implicit none

  integer, parameter :: lxyz=(lx1+2)*(ly1+2)*(lz1+2)

  integer(i8), allocatable :: glo_num(:)
  integer :: nx,ny,nz, ncrnr
  integer :: l
  
 
!     set up direct stiffness summation for each level
  ncrnr = 2**ndim
  call get_vert()

  allocate(glo_num(lxyz*lelt)) 
  do l=1,mg_lmax 
      nx=mg_nh(l)
      ny=mg_nh(l)
      nz=mg_nhz(l)
      call setupds(mg_gsh_handle(l,mg_fld),nx,ny,nz &
      ,nelv,nelgv,vertex,glo_num)
      nx=nx+2
      ny=ny+2
      nz=nz+2
      if( .NOT. if3d) nz=1
      call setupds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz &
      ,nelv,nelgv,vertex,glo_num)
  enddo
end subroutine h1mg_setup_dssum

!----------------------------------------------------------------------
subroutine mg_set_msk(p_msk ,l0)
  use kinds, only : DP
  use hsmg, only : mg_h1_lmax, p_mg_msk, mg_nh, mg_nhz, mg_imask
  use hsmg, only : mg_h1_n, mg_fld
  use tstep, only : ifield, nelfld
  implicit none

  integer :: p_msk, l0

  real(DP), allocatable :: work(:) 
  integer :: l, n, nx, ny, nz, nm 
  l                  = mg_h1_lmax
  p_mg_msk(l,mg_fld) = 0
  n                  = mg_h1_n(l,mg_fld)

  allocate(work(maxval(mg_nh(1:mg_h1_lmax))*maxval(mg_nh(1:mg_h1_lmax))*maxval(mg_nhz(1:mg_h1_lmax))*nelfld(ifield)))

  do l=mg_h1_lmax,1,-1
      nx = mg_nh  (l)
      ny = mg_nh  (l)
      nz = mg_nhz (l)

      p_msk = p_mg_msk(l,mg_fld)

      call h1mg_setup_mask &
      (mg_imask(p_msk),nm,nx,ny,nz,nelfld(ifield),l,work)

      if (l > 1) p_mg_msk(l-1,mg_fld)=p_mg_msk(l,mg_fld)+nm

  enddo

  p_msk = p_mg_msk(l0,mg_fld)

  return
end subroutine mg_set_msk

!----------------------------------------------------------------------
subroutine h1mg_setup_mask(mask,nm,nx,ny,nz,nel,l,w)
  use kinds, only : DP
  use size_m, only : nid
  use input, only : if3d
  implicit none

  integer :: mask(*)        ! Pointer to Dirichlet BCs
  integer :: nx,ny,nz,nel,l
  real(DP) :: w(nx,ny,nz,nel)
        
  integer :: e,count,ptr
  integer :: lbr,rbr,lbs,rbs,lbt,rbt,two
  integer :: nxyz, n, ierrmx, i, ierr, nm
  integer, external :: iglmax
  real(DP) :: zero
  real(DP) :: w_flat(nx*ny*nz)

  zero = 0
  nxyz = nx*ny*nz
  n    = nx*ny*nz*nel

  w = 1._dp   ! Init everything to 1

  ierrmx = 0       ! BC verification
  two    = 2
  do e=1,nel       ! Set dirichlet nodes to zero

      call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,two,ierr)
  !        write(6,6) e,lbr,rbr,lbs,rbs,ierr,nx
  !   6    format(i5,2x,4i3,2x,i2,3x,i5,'  lbr,rbr,lbs')

      if (lbr == 1) call facev(w,e,4,zero,nx,ny,nz)
      if (rbr == 1) call facev(w,e,2,zero,nx,ny,nz)
      if (lbs == 1) call facev(w,e,1,zero,nx,ny,nz)
      if (rbs == 1) call facev(w,e,3,zero,nx,ny,nz)
      if (if3d) then
          if (lbt == 1) call facev(w,e,5,zero,nx,ny,nz)
          if (rbt == 1) call facev(w,e,6,zero,nx,ny,nz)
      endif
      ierrmx = max(ierrmx,ierr)
  enddo

  call hsmg_dsprod(w,l)    ! direct stiffness multiply


!   Prototypical mask layout, nel=5:

!  e=1 ...             10
!    1  2  3  4  5 ... 10 | 11 12 13 14 | 15 | 16 |
!   11 15 16 ...          |  3 p1 p2 p3 |  0 |  0 | ...
!                            ^
!                            |
!                            |_count for e=1


  nm  = 1                  ! store mask
  do e=1,nel

      mask(e) = nel+nm
      count   = 0          ! # Dirchlet points on element e
      ptr     = mask(e)

      w_flat = reshape(w(:,:,:,e), (/nxyz/))
      do i=1,nxyz
          if (w_flat(i) == 0) then
              nm    = nm   +1
              count = count+1
              ptr   = ptr  +1
              mask(ptr) = i + nxyz*(e-1)   ! where I mask on element e
          endif
      enddo


      ptr       = mask(e)
      mask(ptr) = count

      nm        = nm+1     ! bump pointer to hold next count

  enddo

  nm = nel + nm-1 ! Return total number of mask pointers/counters

  ierrmx = iglmax(ierrmx,1)
  if (ierrmx > 0) then
      if (ierr > 0) write(6,*) nid,ierr,' BC FAIL h1'
      call exitti('D INVALID BC FOUND in h1mg_setup_mask$',ierrmx)
  endif

  return
end subroutine h1mg_setup_mask

!----------------------------------------------------------------------
subroutine gxfer_e (g,ng,e)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use geom, only : g1m1, g2m1, g3m1, g4m1, g5m1, g6m1
  use input, only : if3d
  implicit none

  integer :: ng, e
  real(DP) :: g(ng,1)
  integer :: nxyz, i

  nxyz = nx1*ny1*nz1

!   ifdfrm(e) = .true.  ! TOO LATE

  if (if3d) then
      do i=1,nxyz
          g(1,i) = g1m1(i,1,1,e)
          g(2,i) = g2m1(i,1,1,e)
          g(3,i) = g3m1(i,1,1,e)
          g(4,i) = g4m1(i,1,1,e)
          g(5,i) = g5m1(i,1,1,e)
          g(6,i) = g6m1(i,1,1,e)
      enddo
  else
      do i=1,nxyz
          g(1,i) = g1m1(i,1,1,e)
          g(2,i) = g2m1(i,1,1,e)
          g(3,i) = g4m1(i,1,1,e)
      enddo
  endif

  return
end subroutine gxfer_e

!-----------------------------------------------------------------------
subroutine h1mg_setup_schwarz_wt_2(wt,ie,n,work,ifsqrt)
  use kinds, only : DP
  use size_m, only : ndim
  implicit none
  real(DP) :: wt(1),work(1)
  integer :: ie, n
  logical :: ifsqrt

!max    if (ndim == 2) call h1mg_setup_schwarz_wt2d_2(wt,ie,n,work,ifsqrt)
  if (ndim == 3) call h1mg_setup_schwarz_wt3d_2(wt,ie,n,work,ifsqrt)

  return
end subroutine h1mg_setup_schwarz_wt_2

!----------------------------------------------------------------------
subroutine h1mg_setup_schwarz_wt3d_2(wt,ie,n,work,ifsqrt)
  use kinds, only : DP
  use size_m, only : nelv
  implicit none

  logical :: ifsqrt
  integer :: ie, n
  real(DP) :: wt(n,n,4,3,nelv)
  real(DP) :: work(n,n,n)
        
  integer :: i,j,k, ii, ierr

  ierr = 0
  do k=1,n
      do j=1,n
          wt(j,k,1,1,ie)=1.0/work(1,j,k)
          wt(j,k,2,1,ie)=1.0/work(2,j,k)
          wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
          wt(j,k,4,1,ie)=1.0/work(n,j,k)
      enddo
  enddo
  do k=1,n
      do i=1,n
          wt(i,k,1,2,ie)=1.0/work(i,1,k)
          wt(i,k,2,2,ie)=1.0/work(i,2,k)
          wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
          wt(i,k,4,2,ie)=1.0/work(i,n,k)
      enddo
  enddo
  do j=1,n
      do i=1,n
          wt(i,j,1,3,ie)=1.0/work(i,j,1)
          wt(i,j,2,3,ie)=1.0/work(i,j,2)
          wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
          wt(i,j,4,3,ie)=1.0/work(i,j,n)
      enddo
  enddo
  if(ifsqrt) then
      do ii=1,3
          do k=1,4
              do j=1,4
                  do i=1,n
                      wt(i,j,k,ii,ie)=sqrt(wt(i,j,k,ii,ie))
                  enddo
              enddo
          enddo
      enddo
  endif

  return
end subroutine h1mg_setup_schwarz_wt3d_2

!----------------------------------------------------------------------
subroutine h1mg_setup_schwarz_wt_1(wt,l,ifsqrt)
  use kinds, only : DP
  use input, only : if3d
  use hsmg, only : mg_h1_n, p_mg_msk, mg_nh, mg_fld
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: wt(1)
  integer :: l
  logical :: ifsqrt

  integer :: enx,eny,enz,pm, n, ns, i, nx, ny, nz, nxyz, k, ie
  real(DP) :: zero, one, onem
  real(DP), allocatable :: work(:)

  zero =  0
  one  =  1
  onem = -1

  n  = mg_h1_n (l,mg_fld)
  pm = p_mg_msk(l,mg_fld)

  enx=mg_nh(l)+2
  eny=mg_nh(l)+2
  enz=mg_nh(l)+2
  if( .NOT. if3d) enz=1
  ns = enx*eny*enz*nelfld(ifield)
  i  = ns+1

  allocate(work(2*ns))
  work(1:ns) = 0._dp; work(i:i+ns-1) = 1._dp
   
!   Sum overlap region (border excluded)
  call hsmg_extrude(work,0,zero,work(i),0,one ,enx,eny,enz)
  call hsmg_schwarz_dssum(work(i),l)
  call hsmg_extrude(work(i),0,one ,work,0,onem,enx,eny,enz)
  call hsmg_extrude(work(i),2,one, work(i),0,one,enx,eny,enz)

  if( .NOT. if3d) then ! Go back to regular size array
!max        call hsmg_schwarz_toreg2d(mg_work,mg_work(i),mg_nh(l))
  else
      call hsmg_schwarz_toreg3d(work,work(i),mg_nh(l))
  endif

  call hsmg_dssum(work,l)                           ! sum border nodes


  nx = mg_nh(l)
  ny = mg_nh(l)
  nz = mg_nh(l)
  if ( .NOT. if3d) nz=1
  nxyz = nx*ny*nz
  k    = 1
  do ie=1,nelfld(ifield)
  !        call outmat(mg_work(k),nx,ny,'NEW WT',ie)
      call h1mg_setup_schwarz_wt_2(wt,ie,nx,work(k),ifsqrt)
      k = k+nxyz
  enddo
!   stop

  return
end subroutine h1mg_setup_schwarz_wt_1

!----------------------------------------------------------------------
end module hsmg_routines
