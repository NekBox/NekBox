!-----------------------------------------------------------------------
!> \brief THE ROUTINES BELOW ARE THE NEW Helmholtz projectors
!-----------------------------------------------------------------------

module helmholtz
  use kinds, only : DP
  implicit none

  private

  public :: hsolve, approx_space


  type approx_space
    real(DP), allocatable :: projectors(:,:)
    integer :: n_max
    integer :: n_sav
    integer :: next
    real(DP), allocatable :: A_red(:,:) 
  end type approx_space

contains

!> \brief Orthogonalize the rhs wrt previous rhs's for which we already
!! know the soln.
subroutine projh(r,h1,h2,bi,vml,vmk, apx, wl,ws,name4)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use geom, only : voltm1, volvm1
  use tstep, only : istep, ifield, nelfld
  use parallel, only : nid
  implicit none

  real(DP), intent(inout) :: r(*) !>!< residual
  real(DP), intent(in)    :: h1(*) !>!< coefficient of A (stiffness)
  real(DP), intent(in)    :: h2(*) !>!< coefficient of M (mass)
  real(DP), intent(in)    :: vml(*) !>!< multiplicity array
  real(DP), intent(in)    :: vmk(*) !>!< mask array
  real(DP), intent(in)    :: bi(*) !>!< inverse mass matrix
  real(DP), intent(out)   :: wl(*) !>!< large work array (size lx1*ly1*lz1*nelv)
  real(DP), intent(out)   :: ws(*) !>!< small work array (size 2*max vecs)
  type(approx_space), intent(inout) :: apx
  character(4) :: name4

  integer :: nel, ntot, i, j, n10
  real(DP) :: vol, alpha1, alpha2, ratio
  real(DP), external :: glsc23, vlsc3
  real(DP), external :: glsc2, glsc3
  real(DP), allocatable :: evecs(:,:), work(:), ev(:)
  integer :: lwork, ierr
  real(DP), parameter :: one = 1._dp, zero = 0._dp

  if (apx%n_sav == 0) then
    apx%projectors(:,0) = 0._dp
    return
  endif

  nel =nelfld(ifield)
  ntot=nx1*ny1*nz1*nel

  vol = voltm1
  if (nel == nelv) vol = volvm1

!   Diag to see how much reduction in the residual is attained.

  alpha1 = glsc23(r,bi,vml,ntot)
  if (alpha1 > 0) alpha1 = sqrt(alpha1/vol)

!   Update approximation space if dt has changed
!  call updrhsh(apx,h1,h2,vml,vmk,ws,name4)


  allocate(evecs(apx%n_sav, apx%n_sav))
  evecs = apx%A_red(1:apx%n_sav,1:apx%n_sav)

  lwork = 10 * apx%n_sav + 100
  allocate(work(lwork), ev(apx%n_sav))
  call dsyev('V', 'U', apx%n_sav, &
             evecs, apx%n_sav, &
             ev, &
             work, lwork, ierr) 
  if (nid == 0 .and. ierr /= 0) write(*,*) "DSYEV failed", ierr
  wl(1:ntot) = r(1:ntot) * vml(1:ntot)

#if 0
  do i = 1, apx%n_sav
    !ws(i) = glsc3(wl, approx(:,i,1), vml, ntot)
    ws(i) = glsc2(wl, apx%projectors(:,i), ntot)
  enddo
#else
  call dgemv('T', ntot, apx%n_sav, &
             one,  apx%projectors(:,1:apx%n_sav), ntot, &
                   wl, 1, &
             zero, ws, 1)
  call gop(ws, ws(1+apx%n_sav), '+  ', apx%n_sav)
#endif
 
  do i = 1, apx%n_sav
    ev(i) = sum(evecs(:,i) * ws(1:apx%n_sav)) / ev(i)
  enddo

  do i = 1, apx%n_sav
    ws(i) = sum(evecs(i,:) * ev(:))
  enddo

  call dgemv('N', ntot, apx%n_sav, &
             one,  apx%projectors(:,1:apx%n_sav), ntot, &
                   ws, 1, &
             zero, apx%projectors(:,0), 1)

  call axhelm  (wl,apx%projectors(:,0),h1,h2,1,1)
  wl(1:ntot) = wl(1:ntot) * vmk(1:ntot)
  call dssum   (wl)
  r(1:ntot) = r(1:ntot) - wl(1:ntot)

!...............................................................
! Diag.
  alpha2 = glsc23(r,bi,vml,ntot)
  if (alpha2 > 0) alpha2 = sqrt(alpha2/vol)
  ratio  = alpha1/alpha2
  n10=min(10,apx%n_sav)

  if (nid == 0) write(6,10) istep,name4,alpha1,alpha2,ratio,apx%n_sav
  10 format(4X,I7,4x,a4,' alph1n',1p3e12.4,i6)

  if (nid == 0) write(6,11) istep,name4,apx%n_sav,(ws(i),i=1,n10)
  11 format(4X,I7,4x,a4,' halpha',i6,10(1p10e12.4,/,17x))

  return
end subroutine projh

!-----------------------------------------------------------------------
!> \brief Reconstruct the solution to the original problem by adding back
!!     the previous solutions
subroutine gensh(v1,h1,h2,vml,vmk,apx,ws,name4)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use size_m, only : lx1, ly1, lz1
  use parallel, only : nid
  use mesh, only : niterhm
  use tstep, only : nelfld, ifield
  implicit none

  REAL(DP), intent(inout) :: V1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: H1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: H2   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: vmk  (*)
  REAL(DP), intent(in)    :: vml  (*)
  type(approx_space), intent(inout) :: apx
  real(DP), intent(out)   :: ws(:) !>!< workspace?

  character(4), intent(in) :: name4

  real(DP) :: alpha, eps
  real(DP), external :: glsc2, glsc3
  integer :: ntot, ierr, i

  ntot=nx1*ny1*nz1*nelfld(ifield)

!   Reconstruct solution and save current du

  if (apx%n_sav < apx%n_max) then
  
      if (niterhm > 0) then      ! new vector not in space
          apx%n_sav = apx%n_sav+1
          v1(:,:,:,1:nelfld(ifield)) = v1(:,:,:,1:nelfld(ifield))  &
                                     + reshape(apx%projectors(:,0), (/lx1,ly1,lz1,nelfld(ifield)/))
          call copy(apx%projectors(:,apx%n_sav),v1,ntot)
          call hconj(apx,apx%n_sav,h1,h2,vml,vmk,ws,name4,ierr)
          apx%next = mod(apx%n_sav, apx%n_max) + 1
      else
        if (nid == 0) write(*,*) "Freak out!" 
      endif
  else
      v1(:,:,:,1:nelfld(ifield)) = v1(:,:,:,1:nelfld(ifield))  &
                                 + reshape(apx%projectors(:,0), (/lx1,ly1,lz1,nelfld(ifield)/))
      call copy(apx%projectors(:,apx%next),v1,ntot)
      call hconj(apx,apx%next,h1,h2,vml,vmk,ws,name4,ierr)
      apx%next = mod(apx%next, apx%n_max) + 1
  endif

  return
end subroutine gensh

!-----------------------------------------------------------------------
!> \brief Orthonormalize the kth vector against vector set
subroutine hconj(apx,k,h1,h2,vml,vmk,ws,name4,ierr)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nid
  use parallel, only : wdsize
  use tstep, only : istep, ifield, nelfld
  implicit none

  integer, intent(in) :: k
  type(approx_space), intent(inout) :: apx
  real(DP), intent(in) :: h1(*),h2(*),vml(*),vmk(*)
  real(DP), intent(out) :: ws(*)
  character(4) :: name4
  integer, intent(out) :: ierr

  integer :: i, ntot, km1 , nel
  real(DP) :: alpha, ratio, eps, alpham, alph1
  real(DP), external :: glsc2, ddot
  real(DP), parameter :: one = 1._dp, zero = 0._dp

  ierr=0
  nel = nelfld(ifield)
  ntot=nx1*ny1*nz1*nel

  call axhelm  (apx%projectors(:,0),apx%projectors(:,k),h1,h2,1,1)
  apx%projectors(:,0) = apx%projectors(:,0) * vmk(1:ntot)
  call dssum   (apx%projectors(:,0))
  apx%projectors(:,0) = apx%projectors(:,0) * vml(1:ntot)

  call dgemv('T', ntot, apx%n_sav, &
             one,  apx%projectors(:,1:apx%n_sav), ntot, &
                   apx%projectors(:,0), 1, &
             zero, apx%A_red(:,k), 1)
  call gop(apx%A_red(:,k), ws, '+  ', apx%n_sav)

  do i = 1, apx%n_sav
    apx%A_red(k,i) = apx%A_red(i,k)
  enddo

  return
end subroutine hconj

!-----------------------------------------------------------------------
!> \brief Reorthogonalize approx if dt has changed
subroutine updrhsh(apx,h1,h2,vml,vmk,ws,name4)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1
  use input, only : ifvarp, iflomach
  use tstep, only : dt, ifield, nelfld
  implicit none

  integer, parameter :: lt=lx1*ly1*lz1*lelt
  type(approx_space), intent(inout) :: apx
  real(DP), intent(in)    :: h1(1),h2(1),vml(1),vmk(1)
  real(DP), intent(out) :: ws(1)
  character(4) :: name4

  logical :: ifupdate
  logical, save :: ifnewdt = .false.
  integer :: n_sav, l, k, ntot, ierr

  character(4), save :: name_old = 'DMIR'

  real(DP), save :: dtold = 0.0

!   First, we have to decide if the dt has changed.

  ifupdate = .FALSE. 
  if (dt /= dtold) then
      dtold    = dt
      name_old = name4
      ifnewdt  = .TRUE. 
      ifupdate = .TRUE. 
  elseif (ifnewdt) then
      if (name4 == name_old) then
          ifnewdt = .FALSE. 
      else
          ifupdate = .TRUE. 
      endif
  endif
  if (ifvarp(ifield)) ifupdate = .TRUE. 
  if (iflomach)       ifupdate = .TRUE. 

  if (ifupdate) then    ! reorthogonalize
      n_sav = apx%n_sav 
      l     = 1
      do k=1,n_sav
#if 0
      !           Orthogonalize kth vector against {v_1,...,v_k-1}
          if (k /= l) then
              ntot = nx1*ny1*nz1*nelfld(ifield)
              call copy(approx(1,l),approx(1,k),ntot)
          endif
          call hconj(approx,l,h1,h2,vml,vmk,ws,name4,ierr)
          if (ierr == 0) l=l+1
#else
      apx%n_sav = k 
      call hconj(apx,apx%n_sav,h1,h2,vml,vmk,ws,name4,ierr)
#endif 
      enddo
  endif

  return
end subroutine updrhsh

!-----------------------------------------------------------------------
subroutine hmhzpf(name,u,r,h1,h2,mask,mult,imesh,tli,maxit,isd,bi)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1
  use size_m, only : nx1, ny1, nz1, nelv, nelt, ndim
  use ctimer, only : etime1, dnekclock, thmhz
  use fdmh1, only : kfldfdm
  use input, only : param
  implicit none

  CHARACTER(4) ::    NAME
  REAL(DP), intent(out) :: U    (LX1,LY1,LZ1,1) !>!< solution vector
  REAL(DP), intent(in)  :: R    (LX1,LY1,LZ1,1) !>!< right hand side
  REAL(DP), intent(in)  :: H1   (LX1,LY1,LZ1,1) !>!< coefficient of A (stiffness)
  REAL(DP), intent(in)  :: H2   (LX1,LY1,LZ1,1) !>!< coefficient of M (mass)
  REAL(DP), intent(in)  :: MASK (LX1,LY1,LZ1,1) !>!< mask array
  REAL(DP), intent(in)  :: MULT (LX1,LY1,LZ1,1) !>!< multiplicity array
  REAL(DP), intent(in)  :: bi   (LX1,LY1,LZ1,1) !>!< inverse of mass matrix
  real(DP) :: tli
  integer :: imesh, maxit, isd

  integer :: ntot
  real(DP) :: tol

  etime1=dnekclock()

  IF (IMESH == 1) NTOT = NX1*NY1*NZ1*NELV
  IF (IMESH == 2) NTOT = NX1*NY1*NZ1*NELT

  tol = tli
  if (param(22) /= 0) tol = abs(param(22))
  CALL CHKTCG1 (TOL,R,H1,H2,MASK,MULT,IMESH,ISD)


!   Set flags for overlapping Schwarz preconditioner (pff 11/12/98)

  kfldfdm = -1
!   if (name.eq.'TEMP') kfldfdm =  0
!   if (name.eq.'VELX') kfldfdm =  1
!   if (name.eq.'VELY') kfldfdm =  2
!   if (name.eq.'VELZ') kfldfdm =  3
  if (name == 'PRES') kfldfdm =  ndim+1

  call cggo &
  (u,r,h1,h2,mask,mult,imesh,tol,maxit,isd,bi,name)
  thmhz=thmhz+(dnekclock()-etime1)


  return
end subroutine hmhzpf

!-----------------------------------------------------------------------
!> \brief Either std. Helmholtz solve, or a projection + Helmholtz solve
subroutine hsolve(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd &
    ,apx,bi)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, lelv
  use parallel, only : nid
  use input, only : ifflow, param
  use string, only : capit
  use tstep, only : ifield, nelfld, istep
  use geom, only : binvm1
  use poisson, only : spectral_solve
  implicit none

  CHARACTER(4), intent(in) :: NAME !>!< name of field we're solving for
  REAL(DP), intent(out)   :: U    (LX1,LY1,LZ1,lelv) !>!< solution vector
  REAL(DP), intent(inout) :: R    (LX1,LY1,LZ1,lelv) !>!< right hand side
  REAL(DP), intent(in)    :: H1   (LX1,LY1,LZ1,lelv) !>!< coefficient of A (stiffness)
  REAL(DP), intent(in)    :: H2   (LX1,LY1,LZ1,lelv) !>!< coefficient of M (mass)
  REAL(DP), intent(in)    :: vmk  (LX1,LY1,LZ1,lelv) !>!< mask array
  REAL(DP), intent(in)    :: vml  (LX1,LY1,LZ1,lelv) !>!< multiplicity array
  integer,  intent(in)    :: imsh                 !>!< imesh?
  real(DP), intent(in)    :: tol                  !>!< residual tolerance
  integer,  intent(in)    :: maxit                !>!< maximum number of iterations
  integer,  intent(in)    :: isd                  !>!< something to do with axi-symmetric
  type(approx_space), intent(inout) :: apx
  REAL(DP), intent(in)    :: bi   (LX1,LY1,LZ1,*) !>!< inverse of mass matrix

  real(DP), allocatable :: w1(:)
  real(DP), allocatable :: w2(:)

  logical :: ifstdh, spectral_h
  character(4) ::  cname
  integer :: n, nel


  call chcopy(cname,name,4)
  call capit (cname,4)

  ! figure out if we're projecting or not
  ifstdh = .TRUE. 
  if (cname == 'PRES') then
    if (param(95) /= 0 .AND. istep > param(95) .and. param(93) > 0) then
      ifstdh = .FALSE.
    endif
  elseif (cname == 'VELX' .or. cname == 'VELY' .or. cname == 'VELZ') then
    if (param(94) /= 0 .AND. istep > param(94) .and. param(92) > 0) then
      ifstdh = .FALSE. 
    endif
  endif


  if (ifstdh) then

    call hmholtz(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)

  else

      nel = nelfld(ifield)
      n = nx1*ny1*nz1*nel

      call dssum  (r)
      r(:,:,:,1:nel) = r(:,:,:,1:nel) * vmk(:,:,:,1:nel)

      allocate(w2(2+2*apx%n_max))
      allocate(w1(lx1*ly1*lz1*lelt))
      call projh  (r,h1,h2,bi,vml,vmk,apx,w1,w2,name)
      deallocate(w1)

      call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)
      call gensh  (u,h1,h2,vml,vmk,apx,w2,name)

  endif


  return
end subroutine hsolve
!-----------------------------------------------------------------------

end module helmholtz
