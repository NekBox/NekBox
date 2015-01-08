!-----------------------------------------------------------------------
!> \brief THE ROUTINES BELOW ARE THE NEW Helmholtz projectors
!-----------------------------------------------------------------------

module helmholtz
  implicit none

  private

  public :: hsolve

contains

!> \brief Orthogonalize the rhs wrt previous rhs's for which we already
!! know the soln.
subroutine projh(r,h1,h2,bi,vml,vmk,approx,napprox,wl,ws,name4)
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
  real(DP), intent(inout) :: approx(:,0:,:) !>!< approximation space
  integer, intent(inout)  :: napprox(2) !>!< (/ max vecs, current number of vecs /)
  character(4) :: name4

  integer :: n_max, n_sav, nel, ntot, i, n10
  real(DP) :: vol, alpha1, alpha2, ratio
  real(DP), external :: glsc23, vlsc3
  real(DP), external :: glsc2, glsc3

  n_max = napprox(1)
  n_sav = napprox(2)
  if (n_sav == 0) return

  nel =nelfld(ifield)
  ntot=nx1*ny1*nz1*nel

  vol = voltm1
  if (nel == nelv) vol = volvm1

!   Diag to see how much reduction in the residual is attained.

  alpha1 = glsc23(r,bi,vml,ntot)
  if (alpha1 > 0) alpha1 = sqrt(alpha1/vol)

!   Update approximation space if dt has changed
  call updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)


!   Perform Gram-Schmidt for previous soln's

#if 0
  wl(1:ntot) = r(1:ntot)
  do i=1,n_sav
      ws(i) = vlsc3(wl,approx(:,i,2),vml,ntot)
      call gop(ws(i),ws(n_sav+1),'+  ',1)
      wl(1:ntot) = wl(1:ntot) - ws(i) * approx(:,i,2)
  enddo

  approx(:,0,1) = approx(:,1,1) * ws(1)
  do i=2,n_sav
    approx(:,0,1) = approx(:,0,1) + approx(:,i,1) * ws(i)
  enddo
  approx(:,0,2) = r(1:ntot) - wl(1:ntot)
  r(1:ntot) = wl(1:ntot)
#else
  wl(1:ntot) = r(1:ntot) * vml(1:ntot)
  do i=n_sav,1,-1
      !ws(i) = glsc3(wl,approx(:,i,2),vml,ntot)
      ws(i) = glsc2(wl,approx(:,i,2),ntot)
      wl(1:ntot) = wl(1:ntot) - ws(i) * approx(:,i,2)
  enddo

  approx(:,0,1) = 0._dp
  do i=n_sav,1,-1
    approx(:,0,1) = approx(:,0,1) + approx(:,i,1) * ws(i)
    call axhelm  (approx(:,0,2),approx(:,0,1),h1,h2,1,1)
    approx(:,0,2) = approx(:,0,2) * vmk(1:ntot)
    call dssum   (approx(:,0,2))
    wl(1:ntot) = r(1:ntot) - approx(1:ntot,0,2) 

  alpha2 = glsc23(wl,bi,vml,ntot)
  if (alpha2 > 0) alpha2 = sqrt(alpha2/vol)
  ratio  = alpha1/alpha2
  if (nid == 0) write(*,*) "RATIO: ", i, ratio 

  enddo

  call axhelm  (approx(:,0,2),approx(:,0,1),h1,h2,1,1)
  approx(:,0,2) = approx(:,0,2) * vmk(1:ntot)
  call dssum   (approx(:,0,2))
  r(1:ntot) = r(1:ntot) - approx(1:ntot,0,2)
#endif
  !r(1:ntot) = wl(1:ntot) * (1._dp / vml(1:ntot))


!...............................................................
! Diag.
  alpha2 = glsc23(r,bi,vml,ntot)
  if (alpha2 > 0) alpha2 = sqrt(alpha2/vol)
  ratio  = alpha1/alpha2
  n10=min(10,n_sav)

  if (nid == 0) write(6,10) istep,name4,alpha1,alpha2,ratio,n_sav
  10 format(4X,I7,4x,a4,' alph1n',1p3e12.4,i6)

  if (nid == 0) write(6,11) istep,name4,n_sav,(ws(i),i=1,n10)
  11 format(4X,I7,4x,a4,' halpha',i6,10(1p10e12.4,/,17x))

  return
end subroutine projh

!-----------------------------------------------------------------------
!> \brief Reconstruct the solution to the original problem by adding back
!!     the previous solutions
subroutine gensh(v1,h1,h2,vml,vmk,approx,napprox,ws,name4)
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
  real(DP), intent(out)   :: ws(:) !>!< workspace?

  real(DP), intent(inout) :: approx(:,0:,:)
  integer, intent(inout) :: napprox(2)
  character(4), intent(in) :: name4

  real(DP) :: alpha, eps
  real(DP), external :: glsc2, glsc3
  integer :: n_max, n_sav, ntot, ierr

  eps = 1.e-30_dp

  n_max = napprox(1)
  n_sav = napprox(2)
  ntot=nx1*ny1*nz1*nelfld(ifield)

!   Reconstruct solution and save current du

  if (n_sav < n_max) then
  
      if (niterhm > 0) then      ! new vector not in space
          n_sav = n_sav+1
          v1(:,:,:,1:nelfld(ifield)) = v1(:,:,:,1:nelfld(ifield))  &
                                     + reshape(approx(:,0,1), (/lx1,ly1,lz1,nelfld(ifield)/))
          call copy(approx(:,n_sav,1),v1,ntot)

          call axhelm        (approx(:,n_sav,2),approx(:,n_sav,1),h1,h2,1,1)
          approx(:,n_sav,2) = approx(:,n_sav,2) * vmk(1:ntot)
          call dssum         (approx(:,n_sav,2))
          approx(:,n_sav,2) = approx(:,n_sav,2) * vml(1:ntot)
          !alpha = glsc3(approx(:,n_sav,2),approx(:,n_sav,2),vml,ntot)
          alpha = glsc2(approx(:,n_sav,2),approx(:,n_sav,2),ntot)
          if (alpha > eps) approx(:,n_sav,:) = approx(:,n_sav,:) *( 1._dp/ sqrt(alpha))

      else
          n_sav = n_sav+1
          v1(:,:,:,1:nelfld(ifield)) = v1(:,:,:,1:nelfld(ifield))  &
                                     + reshape(approx(:,0,1), (/lx1,ly1,lz1,nelfld(ifield)/))
          call copy(approx(:,n_sav,1),v1,ntot)
          call axhelm        (approx(:,n_sav,2),approx(:,n_sav,1),h1,h2,1,1)
          approx(:,n_sav,2) = approx(:,n_sav,2) * vmk(1:ntot)
          call dssum         (approx(:,n_sav,2))
          approx(:,n_sav,2) = approx(:,n_sav,2) * vml(1:ntot)
          !alpha = glsc3(approx(:,n_sav,2),approx(:,n_sav,2),vml,ntot)
          alpha = glsc2(approx(:,n_sav,2),approx(:,n_sav,2),ntot)
          approx(:,n_sav,:) = approx(:,n_sav,:) *( 1._dp/ sqrt(alpha))

      endif
  else
      n_sav = 1
      v1(:,:,:,1:nelfld(ifield)) = v1(:,:,:,1:nelfld(ifield))  &
                                 + reshape(approx(:,0,1), (/lx1,ly1,lz1,nelfld(ifield)/))
      call copy(approx(:,n_sav,1),v1,ntot)
          call axhelm        (approx(:,n_sav,2),approx(:,n_sav,1),h1,h2,1,1)
          approx(:,n_sav,2) = approx(:,n_sav,2) * vmk(1:ntot)
          call dssum         (approx(:,n_sav,2))
          approx(:,n_sav,2) = approx(:,n_sav,2) * vml(1:ntot)
          !alpha = glsc3(approx(:,n_sav,2),approx(:,n_sav,2),vml,ntot)
          alpha = glsc2(approx(:,n_sav,2),approx(:,n_sav,2),ntot)
          if (alpha > eps) approx(:,n_sav,:) = approx(:,n_sav,:) *( 1._dp/ sqrt(alpha))

  endif

  napprox(2)=n_sav

  return
end subroutine gensh

#if 0
!-----------------------------------------------------------------------
!> \brief Orthonormalize the kth vector against vector set
subroutine hconj(approx,k,h1,h2,vml,vmk,ws,name4,ierr)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nid
  use parallel, only : wdsize
  use tstep, only : istep, ifield, nelfld
  implicit none

  integer, intent(in) :: k
  real(DP), intent(inout) :: approx(:,0:,:)
  real(DP), intent(in) :: h1(*),h2(*),vml(*),vmk(*)
  real(DP), intent(out) :: ws(*)
  character(4) :: name4
  integer, intent(out) :: ierr

  integer :: i, ntot, km1 , nel
  real(DP) :: alpha, ratio, eps, alpham, alph1
  real(DP), external :: glsc2, ddot

  ierr=0
  nel = nelfld(ifield)
  ntot=nx1*ny1*nz1*nel

  call axhelm  (approx(:,0),approx(:,k),h1,h2,1,1)
  approx(:,0) = approx(:,0) * vmk(1:ntot)
  call dssum   (approx(:,0))
  approx(:,0) = approx(:,0) * vml(1:ntot)

!   Compute part of the norm   (Note:  a(0) already scaled by vml)

  alpha = glsc2(approx(:,0),approx(:,k),ntot)
  alph1 = alpha

!   Gram-Schmidt

  km1=k-1
  do i=1,km1
      ws(i) = ddot(ntot, approx(:,0), 1, approx(:,i), 1)
  enddo
  if (km1 > 0) call gop(ws,ws(k),'+  ',km1)

  do i=1,km1
      alpham = -ws(i)
      approx(:,k) = approx(:,k) + alpham * approx(:,i)
      alpha = alpha - ws(i)**2
  enddo

!  .Normalize new element in approximation space

  eps = 1.e-7
  if (wdsize == 8) eps = 1.e-15
  ratio = alpha/alph1

  if (ratio <= 0) then
      ierr=1
      if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
      12 format(I6,1x,a4,' alpha b4 sqrt:',i4,1p2e12.4)
  elseif (ratio <= eps) then
      ierr=2
      if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
  else
      ierr=0
      alpha = 1.0/sqrt(alpha)
      approx(:,k) = alpha * approx(:,k)
  endif

  if (ierr /= 0) then
      call axhelm  (approx(:,0),approx(:,k),h1,h2,1,1)
      approx(:,0) = approx(:,0) * vmk(1:ntot)
      call dssum   (approx(:,0))
      approx(:,0) = approx(:,0) * vml(1:ntot)
  
  !        Compute part of the norm   (Note:  a(0) already scaled by vml)
  
      alpha = glsc2(approx(:,0),approx(:,k),ntot)
      if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
      if (alpha <= 0) then
          ierr=3
          if (nid == 0) write(6,12) istep,name4,k,alpha,alph1
          return
      endif
      alpha = 1.0/sqrt(alpha)
      approx(:,k) = alpha * approx(:,k)
      ierr = 0
  endif

  return
end subroutine hconj
#endif

!-----------------------------------------------------------------------
!> \brief Reorthogonalize approx if dt has changed
subroutine updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1
  use input, only : ifvarp, iflomach
  use tstep, only : dt, ifield, nelfld
  implicit none

  integer, parameter :: lt=lx1*ly1*lz1*lelt
  real(DP), intent(inout) :: approx(:,0:,:)
  real(DP), intent(in)    :: h1(1),h2(1),vml(1),vmk(1)
  real(DP), intent(out) :: ws(1)
  integer, intent(inout) :: napprox(2)
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
      n_sav = napprox(2)
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
      call axhelm(approx(:,k,2),approx(:,k,1),h1,h2,1,1)
#endif 
      enddo
      napprox(2)=min(l,n_sav)
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
    ,approx,napprox,bi)
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
  REAL(DP), intent(inout) :: approx (:,0:,:)        !>!< past solutions for projection
  integer,  intent(inout) :: napprox(2)           !>!< (/ max vecs, current number of vecs /)
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

      allocate(w2(2+2*napprox(1)))
      allocate(w1(lx1*ly1*lz1*lelt))
      call projh  (r,h1,h2,bi,vml,vmk,approx,napprox,w1,w2,name)
      deallocate(w1)

      call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)
      call gensh  (u,h1,h2,vml,vmk,approx,napprox,w2,name)

  endif


  return
end subroutine hsolve
!-----------------------------------------------------------------------

end module helmholtz
