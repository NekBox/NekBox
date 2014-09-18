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
  implicit none

  real(DP), intent(inout) :: r(*) !>!< residual
  real(DP), intent(in)    :: h1(*) !>!< coefficient of A (stiffness)
  real(DP), intent(in)    :: h2(*) !>!< coefficient of M (mass)
  real(DP), intent(in)    :: vml(*) !>!< multiplicity array
  real(DP), intent(in)    :: vmk(*) !>!< mask array
  real(DP), intent(in)    :: bi(*) !>!< inverse mass matrix
  real(DP), intent(out)   :: wl(*) !>!< large work array (size lx1*ly1*lz1*nelv)
  real(DP), intent(out)   :: ws(*) !>!< small work array (size 2*max vecs)
  real(DP), intent(inout) :: approx(:,0:) !>!< approximation space
  integer, intent(inout)  :: napprox(2) !>!< (/ max vecs, current number of vecs /)
  character(4) :: name4

  integer :: n_max, n_sav, nel, ntot, i, n10
  real(DP) :: vol, alpha1, alpha2, ratio
  real(DP), external :: glsc23, vlsc3

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

  do i=1,n_sav
      ws(i) = vlsc3(r,approx(:,i),vml,ntot)
  enddo
  call gop    (ws,ws(n_sav+1),'+  ',n_sav)

  approx(:,0) = approx(:,1) * ws(1)
  do i=2,n_sav
    approx(:,0) = approx(:,0) + approx(:,i) * ws(i)
  enddo

  call axhelm  (wl,approx(:,0),h1,h2,1,1)
  call col2    (wl,vmk,ntot)
  call dssum   (wl)
  r(1:ntot) = r(1:ntot) - wl(1:ntot)

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
  use size_m, only : lx1, ly1, lz1, mxprev
  use mesh, only : niterhm
  use tstep, only : nelfld, ifield
  implicit none

  REAL(DP), intent(inout) :: V1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: H1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: H2   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: vmk  (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)    :: vml  (LX1,LY1,LZ1,*)
  real(DP), intent(out)   :: ws(2+2*mxprev) !>!< workspace?

  real(DP) :: approx(:,0:)
  integer :: napprox(2)
  character(4) :: name4

  integer :: n_max, n_sav, ntot, ierr

  n_max = napprox(1)
  n_sav = napprox(2)
  ntot=nx1*ny1*nz1*nelfld(ifield)

!   Reconstruct solution and save current du

  if (n_sav < n_max) then
  
      if (niterhm > 0) then      ! new vector not in space
          n_sav = n_sav+1
          call copy(approx(:,n_sav),v1,ntot)
          call add2(v1,approx(:,0),ntot)
      !           orthogonalize rhs against previous rhs and normalize
          call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)

      !           if (ierr.ne.0) n_sav = n_sav-1
          if (ierr /= 0) n_sav = 0

      else

          call add2(v1,approx(:,0),ntot)

      endif
  else
      n_sav = 1
      call add2(v1,approx(:,0),ntot)
      call copy(approx(:,n_sav),v1,ntot)
  !        normalize
      call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)
      if (ierr /= 0) n_sav = 0
  endif

  napprox(2)=n_sav

  return
end subroutine gensh

!-----------------------------------------------------------------------
!> \brief Orthonormalize the kth vector against vector set
subroutine hconj(approx,k,h1,h2,vml,vmk,ws,name4,ierr)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nid
  use parallel, only : wdsize
  use tstep, only : istep, ifield, nelfld
  implicit none

  integer, intent(in) :: k
  real(DP), intent(inout) :: approx(:,0:)
  real(DP), intent(in) :: h1(*),h2(*),vml(*),vmk(*)
  real(DP), intent(out) :: ws(*)
  character(4) :: name4
  integer, intent(out) :: ierr

  integer :: i, ntot, km1 
  real(DP) :: alpha, ratio, eps, alpham, alph1
  real(DP), external :: glsc2, vlsc2

  ierr=0
  ntot=nx1*ny1*nz1*nelfld(ifield)

  call axhelm  (approx(:,0),approx(:,k),h1,h2,1,1)
  call col2    (approx(:,0),vmk,ntot)
  call dssum   (approx(:,0))
  call col2    (approx(:,0),vml        ,ntot)

!   Compute part of the norm   (Note:  a(0) already scaled by vml)

  alpha = glsc2(approx(:,0),approx(:,k),ntot)
  alph1 = alpha

!   Gram-Schmidt

  km1=k-1
  do i=1,km1
      ws(i) = vlsc2(approx(:,0),approx(:,i),ntot)
  enddo
  if (km1 > 0) call gop(ws,ws(k),'+  ',km1)

  do i=1,km1
      alpham = -ws(i)
      call add2s2(approx(:,k),approx(:,i),alpham,ntot)
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
      call col2    (approx(:,0),vmk,ntot)
      call dssum   (approx(:,0))
      call col2    (approx(:,0),vml        ,ntot)
  
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

!-----------------------------------------------------------------------
!> \brief Reorthogonalize approx if dt has changed
subroutine updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1
  use input, only : ifvarp, iflomach
  use tstep, only : dt, ifield, nelfld
  implicit none

  integer, parameter :: lt=lx1*ly1*lz1*lelt
  real(DP), intent(inout) :: approx(lt,0:1)
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
      !           Orthogonalize kth vector against {v_1,...,v_k-1}
          if (k /= l) then
              ntot = nx1*ny1*nz1*nelfld(ifield)
              call copy(approx(1,l),approx(1,k),ntot)
          endif
          call hconj(approx,l,h1,h2,vml,vmk,ws,name4,ierr)
          if (ierr == 0) l=l+1
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
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, mxprev
  use input, only : ifflow, param
  use string, only : capit
  use tstep, only : ifield, nelfld, istep
  implicit none

  CHARACTER(4), intent(in) :: NAME !>!< name of field we're solving for
  REAL(DP), intent(out)   :: U    (LX1,LY1,LZ1,*) !>!< solution vector
  REAL(DP), intent(inout) :: R    (LX1,LY1,LZ1,*) !>!< right hand side
  REAL(DP), intent(in)    :: H1   (LX1,LY1,LZ1,*) !>!< coefficient of A (stiffness)
  REAL(DP), intent(in)    :: H2   (LX1,LY1,LZ1,*) !>!< coefficient of M (mass)
  REAL(DP), intent(in)    :: vmk  (LX1,LY1,LZ1,*) !>!< mask array
  REAL(DP), intent(in)    :: vml  (LX1,LY1,LZ1,*) !>!< multiplicity array
  integer,  intent(in)    :: imsh                 !>!< imesh?
  real(DP), intent(in)    :: tol                  !>!< residual tolerance
  integer,  intent(in)    :: maxit                !>!< maximum number of iterations
  integer,  intent(in)    :: isd                  !>!< something to do with axi-symmetric
  REAL(DP), intent(inout) :: approx (:,0:)        !>!< past solutions for projection
  integer,  intent(inout) :: napprox(2)           !>!< (/ max vecs, current number of vecs /)
  REAL(DP), intent(in)    :: bi   (LX1,LY1,LZ1,*) !>!< inverse of mass matrix

  real(DP), allocatable :: w1(:)
  real(DP) :: w2(2+2*mxprev)

  logical :: ifstdh
  character(4) ::  cname
  integer :: n

  call chcopy(cname,name,4)
  call capit (cname,4)

  ifstdh = .TRUE. 

  if ( .NOT. ifflow) ifstdh = .FALSE. 

  if (param(95) /= 0 .AND. istep > param(95)) then
      if (cname == 'PRES') ifstdh = .FALSE. 
  elseif (param(94) /= 0 .AND. istep > param(94)) then
      ifstdh = .FALSE. 
  endif

  if (param(93) == 0) ifstdh = .TRUE. 

  if (ifstdh) then
      call hmholtz(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)
  else
      
      n = nx1*ny1*nz1*nelfld(ifield)

      call dssum  (r)
      call col2   (r,vmk,n)
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
