!=======================================================================
subroutine hmholtz(name,u,rhs,h1,h2,mask,mult,imsh,tli,maxit,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, nelv, nelt, ndim, nid
  use ctimer, only : icalld, thmhz, nhmhz, dnekclock
  use ds, only : dssum
  use fdmh1, only : kfldfdm
  use input, only : ifsplit, param
  use geom, only : binvm1, bintm1
  use tstep, only : istep, nelfld, ifield
  implicit none

  CHARACTER(4)      NAME
  REAL(DP) ::           U    (LX1,LY1,LZ1,*)
  REAL(DP) ::           RHS  (LX1,LY1,LZ1,*)
  REAL(DP) ::           H1   (LX1,LY1,LZ1,*)
  REAL(DP) ::           H2   (LX1,LY1,LZ1,*)
  REAL(DP) ::           MASK (LX1,LY1,LZ1,*)
  REAL(DP) ::           MULT (LX1,LY1,LZ1,*)
  real(DP) :: tli
  integer :: imsh, maxit, isd

  logical :: iffdm
  character(3) :: nam3
  integer :: ntot
  real(DP) :: tol
  real(DP) :: etime

  tol = abs(tli)

  iffdm = .FALSE. 
!   iffdm = .true.
  if (ifsplit) iffdm = .TRUE. 

!max    if (icalld == 0 .AND. iffdm) call set_fdm_prec_h1A

  nhmhz=nhmhz + 1
  etime=dnekclock()
  ntot = nx1*ny1*nz1*nelfld(ifield)
  if (imsh == 1) ntot = nx1*ny1*nz1*nelv
  if (imsh == 2) ntot = nx1*ny1*nz1*nelt

!   Determine which field is being computed for FDM based preconditioner bc's

  call chcopy(nam3,name,3)

  kfldfdm = -1
!   if (nam3.eq.'TEM' ) kfldfdm =  0
!   if (name.eq.'TEM1') kfldfdm =  0  ! hardcode for temp only, for mpaul
!   if (name.eq.'VELX') kfldfdm =  1
!   if (name.eq.'VELY') kfldfdm =  2
!   if (name.eq.'VELZ') kfldfdm =  3
  if (name == 'PRES') kfldfdm =  ndim+1
!   if (.not.iffdm) kfldfdm=-1

  call dssum   (rhs(:,1,1,1))
  rhs(:,:,:,1:nelfld(ifield)) = rhs(:,:,:,1:nelfld(ifield)) * mask(:,:,:,1:nelfld(ifield))
  if (nid == 0 .AND. istep <= 10) &
  write(6,*) param(22),' p22 ',istep,imsh
  if (param(22) == 0 .OR. istep <= 10) &
  call chktcg1 (tol,rhs,h1,h2,mask,mult,imsh,isd)

  if (tli < 0) tol=tli ! caller-specified relative tolerance

  if (imsh == 1) call cggo &
  (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1,name)
  if (imsh == 2) call cggo &
  (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,bintm1,name)


  thmhz=thmhz+(dnekclock()-etime)
  return
end subroutine hmholtz

!=======================================================================
!------------------------------------------------------------------
!> \brief Compute the (Helmholtz) matrix-vector product,
!!  AU = helm1*[A]u + helm2*[B]u, for NEL elements.
!------------------------------------------------------------------
subroutine axhelm (au,u,helm1,helm2,imesh,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1
  use size_m, only : nx1, ny1, nz1, nelt, nelv, ndim
  use ctimer, only : taxhm, naxhm, etime1, dnekclock
  use ctimer, only : axhelm_flop, axhelm_mop
  use geom, only : g4m1, g5m1, g6m1
  use dxyz, only : wddx, wddyt, wddzt
  use input, only : ifaxis
  use geom, only : bm1
  use mesh, only : ifsolv
  implicit none

  real(DP), intent(out) :: AU    (LX1,LY1,LZ1,*) !>!< H u
  real(DP), intent(in)  :: U     (LX1,LY1,LZ1,*) !>!< u
  real(DP), intent(in)  :: HELM1 (LX1,LY1,LZ1,*) !>!< coefficient of stif
  real(DP), intent(in)  :: HELM2 (LX1,LY1,LZ1,*) !>!< coefficient of mass
  integer,  intent(in)  :: imesh !>!< mesh index (v or t)
  integer,  intent(in)  :: isd !>!< axi-symmetric flag of sorts

   ! locals
  REAL(DP) ::           TM1   (LX1,LY1,LZ1)
  REAL(DP) ::           TM2   (LX1,LY1,LZ1)
  REAL(DP) ::           TM3   (LX1,LY1,LZ1)

  integer :: nel, nxy, nyz, nxz, nxyz, ntot
  integer :: e, iz

  nel=nelt
  if (imesh == 1) nel=nelv

  NXY=NX1*NY1
  NYZ=NY1*NZ1
  NXZ=NX1*NZ1
  NXYZ=NX1*NY1*NZ1
  NTOT=NXYZ*NEL

  if (naxhm == 0) taxhm=0.0
  naxhm=naxhm + 1

  etime1=dnekclock()

  if (helm2(1,1,1,1) /= 0._dp) then
    axhelm_mop = axhelm_mop + 6*nx1*ny1*nz1*nel
  else
    axhelm_mop = axhelm_mop + 5*nx1*ny1*nz1*nel
  endif

  do e=1,nel
    ! Fast 3-d mode: constant properties and undeformed element
    axhelm_flop = axhelm_flop + (2*nx1-1)*nx1*nyz
    axhelm_flop = axhelm_flop + 9*nx1*ny1*nz1
    axhelm_flop = axhelm_flop + (2*ny1-1)*nx1*ny1*nz1
    axhelm_flop = axhelm_flop + nxy*(2*nz1-1)*nz1

    call helmholtz(helm1(1,1,1,1), helm2(1,1,1,1), nx1, ny1, nz1, &
                   u(1,1,1,e), au(1,1,1,e), &
                   g4m1(1,1,1,e), g5m1(1,1,1,e), g6m1(1,1,1,e), bm1(1,1,1,e), &
                   tm1, tm2, tm3)
  END DO

  taxhm=taxhm+(dnekclock()-etime1)
  return
end subroutine axhelm

!----------------------------------------------------------------------
!> \brief For undeformed elements, set up appropriate elemental matrices
!!        and geometric factors for fast evaluation of Ax.
!----------------------------------------------------------------------
subroutine sfastax()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelt, ndim
  use dxyz, only : dxm1, dym1, dzm1
  use dxyz, only : wddx, wddyt, wddzt
  use geom, only : g1m1, g2m1, g3m1, g4m1, g5m1, g6m1
  use mesh, only : ifdfrm
  use wz_m, only : wxm1, wym1, wzm1
  implicit none

  logical, save :: IFIRST = .TRUE.
  integer :: nxx, nyy, nzz
  integer :: i, j, ip, ie, ix, iy, iz

  NXX=NX1*NX1
  IF (IFIRST) THEN
      wddx = 0._dp
      DO I=1,NX1
          DO J=1,NX1
              DO IP=1,NX1
                  WDDX(I,J) = WDDX(I,J) + WXM1(IP)*DXM1(IP,I)*DXM1(IP,J)
              enddo
          enddo
      END DO
      NYY=NY1*NY1
      wddyt = 0._dp
      DO I=1,NY1
          DO J=1,NY1
              DO IP=1,NY1
                  WDDYT(J,I) = WDDYT(J,I) + WYM1(IP)*DYM1(IP,I)*DYM1(IP,J)
              enddo
          enddo
      END DO
      NZZ=NZ1*NZ1
      wddzt = 0._dp
      DO I=1,NZ1
          DO J=1,NZ1
              DO IP=1,NZ1
                  WDDZT(J,I) = WDDZT(J,I) + WZM1(IP)*DZM1(IP,I)*DZM1(IP,J)
              enddo
          enddo
      END DO
      IFIRST= .FALSE. 
  ENDIF

  IF (NDIM == 3) THEN
      DO 1001 IE=1,NELT
          IF ( .NOT. IFDFRM(IE)) THEN
              DO IZ=1,NZ1
                  DO IY=1,NY1
                      DO IX=1,NX1
                          G4M1(IX,IY,IZ,IE)=G1M1(IX,IY,IZ,IE)/WXM1(IX)
                          G5M1(IX,IY,IZ,IE)=G2M1(IX,IY,IZ,IE)/WYM1(IY)
                          G6M1(IX,IY,IZ,IE)=G3M1(IX,IY,IZ,IE)/WZM1(IZ)
                      enddo
                  enddo
              END DO
          ENDIF
      1001 END DO
  ELSE
      DO 2001 IE=1,NELT
          IF ( .NOT. IFDFRM(IE)) THEN
              DO IY=1,NY1
                  DO IX=1,NX1
                      G4M1(IX,IY,1,IE)=G1M1(IX,IY,1,IE)/WXM1(IX)
                      G5M1(IX,IY,1,IE)=G2M1(IX,IY,1,IE)/WYM1(IY)
                  enddo
              END DO
          ENDIF
      2001 END DO
  ENDIF
  return
  end subroutine sfastax

!=======================================================================
!> \brief Generate diagonal preconditioner for the Helmholtz operator.
!-------------------------------------------------------------------
subroutine setprec (dpcm1,helm1,helm2,imsh,isd)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1, nelt, nelv, ndim
  use ds, only : dssum
  use dxyz, only : dxtm1, dytm1, dztm1, datm1, dam1
  use geom, only : ifrzer, g1m1, g2m1, g3m1, ym1, jacm1
  use input, only : ifaxis
  use geom, only : bm1
  use mesh, only : ifdfrm
  use wz_m, only : wxm1, wam1
  use ctimer, only : tdpc, ndpc, dnekclock
  implicit none

  REAL(DP), intent(out) ::  DPCM1 (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  ::  HELM1(NX1,NY1,NZ1,*), HELM2(NX1,NY1,NZ1,*)
  integer, intent(in) :: imsh, isd

  REAL(DP) :: YSM1(LY1)

  integer :: nel, ntot, ie, iq, iz, iy, ix, j, i, iel
  real(DP) :: term1, term2
  real(DP) :: etime

  ndpc = ndpc + 1
  etime = dnekclock()

  nel=nelt
  if (imsh == 1) nel=nelv

  ntot = nel*nx1*ny1*nz1

!   The following lines provide a convenient debugging option
!   call rone(dpcm1,ntot)
!   if (ifield.eq.1) call copy(dpcm1,binvm1,ntot)
!   if (ifield.eq.2) call copy(dpcm1,bintm1,ntot)
!   return

  dpcm1(:,:,:,1:nel) = 0._dp
  DO IE=1,NEL

          DO IZ=1,NZ1
              DO IY=1,NY1
                  DO IX=1,NX1
      DO IQ=1,NX1
                      DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                      G1M1(IQ,IY,IZ,IE) * DXTM1(IX,IQ)**2
      END DO
                  enddo
              enddo
          enddo
          DO IZ=1,NZ1
              DO IY=1,NY1
                  DO IX=1,NX1
      DO IQ=1,NY1
                      DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                      G2M1(IX,IQ,IZ,IE) * DYTM1(IY,IQ)**2
      END DO
                  enddo
              enddo
          enddo

              DO IZ=1,NZ1
                  DO IY=1,NY1
                      DO IX=1,NX1
          DO IQ=1,NZ1
                          DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                          G3M1(IX,IY,IQ,IE) * DZTM1(IZ,IQ)**2
          END DO
                      enddo
                  enddo
              enddo
      
  END DO

  dpcm1(:,:,:,1:nel) = dpcm1(:,:,:,1:nel)*helm1(:,:,:,1:nel) + bm1*helm2(:,:,:,1:nel)

  CALL DSSUM (DPCM1(:,1,1,1))
  dpcm1(:,:,:,1:nel) = 1._dp / dpcm1(:,:,:,1:nel)
  tdpc = tdpc + (dnekclock() - etime)

  return
end subroutine setprec

!-------------------------------------------------------------------
!> \brief Check that the tolerances are not too small for the CG-solver.
!! Important when calling the CG-solver (Gauss-Lobatto mesh) with
!! zero Neumann b.c.
!-------------------------------------------------------------------
subroutine chktcg1 (tol,res,h1,h2,mask,mult,imesh,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1
  use size_m, only : nx1, ny1, nz1, nelv, nelt, nid
  use eigen, only : eigaa, eigga
  use input, only : ifprint
  use geom, only : volvm1, voltm1, binvm1, bintm1, bm1
  use ctimer, only : nfoo, tfoo, dnekclock
  use ctimer, only : othr_flop, othr_mop
  implicit none

  real(DP), intent(inout) :: tol
  REAL(DP), intent(in)  :: RES  (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: H1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: H2   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: MASK (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: MULT (LX1,LY1,LZ1,*)
  integer,  intent(in)  :: imesh, isd

  real(DP), allocatable :: W1(:,:,:,:), W2(:,:,:,:)

  real(DP) :: acondno, delta, x, y, diff, eps
  real(DP) :: vol, rinit, rmin, bctest, bcrob, tolmin 
  real(DP) :: bcneu1, bcneu2 
  real(DP) :: etime
  integer :: nl, ntot1
  real(DP), external :: glsc3, glsum

  return

  nfoo = nfoo + 1
  etime = dnekclock()

  IF (EIGAA /= 0.) THEN
      ACONDNO = EIGGA/EIGAA
  ELSE
      ACONDNO = 10.
  ENDIF

!   Single or double precision???

  DELTA = 1.E-9
  X     = 1.+DELTA
  Y     = 1.
  DIFF  = ABS(X-Y)
  IF (DIFF == 0.) EPS = 1.E-6
  IF (DIFF > 0.) EPS = 1.E-13
  ! max
  IF (DIFF > 0.) EPS = 1.E-16

  IF (IMESH == 1) THEN
      NL  = NELV
      VOL = VOLVM1
  ELSEIF (IMESH == 2) THEN
      NL  = NELT
      VOL = VOLTM1
  else
    nl = -1
    vol = -1.0
    write(*,*) "Got somewhere bad"
  ENDIF

  allocate(W1(nx1,ny1,nz1,nl), W2(nx1,ny1,nz1,nl))
  NTOT1 = NX1*NY1*NZ1*NL

  othr_mop = othr_mop + ntot1*5
  othr_flop = othr_flop + ntot1*4
  IF (IMESH == 1) THEN
      w2 = binvm1 * res(:,:,:,1:nl)
      RINIT  = SQRT(GLSC3 (W2,res,MULT,NTOT1)/VOLVM1)
  ELSE
      w2 = bintm1 * res(:,:,:,1:nl)
      RINIT  = SQRT(GLSC3 (W2,res,MULT,NTOT1)/VOLTM1)
  ENDIF
  RMIN   = EPS*RINIT
  IF (TOL < RMIN) THEN
      IF (NID == 0 .AND. IFPRINT) &
      WRITE (6,*) 'New CG1-tolerance (RINIT*epsm) = ',RMIN,TOL
      TOL = RMIN
  ENDIF

  othr_mop = othr_mop + ntot1*4
  othr_flop = othr_flop + ntot1*6
  w1 = 1._dp
  BCNEU1 = GLSC3(W1,MASK,MULT,NTOT1)
  BCNEU2 = GLSC3(W1,W1  ,MULT,NTOT1)
  BCTEST = ABS(BCNEU1-BCNEU2)

  if (BCTEST >= .1) then
    tfoo = tfoo + (dnekclock() - etime)
    return
  endif

  othr_mop = othr_mop + ntot1*4
  othr_flop = othr_flop + ntot1*3
  etime = etime - dnekclock()
  CALL AXHELM (W2,W1,H1,H2,IMESH,ISD)
  etime = etime + dnekclock()
  w2 = w2 * w2 * bm1
  BCROB  = SQRT(GLSUM(W2,NTOT1)/VOL)

  IF (BCROB < (EPS*ACONDNO)) THEN
  !         OTR = GLSC3 (W1,RES,MULT,NTOT1)
      TOLMIN = RINIT*EPS*10.
      IF (TOL < TOLMIN) THEN
          TOL = TOLMIN
          IF (NID == 0 .AND. IFPRINT) &
          WRITE(6,*) 'New CG1-tolerance (Neumann) = ',TOLMIN
      ENDIF
  ENDIF

  tfoo = tfoo + (dnekclock() - etime)

  return
  end subroutine chktcg1

!=======================================================================
!-------------------------------------------------------------------------
!> \brief Solve the Helmholtz equation, H*U = RHS,
!! using preconditioned conjugate gradient iteration.
!! Preconditioner: diag(H).
!------------------------------------------------------------------------



subroutine cggo(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name)
  use kinds, only : DP
  use size_m, only : nid, nx1, ny1, nz1, nelt, nelv
  use size_m, only : lx1, ly1, lz1, lelt
  use ds, only : dssum
  use fdmh1, only : kfldfdm
  use input, only : ifsplit, param, ifprint
  use geom, only : volvm1, voltm1
  use mesh, only : ifsolv, niterhm
  use tstep, only : istep, imesh
#ifdef XSMM
  use STREAM_UPDATE_KERNELS, only : stream_vector_compscale
#endif
  use ctimer, only : ncggo, tcggo, cggo_flop, cggo_mop, dnekclock
  implicit none

  real(DP), intent(out) :: x(lx1*ly1*lz1,*) !>!< solution vector
  real(DP), intent(in)  :: f(lx1*ly1*lz1,*) !>!< residual vector
  real(DP), intent(in)  :: h1(lx1*ly1*lz1,*) !>!< coefficient of A (stiffness)
  real(DP), intent(in)  :: h2(lx1*ly1*lz1,*) !>!< coefficient of M (mass)
  real(DP), intent(in)  :: mask(lx1*ly1*lz1,*) !>!< mask array
  real(DP), intent(in)  :: mult(lx1*ly1*lz1,*) !>!< multiplicity array
  real(DP), intent(in)  :: binv(lx1*ly1*lz1,*) !>!< inverse of mass matrix
  integer,  intent(in)  :: imsh, isd, maxit
  real(DP), intent(in)  :: tin !>!< input tolerance
  character(4) :: name

  logical :: ifmcor,ifprint_hmh

  integer, parameter :: lxyz = lx1*ly1*lz1
  real(DP) :: scalar(2)
  real(DP), allocatable, save :: d (:,:)
  real(DP), allocatable :: r (:,:) , w (:,:) , p (:,:) , z (:,:)

  integer :: i, n, iter, nxyz, nel, niter
  real(DP) :: rho, vol, tol, h2max, skmin
  real(DP) :: rtz1, rtz2, rbn2, rbn0, beta, rho0, alpha, alphm
  real(DP), external :: glmax, glmin, glsum, glsc2, vlsc3, glsc3
  real(DP), save :: h1_prec = 0._dp, h2_prec = 0._dp
  real(DP) :: etime

  if (ifsplit .AND. name == 'PRES' .AND. param(42) == 0) then
      n = nx1*ny1*nz1*nelv
      call copy      (x,f,n)
      call hmh_gmres (x,h1,h2,mult,iter)
      niterhm = iter
      return
  endif
!      write(6,*) ifsplit,name,param(44),' P44 C'

! **  zero out stuff for Lanczos eigenvalue estimator
  ncggo = ncggo + 1

  rho = 0.00

!     Initialization
  NXYZ   = NX1*NY1*NZ1
  NEL    = NELV
  VOL    = VOLVM1
  IF (IMSH == 2) NEL=NELT
  IF (IMSH == 2) VOL=VOLTM1
  n      = NEL*NXYZ

  tol=abs(tin)
  if (param(22) /= 0) tol=abs(param(22))
  if (name == 'PRES' .AND. param(21) /= 0) tol=abs(param(21))
  if (tin < 0)       tol=abs(tin)
  niter = maxit

!     Speed-up for undeformed elements and constant properties.
  if ( .NOT. ifsolv) then
      ifsolv = .TRUE. 
  endif

!     Set up diag preconditioner.

  if (.not. allocated(d)) allocate(d(lxyz,lelt))

  if (kfldfdm < 0) then
    if (h1(1,1) /= h1_prec .or. h2(1,1) /= h2_prec) then
      call setprec(D,h1,h2,imsh,isd)
      h1_prec = h1(1,1); h2_prec = h2(1,1)
    endif
  elseif(param(100) /= 2) then
!max      call set_fdm_prec_h1b(d,h1,h2,nel)
  endif
  etime = dnekclock()

  allocate(r(nxyz,nel))
  cggo_mop = cggo_mop + 3*n
  r = f(:,1:nel)
  x(:,1:nel) = 0._dp

!     Check for non-trivial null-space

  ifmcor = .FALSE. 
  cggo_mop = cggo_mop + 2*n
  h2max = glmax(h2  ,n)
  skmin = glmin(mask,n)
  if (skmin > 0 .AND. h2max == 0) then
    ifmcor = .TRUE. 
    write(*,*) "Oops! ifmcor"
  endif
  if (name == 'PRES') then
    write(*,*) "Oops! name == PRES"
  endif
  if (kfldfdm >= 0) then
    write(*,*) "Oops! kfldfdm"
  endif

  ifprint_hmh = .FALSE. 
  if (nid == 0 .AND. ifprint .AND. param(74) /= 0) ifprint_hmh= .TRUE. 
  if (nid == 0 .AND. istep == 1)                 ifprint_hmh= .TRUE. 

  niterhm = 0

  cggo_mop = cggo_mop + n
  allocate(w(nxyz,nel))
  allocate(p(nxyz,nel))
  iter = 1

  rtz1=1._dp; rtz2 = 1._dp

  cggo_flop = cggo_flop + 8*n
  cggo_mop  = cggo_mop  + 5*n
  scalar(1) = 0._dp; scalar(2) = 0._dp
  do i = 1, nel
#ifdef XSMM
    call stream_vector_compscale(r(:,i), d(:,i), p(:,i), nx1*ny1*nz1)
#else
    p(:,i) = r(:,i) * d(:,i) 
#endif
    scalar(1) = scalar(1) + sum(r(:,i) * d(:,i) * r(:,i) * mult(:,i))
    scalar(2) = scalar(2) + sum(r(:,i) * r(:,i) * mult(:,i) * binv(:,i))
  enddo
  call gop(scalar,w,'+  ',2)
  rtz1=scalar(1)
  rbn2=sqrt(scalar(2)/vol)
  rbn0 = rbn2
  if (param(22) < 0) tol=abs(param(22))*rbn0
  if (tin < 0)       tol=abs(tin)*rbn0
  if (rbn0 < tol) return

  if (ifprint_hmh) &
  write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1,1),h2(1,1)
  
  beta = 0._dp

  etime = etime - dnekclock()
  call axhelm (w,p,h1,h2,imsh,isd)
  etime = etime + dnekclock()
  call dssum  (w(:,1))

  cggo_flop = cggo_flop + 4*n
  cggo_mop  = cggo_mop  + 5*n
  rho0 = rho; rho = 0._dp
  do i = 1, nel
    w(:,i)   = w(:,i) * mask(:,i)
    rho = rho + sum(w(:,i) * p(:,i) * mult(:,i))
  enddo
  allocate(z(nxyz,nel))
  call gop(rho,z,'+  ',1)
  alpha=rtz1/rho
  alphm=-alpha

  do iter=2,niter
    scalar(1) = 0._dp; scalar(2) = 0._dp
    rtz2=rtz1
    cggo_flop = cggo_flop + 10*n
    cggo_mop  = cggo_mop  + 7*n
    do i = 1, nel
      r(:,i) = r(:,i) + alphm * w(:,i)
#ifdef XSMM
      call stream_vector_compscale(r(:,i), d(:,i), z(:,i), nx1*ny1*nz1)
#else
      z(:,i) = r(:,i) * d(:,i) 
#endif
      scalar(1) = scalar(1) + sum(r(:,i) * d(:,i) * r(:,i) * mult(:,i))
      scalar(2) = scalar(2) + sum(r(:,i) * r(:,i) * mult(:,i) * binv(:,i))
    enddo
    call gop(scalar,w,'+  ',2)
    rtz1=scalar(1)
    rbn2=sqrt(scalar(2)/vol)

    if (ifprint_hmh) &
    write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1,1),h2(1,1)


    !   Always take at least one iteration   (for projection) pff 11/23/98
#ifndef TST_WSCAL
    !IF (rbn2 <= TOL .AND. iter > 10) THEN
    IF (rbn2 <= TOL .AND. (iter > 1 .OR. istep <= 5)) THEN
#else
    iter_max = param(150)
    if (name == 'PRES') iter_max = param(151)
    if (iter > iter_max) then
#endif
    !        IF (rbn2.LE.TOL) THEN
        cggo_flop = cggo_flop + 2*n
        cggo_mop  = cggo_mop  + 3*n
        x(:,1:nel) = x(:,1:nel) + alpha * p
        NITER = ITER-1
    !           IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
        if (nid == 0) &
        write(6,3000) istep,name,niter,rbn2,rbn0,tol
        goto 9999
    ENDIF
    
    cggo_flop = cggo_flop + 4*n 
    cggo_mop  = cggo_mop  + 5*n 
    beta = rtz1/rtz2
    do i = 1, nel
      x(:,i) = x(:,i) + alpha * p(:,i)
      p(:,i) = z(:,i) + beta  * p(:,i)
    enddo
    etime = etime - dnekclock()
    call axhelm (w,p,h1,h2,imsh,isd)
    etime = etime + dnekclock()
    call dssum  (w(:,1))    

    cggo_flop = cggo_flop + 4*n
    cggo_mop  = cggo_mop  + 5*n
    rho0 = rho; rho = 0._dp
    do i = 1, nel
      w(:,i)   = w(:,i) * mask(:,i)
      rho = rho + sum(w(:,i) * p(:,i) * mult(:,i))
    enddo
    call gop(rho,z,'+  ',1)
    alpha=rtz1/rho
    alphm=-alpha
  enddo

  niter = iter-1

  if (nid == 0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
  3000 format(4x,i7,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
  3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
  3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
  9999 continue
  niterhm = niter
  ifsolv = .FALSE. 
  tcggo = tcggo + (dnekclock() - etime)
    
  return
end subroutine cggo

!=======================================================================
subroutine set_fdm_prec_h1b(d,h1,h2,nel)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  implicit none

  real(DP), intent(out) :: d (nx1,ny1,nz1,1)
  real(DP), intent(in)  :: h1(nx1,ny1,nz1,1)
  real(DP), intent(in)  :: h2(nx1,ny1,nz1,1)
  integer,  intent(in)  :: nel

  d = 0._dp

  return
end subroutine set_fdm_prec_h1b
!-----------------------------------------------------------------------
