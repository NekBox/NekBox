!> \file plan4.F90
!! \brief Routine plan4 and supporting routines
!!
!! This file contains high-level routines for computing the 
!! right hand sides for the Helmholtz and Poisson solves in the
!! P(N)-P(N) formalism (colocation). 

!-----------------------------------------------------------------------
!> \brief Splitting scheme \cite Tomboulides1989 
!!
!! NOTE: QTL denotes the so called thermal
!!       divergence and has to be provided
!!       by an external subroutine e.g qthermal
subroutine plan4()
  use kinds, only : DP
  use size_m, only : nid 
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : lx2, ly2, lz2, lelv
  use size_m, only : nx1, ny1, nz1, nelv
  use ctimer, only : icalld, tpres, npres, etime1, dnekclock
  use ctimer, only : np4misc, tp4misc
  use geom, only : binvm1, bm1, volvm1
  use helmholtz, only : hsolve, approx_space, init_approx_space
  use soln, only : vx, vy, vz, v1mask, v2mask, v3mask
  use soln, only : vtrans, pmask, vmult, pr
  use tstep, only : imesh, nmxh, tolhv
  use input, only : param
  implicit none

  integer, parameter :: ktot = lx1*ly1*lz1*lelt
  integer :: laxt, v_proj_size

  type(approx_space), save :: p_apx, vx_apx, vy_apx, vz_apx

  real(DP), allocatable :: RES1(:,:,:,:)
  real(DP), allocatable :: RES2(:,:,:,:)   
  real(DP), allocatable :: RES3(:,:,:,:)      
  real(DP), allocatable :: DV1 (:,:,:,:)   
  real(DP), allocatable :: DV2 (:,:,:,:)   
  real(DP), allocatable :: RESPR (:,:,:,:)

  real(DP), allocatable :: h1(:,:,:,:), h2(:,:,:,:)
  real(DP), allocatable :: VEXT(:,:)
  REAL(DP), allocatable :: DPR(:,:,:,:)

  REAL(DP), allocatable :: DVC (:,:,:,:), DFC(:,:,:,:)
  REAL(DP) :: DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
  real(DP) :: tolspl
  real(DP) :: etime
  real(DP), external :: glsum

  integer :: intype, ntot1 

  etime = dnekclock()
  np4misc = np4misc + 1
  v_proj_size = int(param(92))
  laxt = int(param(93))

  if (.not. allocated(p_apx%projectors)) then
    call init_approx_space(p_apx, laxt, ktot)
  endif
  if (.not. allocated(vx_apx%projectors)) then
    call init_approx_space(vx_apx, v_proj_size, ktot)
  endif
  if (.not. allocated(vy_apx%projectors)) then
    call init_approx_space(vy_apx, v_proj_size, ktot)
  endif
  if (.not. allocated(vz_apx%projectors)) then
    call init_approx_space(vz_apx, v_proj_size, ktot)
  endif

  INTYPE = -1
  NTOT1  = NX1*NY1*NZ1*NELV

  ! add user defined divergence to qtl
!max  call add2 (qtl,usrdiv,ntot1)

  allocate(vext(lx1*ly1*lz1*lelv,3))

  ! Time-advance velocity with AB(k)
  CALL V_EXTRAP(vext)

  ! compute explicit contributions (bf{x,y,z}) with BDF(k)
  etime = etime - dnekclock()
  CALL MAKEF()
  etime = etime + dnekclock()

  ! store the velocity field for later
  CALL LAGVEL()

! split viscosity into explicit/implicit part
!    if (ifexplvis) call split_vis

! extrapolate velocity

  ! mask Dirichlet boundaries
  CALL BCDIRVC  (VX,VY,VZ,v1mask,v2mask,v3mask)

!     first, compute pressure
#ifndef NOTIMER
  if (icalld == 0) tpres=0.0
  icalld=icalld+1
  npres=icalld
  etime1=dnekclock()
#endif

  allocate(RESPR(LX2,LY2,LZ2,LELV))
  etime = etime - dnekclock()
  call crespsp  (respr, vext)
  etime = etime + dnekclock()
  deallocate(vext)

  allocate(h1(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv))
  h1(1,1,1,1) = 1._dp / vtrans(1,1,1,1,1)
  h2(1,1,1,1) = 0._dp
  call ctolspl  (tolspl,respr)

  allocate(dpr(lx2,ly2,lz2,lelv))
  etime = etime - dnekclock()
  call hsolve ('PRES', dpr, respr, h1, h2, pmask, vmult &
  ,imesh,tolspl,nmxh,1 &
  ,p_apx,binvm1)
  etime = etime + dnekclock()
  deallocate(respr)
  pr = pr + dpr
  deallocate(dpr)
  call ortho   (pr)
#ifndef NOTIMER
  tpres=tpres+(dnekclock()-etime1)
#endif

!     Compute velocity
  allocate(RES1(lx1,ly1,lz1,lelv), &
           RES2(LX1,LY1,LZ1,LELV), &
           RES3(LX1,LY1,LZ1,LELV)  )

  etime = etime - dnekclock()
  call cresvsp (res1,res2,res3,h1,h2)
  etime = etime + dnekclock()

  allocate(DV1 (LX1,LY1,LZ1,LELV)) 

  !> \note These three calls are task-parallel
  etime = etime - dnekclock()
  call hsolve('VELX', dv1, res1, h1, h2, v1mask, vmult, imesh, tolhv, nmxh, 1, &
              vx_apx, binvm1)
  vx = vx + dv1
  call hsolve('VELY', dv1, res2, h1, h2, v2mask, vmult, imesh, tolhv, nmxh, 2, &
              vy_apx, binvm1)
  vy = vy + dv1
  call hsolve('VELZ', dv1, res3, h1, h2, v3mask, vmult, imesh, tolhv, nmxh, 3, &
              vz_apx, binvm1)
  vz = vz + dv1
  etime = etime + dnekclock()
  deallocate(res1, res2, res3, h1, h2)

!    if (ifexplvis) call redo_split_vis

! Below is just for diagnostics...
  if (param(75) < 1) then

    allocate(DV2 (LX1,LY1,LZ1,LELV)) 
    ! Calculate Divergence norms of new VX,VY,VZ
    allocate(dvc(lx1,ly1,lz1,lelv), dfc(lx1,ly1,lz1,lelv))
    CALL OPDIV   (DVC,VX,VY,VZ)
    CALL DSSUM   (DVC)
    dvc = dvc * binvm1
 
    dv1 = dvc * bm1
    DIV1 = GLSUM (DV1,NTOT1)/VOLVM1
 
    dv2 = dvc * dvc * bm1
    DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
    DIV2 = SQRT  (DIV2)
    !  Calculate Divergence difference norms
    dfc = dvc! - qtl
    dv1 = dfc * bm1
    DIF1 = GLSUM (DV1,NTOT1)/VOLVM1
      
    dv2 = dfc * dfc * bm1
    DIF2 = GLSUM (DV2,NTOT1)/VOLVM1
    DIF2 = SQRT  (DIF2)
 
    dv1 = 0._dp !qtl * bm1
    QTL1 = GLSUM (DV1,NTOT1)/VOLVM1
      
    dv2 = 0._dp !qtl * qtl * bm1
    QTL2 = GLSUM (DV2,NTOT1)/VOLVM1
    QTL2 = SQRT  (QTL2)
 
    IF (NID == 0) THEN
        WRITE(6,'(15X,A,1p2e13.4)') &
        'L1/L2 DIV(V)    :',DIV1,DIV2
        WRITE(6,'(15X,A,1p2e13.4)') &
        'L1/L2 QTL       :',QTL1,QTL2
        WRITE(6,'(15X,A,1p2e13.4)') &
        'L1/L2 DIV(V)-QTL:',DIF1,DIF2
        IF (DIF2 > 0.1) WRITE(6,'(15X,A)') &
        'WARNING: DIV(V)-QTL too large!'
    ENDIF 
  endif

  tp4misc = tp4misc + (dnekclock() - etime)
   
  return
end subroutine plan4

!-----------------------------------------------------------------------
!> \brief Compute start residual/right-hand-side in the pressure
!!
!! \details \f$ R = \vec{F_b} 
!!   - \nabla \times \nabla \times \vec{v} 
!!   - \nabla \left(P + \frac{1}{3} H_1 \nabla \cdot Q \right) \f$
!!
!!   where \f$ F_b \f$ is the boundary force, 
!!   \f$ A \f$ is the stiffness matrix, 
!!   \f$ B \f$ is the mass matrix, 
!!   \f$ \vec{v} \f$ is the velocity, 
!!   \f$ P \f$ is the pressure, and
!!   \f$ \nabla \cdot Q \f$ is the "thermal divergence"
subroutine crespsp (respr, vext)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, lelv
  use size_m, only : nx1, ny1, nz1, nelv, ndim, nx2, ny2, nz2
  use geom, only : rxm2, sxm2, txm2, rym2, sym2, tym2, rzm2, szm2, tzm2
  use geom, only : unx, uny, unz, area
  use input, only : ifaxis, cbc
  use geom, only : bm1, binvm1, jacmi
  use soln, only : bfx, bfy, bfz, pr
  use soln, only : vx, vy, vz
  use soln, only : vtrans, vdiff, vmult!, qtl
  use soln, only : v1mask, v2mask, v3mask
  use tstep, only : imesh, bd, dt, ifield
  use mesh, only : if_ortho
  use dxyz, only : dxtm12, dym12, dzm12
  use ctimer, only : ncrespsp, tcrespsp, dnekclock
  use parallel, only : nid
  implicit none

  REAL(DP), intent(out) :: RESPR (LX2,LY2,LZ2,LELV)
  real(DP), intent(in)  :: VEXT  (LX1*LY1*LZ1*LELV,3)

  real(DP), allocatable, dimension(:,:,:,:) :: TA1, TA2, TA3
  real(DP), allocatable, dimension(:,:,:,:) :: WA1, WA2, WA3
  real(DP), allocatable, dimension(:,:,:,:) :: W1,  W2

  CHARACTER(3) :: CB
 
  integer :: nxyz1, ntot1, nfaces       
  integer :: n, ifc, iel, e, iz
  integer :: nyz2, nxy1
  real(DP) :: scale, dtbd
  real(DP) :: etime
  real(DP), allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)
  logical :: any_face

  ncrespsp = ncrespsp + 1
  etime = dnekclock()

  allocate(TA1 (LX1,LY1,LZ1,LELV) &
  , TA2 (LX1,LY1,LZ1,LELV) &
  , TA3 (LX1,LY1,LZ1,LELV) )
  allocate(W1 (LX1,LY1,LZ1,LELV) &
  , W2  (LX1,LY1,LZ1,LELV) )

 
  NXYZ1  = NX1*NY1*NZ1
  NTOT1  = NXYZ1*NELV
  NFACES = 2*NDIM

!   -mu*curl(curl(v))
  call op_curl (ta1,ta2,ta3,vext(1,1),vext(1,2),vext(1,3), &
   .TRUE. ,w1,w2)
  if(IFAXIS) then
!max      CALL COL2 (TA2, OMASK,NTOT1)
!max      CALL COL2 (TA3, OMASK,NTOT1)
  endif

  allocate( WA1 (LX1,LY1,LZ1,LELV) &
  , WA2 (LX1,LY1,LZ1,LELV) &
  , WA3 (LX1,LY1,LZ1,LELV) )
  call op_curl  (wa1,wa2,wa3,ta1,ta2,ta3, .TRUE. ,w1,w2)
  deallocate(w1, w2)

  allocate(tmp1(lx1, ly1, lz1), tmp2(lx1,ly1,lz1), tmp3(lx1,ly1,lz1))

  if(IFAXIS) then
!max      CALL COL2  (WA2, OMASK,NTOT1)
!max      CALL COL2  (WA3, OMASK,NTOT1)
  endif
  !wa1 = wa1 * bm1; wa2 = wa2 * bm1; wa3 = wa3 * bm1

  !call opgrad   (ta1,ta2,ta3,QTL)
  !ta1 = 0._dp; ta2 = 0._dp; ta3 = 0._dp
  if(IFAXIS) then
!max      CALL COL2  (ta2, OMASK,ntot1)
!max      CALL COL2  (ta3, OMASK,ntot1)
  endif

!   add old pressure term because we solve for delta p
  ta2 = 0._dp
  ta1 = 1._dp / vtrans(:,:,:,:,1)

  etime = etime - dnekclock()
  CALL AXHELM  (RESPR,PR,TA1,TA2,IMESH,1)
  etime = etime + dnekclock()

  !   add explicit (NONLINEAR) terms
  n = nx1*ny1*nz1*nelv
  do iel = 1, nelv
    !tmp1 = vdiff(:,:,:,iel,1) *ta1(:,:,:,iel)
    tmp1 = bm1(:,:,:,iel) * vdiff(:,:,:,iel,1) *ta1(:,:,:,iel)
    ta3(:,:,:,iel) = binvm1(:,:,:,iel) * (bfz(:,:,:,iel)*ta1(:,:,:,iel)-tmp1*wa3(:,:,:,iel))
    ta2(:,:,:,iel) = binvm1(:,:,:,iel) * (bfy(:,:,:,iel)*ta1(:,:,:,iel)-tmp1*wa2(:,:,:,iel))
    ta1(:,:,:,iel) = binvm1(:,:,:,iel) * (bfx(:,:,:,iel)*ta1(:,:,:,iel)-tmp1*wa1(:,:,:,iel))
  enddo

  call opdssum (ta1,ta2,ta3)

  if (if_ortho) then
    !deallocate(wa1,wa2,wa3)

    do e = 1, nelv
      call div_diag(1.0_dp, -1.0_dp, nx1, ny1, nz1, &
        jacmi(:,:,:,e) * bm1(:,:,:,e),  &
        ta1(:,:,:,e), rxm2(:,:,:,e), &
        ta2(:,:,:,e), sym2(:,:,:,e), &
        ta3(:,:,:,e), tzm2(:,:,:,e), & 
        respr(:,:,:,e), tmp2, tmp3)
    enddo
  else
    call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
    call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
    call cdtp    (wa3,ta3,rzm2,szm2,tzm2,1)
    respr = -respr + wa1 + wa2 + wa3
    !deallocate(wa1,wa2,wa3)
  endif
  deallocate(wa1,wa2,wa3)

!   add thermal divergence
  dtbd = BD(1)/DT
  respr = respr !+ qtl * bm1 * dtbd
   
!   surface terms
  DO IEL=1,NELV
      DO IFC=1,NFACES
          CB = CBC(IFC,IEL,IFIELD)
          any_face = .false.
          IF (CB(1:1) == 'V' .OR. CB(1:1) == 'v') THEN
              any_face = .true.
              tmp1 = 0._dp; tmp2 = 0._dp
              IF (NDIM == 3) tmp3 = 0._dp
              CALL FACCL3 &
              (tmp1,VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
              CALL FACCL3 &
              (tmp2,VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
              IF (NDIM == 3) &
              CALL FACCL3 &
              (tmp3,VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
          else if (cb(1:3) == 'SYM') then
            any_face = .true.
            tmp1 = 0._dp; tmp2 = 0._dp
            IF (NDIM == 3) tmp3 = 0._dp

            CALL FACCL3 (tmp1,ta1(:,:,:,iel),UNX(1,1,IFC,IEL),IFC)
            CALL FACCL3 (tmp2,ta2(:,:,:,iel),UNY(1,1,IFC,IEL),IFC)
            IF (NDIM == 3) then
              CALL FACCL3 (tmp3,ta3(:,:,:,iel),UNZ(1,1,IFC,IEL),IFC)
            endif
          endif

          if (any_face) then
            IF (NDIM == 3) then 
              tmp1 = tmp1 + tmp2 + tmp3
            else
              tmp1 = tmp1 + tmp2
            endif
            CALL FACCL2 (tmp1,AREA(1,1,IFC,IEL),IFC)
          endif

          if (CB(1:1) == 'V' .OR. CB(1:1) == 'v') then
            respr(:,:,:,iel) = respr(:,:,:,iel) - dtbd*tmp1
          else if (cb(1:3) == 'SYM') then
            respr(:,:,:,iel) = respr(:,:,:,iel) - tmp1 
          endif
      END DO
      !call dssum(ta1) ! maybe this should be here for SYM? (maxhutch)
  END DO
  deallocate(ta1, ta2, ta3)

!   Assure that the residual is orthogonal to (1,1,...,1)T
!   (only if all Dirichlet b.c.)
  CALL ORTHO (RESPR)
  tcrespsp = tcrespsp + (dnekclock() - etime)

  return
end subroutine crespsp

!----------------------------------------------------------------------
!> \brief Compute the residual for the velocity
subroutine cresvsp (resv1,resv2,resv3,h1,h2)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv
  use size_m, only : lx1, ly1, lz1, lelv
  use input, only : ifaxis
  use soln, only : vx, vy, vz, pr, bfx, bfy, bfz !, qtl
  use ctimer, only : ncresvsp, tcresvsp, dnekclock
  implicit none

  real(DP), intent(out) :: resv1(lx1,ly1,lz1,lelv) 
  real(DP), intent(out) :: resv2(lx1,ly1,lz1,lelv) 
  real(DP), intent(out) :: resv3(lx1,ly1,lz1,lelv) 
  real(DP), intent(out)  :: h1   (lx1,ly1,lz1,lelv) 
  real(DP), intent(out)  :: h2   (lx1,ly1,lz1,lelv)

  real(DP), allocatable, dimension(:,:,:,:) :: TA1, TA2, TA3, TA4

  integer :: ntot, intype
  real(DP) :: scale
  real(DP) :: etime

  etime = dnekclock()
  ncresvsp = ncresvsp + 1

  NTOT = NX1*NY1*NZ1*NELV
  INTYPE = -1

  CALL SETHLM  (H1,H2,INTYPE)

  etime = etime - dnekclock()
  ! just three axhelm calls
  CALL OPHX    (RESV1,RESV2,RESV3,VX,VY,VZ,H1,H2)
  etime = etime + dnekclock()

  scale = -1./3.
  allocate(TA1(LX1,LY1,LZ1,LELV) &
  ,             TA2   (LX1,LY1,LZ1,LELV) &
  ,             TA3   (LX1,LY1,LZ1,LELV) &
  ,             TA4   (LX1,LY1,LZ1,LELV) )
  ta4 = pr !+ scale*(vdiff(:,:,:,:,1) * qtl)

  call opgrad  (ta1,ta2,ta3,TA4)
  deallocate(ta4)

  if(IFAXIS) then
!max      CALL COL2 (TA2, OMASK,NTOT)
!max      CALL COL2 (TA3, OMASK,NTOT)
  endif

  resv1 = -resv1 + bfx - ta1
  resv2 = -resv2 + bfy - ta2
  resv3 = -resv3 + bfz - ta3

  tcresvsp = tcresvsp + (dnekclock() - etime)

  return
end subroutine cresvsp

!-----------------------------------------------------------------------
!> \brief Extrapolate the velocity forward in time with AB(k)
!!
!! \details This is the first half of (6.5.8) in HOMfIFF \cite HOMfIFF :
!! \f$ \hat{v} = \sum_{j=1}^k \beta_{k-j} v^{n+1-j} + ... \f$
subroutine v_extrap(vext)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv
  use soln, only : vx, vy, vz, vxlag, vylag, vzlag
  use tstep, only : ab, nab
  implicit none

  !> velocity (3 components) extrapolated forward in time
  real(DP), intent(out) :: vext(lx1,ly1,lz1,lelv,3) 
  real(DP) :: AB0, AB1, AB2

  AB0 = AB(1)
  AB1 = AB(2)
  AB2 = AB(3)

  if(nab == 3) then
             vext(:,:,:,:,1) = ab0*vx + ab1*vxlag(:,:,:,:,1) + ab2*vxlag(:,:,:,:,2)
             vext(:,:,:,:,2) = ab0*vy + ab1*vylag(:,:,:,:,1) + ab2*vylag(:,:,:,:,2)
             vext(:,:,:,:,3) = ab0*vz + ab1*vzlag(:,:,:,:,1) + ab2*vzlag(:,:,:,:,2)
  else
             vext(:,:,:,:,1) = ab0*vx + ab1*vxlag(:,:,:,:,1)
             vext(:,:,:,:,2) = ab0*vy + ab1*vylag(:,:,:,:,1)
             vext(:,:,:,:,3) = ab0*vz + ab1*vzlag(:,:,:,:,1)
  endif

  return
end subroutine v_extrap

!-----------------------------------------------------------------------
