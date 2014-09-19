!> \file plan4.F90
!! \brief Routine plan4 and supporting routines
!!
!! This file contains high-level routines for computing the 
!! right hand sides for the Helmholtz and Poisson solves in the
!! P(N)-P(N) formalism (colocation). 

!-----------------------------------------------------------------------
!> \brief Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
!!
!! NOTE: QTL denotes the so called thermal
!!       divergence and has to be provided
!!       by an external subroutine e.g qthermal
subroutine plan4
  use kinds, only : DP
  use size_m, only : nid, mxprev
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : lx2, ly2, lz2, lelv
  use size_m, only : nx1, ny1, nz1, nelv
  use ctimer, only : icalld, tpres, npres, etime1, dnekclock
  use geom, only : binvm1, bm1, volvm1
  use helmholtz, only : hsolve
  use soln, only : qtl, vx, vy, vz, v1mask, v2mask, v3mask
  use soln, only : vtrans, pmask, vmult, pr
  use tstep, only : imesh, nmxh, tolhv
  implicit none

  integer, parameter :: ktot = lx1*ly1*lz1*lelt
  integer, parameter :: laxt = mxprev

  integer, save :: napprox(2)
  real(DP), allocatable, save :: approx(:,:)

  real(DP), allocatable :: RES1(:,:,:,:)
  real(DP), allocatable :: RES2(:,:,:,:)   
  real(DP), allocatable :: RES3(:,:,:,:)      
  real(DP), allocatable :: DV1 (:,:,:,:)   
  real(DP), allocatable :: DV2 (:,:,:,:)   
  real(DP), allocatable :: DV3 (:,:,:,:)   
  real(DP), allocatable :: RESPR (:,:,:,:)

  real(DP), allocatable :: h1(:,:,:,:), h2(:,:,:,:)
  real(DP), allocatable :: VEXT(:,:)
  REAL(DP), allocatable :: DPR(:,:,:,:)

  REAL(DP), allocatable :: DVC (:,:,:,:), DFC(:,:,:,:)
  REAL(DP) :: DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
  real(DP) :: tolspl
  real(DP), external :: glsum

  integer :: intype, ntot1 

  if (.not. allocated(approx)) then
    allocate(approx(ktot,0:laxt))
    approx = 0._dp
  endif

  INTYPE = -1
  NTOT1  = NX1*NY1*NZ1*NELV

  ! add user defined divergence to qtl
!max  call add2 (qtl,usrdiv,ntot1)

  allocate(vext(lx1*ly1*lz1*lelv,3))

  ! Time-advance velocity with AB(k)
  CALL V_EXTRAP(vext)

  ! compute explicit contributions (bf{x,y,z}) with BDF(k)
  CALL MAKEF()

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
  call crespsp  (respr, vext)
  deallocate(vext)

  allocate(h1(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv))
  h1 = 1_dp / vtrans(:,:,:,:,1)
  call rzero    (h2,ntot1)
  call ctolspl  (tolspl,respr)

  allocate(dpr(lx2,ly2,lz2,lelv))
  napprox(1) = laxt
  call hsolve ('PRES', dpr, respr, h1, h2, pmask, vmult &
  ,imesh,tolspl,nmxh,1 &
  ,approx,napprox,binvm1)

  deallocate(respr)
  call add2    (pr,dpr,ntot1)
  deallocate(dpr)
  call ortho   (pr)
#ifndef NOTIMER
  tpres=tpres+(dnekclock()-etime1)
#endif

!     Compute velocity
  allocate(RES1(lx1,ly1,lz1,lelv), &
           RES2(LX1,LY1,LZ1,LELV), &
           RES3(LX1,LY1,LZ1,LELV)  )

  call cresvsp (res1,res2,res3,h1,h2)

  allocate(DV1 (LX1,LY1,LZ1,LELV), &
           DV2 (LX1,LY1,LZ1,LELV), &
           DV3 (LX1,LY1,LZ1,LELV)  )
  call ophinv  (dv1,dv2,dv3,res1,res2,res3,h1,h2,tolhv,nmxh)
  deallocate(res1, res2, res3, h1, h2)

  vx = vx + dv1; vy = vy + dv2; vz = vz + dv3

!    if (ifexplvis) call redo_split_vis

! Below is just for diagnostics...

!     Calculate Divergence norms of new VX,VY,VZ
  allocate(dvc(lx1,ly1,lz1,lelv), dfc(lx1,ly1,lz1,lelv))
  CALL OPDIV   (DVC,VX,VY,VZ)
  CALL DSSUM   (DVC)
  dvc = dvc * binvm1

  dv1 = dvc * bm1
  DIV1 = GLSUM (DV1,NTOT1)/VOLVM1

  dv2 = dvc * dvc * bm1
  DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
  DIV2 = SQRT  (DIV2)
!     Calculate Divergence difference norms
  dfc = dvc - qtl
  dv1 = dfc * bm1
  DIF1 = GLSUM (DV1,NTOT1)/VOLVM1
    
  dv2 = dfc * dfc * bm1
  DIF2 = GLSUM (DV2,NTOT1)/VOLVM1
  DIF2 = SQRT  (DIF2)

  dv1 = qtl * bm1
  QTL1 = GLSUM (DV1,NTOT1)/VOLVM1
    
  dv2 = qtl * qtl * bm1
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
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use geom, only : rxm2, sxm2, txm2, rym2, sym2, tym2, rzm2, szm2, tzm2
  use geom, only : area
  use input, only : ifaxis, if3d, cbc
  use geom, only : bm1, binvm1
  use soln, only : bfx, bfy, bfz, pr
  use soln, only : vtrans, vdiff, qtl
  use tstep, only : imesh, bd, dt, ifield
  implicit none

  REAL(DP), intent(out) :: RESPR (LX2,LY2,LZ2,LELV)
  real(DP), intent(in)  :: VEXT  (LX1*LY1*LZ1*LELV,3)

  real(DP), allocatable, dimension(:,:,:,:) :: TA1, TA2, TA3
  real(DP), allocatable, dimension(:,:,:,:) :: WA1, WA2, WA3
  real(DP), allocatable, dimension(:,:,:,:) :: W1,  W2

  CHARACTER(3) :: CB
 
  integer :: nxyz1, ntot1, nfaces       
  integer :: n, ifc, iel
  real(DP) :: scale, dtbd

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
  deallocate(w2)

  if(IFAXIS) then
!max      CALL COL2  (WA2, OMASK,NTOT1)
!max      CALL COL2  (WA3, OMASK,NTOT1)
  endif
  wa1 = wa1 * bm1; wa2 = wa2 * bm1; wa3 = wa3 * bm1

  call opgrad   (ta1,ta2,ta3,QTL)
  if(IFAXIS) then
!max      CALL COL2  (ta2, OMASK,ntot1)
!max      CALL COL2  (ta3, OMASK,ntot1)
  endif

  scale = -4./3.
  w1 = vdiff(:,:,:,:,1) / vtrans(:,:,:,:,1)
  wa1 = w1*(wa1 + scale*ta1)
  wa2 = w1*(wa2 + scale*ta2)
  wa3 = w1*(wa3 + scale*ta3)
  deallocate(w1)

!   add old pressure term because we solve for delta p
  ta1 = 1_dp / vtrans(:,:,:,:,1)
  CALL RZERO   (TA2,NTOT1)
  CALL AXHELM  (RESPR,PR,TA1,TA2,IMESH,1)
  CALL CHSIGN  (RESPR,NTOT1)

!   add explicit (NONLINEAR) terms
  n = nx1*ny1*nz1*nelv
  ta1 = bfx/vtrans(:,:,:,:,1)-wa1
  ta2 = bfy/vtrans(:,:,:,:,1)-wa2
  ta3 = bfz/vtrans(:,:,:,:,1)-wa3

  call opdssum (ta1,ta2,ta3)
  ta1 = ta1*binvm1
  ta2 = ta2*binvm1
  ta3 = ta3*binvm1

  if (if3d) then
      call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
      call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
      call cdtp    (wa3,ta3,rzm2,szm2,tzm2,1)
      respr = respr + wa1 + wa2 + wa3
  else
      call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
      call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
      respr = respr + wa1 + wa2 
  endif
  deallocate(wa1,wa2,wa3)

!   add thermal divergence
  dtbd = BD(1)/DT
  call admcol3(respr,QTL,bm1,dtbd,ntot1)
   
!   surface terms
  DO IFC=1,NFACES
      ta1 = 0._dp; ta2 = 0._dp
      IF (NDIM == 3) ta3 = 0._dp
      DO IEL=1,NELV
          CB = CBC(IFC,IEL,IFIELD)
          IF (CB(1:1) == 'V' .OR. CB(1:1) == 'v') THEN
            write(*,*) "Oops: cb = v"
#if 0
              CALL FACCL3 &
              (TA1(1,IEL),VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
              CALL FACCL3 &
              (TA2(1,IEL),VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
              IF (NDIM == 3) &
              CALL FACCL3 &
              (TA3(1,IEL),VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
#endif
          ENDIF
          CALL ADD2   (TA1(:,:,:,IEL),TA2(:,:,:,IEL),NXYZ1)
          IF (NDIM == 3) &
          CALL ADD2   (TA1(:,:,:,IEL),TA3(:,:,:,IEL),NXYZ1)
          CALL FACCL2 (TA1(:,:,:,IEL),AREA(1,1,IFC,IEL),IFC)
      END DO
      respr = respr - dtbd*ta1
  END DO
  deallocate(ta1, ta2, ta3)

!   Assure that the residual is orthogonal to (1,1,...,1)T
!   (only if all Dirichlet b.c.)
  CALL ORTHO (RESPR)

  return
end subroutine crespsp

!----------------------------------------------------------------------
!> \brief Compute the residual for the velocity
subroutine cresvsp (resv1,resv2,resv3,h1,h2)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv
  use size_m, only : lx1, ly1, lz1, lelv
  use input, only : ifaxis
  use soln, only : vx, vy, vz, vdiff, qtl, pr, bfx, bfy, bfz
  implicit none

  real(DP), intent(out) :: resv1(lx1,ly1,lz1,lelv) 
  real(DP), intent(out) :: resv2(lx1,ly1,lz1,lelv) 
  real(DP), intent(out) :: resv3(lx1,ly1,lz1,lelv) 
  real(DP), intent(out)  :: h1   (lx1,ly1,lz1,lelv) 
  real(DP), intent(out)  :: h2   (lx1,ly1,lz1,lelv)

  real(DP), allocatable, dimension(:,:,:,:) :: TA1, TA2, TA3, TA4

  integer :: ntot, intype
  real(DP) :: scale

  NTOT = NX1*NY1*NZ1*NELV
  INTYPE = -1

  CALL SETHLM  (H1,H2,INTYPE)

  CALL OPHX    (RESV1,RESV2,RESV3,VX,VY,VZ,H1,H2)
  CALL OPCHSGN (RESV1,RESV2,RESV3)

  scale = -1./3.
  allocate(TA1(LX1,LY1,LZ1,LELV) &
  ,             TA2   (LX1,LY1,LZ1,LELV) &
  ,             TA3   (LX1,LY1,LZ1,LELV) &
  ,             TA4   (LX1,LY1,LZ1,LELV) )
  ta4 = scale*(vdiff(:,:,:,:,1) * qtl) + pr 

  call opgrad  (ta1,ta2,ta3,TA4)
  deallocate(ta4)

  if(IFAXIS) then
!max      CALL COL2 (TA2, OMASK,NTOT)
!max      CALL COL2 (TA3, OMASK,NTOT)
  endif

  resv1 = resv1 + bfx - ta1
  resv2 = resv2 + bfy - ta2
  resv3 = resv3 + bfz - ta3

  return
end subroutine cresvsp

!-----------------------------------------------------------------------
!> \brief Compute curl of U.
!!
!! \f$ (w_1, w_2, w_3) = \nabla \times (u_1, u_2, u_3) \f$
subroutine op_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv
  use geom, only : rxm1, rym1, rzm1, sxm1, sym1, szm1, txm1, tym1, tzm1
  use geom, only : jacm1, bm1, binvm1
  use input, only : if3d, ifaxis, ifcyclic
  use tstep, only : ifield
  implicit none


  real(DP), intent(out) :: w1(lx1,ly1,lz1,lelv) !>!< 1st component of curl U
  real(DP), intent(out) :: w2(lx1,ly1,lz1,lelv) !>!< 2nd component of curl U
  real(DP), intent(out) :: w3(lx1,ly1,lz1,lelv) !>!< 3rd component of curl U
  real(DP), intent(in)  :: u1(*) !>!< 1st component of U
  real(DP), intent(in)  :: u2(*) !>!< 2nd component of U
  real(DP), intent(in)  :: u3(*) !>!< 3rd component of U
  real(DP), intent(out) :: work1(lx1,ly1,lz1,*) !>!< work array
  real(DP), intent(out) :: work2(lx1,ly1,lz1,*) !>!< work array
  logical, intent(in)   :: ifavg !>!< Average at boundary? 

  integer :: ntot, nxyz, ifielt

  ntot  = nx1*ny1*nz1*nelv
  nxyz  = nx1*ny1*nz1
!   work1=dw/dy ; work2=dv/dz
  call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2)
  if (if3d) then
      call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
      w1 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)
  else
      call copy(w1,work1,ntot)

      if(ifaxis) then
        write(*,*) "Oops: ifaxis"
#if 0
          call copy (ta,u3,ntot)
          do iel = 1,nelv
              if(IFRZER(iel)) then
                  call rzero (ta(1,1,1,iel),nx1)
                  call MXM   (ta(1,1,1,iel),nx1,DATM1,ny1,duax,1)
                  call copy  (ta(1,1,1,iel),duax,nx1)
              endif
              call col2    (ta(1,1,1,iel),yinvm1(1,1,1,iel),nxyz)
          enddo
          call add2      (w1,ta,ntot)
#endif
      endif
  endif
!   work1=du/dz ; work2=dw/dx
  if (if3d) then
      call dudxyz(work1,u1,rzm1,szm1,tzm1,jacm1,1,3)
      call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
      w2 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)
  else
      call rzero (work1,ntot)
      call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
      w2 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)
  endif
!   work1=dv/dx ; work2=du/dy
  call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1)
  call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
  w3 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)

!  Avg at bndry

!   if (ifavg) then
  if (ifavg .AND. .NOT. ifcyclic) then

      ifielt = ifield
      ifield = 1
             
      w1 = w1 * bm1; w2 = w2 * bm1; w3 = w3 * bm1
      call opdssum (w1,w2,w3)
      w1 = w1 * binvm1; w2 = w2 * binvm1; w3 = w3 * binvm1

      ifield = ifielt

  endif

  return
end subroutine op_curl

!-----------------------------------------------------------------------
!> \brief Extrapolate the velocity forward in time with AB(k)
!!
!! \details This is the first half of (6.5.8) in HOMfIFF:
!! \f$ \hat{v} = \sum_{j=1}^k \beta_{k-j} v^{n+1-j} + ... \f$
subroutine v_extrap(vext)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv
  use input, only : if3d
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
    if(if3d) vext(:,:,:,:,3) = ab0*vz + ab1*vzlag(:,:,:,:,1) + ab2*vzlag(:,:,:,:,2)
  else
             vext(:,:,:,:,1) = ab0*vx + ab1*vxlag(:,:,:,:,1)
             vext(:,:,:,:,2) = ab0*vy + ab1*vylag(:,:,:,:,1)
    if(if3d) vext(:,:,:,:,3) = ab0*vz + ab1*vzlag(:,:,:,:,1)
  endif

  return
end subroutine v_extrap

!-----------------------------------------------------------------------
