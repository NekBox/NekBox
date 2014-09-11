!> \file plan4.F90
!! \brief Routine plan4 and supporting routines

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
  use soln, only : qtl, usrdiv, vx, vy, vz, v1mask, v2mask, v3mask
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
  CALL V_EXTRAP(vext)

! compute explicit contributions bfx,bfy,bfz
  CALL MAKEF()
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
  call invers2  (h1,vtrans,ntot1)
  call rzero    (h2,ntot1)
  call ctolspl  (tolspl,respr)

  allocate(dpr(lx2,ly2,lz2,lelv))
  napprox(1) = laxt
  call hsolve   ('PRES',dpr,respr,h1,h2 &
  ,pmask,vmult &
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
  CALL DSSUM   (DVC,NX1,NY1,NZ1)
  CALL COL2    (DVC,BINVM1,NTOT1)

  CALL COL3    (DV1,DVC,BM1,NTOT1)
  DIV1 = GLSUM (DV1,NTOT1)/VOLVM1

  CALL COL3    (DV2,DVC,DVC,NTOT1)
  CALL COL2    (DV2,BM1   ,NTOT1)
  DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
  DIV2 = SQRT  (DIV2)
!     Calculate Divergence difference norms
  CALL SUB3    (DFC,DVC,QTL,NTOT1)
  CALL COL3    (DV1,DFC,BM1,NTOT1)
  DIF1 = GLSUM (DV1,NTOT1)/VOLVM1
    
  CALL COL3    (DV2,DFC,DFC,NTOT1)
  CALL COL2    (DV2,BM1   ,NTOT1)
  DIF2 = GLSUM (DV2,NTOT1)/VOLVM1
  DIF2 = SQRT  (DIF2)

  CALL COL3    (DV1,QTL,BM1,NTOT1)
  QTL1 = GLSUM (DV1,NTOT1)/VOLVM1
    
  CALL COL3    (DV2,QTL,QTL,NTOT1)
  CALL COL2    (DV2,BM1   ,NTOT1)
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
subroutine crespsp (respr, vext)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, lelv
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use geom, only : rxm2, sxm2, txm2, rym2, sym2, tym2, rzm2, szm2, tzm2
  use geom, only : area
  use input, only : ifaxis, if3d, cbc
  use geom, only : bm1, binvm1
  use soln, only : vx, vy, vz, bfx, bfy, bfz, pr
  use soln, only : vtrans, vdiff, qtl, omask
  use tstep, only : imesh, bd, dt, ifield
  implicit none

  REAL(DP) :: RESPR (LX2*LY2*LZ2*LELV)
  real(DP) :: VEXT  (LX1*LY1*LZ1*LELV,3)

  real(DP), allocatable :: TA1 (:,:), TA2(:,:), TA3(:,:)
  real(DP), allocatable :: WA1(:),    WA2(:),   WA3(:)
  real(DP), allocatable ::  W1(:),     W2(:)

  CHARACTER CB*3

  integer :: nxyz1, ntot1, nfaces       
  integer :: i, n, ifc, iel
  real(DP) :: scale, dtbd

  allocate(TA1 (LX1*LY1*LZ1,LELV) &
  , TA2 (LX1*LY1*LZ1,LELV) &
  , TA3 (LX1*LY1*LZ1,LELV) )
  allocate(W1 (LX1*LY1*LZ1*LELV) &
  , W2  (LX1*LY1*LZ1*LELV) )
 
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

  allocate( WA1 (LX1*LY1*LZ1*LELV) &
  , WA2 (LX1*LY1*LZ1*LELV) &
  , WA3 (LX1*LY1*LZ1*LELV) )
  call op_curl  (wa1,wa2,wa3,ta1,ta2,ta3, .TRUE. ,w1,w2)
  deallocate(w2)

  if(IFAXIS) then
!max      CALL COL2  (WA2, OMASK,NTOT1)
!max      CALL COL2  (WA3, OMASK,NTOT1)
  endif
  call opcolv   (wa1,wa2,wa3,bm1)

  call opgrad   (ta1,ta2,ta3,QTL)
  if(IFAXIS) then
!max      CALL COL2  (ta2, OMASK,ntot1)
!max      CALL COL2  (ta3, OMASK,ntot1)
  endif

  scale = -4./3.
  call opadd2cm (wa1,wa2,wa3,ta1,ta2,ta3,scale)
  call invcol3  (w1,vdiff,vtrans,ntot1)
  call opcolv   (wa1,wa2,wa3,w1)
  deallocate(w1)

!   add old pressure term because we solve for delta p
  CALL INVERS2 (TA1,VTRANS,NTOT1)
  CALL RZERO   (TA2,NTOT1)
  CALL AXHELM  (RESPR,PR,TA1,TA2,IMESH,1)
  CALL CHSIGN  (RESPR,NTOT1)

!   add explicit (NONLINEAR) terms
  n = nx1*ny1*nz1*nelv
  do i=1,n
      ta1(i,1) = bfx(i,1,1,1)/vtrans(i,1,1,1,1)-wa1(i)
      ta2(i,1) = bfy(i,1,1,1)/vtrans(i,1,1,1,1)-wa2(i)
      ta3(i,1) = bfz(i,1,1,1)/vtrans(i,1,1,1,1)-wa3(i)
  enddo
  call opdssum (ta1,ta2,ta3)
  do i=1,n
      ta1(i,1) = ta1(i,1)*binvm1(i,1,1,1)
      ta2(i,1) = ta2(i,1)*binvm1(i,1,1,1)
      ta3(i,1) = ta3(i,1)*binvm1(i,1,1,1)
  enddo
  if (if3d) then
      call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
      call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
      call cdtp    (wa3,ta3,rzm2,szm2,tzm2,1)
      do i=1,n
          respr(i) = respr(i)+wa1(i)+wa2(i)+wa3(i)
      enddo
  else
      call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
      call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
      do i=1,n
          respr(i) = respr(i)+wa1(i)+wa2(i)
      enddo
  endif
  deallocate(wa1,wa2,wa3)

!   add thermal divergence
  dtbd = BD(1)/DT
  call admcol3(respr,QTL,bm1,dtbd,ntot1)
   
!   surface terms
  DO IFC=1,NFACES
      CALL RZERO  (TA1,NTOT1)
      CALL RZERO  (TA2,NTOT1)
      IF (NDIM == 3) &
      CALL RZERO  (TA3,NTOT1)
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
          CALL ADD2   (TA1(1,IEL),TA2(1,IEL),NXYZ1)
          IF (NDIM == 3) &
          CALL ADD2   (TA1(1,IEL),TA3(1,IEL),NXYZ1)
          CALL FACCL2 (TA1(1,IEL),AREA(1,1,IFC,IEL),IFC)
      END DO
      CALL CMULT(TA1,dtbd,NTOT1)
      CALL SUB2 (RESPR,TA1,NTOT1)
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
  use soln, only : vx, vy, vz, vdiff, qtl, pr, omask, bfx, bfy, bfz
  implicit none

  real(DP) :: resv1(lx1,ly1,lz1,lelv) &
  , resv2(lx1,ly1,lz1,lelv) &
  , resv3(lx1,ly1,lz1,lelv) &
  , h1   (lx1,ly1,lz1,lelv) &
  , h2   (lx1,ly1,lz1,lelv)

  real(DP), allocatable :: TA1(:,:,:,:), TA2(:,:,:,:) &
                         , TA3(:,:,:,:), TA4(:,:,:,:)

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
subroutine op_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv
  use dxyz, only : datm1
  use geom, only : rxm1, rym1, rzm1, sxm1, sym1, szm1, txm1, tym1, tzm1
  use geom, only : jacm1, ifrzer
  use input, only : if3d, ifaxis, ifcyclic
  use geom, only : yinvm1, bm1, binvm1
  use tstep, only : ifield
  implicit none

  real(DP) :: w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
  logical :: ifavg

  real :: duax(lx1)!, ta(lx1,ly1,lz1,lelv)
  integer :: ntot, nxyz, iel, ifielt

  ntot  = nx1*ny1*nz1*nelv
  nxyz  = nx1*ny1*nz1
!   work1=dw/dy ; work2=dv/dz
  call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2)
  if (if3d) then
      call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
      call sub3(w1,work1,work2,ntot)
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
      call sub3(w2,work1,work2,ntot)
  else
      call rzero (work1,ntot)
      call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
      call sub3(w2,work1,work2,ntot)
  endif
!   work1=dv/dx ; work2=du/dy
  call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1)
  call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
  call sub3(w3,work1,work2,ntot)

!  Avg at bndry

!   if (ifavg) then
  if (ifavg .AND. .NOT. ifcyclic) then

      ifielt = ifield
      ifield = 1
             
      call opcolv  (w1,w2,w3,bm1)
      call opdssum (w1,w2,w3)
      call opcolv  (w1,w2,w3,binvm1)

      ifield = ifielt

  endif

  return
end subroutine op_curl

!-----------------------------------------------------------------------
subroutine opadd2cm (a1,a2,a3,b1,b2,b3,c)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  implicit none

  REAL(DP) :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1),C
  integer :: ntot1, i

  NTOT1=NX1*NY1*NZ1*NELV
  if (ndim == 3) then
      do i=1,ntot1
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
          a3(i) = a3(i) + b3(i)*c
      enddo
  else
      do i=1,ntot1
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
      enddo
  endif
  return
end subroutine opadd2cm

!-----------------------------------------------------------------------
!> \brief extrapolate velocity
subroutine v_extrap(vext)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv
  use input, only : if3d
  use soln, only : vx, vy, vz, vxlag, vylag, vzlag
  use tstep, only : ab, nab
  implicit none
       
  real(DP) :: vext(lx1*ly1*lz1*lelv,3)

  integer :: ntot
  real(DP) :: AB0, AB1, AB2

  NTOT = NX1*NY1*NZ1*NELV

  AB0 = AB(1)
  AB1 = AB(2)
  AB2 = AB(3)

!   call copy(vext(1,1),vx,ntot)
!   call copy(vext(1,2),vy,ntot)
!   call copy(vext(1,3),vz,ntot)
!   return

  call add3s2(vext(1,1),vx,vxlag,ab0,ab1,ntot)
  call add3s2(vext(1,2),vy,vylag,ab0,ab1,ntot)
  if(if3d) call add3s2(vext(1,3),vz,vzlag,ab0,ab1,ntot)

  if(nab == 3) then
      call add2s2(vext(1,1),vxlag(1,1,1,1,2),ab2,ntot)
      call add2s2(vext(1,2),vylag(1,1,1,1,2),ab2,ntot)
      if(if3d) call add2s2(vext(1,3),vzlag(1,1,1,1,2),ab2,ntot)
  endif

  return
end subroutine v_extrap

!-----------------------------------------------------------------------
