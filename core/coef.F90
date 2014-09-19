!-----------------------------------------------------------------
!> \brief Generate derivative and interpolation operators.
!!  GENERATE
!!         - DERIVATIVE OPERATORS
!!         - INTERPOLATION OPERATORS
!!         - WEIGHTS
!!         - COLLOCATION POINTS
!!  ASSOCIATED WITH THE
!!         - GAUSS-LOBATTO LEGENDRE MESH (SUFFIX M1/M2/M3)
!!         - GAUSS LEGENDRE         MESH (SUFFIX M2)
!!         - GAUSS-LOBATTO JACOBI   MESH (SUFFIX M1/M2/M3)
!-----------------------------------------------------------------
subroutine genwz
  use size_m, only : ndim, nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3
  use dxyz, only : dxm1, dxtm1, dym1, dytm1, dzm1, dztm1, dxm3, dxtm3
  use dxyz, only : dym3, dytm3, dzm3, dztm3, dxm12, dxtm12, dym12, dytm12
  use dxyz, only : dzm12, dztm12
  use input, only : ifsplit
  use ixyz, only : ixm12, ixtm12, iym12, iytm12, izm12, izm12
  use ixyz, only : iztm12, ixm21, ixtm21, iym21, iytm21, izm21, iztm21
  use ixyz, only : ixm13, ixtm13, iym13, iytm13, izm13, iztm13
  use ixyz, only : ixm31, ixtm31, iym31, iytm31, izm31, iztm31
  use wz_m, only : zgm1, zgm2, zgm3
  use wz_m, only : wxm1, wym1, wzm1, wxm2, wym2, wzm2, wxm3, wym3, wzm3
  use wz_m, only : w3m1, w3m2, w3m3
  implicit none

  integer :: ix, iy, iz

  IF (NDIM == 2) THEN
    write (*,*) "Oops: ndim == 2"
#if 0    
  !***  Two-dimensional case  **********************
  
  
  !     Gauss-Lobatto Legendre mesh (suffix M1)
  !     Generate collocation points and weights
  
      CALL ZWGLL (ZGM1(1,1),WXM1,NX1)
      CALL ZWGLL (ZGM1(1,2),WYM1,NY1)
      ZGM1(NZ1,3) = 0.
      WZM1(NZ1)   = 1.
      DO 100 IY=1,NY1
          DO 100 IX=1,NX1
              W3M1(IX,IY,1)=WXM1(IX)*WYM1(IY)
      100 END DO
  
  !     Compute derivative matrices
  
      CALL DGLL (DXM1,DXTM1,ZGM1(1,1),NX1,NX1)
      CALL DGLL (DYM1,DYTM1,ZGM1(1,2),NY1,NY1)
      CALL RZERO (DZM1 ,NZ1*NZ1)
      CALL RZERO (DZTM1,NZ1*NZ1)
  
  !     Gauss Legendre mesh (suffix M2)
  !     Generate collocation points and weights
  
      IF(IFSPLIT)THEN
          CALL ZWGLL (ZGM2(1,1),WXM2,NX2)
          CALL ZWGLL (ZGM2(1,2),WYM2,NY2)
      ELSE
          CALL ZWGL  (ZGM2(1,1),WXM2,NX2)
          CALL ZWGL  (ZGM2(1,2),WYM2,NY2)
      ENDIF
      ZGM2(NZ2,3) = 0.
      WZM2(NZ2)   = 1.
      DO 200 IY=1,NY2
          DO 200 IX=1,NX2
              W3M2(IX,IY,1)=WXM2(IX)*WYM2(IY)
      200 END DO
  
  !     Gauss-Lobatto Legendre mesh (suffix M3).
  !     Generate collocation points and weights.
  
      CALL ZWGLL (ZGM3(1,1),WXM3,NX3)
      CALL ZWGLL (ZGM3(1,2),WYM3,NY3)
      ZGM3(NZ3,3) = 0.
      WZM3(NZ3)   = 1.
      DO 300 IY=1,NY3
          DO 300 IX=1,NX3
              W3M3(IX,IY,1)=WXM3(IX)*WYM3(IY)
      300 END DO
  
  !     Compute derivative matrices
  
      CALL DGLL (DXM3,DXTM3,ZGM3(1,1),NX3,NX3)
      CALL DGLL (DYM3,DYTM3,ZGM3(1,2),NY3,NY3)
      CALL RZERO (DZM3 ,NZ3*NZ3)
      CALL RZERO (DZTM3,NZ3*NZ3)
  
  !     Generate interpolation operators for the staggered mesh
  
      CALL IGLLM (IXM12,IXTM12,ZGM1(1,1),ZGM2(1,1),NX1,NX2,NX1,NX2)
      CALL IGLLM (IYM12,IYTM12,ZGM1(1,2),ZGM2(1,2),NY1,NY2,NY1,NY2)
      IZM12 (NZ2,NZ1) = 1.
      IZTM12(NZ1,NZ2) = 1.
  
  !     NOTE: The splitting scheme has only one mesh!!!!!
  
      IF (IFSPLIT) THEN
          CALL IGLLM (IXM21,IXTM21,ZGM1(1,1),ZGM2(1,1),NX1,NX2,NX1,NX2)
          CALL IGLLM (IYM21,IYTM21,ZGM1(1,2),ZGM2(1,2),NY1,NY2,NY1,NY2)
      ELSE
        write(*,*) "Oops: ifsplit false"
#if 0
          CALL IGLM  (IXM21,IXTM21,ZGM2(1,1),ZGM1(1,1),NX2,NX1,NX2,NX1)
          CALL IGLM  (IYM21,IYTM21,ZGM2(1,2),ZGM1(1,2),NY2,NY1,NY2,NY1)
#endif
      ENDIF
      IZM21 (NZ1,NZ2) = 1.
      IZTM21(NZ2,NZ1) = 1.
  
  !     Compute derivative operators for the staggered mesh
  
      IF(IFSPLIT)THEN
          CALL COPY (DXM12, DXM1, NX1*NX2)
          CALL COPY (DXTM12,DXTM1,NX1*NX2)
          CALL COPY (DYM12, DYM1, NY1*NY2)
          CALL COPY (DYTM12,DYTM1,NY1*NY2)
          CALL COPY (DZM12, DZM1, NZ1*NZ2)
          CALL COPY (DZTM12,DZTM1,NZ1*NZ2)
      ELSE
          CALL DGLLGL (DXM12,DXTM12,ZGM1(1,1),ZGM2(1,1),IXM12, &
          NX1,NX2,NX1,NX2)
          CALL DGLLGL (DYM12,DYTM12,ZGM1(1,2),ZGM2(1,2),IYM12, &
          NY1,NY2,NY1,NY2)
          DZM12 (NZ2,NZ1) = 0.
          DZTM12(NZ2,NZ1) = 0.
      ENDIF
  
  !     Compute interpolation operators for the geometry mesh M3.
  
      CALL IGLLM (IXM13,IXTM13,ZGM1(1,1),ZGM3(1,1),NX1,NX3,NX1,NX3)
      CALL IGLLM (IYM13,IYTM13,ZGM1(1,2),ZGM3(1,2),NY1,NY3,NY1,NY3)
      CALL IGLLM (IXM31,IXTM31,ZGM3(1,1),ZGM1(1,1),NX3,NX1,NX3,NX1)
      CALL IGLLM (IYM31,IYTM31,ZGM3(1,2),ZGM1(1,2),NY3,NY1,NY3,NY1)
      IZM13 (NZ3,NZ1) = 1.
      IZTM13(NZ1,NZ3) = 1.
      IZM31 (NZ1,NZ3) = 1.
      IZTM31(NZ3,NZ1) = 1.
  
  
      IF (IFAXIS) THEN
      
      !     Special treatment for the axisymmetric case
      !     Generate additional points, weights, derivative operators and
      !     interpolation operators required for elements close to the axis.
      
      
      !     Gauss-Lobatto Jacobi mesh (suffix M1).
      !     Generate collocation points and weights (alpha=0, beta=1).
      
          ALPHA = 0.
          BETA  = 1.
          CALL ZWGLJ (ZAM1,WAM1,NY1,ALPHA,BETA)
          DO 400 IY=1,NY1
              DO 400 IX=1,NX1
                  W2AM1(IX,IY)=WXM1(IX)*WAM1(IY)
                  W2CM1(IX,IY)=WXM1(IX)*WYM1(IY)
          400 END DO
      
      !     Compute derivative matrices
      
          CALL COPY (DCM1,DYM1,NY1*NY1)
          CALL COPY (DCTM1,DYTM1,NY1*NY1)
          CALL DGLJ (DAM1,DATM1,ZAM1,NY1,NY1,ALPHA,BETA)
      
      !     Gauss Jacobi mesh (suffix M2)
      !     Generate collocation points and weights
      
          IF(IFSPLIT)THEN
              CALL ZWGLJ (ZAM2,WAM2,NY2,ALPHA,BETA)
          ELSE
              CALL ZWGJ  (ZAM2,WAM2,NY2,ALPHA,BETA)
          ENDIF
          DO 500 IY=1,NY2
              DO 500 IX=1,NX2
                  W2CM2(IX,IY)=WXM2(IX)*WYM2(IY)
                  W2AM2(IX,IY)=WXM2(IX)*WAM2(IY)
          500 END DO
      
      !     Gauss-Lobatto Jacobi mesh (suffix M3).
      !     Generate collocation points and weights.
      
          CALL ZWGLJ (ZAM3,WAM3,NY3,ALPHA,BETA)
          DO 600 IY=1,NY3
              DO 600 IX=1,NX3
                  W2CM3(IX,IY)=WXM3(IX)*WYM3(IY)
                  W2AM3(IX,IY)=WXM3(IX)*WAM3(IY)
          600 END DO
      
      !     Compute derivative matrices
      
          CALL COPY (DCM3,DYM3,NY3*NY3)
          CALL COPY (DCTM3,DYTM3,NY3*NY3)
          CALL DGLJ (DAM3,DATM3,ZAM3,NY3,NY3,ALPHA,BETA)
      
      !     Generate interpolation operators for the staggered mesh
      
          CALL COPY  (ICM12,IYM12,NY2*NY1)
          CALL COPY  (ICTM12,IYTM12,NY1*NY2)
          CALL IGLJM (IAM12,IATM12,ZAM1,ZAM2,NY1,NY2,NY1,NY2,ALPHA,BETA)
          CALL COPY  (ICM21,IYM21,NY1*NY2)
          CALL COPY  (ICTM21,IYTM21,NY2*NY1)
          IF (IFSPLIT) THEN
              CALL IGLJM (IAM21,IATM21,ZAM2,ZAM1,NY1,NY2,NY1,NY2,ALPHA,BETA)
          ELSE
              CALL IGJM  (IAM21,IATM21,ZAM2,ZAM1,NY2,NY1,NY2,NY1,ALPHA,BETA)
          ENDIF
      
      !     Compute derivative operators for the staggered mesh
      
          CALL COPY  (DCM12,DYM12,NY2*NY1)
          CALL COPY  (DCTM12,DYTM12,NY1*NY2)
          IF(IFSPLIT)THEN
              CALL COPY (DAM12, DAM1, NY1*NY2)
              CALL COPY (DATM12,DATM1,NY1*NY2)
          ELSE
              CALL DGLJGJ (DAM12,DATM12,ZAM1,ZAM2,IAM12, &
              NY1,NY2,NY1,NY2,ALPHA,BETA)
          ENDIF
      
      !     Compute interpolation operators for the geometry mesh M3.
      
          CALL COPY  (ICM13,IYM13,NY3*NY1)
          CALL COPY  (ICTM13,IYTM13,NY1*NY3)
          CALL IGLJM (IAM13,IATM13,ZAM1,ZAM3,NY1,NY3,NY1,NY3,ALPHA,BETA)
          CALL COPY  (ICM31,IYM31,NY1*NY3)
          CALL COPY  (ICTM31,IYTM31,NY3*NY1)
          CALL IGLJM (IAM31,IATM31,ZAM3,ZAM1,NY3,NY1,NY3,NY1,ALPHA,BETA)
      
      !     Compute interpolation operators between Gauss-Lobatto Jacobi
      !     and Gauss-Lobatto Legendre (to be used in PREPOST).
      
          CALL IGLJM(IAJL1,IATJL1,ZAM1,ZGM1(1,2),NY1,NY1,NY1,NY1,ALPHA,BETA)
          IF (IFSPLIT) THEN
              CALL IGLJM(IAJL2,IATJL2,ZAM2,ZGM2(1,2),NY2,NY2,NY2,NY2,ALPHA,BETA)
          ELSE
              CALL IGJM (IAJL2,IATJL2,ZAM2,ZGM2(1,2),NY2,NY2,NY2,NY2,ALPHA,BETA)
          ENDIF

          CALL INVMT(IAJL1 ,IALJ1 ,TMP ,NY1)
          CALL INVMT(IATJL1,IATLJ1,TMPT,NY1)
          CALL MXM (IATJL1,NY1,IATLJ1,NY1,TMPT,NY1)
          CALL MXM (IAJL1 ,NY1,IALJ1 ,NY1,TMP ,NY1)

      
      !     Compute interpolation operators between Gauss-Lobatto Legendre
      !     and Gauss-Lobatto Jacobi (to be used in subr. genxyz IN postpre).
      
      
      !     This call is not right, and these arrays are not used. 3/27/02. pff
      !     CALL IGLLM(IALJ3,IATLJ3,ZGM3(1,2),ZAM3,NY3,NY3,NY3,NY3,ALPHA,BETA)
          CALL IGLJM(IALJ3,IATLJ3,ZGM3(1,2),ZAM3,NY3,NY3,NY3,NY3,ALPHA,BETA)
      
      ENDIF
#endif    
    
  ELSE
  
  !***  Three-dimensional case ************************************
  
  
  !     Gauss-Lobatto Legendre mesh (suffix M1)
  !     Generate collocation points and weights
  
      CALL ZWGLL (ZGM1(1,1),WXM1,NX1)
      CALL ZWGLL (ZGM1(1,2),WYM1,NY1)
      CALL ZWGLL (ZGM1(1,3),WZM1,NZ1)
      DO IZ=1,NZ1
          DO IY=1,NY1
              DO IX=1,NX1
                  W3M1(IX,IY,IZ)=WXM1(IX)*WYM1(IY)*WZM1(IZ)
              enddo
          enddo
      END DO
  
  !     Compute derivative matrices
  
      CALL DGLL (DXM1,DXTM1,ZGM1(1,1),NX1,NX1)
      CALL DGLL (DYM1,DYTM1,ZGM1(1,2),NY1,NY1)
      CALL DGLL (DZM1,DZTM1,ZGM1(1,3),NZ1,NZ1)
  
  !     Gauss Legendre mesh (suffix M2)
  !     Generate collocation points and weights
  
      IF(IFSPLIT)THEN
          CALL ZWGLL (ZGM2(1,1),WXM2,NX2)
          CALL ZWGLL (ZGM2(1,2),WYM2,NY2)
          CALL ZWGLL (ZGM2(1,3),WZM2,NZ2)
      ELSE
          CALL ZWGL  (ZGM2(1,1),WXM2,NX2)
          CALL ZWGL  (ZGM2(1,2),WYM2,NY2)
          CALL ZWGL  (ZGM2(1,3),WZM2,NZ2)
      ENDIF
      DO IZ=1,NZ2
          DO IY=1,NY2
              DO IX=1,NX2
                  W3M2(IX,IY,IZ)=WXM2(IX)*WYM2(IY)*WZM2(IZ)
              enddo
          enddo
      END DO
  
  !     Gauss-Loabtto Legendre mesh (suffix M3).
  !     Generate collocation points and weights.
  
      CALL ZWGLL (ZGM3(1,1),WXM3,NX3)
      CALL ZWGLL (ZGM3(1,2),WYM3,NY3)
      CALL ZWGLL (ZGM3(1,3),WZM3,NZ3)
      DO IZ=1,NZ3
          DO IY=1,NY3
              DO IX=1,NX3
                  W3M3(IX,IY,IZ)=WXM3(IX)*WYM3(IY)*WZM3(IZ)
              enddo
          enddo
      END DO
  
  !     Compute derivative matrices
  
      CALL DGLL (DXM3,DXTM3,ZGM3(1,1),NX3,NX3)
      CALL DGLL (DYM3,DYTM3,ZGM3(1,2),NY3,NY3)
      CALL DGLL (DZM3,DZTM3,ZGM3(1,3),NZ3,NZ3)
  
  !     Generate interpolation operators for the staggered mesh
  
      CALL IGLLM (IXM12,IXTM12,ZGM1(1,1),ZGM2(1,1),NX1,NX2,NX1,NX2)
      CALL IGLLM (IYM12,IYTM12,ZGM1(1,2),ZGM2(1,2),NY1,NY2,NY1,NY2)
      CALL IGLLM (IZM12,IZTM12,ZGM1(1,3),ZGM2(1,3),NZ1,NZ2,NZ1,NZ2)
  
  !     NOTE: The splitting scheme has only one mesh!!!!!
  
      IF (IFSPLIT) THEN
          CALL IGLLM (IXM21,IXTM21,ZGM1(1,1),ZGM2(1,1),NX1,NX2,NX1,NX2)
          CALL IGLLM (IYM21,IYTM21,ZGM1(1,2),ZGM2(1,2),NY1,NY2,NY1,NY2)
          CALL IGLLM (IZM21,IZTM21,ZGM1(1,3),ZGM2(1,3),NZ1,NZ2,NZ1,NZ2)
      ELSE
          write(*,*) "Oops: ifsplit is false"
#if 0
          CALL IGLM  (IXM21,IXTM21,ZGM2(1,1),ZGM1(1,1),NX2,NX1,NX2,NX1)
          CALL IGLM  (IYM21,IYTM21,ZGM2(1,2),ZGM1(1,2),NY2,NY1,NY2,NY1)
          CALL IGLM  (IZM21,IZTM21,ZGM2(1,3),ZGM1(1,3),NZ2,NZ1,NZ2,NZ1)
#endif
      ENDIF
  
  !     Compute derivative operators for the staggered mesh
  
      IF(IFSPLIT)THEN
          CALL COPY (DXM12, DXM1, NX1*NX2)
          CALL COPY (DXTM12,DXTM1,NX1*NX2)
          CALL COPY (DYM12, DYM1, NY1*NY2)
          CALL COPY (DYTM12,DYTM1,NY1*NY2)
          CALL COPY (DZM12, DZM1, NZ1*NZ2)
          CALL COPY (DZTM12,DZTM1,NZ1*NZ2)
      ELSE
          CALL DGLLGL (DXM12,DXTM12,ZGM1(1,1),ZGM2(1,1),IXM12, &
          NX1,NX2,NX1,NX2)
          CALL DGLLGL (DYM12,DYTM12,ZGM1(1,2),ZGM2(1,2),IYM12, &
          NY1,NY2,NY1,NY2)
          CALL DGLLGL (DZM12,DZTM12,ZGM1(1,3),ZGM2(1,3),IZM12, &
          NZ1,NZ2,NZ1,NZ2)
      ENDIF
  
  !     Compute interpolation operators for the geometry mesh M3.
  
      CALL IGLLM (IXM13,IXTM13,ZGM1(1,1),ZGM3(1,1),NX1,NX3,NX1,NX3)
      CALL IGLLM (IYM13,IYTM13,ZGM1(1,2),ZGM3(1,2),NY1,NY3,NY1,NY3)
      CALL IGLLM (IZM13,IZTM13,ZGM1(1,3),ZGM3(1,3),NZ1,NZ3,NZ1,NZ3)
      CALL IGLLM (IXM31,IXTM31,ZGM3(1,1),ZGM1(1,1),NX3,NX1,NX3,NX1)
      CALL IGLLM (IYM31,IYTM31,ZGM3(1,2),ZGM1(1,2),NY3,NY1,NY3,NY1)
      CALL IGLLM (IZM31,IZTM31,ZGM3(1,3),ZGM1(1,3),NZ3,NZ1,NZ3,NZ1)
  
  ENDIF

  RETURN
end subroutine genwz

!-----------------------------------------------------------------------
!> \brief Routine to generate all elemental geometric data for mesh 1.
!!  Velocity formulation : global-to-local mapping based on mesh 3
!!  Stress   formulation : global-to-local mapping based on mesh 1
!-----------------------------------------------------------------------
subroutine geom1 ()!xm3,ym3,zm3)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use geom, only : ifgmsh3
  use tstep, only : istep
  implicit none

  real(DP), allocatable, dimension(:,:,:,:) :: &
    XRM1, YRM1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1 

  allocate(XRM1(LX1,LY1,LZ1,LELT) &
  ,           YRM1(LX1,LY1,LZ1,LELT) &
  ,           XSM1(LX1,LY1,LZ1,LELT) &
  ,           YSM1(LX1,LY1,LZ1,LELT) &
  ,           XTM1(LX1,LY1,LZ1,LELT) &
  ,           YTM1(LX1,LY1,LZ1,LELT) &
  ,           ZRM1(LX1,LY1,LZ1,LELT) &
  ,           ZSM1(LX1,LY1,LZ1,LELT) &
  ,           ZTM1(LX1,LY1,LZ1,LELT) )


  IF (IFGMSH3 .AND. ISTEP == 0) THEN
    write(*,*) "Oops: IFGMSH3"
!max        CALL GLMAPM3 (XM3,YM3,ZM3)
  ELSE
      CALL GLMAPM1(XRM1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  ENDIF

  CALL GEODAT1(XRM1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)

  RETURN
end subroutine geom1

!-----------------------------------------------------------------------
!> \brief Routine to generate mapping data based on mesh 1
!!        (Gauss-Legendre Lobatto meshes).
!!
!!         XRM1,  YRM1,  ZRM1   -   dx/dr, dy/dr, dz/dr
!!         XSM1,  YSM1,  ZSM1   -   dx/ds, dy/ds, dz/ds
!!         XTM1,  YTM1,  ZTM1   -   dx/dt, dy/dt, dz/dt
!!         RXM1,  RYM1,  RZM1   -   dr/dx, dr/dy, dr/dz
!!         SXM1,  SYM1,  SZM1   -   ds/dx, ds/dy, ds/dz
!!         TXM1,  TYM1,  TZM1   -   dt/dx, dt/dy, dt/dz
!!         JACM1                -   Jacobian
!-----------------------------------------------------------------------
subroutine glmapm1(XRM1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  use kinds, only : DP
  use size_m, only : ndim, nid
  use size_m, only : nx1, ny1, nz1, nelt
  use size_m, only : lx1, ly1, lz1, lelt
  use geom, only : jacm1, jacmi
  use geom, only : rxm1, rym1, rzm1, sxm1, sym1, szm1, txm1, tym1, tzm1
  use geom, only : xm1, ym1, zm1
  use input, only : ifaxis, ifxyo, ifvo, ifpo, ifto, param
  use soln, only : vx, vy, vz, pr, t
  implicit none

!     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3
!           share the same array structure in Scratch Common /SCRNS/.
!  real(DP) :: xrm1, yrm1, zrm1
!  real(DP) :: xsm1, ysm1
!  real(DP) :: xtm1, ytm1

  real(DP) :: XRM1(LX1,LY1,LZ1,LELT) &
  ,             YRM1(LX1,LY1,LZ1,LELT) &
  ,             XSM1(LX1,LY1,LZ1,LELT) &
  ,             YSM1(LX1,LY1,LZ1,LELT) &
  ,             XTM1(LX1,LY1,LZ1,LELT) &
  ,             YTM1(LX1,LY1,LZ1,LELT) &
  ,             ZRM1(LX1,LY1,LZ1,LELT)
  real(DP) :: ZSM1(LX1,LY1,LZ1,LELT), ZTM1(LX1,LY1,LZ1,LELT)

  integer :: nxy1, nyz1, nxyz1, ntot1
  integer :: ie, ierr, kerr
  integer, external :: iglsum

  NXY1  = NX1*NY1
  NYZ1  = NY1*NZ1
  NXYZ1 = NX1*NY1*NZ1
  NTOT1 = NXYZ1*NELT

  CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1, &
  IFAXIS)

  IF (NDIM == 2) THEN
#if 0
      CALL RZERO   (JACM1,NTOT1)
      CALL ADDCOL3 (JACM1,XRM1,YSM1,NTOT1)
      CALL SUBCOL3 (JACM1,XSM1,YRM1,NTOT1)
      CALL COPY    (RXM1,YSM1,NTOT1)
      CALL COPY    (RYM1,XSM1,NTOT1)
      CALL CHSIGN  (RYM1,NTOT1)
      CALL COPY    (SXM1,YRM1,NTOT1)
      CALL CHSIGN  (SXM1,NTOT1)
      CALL COPY    (SYM1,XRM1,NTOT1)
      CALL RZERO   (RZM1,NTOT1)
      CALL RZERO   (SZM1,NTOT1)
      CALL RONE    (TZM1,NTOT1)
#endif
  ELSE
    jacm1 = xrm1*ysm1*ztm1 + xtm1*yrm1*zsm1 + xsm1*ytm1*zrm1 &
          - xrm1*ytm1*zsm1 - xsm1*yrm1*ztm1 - xtm1*ysm1*zrm1

    rxm1 = ysm1*ztm1 - ytm1*zsm1
    rym1 = xtm1*zsm1 - xsm1*ztm1
    rzm1 = xsm1*ytm1 - xtm1*ysm1 

    sxm1 = ytm1*zrm1 - yrm1*ztm1
    sym1 = xrm1*ztm1 - xtm1*zrm1
    szm1 = xtm1*yrm1 - xrm1*ytm1

    txm1 = yrm1*zsm1 - ysm1*zrm1
    tym1 = xsm1*zrm1 - xrm1*zsm1
    tzm1 = xrm1*ysm1 - xsm1*yrm1
  ENDIF

  kerr = 0
  DO 500 ie=1,NELT
      CALL CHKJAC(JACM1(1,1,1,ie),NXYZ1,ie,xm1,ym1,zm1,ierr)
      if (ierr /= 0) kerr = kerr+1
  500 END DO
  kerr = iglsum(kerr,1)
  if (kerr > 0) then
      ifxyo = .TRUE. 
      ifvo  = .FALSE. 
      ifpo  = .FALSE. 
      ifto  = .FALSE. 
      param(66) = 4
      call outpost(vx,vy,vz,pr,t,'xyz')
      if (nid == 0) write(6,*) 'Jac error 1, setting p66=4, ifxyo=t'
      call exitt
  endif

  jacmi = 1_dp / jacm1

  RETURN
end subroutine glmapm1

!-----------------------------------------------------------------------
!> \brief Routine to generate elemental geometric matrices on mesh 1
!!  (Gauss-Legendre Lobatto mesh).
!-----------------------------------------------------------------------
subroutine geodat1(XRM1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelt, ndim
  use geom, only : ifgmsh3, jacm1, ifrzer, ym1, bm1
  use geom, only : g1m1, rxm1, rym1, rzm1, g2m1, sxm1, sym1, szm1
  use geom, only : g3m1, txm1, tym1, tzm1, g4m1, g5m1, g6m1
  use input, only : ifaxis
  use tstep, only : istep
  use wz_m, only : zam1, w3m1
  implicit none

!   Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3
!         share the same array structure in Scratch Common /SCRNS/.

  real(DP) ::  XRM1(LX1,LY1,LZ1,LELT) &
  ,             YRM1(LX1,LY1,LZ1,LELT) &
  ,             XSM1(LX1,LY1,LZ1,LELT) &
  ,             YSM1(LX1,LY1,LZ1,LELT) &
  ,             XTM1(LX1,LY1,LZ1,LELT) &
  ,             YTM1(LX1,LY1,LZ1,LELT) &
  ,             ZRM1(LX1,LY1,LZ1,LELT)
  real(DP) ::  ZSM1(LX1,LY1,LZ1,LELT) &
  ,             ZTM1(LX1,LY1,LZ1,LELT) 
  real(DP), allocatable ::  WJ(:,:,:,:)

  integer :: nxyz1, ntot1, iel, j, i

  NXYZ1 = NX1*NY1*NZ1
  NTOT1 = NXYZ1*NELT


  IF (IFGMSH3 .AND. ISTEP == 0) then
    write(*,*) "Oops: IFGMSH3"
!max      CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1, &
!                IFAXIS)
  endif

  allocate(WJ(LX1,LY1,LZ1,LELT))
  IF ( .NOT. IFAXIS) THEN
      wj = 1_dp / jacm1
  ELSE
      DO 500 IEL=1,NELT
          IF (IFRZER(IEL)) THEN
              DO J=1,NY1
                  DO I=1,NX1
                      IF (J > 1) THEN
                          WJ(I,J,1,IEL) = YM1(I,J,1,IEL)/ &
                          (JACM1(I,J,1,IEL)*(1.+ZAM1(J)))
                      ELSE
                          WJ(I,J,1,IEL) = YSM1(I,J,1,IEL)/JACM1(I,J,1,IEL)
                      ENDIF
                  enddo
              END DO
          ELSE
              wj(:,:,:,iel) = ym1(:,:,:,iel) / jacm1(:,:,:,iel)
          ENDIF
      500 END DO
  ENDIF

!   Compute geometric factors for integrated del-squared operator.

  IF (NDIM == 2) THEN
      write(*,*) "Whoops! geodat1"
#if 0
      CALL VDOT2 (G1M1,RXM1,RYM1,RXM1,RYM1,NTOT1)
      CALL VDOT2 (G2M1,SXM1,SYM1,SXM1,SYM1,NTOT1)
      CALL VDOT2 (G4M1,RXM1,RYM1,SXM1,SYM1,NTOT1)
      CALL COL2  (G1M1,WJ,NTOT1)
      CALL COL2  (G2M1,WJ,NTOT1)
      CALL COL2  (G4M1,WJ,NTOT1)
      CALL RZERO (G3M1,NTOT1)
      CALL RZERO (G5M1,NTOT1)
      CALL RZERO (G6M1,NTOT1)
#endif
  ELSE
      CALL VDOT3 (G1M1,RXM1,RYM1,RZM1,RXM1,RYM1,RZM1,NTOT1)
      CALL VDOT3 (G2M1,SXM1,SYM1,SZM1,SXM1,SYM1,SZM1,NTOT1)
      CALL VDOT3 (G3M1,TXM1,TYM1,TZM1,TXM1,TYM1,TZM1,NTOT1)
      g1m1 = g1m1 * wj
      g2m1 = g2m1 * wj
      g3m1 = g3m1 * wj
#if 1
      CALL VDOT3 (G4M1,RXM1,RYM1,RZM1,SXM1,SYM1,SZM1,NTOT1)
      CALL VDOT3 (G5M1,RXM1,RYM1,RZM1,TXM1,TYM1,TZM1,NTOT1)
      CALL VDOT3 (G6M1,SXM1,SYM1,SZM1,TXM1,TYM1,TZM1,NTOT1)
      g4m1 = g4m1 * wj
      g5m1 = g5m1 * wj
      g6m1 = g6m1 * wj
#endif
  ENDIF
  deallocate(wj)

!   Multiply the geometric factors GiM1,i=1,5 with the
!   weights on mesh M1.

  DO IEL=1,NELT
!max      IF (IFAXIS) CALL SETAXW1 ( IFRZER(IEL) )
      g1m1(:,:,:,iel) = g1m1(:,:,:,iel) * w3m1
      g2m1(:,:,:,iel) = g2m1(:,:,:,iel) * w3m1
  !            CALL COL2 (G4M1(1,1,1,IEL),W3M1,NXYZ1)
      IF (NDIM == 3) THEN
        g3m1(:,:,:,iel) = g3m1(:,:,:,iel) * w3m1
      !            CALL COL2 (G5M1(1,1,1,IEL),W3M1,NXYZ1)
      !            CALL COL2 (G6M1(1,1,1,IEL),W3M1,NXYZ1)
      ENDIF
  END DO

!   Compute the mass matrix on mesh M1.

  DO IEL=1,NELT
!max      IF (IFAXIS) CALL SETAXW1 ( IFRZER(IEL) )
      bm1(:,:,:,iel) = jacm1(:,:,:,iel) * w3m1
      IF (IFAXIS) THEN
        write(*,*) "Oops: ifaxis"
#if 0
          CALL COL3(BAXM1(1,1,1,IEL),JACM1(1,1,1,IEL),W3M1,NXYZ1)
          IF (IFRZER(IEL)) THEN
              DO 600 J=1,NY1
                  IF (J > 1) THEN
                      DO 610 I=1,NX1
                          BM1(I,J,1,IEL) = BM1(I,J,1,IEL)*YM1(I,J,1,IEL) &
                          /(1.+ZAM1(J))
                          BAXM1(I,J,1,IEL)=BAXM1(I,J,1,IEL)/(1.+ZAM1(J))
                      610 END DO
                  ELSE
                      DO 620 I=1,NX1
                          BM1(I,J,1,IEL) = BM1(I,J,1,IEL)*YSM1(I,J,1,IEL)
                          BAXM1(I,J,1,IEL)=BAXM1(I,J,1,IEL)
                      620 END DO
                  ENDIF
              600 END DO
          ELSE
              CALL COL2 (BM1(1,1,1,IEL),YM1(1,1,1,IEL),NXYZ1)
          ENDIF
#endif
      ENDIF
  
  END DO

  IF(IFAXIS) THEN
    write(*,*) "Oops: ifaxis"
#if 0
      DO IEL=1,NELT
          IF(IFRZER(IEL)) THEN
              DO J=1,NY1
                  DO I=1,NX1
                      IF(J == 1) THEN
                          YINVM1(I,J,1,IEL)=1.0D0/YSM1(I,J,1,IEL)
                      ELSE
                          YINVM1(I,J,1,IEL)=1.0D0/YM1 (I,J,1,IEL)
                      ENDIF
                  ENDDO
              ENDDO
          ELSE
              yinvm1(:,:,:,iel) = 1_dp / ym(:,:,:,iel)
          ENDIF
      ENDDO
#endif
  ENDIF

!   Compute normals, tangents, and areas on elemental surfaces

  CALL SETAREA(XRM1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)

  RETURN
end subroutine geodat1

!-------------------------------------------------------------------
!>\brief Routine to generate all elemental geometric data for mesh 2
!!  (Gauss-Legendre mesh).
!!      RXM2,  RYM2,  RZM2   -   dr/dx, dr/dy, dr/dz
!!      SXM2,  SYM2,  SZM2   -   ds/dx, ds/dy, ds/dz
!!      TXM2,  TYM2,  TZM2   -   dt/dx, dt/dy, dt/dz
!!      JACM2                -   Jacobian
!!      BM2                  -   Mass matrix
!------------------------------------------------------------------
subroutine geom2
  use size_m, only : nx2, ny2, nz2, nelv
  use geom, only : rxm2, rxm1, rym2, rym1, rzm2, rzm1
  use geom, only : sxm2, sxm1, sym2, sym1, szm2, szm1
  use geom, only : txm2, txm1, tym2, tym1, tzm2, tzm1
  use geom, only : jacm2, jacm1, xm2, xm1, ym2, ym1, zm2, zm1
  use input, only : ifsplit
  use geom, only : bm2, bm1
  implicit none

  integer :: nxyz2, ntot2

  NXYZ2 = NX2*NY2*NZ2
  NTOT2 = NXYZ2*NELV

  IF (IFSPLIT) THEN
    ! Mesh 1 and 2 are identical
    rxm2 => rxm1 
    rym2 => rym1 
    rzm2 => rzm1 

    sxm2 => sxm1 
    sym2 => sym1 
    szm2 => szm1 

    txm2 => txm1 
    tym2 => tym1 
    tzm2 => tzm1 

    jacm2 => jacm1
    bm2 => bm1

    xm2 => xm1
    ym2 => ym1
    zm2 => zm1

  ELSE
    write(*,*) "Oops: ifsplit" 
#if 0
  !     Consistent approximation spaces (UZAWA)
  
      IF (NDIM == 2) THEN
          CALL RZERO (RZM2,NTOT2)
          CALL RZERO (SZM2,NTOT2)
          CALL RONE  (TZM2,NTOT2)
      ENDIF
  
      DO 1000 IEL=1,NELV
      
      !        Mapping from mesh M1 to mesh M2
      
          CALL MAP12 (RXM2(1,1,1,IEL),RXM1(1,1,1,IEL),IEL)
          CALL MAP12 (RYM2(1,1,1,IEL),RYM1(1,1,1,IEL),IEL)
          CALL MAP12 (SXM2(1,1,1,IEL),SXM1(1,1,1,IEL),IEL)
          CALL MAP12 (SYM2(1,1,1,IEL),SYM1(1,1,1,IEL),IEL)
          IF (NDIM == 3) THEN
              CALL MAP12 (RZM2(1,1,1,IEL),RZM1(1,1,1,IEL),IEL)
              CALL MAP12 (SZM2(1,1,1,IEL),SZM1(1,1,1,IEL),IEL)
              CALL MAP12 (TXM2(1,1,1,IEL),TXM1(1,1,1,IEL),IEL)
              CALL MAP12 (TYM2(1,1,1,IEL),TYM1(1,1,1,IEL),IEL)
              CALL MAP12 (TZM2(1,1,1,IEL),TZM1(1,1,1,IEL),IEL)
          ENDIF
          CALL MAP12 (JACM2(1,1,1,IEL),JACM1(1,1,1,IEL),IEL)
      
          CALL MAP12 (XM2(1,1,1,IEL),XM1(1,1,1,IEL),IEL)
          CALL MAP12 (YM2(1,1,1,IEL),YM1(1,1,1,IEL),IEL)
          CALL MAP12 (ZM2(1,1,1,IEL),ZM1(1,1,1,IEL),IEL)
      
      !        Compute the mass matrix on mesh M2.
      
          IF (IFAXIS) CALL SETAXW2 ( IFRZER(IEL) )
          CALL COL3 (BM2(1,1,1,IEL),W3M2,JACM2(1,1,1,IEL),NXYZ2)
      
          IF (IFAXIS .AND. IFRZER(IEL)) THEN
              DO 300 J=1,NY2
                  DO 300 I=1,NX2
                      BM2(I,J,1,IEL) = BM2(I,J,1,IEL)*YM2(I,J,1,IEL) &
                      /(1.+ZAM2(J))
              300 END DO
          ELSEIF (IFAXIS .AND. ( .NOT. IFRZER(IEL))) THEN
              CALL COL2 (BM2(1,1,1,IEL),YM2(1,1,1,IEL),NXYZ2)
          ENDIF
      1000 END DO
#endif
  ENDIF

!   Compute inverse of mesh 2 mass matrix, pff 3/5/92
!max  CALL INVERS2(BM2INV,BM2,NTOT2)

  RETURN
end subroutine geom2

!-----------------------------------------------------------------------
!> \brief Compute global-to-local derivatives on mesh 1.
!-----------------------------------------------------------------------
subroutine xyzrst (xrm1,yrm1,zrm1,xsm1,ysm1,zsm1, XTM1,YTM1,ZTM1,IFAXIS)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, nelt, ndim
  use dxyz, only : dxm1, dytm1, dztm1
  use geom, only : xm1, ym1, zm1
  implicit none

  real(DP) :: &
    XRM1(LX1,LY1,LZ1,*),YRM1(LX1,LY1,LZ1,*) &
  , ZRM1(LX1,LY1,LZ1,*),XSM1(LX1,LY1,LZ1,*) &
  , YSM1(LX1,LY1,LZ1,*),ZSM1(LX1,LY1,LZ1,*) &
  , XTM1(LX1,LY1,LZ1,*),YTM1(LX1,LY1,LZ1,*) &
  , ZTM1(LX1,LY1,LZ1,*)
  LOGICAL :: IFAXIS

  integer :: nxy1, nyz1, iel, iz

  NXY1=NX1*NY1
  NYZ1=NY1*NZ1

  DO 100 IEL=1,NELT
    
      IF (IFAXIS) then
        write(*,*) "Oops: ifaxis"
        !CALL SETAXDY ( IFRZER(IEL) )
      endif
    
      CALL MXM (DXM1,NX1,XM1(1,1,1,IEL),NX1,XRM1(1,1,1,IEL),NYZ1)
      CALL MXM (DXM1,NX1,YM1(1,1,1,IEL),NX1,YRM1(1,1,1,IEL),NYZ1)
      CALL MXM (DXM1,NX1,ZM1(1,1,1,IEL),NX1,ZRM1(1,1,1,IEL),NYZ1)
  
      DO 10 IZ=1,NZ1
          CALL MXM (XM1(1,1,IZ,IEL),NX1,DYTM1,NY1,XSM1(1,1,IZ,IEL),NY1)
          CALL MXM (YM1(1,1,IZ,IEL),NX1,DYTM1,NY1,YSM1(1,1,IZ,IEL),NY1)
          CALL MXM (ZM1(1,1,IZ,IEL),NX1,DYTM1,NY1,ZSM1(1,1,IZ,IEL),NY1)
      10 END DO
  
      IF (NDIM == 3) THEN
          CALL MXM (XM1(1,1,1,IEL),NXY1,DZTM1,NZ1,XTM1(1,1,1,IEL),NZ1)
          CALL MXM (YM1(1,1,1,IEL),NXY1,DZTM1,NZ1,YTM1(1,1,1,IEL),NZ1)
          CALL MXM (ZM1(1,1,1,IEL),NXY1,DZTM1,NZ1,ZTM1(1,1,1,IEL),NZ1)
      ELSE
          CALL RZERO (XTM1(1,1,1,IEL),NXY1)
          CALL RZERO (YTM1(1,1,1,IEL),NXY1)
          CALL RONE  (ZTM1(1,1,1,IEL),NXY1)
      ENDIF
  
  100 END DO

  RETURN
end subroutine xyzrst

!> \brief Check the array JAC for a change in sign.
subroutine chkjac(jac,n,iel,X,Y,Z,IERR)
  use kinds, only : DP
  use size_m
  use parallel
  implicit none

  integer :: n, iel, ierr
  REAL(DP) :: JAC(N),x(1),y(1),z(1)
  real(DP) :: sign
  integer :: i, ieg

  ierr = 1
  SIGN = JAC(1)
  DO 100 I=2,N
      IF (SIGN*JAC(I) <= 0.0) THEN
          ieg = lglel(iel)
          WRITE(6,101) nid,I,ieg
          write(6,*) jac(i-1),jac(i)
          if (ndim == 3) then
              write(6,7) nid,x(i-1),y(i-1),z(i-1)
              write(6,7) nid,x(i),y(i),z(i)
          else
              write(6,7) nid,x(i-1),y(i-1)
              write(6,7) nid,x(i),y(i)
          endif
          7 format(i5,' xyz:',1p3e14.5)
      !           if (np.eq.1) call out_xyz_el(x,y,z,iel)
      !           ierr=0
          return
      ENDIF
  100 END DO
  101 FORMAT(//,i5,2x &
  ,'ERROR:  Vanishing Jacobian near',i7,'th node of element' &
  ,I10,'.')


  ierr = 0
  RETURN
end subroutine chkjac

!-----------------------------------------------------------------------
!> \brief Compute the volume based on mesh M1 and mesh M2
subroutine volume
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, nelt, nx2, ny2, nz2, nfield, nid
  use esolv, only : volel
  use input, only : ifmvbd, ifmhd, iftmsh
  use geom, only : volvm1, volvm2, voltm1, voltm2, bm1, bm2
  use tstep, only : volfld
  implicit none

  real(DP), external :: glsum, vlsum
  integer :: e, mfield, nfldt, ifld, nxyz

  volvm1=glsum(bm1,nx1*ny1*nz1*nelv)
  volvm2=glsum(bm2,nx2*ny2*nz2*nelv)
  voltm1=glsum(bm1,nx1*ny1*nz1*nelt)
  voltm2=glsum(bm2,nx2*ny2*nz2*nelt)
  mfield=1
  if (ifmvbd) mfield=0
  nfldt = nfield
  if (ifmhd) nfldt = nfield+1

  do ifld=mfield,nfldt
      if (iftmsh(ifld)) then
          volfld(ifld) = voltm1
      else
          volfld(ifld) = volvm1
      endif
  enddo

  if (nid == 0) write(6,*) 'vol_t,vol_v:',voltm1,volvm1


  nxyz = nx1*ny1*nz1
  do e=1,nelt
      volel(e) = vlsum(bm1(1,1,1,e),nxyz)
  enddo

  return
end subroutine volume

!-----------------------------------------------------------------------
!> \brief Compute surface data: areas, normals and tangents
subroutine setarea(xrm1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, nz1, nelt, ndim
  use geom, only : area, unx, uny, unz
  implicit none

  real(DP) :: XRM1(LX1,LY1,LZ1,LELT) &
  ,           YRM1(LX1,LY1,LZ1,LELT) &
  ,           XSM1(LX1,LY1,LZ1,LELT) &
  ,           YSM1(LX1,LY1,LZ1,LELT) &
  ,           XTM1(LX1,LY1,LZ1,LELT) &
  ,           YTM1(LX1,LY1,LZ1,LELT) &
  ,           ZRM1(LX1,LY1,LZ1,LELT) &
  ,           ZSM1(LX1,LY1,LZ1,LELT) &
  ,           ZTM1(LX1,LY1,LZ1,LELT)

  integer :: nsrf

  NSRF  = 6*NX1*NZ1*NELT

  area = 0_dp
  unx  = 0_dp; uny = 0_dp; unz = 0_dp
!  CALL RZERO3 (T1X,T1Y,T1Z,NSRF)
!  CALL RZERO3 (T2X,T2Y,T2Z,NSRF)

  IF (NDIM == 2) THEN
!max        CALL AREA2
  ELSE
      CALL AREA3(xrm1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  ENDIF

  RETURN
end subroutine setarea

!--------------------------------------------------------------------
!> \brief Compute areas, normals and tangents (3D geom.)
!--------------------------------------------------------------------
subroutine area3(xrm1, yrm1, zrm1, xsm1, ysm1, zsm1, xtm1, ytm1, ztm1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelt, ndim
  use geom, only : area, unx, uny, unz
  use wz_m, only : wxm1, wym1, wzm1
  implicit none

!   Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3
!         share the same array structure in Scratch Common /SCRNS/.

  real(DP) :: XRM1(LX1,LY1,LZ1,LELT) &
  ,             YRM1(LX1,LY1,LZ1,LELT) &
  ,             XSM1(LX1,LY1,LZ1,LELT) &
  ,             YSM1(LX1,LY1,LZ1,LELT) &
  ,             XTM1(LX1,LY1,LZ1,LELT) &
  ,             YTM1(LX1,LY1,LZ1,LELT) &
  ,             ZRM1(LX1,LY1,LZ1,LELT)
  real(DP) ::  ZSM1(LX1,LY1,LZ1,LELT) &
  ,             ZTM1(LX1,LY1,LZ1,LELT)

  real(DP), allocatable :: A(:,:,:,:), B(:,:,:,:), C(:,:,:,:), dot(:,:,:,:)
  integer :: nxy1, nface, ntot, nsrf, iel, iz, iy, ix
  real(DP) :: weight

  NXY1  = NX1*NY1
  NFACE = 2*NDIM
  NTOT  = NX1*NY1*NZ1*NELT
  NSRF  = 6*NX1*NY1*NELT

!      "R"

  allocate(A(lx1,ly1,lz1,lelt), B(lx1,ly1,lz1,lelt), C(lx1,ly1,lz1,lelt))
  allocate(dot(lx1,ly1,lz1,lelt))
  CALL VCROSS(A,B,C,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,NTOT)
  CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)

  DO IEL=1,NELT
      DO IZ=1,NZ1
          DO IY=1,NY1
              WEIGHT = WYM1(IY)*WZM1(IZ)
              AREA(IY,IZ,2,IEL) = SQRT(DOT(NX1,IY,IZ,IEL))*WEIGHT
              AREA(IY,IZ,4,IEL) = SQRT(DOT(  1,IY,IZ,IEL))*WEIGHT
              UNX (IY,IZ,4,IEL) = -A(  1,IY,IZ,IEL)
              UNX (IY,IZ,2,IEL) =  A(NX1,IY,IZ,IEL)
              UNY (IY,IZ,4,IEL) = -B(  1,IY,IZ,IEL)
              UNY (IY,IZ,2,IEL) =  B(NX1,IY,IZ,IEL)
              UNZ (IY,IZ,4,IEL) = -C(  1,IY,IZ,IEL)
              UNZ (IY,IZ,2,IEL) =  C(NX1,IY,IZ,IEL)
          enddo
      enddo
  END DO

!      "S"

  CALL VCROSS(A,B,C,XRM1,YRM1,ZRM1,XTM1,YTM1,ZTM1,NTOT)
  CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)
  DO IEL=1,NELT
      DO IZ=1,NZ1
          DO IX=1,NX1
              WEIGHT=WXM1(IX)*WZM1(IZ)
              AREA(IX,IZ,1,IEL) = SQRT(DOT(IX,  1,IZ,IEL))*WEIGHT
              AREA(IX,IZ,3,IEL) = SQRT(DOT(IX,NY1,IZ,IEL))*WEIGHT
              UNX (IX,IZ,1,IEL) =  A(IX,  1,IZ,IEL)
              UNX (IX,IZ,3,IEL) = -A(IX,NY1,IZ,IEL)
              UNY (IX,IZ,1,IEL) =  B(IX,  1,IZ,IEL)
              UNY (IX,IZ,3,IEL) = -B(IX,NY1,IZ,IEL)
              UNZ (IX,IZ,1,IEL) =  C(IX,  1,IZ,IEL)
              UNZ (IX,IZ,3,IEL) = -C(IX,NY1,IZ,IEL)
          enddo
      enddo
  END DO

!      "T"

  CALL VCROSS(A,B,C,XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,NTOT)
  CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)
  DO IEL=1,NELT
      DO IX=1,NX1
          DO IY=1,NY1
              WEIGHT=WXM1(IX)*WYM1(IY)
              AREA(IX,IY,5,IEL) = SQRT(DOT(IX,IY,  1,IEL))*WEIGHT
              AREA(IX,IY,6,IEL) = SQRT(DOT(IX,IY,NZ1,IEL))*WEIGHT
              UNX (IX,IY,5,IEL) = -A(IX,IY,  1,IEL)
              UNX (IX,IY,6,IEL) =  A(IX,IY,NZ1,IEL)
              UNY (IX,IY,5,IEL) = -B(IX,IY,  1,IEL)
              UNY (IX,IY,6,IEL) =  B(IX,IY,NZ1,IEL)
              UNZ (IX,IY,5,IEL) = -C(IX,IY,  1,IEL)
              UNZ (IX,IY,6,IEL) =  C(IX,IY,NZ1,IEL)
          enddo
      enddo
  END DO

  CALL UNITVEC (UNX,UNY,UNZ,NSRF)

!   COMPUTE UNIT TANGENT T1
#if 0
  DO 600 IEL=1,NELT
      DO 600 IFC=1,NFACE
          IF (IFC == 1 .OR. IFC == 6) THEN
              CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL), &
              T1Z(1,1,IFC,IEL), &
              XRM1(1,1,1,IEL),YRM1(1,1,1,IEL), &
              ZRM1(1,1,1,IEL),IFC,0)
          ELSEIF (IFC == 2 .OR. IFC == 5) THEN
              CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL), &
              T1Z(1,1,IFC,IEL), &
              XSM1(1,1,1,IEL),YSM1(1,1,1,IEL), &
              ZSM1(1,1,1,IEL),IFC,0)
          ELSE
              CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL), &
              T1Z(1,1,IFC,IEL), &
              XTM1(1,1,1,IEL),YTM1(1,1,1,IEL), &
              ZTM1(1,1,1,IEL),IFC,0)
          ENDIF
  600 END DO

  CALL UNITVEC (T1X,T1Y,T1Z,NSRF)

!   COMPUTE UNIT TANGENT T2  ( T2 = Normal X T1 )
  DO 700 IEL=1,NELT
      DO 700 IFC=1,NFACE
          CALL VCROSS (T2X(1,1,IFC,IEL),T2Y(1,1,IFC,IEL), &
          T2Z(1,1,IFC,IEL), &
          UNX(1,1,IFC,IEL),UNY(1,1,IFC,IEL), &
          UNZ(1,1,IFC,IEL), &
          T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL), &
          T1Z(1,1,IFC,IEL),NXY1)
  700 END DO
#endif

  RETURN
end subroutine area3

!--------------------------------------------------------------------
!> \brief Lag the mass matrix (matrices)
!--------------------------------------------------------------------
subroutine lagmass()
  use size_m, only : nx1, ny1, nz1, nelt
  use geom, only : bm1, bm1lag
  use tstep, only : nbdinp
  implicit none

  integer :: ntot1, ilag
  NTOT1 = NX1*NY1*NZ1*NELT
  DO ILAG=NBDINP-1,2,-1
      CALL COPY (BM1LAG(1,1,1,1,ILAG),BM1LAG(1,1,1,1,ILAG-1),NTOT1)
  END DO
  CALL COPY (BM1LAG(1,1,1,1,1),BM1,NTOT1)

  RETURN
end subroutine lagmass

!--------------------------------------------------------------------
!> \brief Invert the mass matrix.
!! 1)  Copy BM1 to BINVM1
!! 2)  Perform direct stiffness summation on BINVM1
!! 3)  Compute BINVM1 = 1/BINVM1
!! 4)  Two inverse mass matrices required because of difference
!!     in DSSUM routine for IMESH=1 and IMESH=2.
!--------------------------------------------------------------------
subroutine setinvm()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, nelt, nfield
  use input, only : ifflow, ifheat, iftmsh
  use geom, only : bm1, binvm1, bintm1
  use tstep, only : ifield
  implicit none

  logical :: any_iftmsh
  integer :: nxyz1, ifld, store_field, ntot
  nxyz1  = nx1*ny1*nz1

  store_field = ifield

  IF (IFFLOW) THEN ! Velocity mass matrix
      IFIELD = 1
      NTOT   = NXYZ1*NELV
      CALL COPY    (BINVM1,BM1,NTOT)
      CALL DSSUM   (BINVM1)
      binvm1 = 1._dp / binvm1 
  ENDIF

  any_iftmsh = .false.
  do ifld = 1,nfield
    if (iftmsh(ifld)) any_iftmsh = .true.
  enddo
   
  IF (IFHEAT .and. any_iftmsh) THEN ! Temperature mass matrix
      IFIELD = 2
      NTOT   = NXYZ1*NELT
      CALL COPY    (BINTM1,BM1,NTOT)
      CALL DSSUM   (BINTM1)
      bintm1 = 1._dp / bintm1
  ENDIF

  ifield = store_field

  return
end subroutine setinvm
!-----------------------------------------------------------------------
