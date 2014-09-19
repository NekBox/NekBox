!=======================================================================
subroutine hmholtz(name,u,rhs,h1,h2,mask,mult,imsh,tli,maxit,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, nelv, nelt, ndim, nid
  use ctimer, only : icalld, thmhz, nhmhz, etime1, dnekclock
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

  tol = abs(tli)

  if (icalld == 0) thmhz=0.0

  iffdm = .FALSE. 
!   iffdm = .true.
  if (ifsplit) iffdm = .TRUE. 

!max    if (icalld == 0 .AND. iffdm) call set_fdm_prec_h1A

  icalld=icalld+1
  nhmhz=icalld
  etime1=dnekclock()
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

  call dssum   (rhs)
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


  thmhz=thmhz+(dnekclock()-etime1)
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
  use ctimer, only : icalld, taxhm, naxhm, etime1, dnekclock
  use geom, only : g4m1, g5m1, g6m1
  use dxyz, only : wddx, wddyt, wddzt
  use input, only : ifaxis
  use geom, only : bm1
  use mesh, only : iffast, ifsolv
  implicit none

  integer, intent(in) :: imesh, isd
  REAL(DP), intent(out) :: AU    (LX1,LY1,LZ1,*)
  real(DP), intent(in)  :: U     (LX1,LY1,LZ1,*) 
  real(DP), intent(in)  :: HELM1 (LX1,LY1,LZ1,*) 
  real(DP), intent(in)  :: HELM2 (LX1,LY1,LZ1,*)

  ! locals
  REAL(DP) ::           TM1   (LX1,LY1,LZ1)
  REAL(DP) ::           TM2   (LX1,LY1,LZ1)
  REAL(DP) ::           TM3   (LX1,LY1,LZ1)

  integer :: nel, nxy, nyz, nxz, nxyz, ntot
  real(DP) :: h1 
  integer :: e, iz

  nel=nelt
  if (imesh == 1) nel=nelv

  NXY=NX1*NY1
  NYZ=NY1*NZ1
  NXZ=NX1*NZ1
  NXYZ=NX1*NY1*NZ1
  NTOT=NXYZ*NEL

  if (icalld == 0) taxhm=0.0
  icalld=icalld+1
  naxhm=icalld
  etime1=dnekclock()

  IF ( .NOT. IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
  au(:,:,:,1:nel) = 0._dp

  do 100 e=1,nel
  
!      if (ifaxis) call setaxdy ( ifrzer(e) )
  
      IF (NDIM == 2) THEN
          write(*,*) "Whoops! axhelm"
#if 0
      
      !       2-d case ...............
      
          if (iffast(e)) then
          
          !          Fast 2-d mode: constant properties and undeformed element
          
              h1 = helm1(1,1,1,e)
              call mxm   (wddx,nx1,u(1,1,1,e),nx1,tm1,nyz)
              call mxm   (u(1,1,1,e),nx1,wddyt,ny1,tm2,ny1)
              call col2  (tm1,g4m1(1,1,1,e),nxyz)
              call col2  (tm2,g5m1(1,1,1,e),nxyz)
              call add3  (au(1,1,1,e),tm1,tm2,nxyz)
              call cmult (au(1,1,1,e),h1,nxyz)
          
          else
          
          !          General case, speed-up for undeformed elements
          
              call mxm  (dxm1,nx1,u(1,1,1,e),nx1,dudr,nyz)
              call mxm  (u(1,1,1,e),nx1,dytm1,ny1,duds,ny1)
              call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
              call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
              if (ifdfrm(e)) then
                  call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
                  call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
              endif
              call col2 (tmp1,helm1(1,1,1,e),nxyz)
              call col2 (tmp2,helm1(1,1,1,e),nxyz)
              call mxm  (dxtm1,nx1,tmp1,nx1,tm1,nyz)
              call mxm  (tmp2,nx1,dym1,ny1,tm2,ny1)
              call add2 (au(1,1,1,e),tm1,nxyz)
              call add2 (au(1,1,1,e),tm2,nxyz)

          endif
      
#endif
      else
      
      !       3-d case ...............
      
          if (iffast(e)) then
          
          !          Fast 3-d mode: constant properties and undeformed element
          
              h1 = helm1(1,1,1,e)

              call mxm   (wddx,nx1,u(1,1,1,e),nx1,tm1,nyz)
              do 5 iz=1,nz1
                  call mxm   (u(1,1,iz,e),nx1,wddyt,ny1,tm2(1,1,iz),ny1)
              5 END DO
              call mxm   (u(1,1,1,e),nxy,wddzt,nz1,tm3,nz1)
#if 0
              call col2  (tm1,g4m1(1,1,1,e),nxyz)
              call col2  (tm2,g5m1(1,1,1,e),nxyz)
              call col2  (tm3,g6m1(1,1,1,e),nxyz)

              call add3  (au(1,1,1,e),tm1,tm2,nxyz)
              call add2  (au(1,1,1,e),tm3,nxyz)

              call cmult (au(1,1,1,e),h1,nxyz)
#else
                 
              au(:,:,:,e) = &
              h1 * ( &
              au(:,:,:,e) + tm1*g4m1(:,:,:,e) &
              + tm2*g5m1(:,:,:,e) + tm3*g6m1(:,:,:,e) &
              )
#endif

            
          else
              write(*,*) "Woops! axhelm"
#if 0
          
          !          General case, speed-up for undeformed elements
          
              call mxm(dxm1,nx1,u(1,1,1,e),nx1,dudr,nyz)
              do 10 iz=1,nz1
                  call mxm(u(1,1,iz,e),nx1,dytm1,ny1,duds(1,1,iz),ny1)
              10 END DO
              call mxm     (u(1,1,1,e),nxy,dztm1,nz1,dudt,nz1)
              call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
              call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
              call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
              if (ifdfrm(e)) then
                  call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
                  call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
                  call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
                  call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
                  call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
                  call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
              endif
              call col2 (tmp1,helm1(1,1,1,e),nxyz)
              call col2 (tmp2,helm1(1,1,1,e),nxyz)
              call col2 (tmp3,helm1(1,1,1,e),nxyz)
              call mxm  (dxtm1,nx1,tmp1,nx1,tm1,nyz)
              do 20 iz=1,nz1
                  call mxm(tmp2(1,1,iz),nx1,dym1,ny1,tm2(1,1,iz),ny1)
              20 END DO
              call mxm  (tmp3,nxy,dzm1,nz1,tm3,nz1)
              call add2 (au(1,1,1,e),tm1,nxyz)
              call add2 (au(1,1,1,e),tm2,nxyz)
              call add2 (au(1,1,1,e),tm3,nxyz)
            
#endif
          endif
      
      endif
  
  100 END DO

  au(:,:,:,1:nel) = au(:,:,:,1:nel) + helm2(:,:,:,1:nel)*bm1(:,:,:,1:nel)*u(:,:,:,1:nel)
!  call addcol4 (au,helm2,bm1,u,ntot)

!   If axisymmetric, add a diagonal term in the radial direction (ISD=2)

  if (ifaxis .AND. (isd == 2)) then
      write(*,*) "Whoops! axhelm 3"
#if 0
      do 200 e=1,nel
      
          if (ifrzer(e)) then
              call mxm(u  (1,1,1,e),nx1,datm1,ny1,duax,1)
              call mxm(ym1(1,1,1,e),nx1,datm1,ny1,ysm1,1)
          endif
      
          do 190 j=1,ny1
              do 190 i=1,nx1
              !               if (ym1(i,j,1,e).ne.0.) then
                  if (ifrzer(e)) then
                      term1 = 0.0
                      if(j /= 1) &
                      term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                      term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i) &
                      *jacm1(i,1,1,e)/ysm1(i)
                  else
                      term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                      term2 = 0.
                  endif
                  au(i,j,1,e) = au(i,j,1,e) &
                  + helm1(i,j,1,e)*(term1+term2)
              !               endif
          190 END DO
      200 END DO
#endif
  endif

  taxhm=taxhm+(dnekclock()-etime1)
  return
end subroutine axhelm

!=======================================================================
!> \brief Set logicals for fast evaluation of A*x
!-------------------------------------------------------------------
subroutine setfast (helm1,helm2,imesh)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelt, nelv
  use input, only : ifaxis, ifmodel
  use mesh, only : iffast, ifdfrm
  implicit none

  integer, intent(in) :: imesh
  REAL(DP), intent(in) :: HELM1(NX1,NY1,NZ1,*), HELM2(NX1,NY1,NZ1,*)
 
  integer :: nel, nxyz, ntot, ie 
  real(DP) :: delta, x, y, diff, epsm, h1min, h1max, testh1
  real(DP) :: den
  real(DP), external :: vlmin, vlmax, vlamax

  nel = -1
  IF (IMESH == 1) NEL=NELV
  IF (IMESH == 2) NEL=NELT
  NXYZ = NX1*NY1*NZ1
  NTOT = NXYZ*NEL

  DELTA = 1.E-9
  X    = 1.+DELTA
  Y    = 1.
  DIFF = ABS(X-Y)
  IF (DIFF == 0.0) EPSM = 1.E-6
  IF (DIFF > 0.0) EPSM = 1.E-13

  DO 100 ie=1,NEL
      IFFAST(ie) = .FALSE. 
      IF (IFDFRM(ie) .OR. IFAXIS .OR. IFMODEL ) THEN
          IFFAST(ie) = .FALSE. 
      !            write(*,*) IFDFRM(ie), IFAXIS, IFMODEL
      ELSE
          H1MIN  = VLMIN(HELM1(1,1,1,ie),NXYZ)
          H1MAX  = VLMAX(HELM1(1,1,1,ie),NXYZ)
          den    = abs(h1max)+abs(h1min)
          if (den > 0) then
              TESTH1 = ABS((H1MAX-H1MIN)/(H1MAX+H1MIN))
              IF (TESTH1 < EPSM) IFFAST(ie) = .TRUE. 
          else
              iffast(ie) = .TRUE. 
          endif
      ENDIF
  100 END DO

  return
end subroutine setfast

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
  use dxyz, only : dxtm1, dytm1, dztm1, datm1, dam1
  use geom, only : ifrzer, g1m1, g2m1, g3m1, ym1, jacm1
  use input, only : ifaxis
  use geom, only : bm1
  use mesh, only : ifdfrm
  use wz_m, only : wxm1, wam1
  implicit none

  REAL(DP), intent(out) ::  DPCM1 (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  ::  HELM1(NX1,NY1,NZ1,*), HELM2(NX1,NY1,NZ1,*)
  integer, intent(in) :: imsh, isd

  REAL(DP) :: YSM1(LY1)

  integer :: nel, ntot, ie, iq, iz, iy, ix, j, i, iel
  real(DP) :: term1, term2

  nel=nelt
  if (imsh == 1) nel=nelv

  ntot = nel*nx1*ny1*nz1

!   The following lines provide a convenient debugging option
!   call rone(dpcm1,ntot)
!   if (ifield.eq.1) call copy(dpcm1,binvm1,ntot)
!   if (ifield.eq.2) call copy(dpcm1,bintm1,ntot)
!   return

  dpcm1(:,:,:,1:nel) = 0._dp
  DO 1000 IE=1,NEL

!      IF (IFAXIS) CALL SETAXDY ( IFRZER(IE) )

      DO IQ=1,NX1
          DO IZ=1,NZ1
              DO IY=1,NY1
                  DO IX=1,NX1
                      DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                      G1M1(IQ,IY,IZ,IE) * DXTM1(IX,IQ)**2
                  enddo
              enddo
          enddo
      END DO
      DO IQ=1,NY1
          DO IZ=1,NZ1
              DO IY=1,NY1
                  DO IX=1,NX1
                      DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                      G2M1(IX,IQ,IZ,IE) * DYTM1(IY,IQ)**2
                  enddo
              enddo
          enddo
      END DO

      IF (NDIM == 3) THEN
          DO IQ=1,NZ1
              DO IZ=1,NZ1
                  DO IY=1,NY1
                      DO IX=1,NX1
                          DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + &
                          G3M1(IX,IY,IQ,IE) * DZTM1(IZ,IQ)**2
                      enddo
                  enddo
              enddo
          END DO
      
      !       Add cross terms if element is deformed.
      
          IF (IFDFRM(IE)) THEN
              write(*,*) "Whoops!"
#if 0
              DO 600 IY=1,NY1
                  DO 600 IZ=1,NZ1
                      DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE) &
                      + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY) &
                      + G5M1(1,IY,IZ,IE) * DXTM1(1,1)*DZTM1(IZ,IZ)
                      DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE) &
                      + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY) &
                      + G5M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DZTM1(IZ,IZ)
              600 END DO
              DO 700 IX=1,NX1
                  DO 700 IZ=1,NZ1
                      DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE) &
                      + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX) &
                      + G6M1(IX,1,IZ,IE) * DYTM1(1,1)*DZTM1(IZ,IZ)
                      DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE) &
                      + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX) &
                      + G6M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DZTM1(IZ,IZ)
              700 END DO
              DO 800 IX=1,NX1
                  DO 800 IY=1,NY1
                      DPCM1(IX,IY,1,IE) = DPCM1(IX,IY,1,IE) &
                      + G5M1(IX,IY,1,IE) * DZTM1(1,1)*DXTM1(IX,IX) &
                      + G6M1(IX,IY,1,IE) * DZTM1(1,1)*DYTM1(IY,IY)
                      DPCM1(IX,IY,NZ1,IE) = DPCM1(IX,IY,NZ1,IE) &
                      + G5M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DXTM1(IX,IX) &
                      + G6M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DYTM1(IY,IY)
              800 END DO
#endif
          ENDIF
      ELSE
      
          IF (IFDFRM(IE)) THEN
              write(*,*) "Whoops!"
#if 0
              IZ=1
              DO 602 IY=1,NY1
                  DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE) &
                  + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY)
                  DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE) &
                  + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY)
              602 END DO
              DO 702 IX=1,NX1
                  DO 702 IZ=1,NZ1
                      DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE) &
                      + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX)
                      DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE) &
                      + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX)
              702 END DO
#endif
          ENDIF
      ENDIF
  1000 END DO

  dpcm1(:,:,:,1:nel) = dpcm1(:,:,:,1:nel)*helm1(:,:,:,1:nel) + bm1*helm2(:,:,:,1:nel)

!   If axisymmetric, add a diagonal term in the radial direction (ISD=2)

  IF (IFAXIS .AND. (ISD == 2)) THEN
      DO 1200 IEL=1,NEL
      
          IF (IFRZER(IEL)) THEN
              CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
          ENDIF
      
          DO J=1,NY1
              DO I=1,NX1
                  IF (YM1(I,J,1,IEL) /= 0.) THEN
                      TERM1 = BM1(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                      IF (IFRZER(IEL)) THEN
                          TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J) &
                          *JACM1(I,1,1,IEL)/YSM1(I)
                      ELSE
                          TERM2 = 0.
                      ENDIF
                      DPCM1(I,J,1,IEL) = DPCM1(I,J,1,IEL) &
                      + HELM1(I,J,1,IEL)*(TERM1+TERM2)
                  ENDIF
              enddo
          END DO
      1200 END DO
  ENDIF

  CALL DSSUM (DPCM1)
  dpcm1(:,:,:,1:nel) = 1._dp / dpcm1(:,:,:,1:nel)

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
  implicit none

  real(DP), intent(out) :: tol
  REAL(DP), intent(in)  :: RES  (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: H1   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: H2   (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: MASK (LX1,LY1,LZ1,*)
  REAL(DP), intent(in)  :: MULT (LX1,LY1,LZ1,*)
  integer,  intent(in)  :: imesh, isd

  real(DP), allocatable :: W1(:,:,:,:), W2(:,:,:,:)

  real(DP) :: acondno, delta, x, y, diff, eps
  real(DP) :: vol, rinit, rmin, bcneu1, bcneu2, bctest, bcrob, tolmin
  integer :: nl, ntot1
  real(DP), external :: glsc3, glsum


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

  CALL RONE (W1,NTOT1)
  BCNEU1 = GLSC3(W1,MASK,MULT,NTOT1)
  BCNEU2 = GLSC3(W1,W1  ,MULT,NTOT1)
  BCTEST = ABS(BCNEU1-BCNEU2)

  CALL AXHELM (W2,W1,H1,H2,IMESH,ISD)
  w2 = w2 * w2 * bm1
  BCROB  = SQRT(GLSUM(W2,NTOT1)/VOL)

  IF ((BCTEST < .1) .AND. (BCROB < (EPS*ACONDNO))) THEN
  !         OTR = GLSC3 (W1,RES,MULT,NTOT1)
      TOLMIN = RINIT*EPS*10.
      IF (TOL < TOLMIN) THEN
          TOL = TOLMIN
          IF (NID == 0 .AND. IFPRINT) &
          WRITE(6,*) 'New CG1-tolerance (Neumann) = ',TOLMIN
      ENDIF
  ENDIF

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
  use fdmh1, only : kfldfdm
  use input, only : ifsplit, param, ifprint
  use geom, only : volvm1, voltm1, bm1
  use mesh, only : ifsolv, niterhm
  use tstep, only : istep, imesh
  implicit none

  real(DP), intent(out) :: x(*) !>!< solution vector
  real(DP), intent(in)  :: f(*) !>!< residual vector
  real(DP), intent(in)  :: h1(*) !>!< coefficient of A (stiffness)
  real(DP), intent(in)  :: h2(*) !>!< coefficient of M (mass)
  real(DP), intent(in)  :: mask(*) !>!< mask array
  real(DP), intent(in)  :: mult(*) !>!< multiplicity array
  real(DP), intent(in)  :: binv(*) !>!< inverse of mass matrix
  integer,  intent(in)  :: imsh, isd, maxit
  real(DP), intent(in)  :: tin !>!< input tolerance
  character(4) :: name

  logical :: ifmcor,ifprint_hmh

  integer, parameter :: lg=lx1*ly1*lz1*lelt
  real(DP) :: scalar(2)
  real(DP), allocatable :: d (:)
  real(DP), allocatable :: r (:) , w (:) , p (:) , z (:)

  integer, parameter :: maxcg=900
  real(DP), allocatable :: diagt(:), upper(:)

  integer :: n, iter, nxyz, nel, niter, krylov
  real(DP) :: rho, vol, tol, h2max, skmin, smean, rmean
  real(DP) :: rtz1, rtz2, rbn2, rbn0, beta, rho0, alpha, alphm
  real(DP), external :: glmax, glmin, glsum, glsc2, vlsc3, vlsc32, glsc3


  if (ifsplit .AND. name == 'PRES' .AND. param(42) == 0) then
      n = nx1*ny1*nz1*nelv
      call copy      (x,f,n)
      call hmh_gmres (x,h1,h2,mult,iter)
      niterhm = iter
      return
  endif
!      write(6,*) ifsplit,name,param(44),' P44 C'

! **  zero out stuff for Lanczos eigenvalue estimator
  allocate(diagt(maxcg),upper(maxcg))
  diagt = 0_dp; upper = 0_dp
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
  niter = min(maxit,maxcg)

!     Speed-up for undeformed elements and constant properties.
  if ( .NOT. ifsolv) then
      call setfast(h1,h2,imesh)
      ifsolv = .TRUE. 
  endif

!     Set up diag preconditioner.

  allocate(d(lg)); d = 0_dp
  if (kfldfdm < 0) then
      call setprec(D,h1,h2,imsh,isd)
  elseif(param(100) /= 2) then
!max      call set_fdm_prec_h1b(d,h1,h2,nel)
  endif

  allocate(r(lg))
  r(1:n) = f(1:n)
  x(1:n) = 0._dp

!     Check for non-trivial null-space

  ifmcor = .FALSE. 
  h2max = glmax(h2  ,n)
  skmin = glmin(mask,n)
  if (skmin > 0 .AND. h2max == 0) ifmcor = .TRUE. 

  smean = 0.
  if (name == 'PRES') then
  !        call ortho (r)           ! Commented out March 15, 2011,pff
  elseif (ifmcor) then

      smean = -1./glsum(bm1,n) ! Modified 5/4/12 pff
      rmean = smean*glsc2(r,mult,n)
      call copy(x,bm1,n)
      call dssum(x)
      r = r + rmean * x(1:n)
      x (1:n) = 0._dp
  endif

  krylov = 0
  rtz1=1.0
  niterhm = 0

  allocate(z(lg)); z = 0_dp
  allocate(w(lg)); w = 0_dp
  allocate(p(lg)); p = 0_dp
  do iter=1,niter
  
      if (kfldfdm < 0) then  ! Jacobi Preconditioner
      !           call copy(z,r,n)
          z = r * d
      else                                       ! Schwarz Preconditioner
      write (*,*) "Oops: kfldfdm"
#if 0
          if (name == 'PRES' .AND. param(100) == 2) then
              call h1_overlap_2(z,r,mask)
              call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
              call add2         (z,w,n)
          else
              call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
              if (name == 'PRES') then
                  call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
                  call add2         (z,w,n)
              endif
          endif
#endif
      endif
  
      if (name == 'PRES') then
          call ortho (z)
      elseif (ifmcor) then
          rmean = smean*glsc2(z,bm1,n)
          z = z + rmean
      endif
  !        write(6,*) rmean,ifmcor,' ifmcor'
  
      rtz2=rtz1
      scalar(1)=vlsc3 (z,r,mult,n)
      scalar(2)=vlsc32(r,mult,binv,n)
      call gop(scalar,w,'+  ',2)
      rtz1=scalar(1)
      rbn2=sqrt(scalar(2)/vol)
      if (iter == 1) rbn0 = rbn2
      if (param(22) < 0) tol=abs(param(22))*rbn0
      if (tin < 0)       tol=abs(tin)*rbn0

      ifprint_hmh = .FALSE. 
      if (nid == 0 .AND. ifprint .AND. param(74) /= 0) ifprint_hmh= .TRUE. 
      if (nid == 0 .AND. istep == 1)                 ifprint_hmh= .TRUE. 

      if (ifprint_hmh) &
      write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1),h2(1)


  !        Always take at least one iteration   (for projection) pff 11/23/98
#ifndef TST_WSCAL
      IF (rbn2 <= TOL .AND. (iter > 1 .OR. istep <= 5)) THEN
#else
          iter_max = param(150)
          if (name == 'PRES') iter_max = param(151)
          if (iter > iter_max) then
#endif
          !        IF (rbn2.LE.TOL) THEN
              NITER = ITER-1
          !           IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
              if (nid == 0) &
              write(6,3000) istep,name,niter,rbn2,rbn0,tol
              goto 9999
          ENDIF
      
          beta = rtz1/rtz2
          if (iter == 1) beta=0.0
          p = beta * p + z
          call axhelm (w,p,h1,h2,imsh,isd)
          call dssum  (w)
          w = w * mask(1:n)
      
          rho0 = rho
          rho  = glsc3(w,p,mult,n)
          alpha=rtz1/rho
          alphm=-alpha
          x(1:n) = x(1:n) + alpha * p
          r = r + alphm * w
      
      !        Generate tridiagonal matrix for Lanczos scheme
          if (iter == 1) then
              krylov = krylov+1
              diagt(iter) = rho/rtz1
          elseif (iter <= maxcg) then
              krylov = krylov+1
              diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
              upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
          endif
  enddo

  niter = iter-1

  if (nid == 0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
  3000 format(4x,i7,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
  3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
  3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
  9999 continue
  niterhm = niter
  ifsolv = .FALSE. 
    
    
    !     Call eigenvalue routine for Lanczos scheme:
    !          two work arrays are req'd if you want to save "diag & upper"
    
    !     if (iter.ge.3) then
    !        niter = iter-1
    !        call calc (diagt,upper,w,z,krylov,dmax,dmin)
    !        cond = dmax/dmin
    !        if (nid.eq.0) write(6,6) istep,cond,dmin,dmax,' lambda'
    !     endif
    !   6 format(i9,1p3e12.4,4x,a7)
    
    !     if (n.gt.0) write(6,*) 'quit in cggo'
    !     if (n.gt.0) call exitt
    !     call exitt
  return
end subroutine cggo

!=======================================================================
real(DP) function vlsc32(r,b,m,n)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  real(DP), intent(in) :: r(n),b(n),m(n)

  integer :: i
  vlsc32 = 0.
  do i=1,n
      vlsc32 = vlsc32 + b(i)*m(i)*r(i)*r(i)
  enddo
  return
end function vlsc32

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

#if 0
!  use kinds, only : DP, DP_eps
!  use size_m, only : nx1, ny1, nz1
!  use fdmh1, only : ktype, elsize, dd, ifbhalf, kfldfdm
!  use geom, only : ym1
!  use input, only : if3d, ifaxis


  integer :: nxyz, i1, i2, ie, i3, k1, k2, k3
  real(DP) :: h1b, h2b, vol, vl1, vl2, vl3, den
  real(DP), external :: vlsum, vlsc2


!   Set up diagonal for FDM for each spectral element
  nxyz = nx1*ny1*nz1
  if (if3d) then
      do ie=1,nel
          h1b = vlsum(h1(1,1,1,ie),nxyz)/nxyz
          h2b = vlsum(h2(1,1,1,ie),nxyz)/nxyz
          k1 = ktype(ie,1,kfldfdm)
          k2 = ktype(ie,2,kfldfdm)
          k3 = ktype(ie,3,kfldfdm)
          vol = elsize(1,ie)*elsize(2,ie)*elsize(3,ie)
          vl1 = elsize(2,ie)*elsize(3,ie)/(elsize(1,ie)+DP_eps)
          vl2 = elsize(1,ie)*elsize(3,ie)/(elsize(2,ie)+DP_eps)
          vl3 = elsize(1,ie)*elsize(2,ie)/(elsize(3,ie)+DP_eps)
          do i3=1,nz1
              do i2=1,ny1
                  do i1=1,nx1
                      den = h1b*(vl1*dd(i1,k1) + vl2*dd(i2,k2) + vl3*dd(i3,k3)) &
                      + h2b*vol
                      if (ifbhalf) den = den/vol
                      if (den /= 0) then
                          d(i1,i2,i3,ie) = 1./den
                      else
                          d(i1,i2,i3,ie) = 0.
                      
                      !                 write(6,3) 'd=0:'
                      !    $                 ,h1(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2),dd(i3,k3)
                      !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2,k3
                          3 format(a4,1p4e12.4,8i8)
                      
                      endif
                  enddo
              enddo
          enddo
      enddo
  else
      do ie=1,nel
          if (ifaxis) then
              h1b = vlsc2(h1(1,1,1,ie),ym1(1,1,1,ie),nxyz)/nxyz
              h2b = vlsc2(h2(1,1,1,ie),ym1(1,1,1,ie),nxyz)/nxyz
          else
              h1b = vlsum(h1(1,1,1,ie),nxyz)/nxyz
              h2b = vlsum(h2(1,1,1,ie),nxyz)/nxyz
          endif
          k1 = ktype(ie,1,kfldfdm)
          k2 = ktype(ie,2,kfldfdm)
          vol = elsize(1,ie)*elsize(2,ie)
          vl1 = elsize(2,ie)/elsize(1,ie)
          vl2 = elsize(1,ie)/elsize(2,ie)
          i3=1
          do i2=1,ny1
              do i1=1,nx1
                  den = h1b*( vl1*dd(i1,k1) + vl2*dd(i2,k2) ) &
                  + h2b*vol
                  if (ifbhalf) den = den/vol
                  if (den /= 0) then
                      d(i1,i2,i3,ie) = 1./den
                  !                 write(6,3) 'dn0:'
                  !    $                 ,d(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2)
                  !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2
                  else
                      d(i1,i2,i3,ie) = 0.
                  !                 write(6,3) 'd=0:'
                  !    $                 ,h1(i1,i2,i3,ie),dd(i1,k1),dd(i2,k2)
                  !    $                 ,i1,i2,i3,ie,kfldfdm,k1,k2
                      2 format(a4,1p3e12.4,8i8)
                  endif
              !           write(6,1) ie,i1,i2,k1,k2,'d:',d(i1,i2,i3,ie),vol,vl1,vl2
              !   1       format(5i3,2x,a2,1p4e12.4)
              enddo
          enddo
      enddo
  endif
#endif

  return
end subroutine set_fdm_prec_h1b

!-----------------------------------------------------------------------
!> \brief Solve the generalized eigenvalue problem  A x = lam B x
!!
!! A -- symm.
!! B -- symm., pos. definite
!!
!! "SIZE" is included here only to deduce WDSIZE, the working
!! precision, in bytes, so as to know whether dsygv or ssygv
!! should be called.
subroutine generalev(a,b,lam,n,w)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nid
  use parallel, only : ifdblas
  implicit none

  integer, intent(in) :: n
  real(DP), intent(inout) :: a(n,n), b(n,n)
  real(DP), intent(out)   :: lam(n)
  real(DP), intent(in) :: w(n,n) !>!< unused

  real(DP) :: aa(100),bb(100)

  integer, parameter :: lbw=4*lx1*ly1*lz1*lelv
  real(DP), allocatable :: bw(:)

  integer :: info, ninf, lw

  allocate(bw(lbw))

  lw = n*n
!   write(6,*) 'in generalev, =',info,n,ninf

  call copy(aa,a,100)
  call copy(bb,b,100)

  if (ifdblas) then
      call dsygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
  else
      call ssygv(1,'V','U',n,a,n,b,n,lam,bw,lbw,info)
  endif

  if (info /= 0) then
      if (nid == 0) then
          call outmat2(aa ,n,n,n,'aa  ')
          call outmat2(bb ,n,n,n,'bb  ')
          call outmat2(a  ,n,n,n,'Aeig')
          call outmat2(lam,1,n,n,'Deig')
      endif

      ninf = n-info
      write(6,*) 'Error in generalev, info=',info,n,ninf
      call exitt
  endif

  return
end subroutine generalev

!-----------------------------------------------------------------------
subroutine outmat2(a,m,n,k,name)
  use size_m, only : nid
  use kinds, only : DP
  implicit none
  integer :: m,n,k
  real(DP) :: a(m,n)
  character(4) :: name
  integer :: n2, i, j

  n2 = min(n,8)
  write(6,2) nid,name,m,n,k
  do i=1,m
      write(6,1) nid,name,(a(i,j),j=1,n2)
  enddo
!  1 format(i3,1x,a4,16f6.2)
  1 format(i3,1x,a4,1p8e14.5)
  2 format(/,'Matrix: ',i3,1x,a4,3i8)
  return
end subroutine outmat2

!-----------------------------------------------------------------------
