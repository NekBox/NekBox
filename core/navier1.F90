!-----------------------------------------------------------------------
!> \brief Compute the pressure tolerance
subroutine ctolspl (tolspl,respr)
  use kinds, only : DP
  use size_m, only : lelv, lx2, ly2, lz2
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use geom, only : binvm1, volvm1
  use tstep, only : dt, tolpdf, tolps, prelax
  use soln, only : pmask, vmult
  implicit none

  real(DP) :: tolspl
  REAL(DP) :: RESPR (LX2,LY2,LZ2,LELV)

  real(DP), allocatable :: WORK(:,:,:,:)
  integer :: ntot1
  real(DP) :: rinit, tolmin, tolold
  real(DP), external :: glsc23

  allocate(WORK(lx2,ly2,lz2,lelv))
  NTOT1 = NX1*NY1*NZ1*NELV
  work = respr
  call dssum(work)
  work = work * pmask
  rinit = glsc23(work, binvm1, vmult, ntot1)
  RINIT  = SQRT (rinit/VOLVM1)
  IF (TOLPDF > 0.) THEN
      TOLSPL = TOLPDF
      TOLMIN = TOLPDF
  ELSE
      TOLSPL = TOLPS/DT
      TOLMIN = RINIT*PRELAX
  ENDIF
  IF (TOLSPL < TOLMIN) THEN
      TOLOLD = TOLSPL
      TOLSPL = TOLMIN
      IF (NID == 0) &
      WRITE (6,*) 'Relax the pressure tolerance ',TOLSPL,TOLOLD
  ENDIF

  return
end subroutine ctolspl

!------------------------------------------------------------------------
!> \brief   Orthogonalize the residual in the pressure solver with respect
!! to (1,1,...,1)T  (only if all Dirichlet b.c.).
subroutine ortho (respr)
  use kinds, only : DP, i8, i4
  use size_m, only : lx2, ly2, lz2, lelv, nx2, ny2, nz2, nelv
  use geom, only : ifvcor, ifbcor
  use input, only : ifldmhd
  use parallel, only : nelgv
  use tstep, only : ifield
  implicit none

  real(DP) :: respr (lx2,ly2,lz2,lelv)
  integer(i8) :: ntotg,nxyz2
  integer :: ntot
  real(DP) :: rlam
  real(DP), external :: glsum

  nxyz2 = nx2*ny2*nz2
  ntot  = int(nxyz2*nelv, kind=i4)
  ntotg = nxyz2*nelgv

  if (ifield == 1) then
      if (ifvcor) then
          rlam  = glsum (respr,ntot)/ntotg
          respr = respr - rlam
      endif
  elseif (ifield == ifldmhd) then
      if (ifbcor) then
          rlam = glsum (respr,ntot)/ntotg
          respr = respr - rlam
      endif
  else
      call exitti('ortho: unaccounted ifield = $',ifield)
  endif

  return
end subroutine ortho

!---------------------------------------------------------------------
!> \brief Compute OUTi = Di*INP, i=1,2,3.
!! the gradient of the scalar field INP.
!! Note: OUTi is defined on the pressure mesh !!!
!---------------------------------------------------------------------
subroutine opgrad (out1,out2,out3,inp)
  use kinds, only : DP
  use size_m, only : nx2, ny2, nz2, ndim, nelv
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, nx1, ny1, nz1
  use geom, only : rxm2, rym2, rzm2, sxm2, sym2, szm2, txm2, tym2, tzm2
  use geom, only : bm1, jacmi
  use input, only : ifsplit, ifaxis
  use mesh, only : if_ortho
  use dxyz, only : dxm12, dytm12, dztm12
  use ixyz, only : ixm12, iytm12, iztm12
  implicit none

  REAL(DP) :: OUT1 (LX2,LY2,LZ2,nelv)
  REAL(DP) :: OUT2 (LX2,LY2,LZ2,nelv)
  REAL(DP) :: OUT3 (LX2,LY2,LZ2,nelv)
  REAL(DP) :: INP  (LX1,LY1,LZ1,nelv)

  real(DP) ::  ta1 (lx1*ly1*lz1) &
  ,             ta2 (lx1*ly1*lz1) &
  ,             ta3 (lx1*ly1*lz1)

  integer :: iflg, ntot2, i1, i2, e, iz, n1, n2, nxy2, nyz1
  iflg = 0

  if (ifsplit .AND. .NOT. ifaxis) then
      call wgradm1(out1,out2,out3,inp,nelv) ! weak grad on FLUID mesh
      return
  endif
  write(*,*) "Oops! opgrad"


  return
end subroutine opgrad


!---------------------------------------------------------------------
!> \brief Compute D*X .
!!  X    : input variable, defined on M1
!!  DX   : output variable, defined on M2 (note: D is rectangular)
!!  RM2 : RXM2, RYM2 or RZM2
!!  SM2 : SXM2, SYM2 or SZM2
!!  TM2 : TXM2, TYM2 or TZM2
!!  ISD : spatial direction (x=1,y=2,z=3)
!!  IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)
!---------------------------------------------------------------------
subroutine multd (dx,x,rm2,sm2,tm2,isd,iflg)
  use kinds, only : DP
  use size_m, only : lx2, ly2, lz2, lx1, ly1, lz1, lelv
  use size_m, only : nx1, ny1, nz1, nx2, ny2, nz2, nelv, ndim
  use ctimer, only : icalld, tmltd, nmltd, etime1, dnekclock
  use dxyz, only : dytm12, datm12, dctm12, dxm12, dztm12
  use geom, only : ifrzer, jacm1
  use input, only : ifaxis, ifsplit
  use ixyz, only : iytm12, iatm12, ictm12, iztm12, ixm12
  use geom, only : bm1
  use wz_m, only : w3m2, w2am2, w2cm2
  implicit none

  integer :: isd, iflg
  real(DP) ::           dx   (lx2,ly2,lz2,lelv)
  real(DP) ::           x    (lx1,ly1,lz1,lelv)
  real(DP) ::           rm2  (lx2,ly2,lz2,lelv)
  real(DP) ::           sm2  (lx2,ly2,lz2,lelv)
  real(DP) ::           tm2  (lx2,ly2,lz2,lelv)

  real(DP) ::  ta1 (lx1*ly1*lz1) &
  ,             ta2 (lx1*ly1*lz1) &
  ,             ta3 (lx1*ly1*lz1)

  integer :: e
  integer :: nxyz1, nxy2, nxyz2, n1, n2, ny12, i1, i2, iz, nyz1

#ifndef NOTIMER
  if (icalld == 0) tmltd=0.0
  icalld=icalld+1
  nmltd=icalld
  etime1=dnekclock()
#endif

  nyz1  = ny1*nz1
  nxy2  = nx2*ny2
  nxyz1 = nx1*ny1*nz1
  nxyz2 = nx2*ny2*nz2

  n1    = nx2*ny1
  n2    = nx2*ny2

  do e=1,nelv


          call mxm (dxm12,nx2,x(:,:,:,e),nx1,ta1,nyz1)
          i1=1
          i2=1
          do iz=1,nz1
              call mxm (ta1(i1),nx2,iytm12,ny1,ta2(i2),ny2)
              i1=i1+n1
              i2=i2+n2
          enddo
          call mxm  (ta2,nxy2,iztm12,nz1,dx(:,:,:,e),nz2)
          dx(:,:,:,e) = dx(:,:,:,e) * rm2(:,:,:,e)

          call mxm  (ixm12,nx2,x(:,:,:,e),nx1,ta3,nyz1) ! reuse ta3 below
          i1=1
          i2=1
          do iz=1,nz1
              call mxm (ta3(i1),nx2,dytm12,ny1,ta2(i2),ny2)
              i1=i1+n1
              i2=i2+n2
          enddo
          call mxm     (ta2,nxy2,iztm12,nz1,ta1,nz2)
          dx(:,:,:,e) = dx(:,:,:,e) + reshape(ta1,(/lx2,ly2,lz2/)) * sm2(:,:,:,e)

      !            call mxm (ixm12,nx2,x(1,e),nx1,ta1,nyz1) ! reuse ta3 from above
          i1=1
          i2=1
          do iz=1,nz1
              call mxm (ta3(i1),nx2,iytm12,ny1,ta2(i2),ny2)
              i1=i1+n1
              i2=i2+n2
          enddo
          call mxm (ta2,nxy2,dztm12,nz1,ta3,nz2)
          dx(:,:,:,e) = dx(:,:,:,e) + reshape(ta3,(/lx2,ly2,lz2/)) * tm2(:,:,:,e)
      !           endif
  
  !        Collocate with the weights on the pressure mesh
          dx(:,:,:,e) = dx(:,:,:,e) * bm1(:,:,:,e)
          dx(:,:,:,e) = dx(:,:,:,e) / jacm1(:,:,:,e)
  
  enddo

#ifndef NOTIMER
  tmltd=tmltd+(dnekclock()-etime1)
#endif
  return
end subroutine multd

!----------------------------------------------------------------------
!> \brief OUT = (H1*A+H2*B) * INP
!----------------------------------------------------------------------
subroutine ophx (out1,out2,out3,inp1,inp2,inp3,h1,h2)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, ndim
  use input, only : ifstrs
  implicit none

  REAL(DP) :: OUT1 (LX1,LY1,LZ1,1)
  REAL(DP) :: OUT2 (LX1,LY1,LZ1,1)
  REAL(DP) :: OUT3 (LX1,LY1,LZ1,1)
  REAL(DP) :: INP1 (LX1,LY1,LZ1,1)
  REAL(DP) :: INP2 (LX1,LY1,LZ1,1)
  REAL(DP) :: INP3 (LX1,LY1,LZ1,1)
  REAL(DP) :: H1   (LX1,LY1,LZ1,1)
  REAL(DP) :: H2   (LX1,LY1,LZ1,1)

  integer :: imesh
  IMESH = 1

  IF (IFSTRS) THEN
#if 0
      MATMOD = 0
      CALL AXHMSF (OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2,MATMOD)
#endif
  ELSE
      CALL AXHELM (OUT1,INP1,H1,H2,IMESH,1)
      CALL AXHELM (OUT2,INP2,H1,H2,IMESH,2)
      IF (NDIM == 3) &
      CALL AXHELM (OUT3,INP3,H1,H2,IMESH,3)
  ENDIF

  return
end subroutine ophx

!--------------------------------------------------------------
!> \brief Compute some derviatives?
!!   DU   - dU/dx or dU/dy or dU/dz
!!   U    - a field variable defined on mesh 1
!!   RM1  - dr/dx or dr/dy or dr/dz
!!   SM1  - ds/dx or ds/dy or ds/dz
!!   TM1  - dt/dx or dt/dy or dt/dz
!!   JM1  - the Jacobian
!!   IMESH - topology: velocity (1) or temperature (2) mesh
!--------------------------------------------------------------
subroutine dudxyz (du,u,rm1,sm1,tm1,jm1,imsh,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1
  use size_m, only : nx1, ny1, nz1, nelv, nelt, ndim
  use dxyz, only : dxm1, dytm1, dztm1
  use geom, only : jacmi
  implicit none

  integer :: imsh, isd
  REAL(DP) ::  DU  (LX1,LY1,LZ1,*)
  REAL(DP) ::  U   (LX1,LY1,LZ1,*)
  REAL(DP) ::  RM1 (LX1,LY1,LZ1,*)
  REAL(DP) ::  SM1 (LX1,LY1,LZ1,*)
  REAL(DP) ::  TM1 (LX1,LY1,LZ1,*)
  REAL(DP) ::  JM1 (LX1,LY1,LZ1,*)

  REAL(DP) ::  DRST(LX1,LY1,LZ1)
  integer :: nel, nxy1, nyz1, nxyz1, ntot, iel, iz

  nel = -1
  IF (imsh == 1) NEL = NELV
  IF (imsh == 2) NEL = NELT
  NXY1  = NX1*NY1
  NYZ1  = NY1*NZ1
  NXYZ1 = NX1*NY1*NZ1
  NTOT  = NXYZ1*NEL

  !> \todo why this loop?
  DO 1000 IEL=1,NEL
  
!max      IF (IFAXIS) CALL SETAXDY (IFRZER(IEL) )
  
      IF (NDIM == 2) THEN
          CALL MXM     (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
          du(:,:,:,iel) = du(:,:,:,iel) * rm1(:,:,:,iel)
          CALL MXM     (U(1,1,1,IEL),NX1,DYTM1,NY1,DRST,NY1)
          du(:,:,:,iel) = du(:,:,:,iel) + drst * sm1(:,:,:,iel)
      ELSE
          CALL MXM   (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
          du(:,:,:,iel) = du(:,:,:,iel) * rm1(:,:,:,iel)
          DO 20 IZ=1,NZ1
              CALL MXM  (U(1,1,IZ,IEL),NX1,DYTM1,NY1,DRST(1,1,IZ),NY1)
          20 END DO
          du(:,:,:,iel) = du(:,:,:,iel) + drst * sm1(:,:,:,iel)
          CALL MXM     (U(1,1,1,IEL),NXY1,DZTM1,NZ1,DRST,NZ1)
          du(:,:,:,iel) = du(:,:,:,iel) + drst * tm1(:,:,:,iel)
      ENDIF
  
  1000 END DO

!    CALL INVCOL2 (DU,JM1,NTOT)
  du(:,:,:,1:nel) = du(:,:,:,1:nel) * jacmi

  return
end subroutine dudxyz

!---------------------------------------------------------------------
!> \brief Compute and add: (1) user specified forcing function (FX,FY,FZ)
!!                  (2) driving force due to natural convection
!!                  (3) convection term
!! !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
!!          current time step is completed.
!----------------------------------------------------------------------
subroutine makef
  use kinds, only : DP
  use input, only : ifnav, ifchar, iftran
  use ctimer, only : tmakef, nmakef, dnekclock
  implicit none
  real(DP) :: etime

  nmakef = nmakef + 1
  etime = dnekclock()
  call makeuf
!  if (ifnatc)                               call natconv
!  if (ifexplvis .AND. ifsplit)                call explstrs
  etime = etime - dnekclock()
  if (ifnav .AND. ( .NOT. ifchar))              call advab
  etime = etime + dnekclock()
!  if (ifmvbd)                               call admeshv
  if (iftran) then
    call makeabf
  endif
  if ((iftran .AND. .NOT. ifchar) .OR. &
  (iftran .AND. .NOT. ifnav .AND. ifchar))   call makebdf
!max    if (ifnav .AND. ifchar .AND. ( .NOT. ifmvbd))   call advchar
#if 0
  if (ifmodel)                              call twallsh
#endif
  tmakef = tmakef + (dnekclock() - etime)
  return
end subroutine makef

!---------------------------------------------------------------------
!> \brief Compute and add: (1) user specified forcing function (FX,FY,FZ)
!----------------------------------------------------------------------
subroutine makeuf
  use size_m, only : nelv, nx1, ny1, nz1
  use geom, only : bm1
  use soln, only : bfx, bfy, bfz
  use tstep, only : time, dt
  use ctimer, only : othr_flop, othr_mop
  implicit none

  integer :: i
  TIME = TIME-DT
  othr_flop = othr_flop + 3*nx1*ny1*nz1*nelv
  othr_mop  = othr_mop  + 6*nx1*ny1*nz1*nelv
  CALL NEKUF   (BFX,BFY,BFZ)
  do i = 1, nelv
    bfx(:,:,:,i) = bfx(:,:,:,i) * bm1(:,:,:,i)
    bfy(:,:,:,i) = bfy(:,:,:,i) * bm1(:,:,:,i)
    bfz(:,:,:,i) = bfz(:,:,:,i) * bm1(:,:,:,i)
  enddo
  TIME = TIME+DT

  return
end subroutine makeuf

subroutine nekuf (f1,f2,f3)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv
  use nekuse, only : ffx, ffy, ffz
  use parallel, only : lglel
  implicit none

  REAL(DP) :: F1 (LX1,LY1,LZ1,LELV)
  REAL(DP) :: F2 (LX1,LY1,LZ1,LELV)
  REAL(DP) :: F3 (LX1,LY1,LZ1,LELV)

  integer :: i, j, k, ielg, iel
  DO IEL=1,NELV
      ielg = lglel(iel)
      DO K=1,NZ1
          DO J=1,NY1
              DO I=1,NX1
                  CALL NEKASGN (I,J,K,IEL)
                  CALL USERF   (I,J,K,IELG)
                  F1(I,J,K,IEL) = FFX
                  F2(I,J,K,IEL) = FFY
                  F3(I,J,K,IEL) = FFZ
              enddo
          enddo
      enddo
  END DO

  return
end subroutine nekuf

!---------------------------------------------------------------
!> \brief Eulerian scheme, add convection term to forcing function
!! at current time step.
!---------------------------------------------------------------
subroutine advab()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use geom, only : bm1
  use soln, only : vx, vy, vz, bfx, bfy, bfz
  use ctimer, only : othr_flop, othr_mop
  implicit none

  real(DP), allocatable:: TA (:,:,:,:)

  integer :: ntot1
  NTOT1 = NX1*NY1*NZ1*NELV

  othr_flop = othr_flop + 6*ntot1
  othr_mop  = othr_mop + 12*ntot1

  allocate(TA(nx1,ny1,nz1,nelv)) 
  CALL CONVOP  (TA,VX)
  bfx = bfx - bm1*ta
  CALL CONVOP  (TA,VY)
  bfy = bfy - bm1*ta
  IF (NDIM /= 2) THEN
      CALL CONVOP  (TA,VZ)
      bfz = bfz - bm1*ta
  ENDIF

  return
end subroutine advab

!-----------------------------------------------------------------------
!> \brief Add contributions to F from lagged BD terms.
subroutine makebdf()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv
  use geom, only : ifgeom
  use geom, only : bm1, bm1lag
  use soln, only : vx, vy, vz, vxlag, vylag, vzlag, bfx, bfy, bfz
  use soln, only : vtrans
  use tstep, only : dt, ifield, bd, nbd
  use ctimer, only : othr_flop, othr_mop
  implicit none

  real(DP), allocatable :: H2 (:,:,:)

  integer :: ilag, i

  allocate(H2(nx1,ny1,nz1))

  othr_flop = othr_flop + (12*nbd+1)*nelv*nx1*ny1*nz1
  othr_mop  = othr_mop  + (10*nbd+1)*nelv*nx1*ny1*nz1 
  do i = 1, nelv 
    h2 = (1./DT) * vtrans(:,:,:,i,ifield) 
    DO ILAG=2,NBD
      IF (IFGEOM) THEN
        bfx(:,:,:,i) = bfx(:,:,:,i) + h2 * (vxlag(:,:,:,i,ilag-1) * bm1lag(:,:,:,i,ilag-1) * bd(ilag+1))
        bfy(:,:,:,i) = bfy(:,:,:,i) + h2 * (vylag(:,:,:,i,ilag-1) * bm1lag(:,:,:,i,ilag-1) * bd(ilag+1))
        bfz(:,:,:,i) = bfz(:,:,:,i) + h2 * (vzlag(:,:,:,i,ilag-1) * bm1lag(:,:,:,i,ilag-1) * bd(ilag+1))
      ELSE
        bfx(:,:,:,i) = bfx(:,:,:,i) + h2 * (vxlag(:,:,:,i,ilag-1) * bm1(:,:,:,i)    * bd(ilag+1))
        bfy(:,:,:,i) = bfy(:,:,:,i) + h2 * (vylag(:,:,:,i,ilag-1) * bm1(:,:,:,i)    * bd(ilag+1))
        bfz(:,:,:,i) = bfz(:,:,:,i) + h2 * (vzlag(:,:,:,i,ilag-1) * bm1(:,:,:,i)    * bd(ilag+1))
      ENDIF
    END DO
    bfx(:,:,:,i) = bfx(:,:,:,i) + h2*vx(:,:,:,i)*bm1(:,:,:,i)*bd(2)
    bfy(:,:,:,i) = bfy(:,:,:,i) + h2*vy(:,:,:,i)*bm1(:,:,:,i)*bd(2)
    bfz(:,:,:,i) = bfz(:,:,:,i) + h2*vz(:,:,:,i)*bm1(:,:,:,i)*bd(2)
  enddo

  return
end subroutine makebdf

!-----------------------------------------------------------------------
!> \brief Sum up contributions to kth order extrapolation scheme.
!! NOTE: rho^{n+1} should multiply all the Sum_q{beta_q} term
!!       if rho is not constant!
!-----------------------------------------------------------------------
subroutine makeabf
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use soln, only : abx1, aby1, abz1, abx2, aby2, abz2, bfx, bfy, bfz, vtrans
  use tstep, only : ab
  use ctimer, only : othr_flop, othr_mop
  implicit none

  real(DP), allocatable :: TA (:,:,:,:) 

  integer :: ntot1, iel
  real(DP) :: ab0, ab1, ab2

  allocate(TA(nx1,ny1,nz1,nelv))

  NTOT1 = NX1*NY1*NZ1*NELV

  AB0 = AB(1)
  AB1 = AB(2)
  AB2 = AB(3)

#if 0
! 11*ntot mops
! 6*ntot flops
  ta = ab1 * abx1 + ab2 * abx2
  CALL COPY   (ABX2,ABX1,NTOT1)
  CALL COPY   (ABX1,BFX,NTOT1)
  bfx = (ab0 * bfx + ta) * vtrans(:,:,:,:,1) ! multiply by density

  ta = ab1 * aby1 + ab2 * aby2
  CALL COPY   (ABY2,ABY1,NTOT1)
  CALL COPY   (ABY1,BFY,NTOT1)
  bfy = (ab0 * bfy + ta) * vtrans(:,:,:,:,1) ! multiply by density

  IF (NDIM == 3) THEN
    ta = ab1 * abz1 + ab2 * abz2
    CALL COPY   (ABZ2,ABZ1,NTOT1)
    CALL COPY   (ABZ1,BFZ,NTOT1)
    bfz = (ab0 * bfz + ta) * vtrans(:,:,:,:,1) ! multiply by density
  ENDIF
#else
! 7*ntot mops
! 6*ntot flops

  othr_flop = othr_flop + 18*ntot1
  othr_mop  = othr_mop + 21*ntot1

  do iel = 1, nelv
    ta(:,:,:,1) = ab1 * abx1(:,:,:,iel) + ab2 * abx2(:,:,:,iel)
    abx2(:,:,:,iel) = abx1(:,:,:,iel)
    abx1(:,:,:,iel) = bfx(:,:,:,iel)
    bfx(:,:,:,iel) = (ab0*bfx(:,:,:,iel) + ta(:,:,:,1)) * vtrans(:,:,:,iel,1)
  enddo

  do iel = 1, nelv
    ta(:,:,:,1) = ab1 * aby1(:,:,:,iel) + ab2 * aby2(:,:,:,iel)
    aby2(:,:,:,iel) = aby1(:,:,:,iel)
    aby1(:,:,:,iel) = bfy(:,:,:,iel)
    bfy(:,:,:,iel) = (ab0*bfy(:,:,:,iel) + ta(:,:,:,1)) * vtrans(:,:,:,iel,1)
  enddo

  do iel = 1, nelv
    ta(:,:,:,1) = ab1 * abz1(:,:,:,iel) + ab2 * abz2(:,:,:,iel)
    abz2(:,:,:,iel) = abz1(:,:,:,iel)
    abz1(:,:,:,iel) = bfz(:,:,:,iel)
    bfz(:,:,:,iel) = (ab0*bfz(:,:,:,iel) + ta(:,:,:,1)) * vtrans(:,:,:,iel,1)
  enddo
#endif

  return
end subroutine makeabf

!-----------------------------------------------------------------------
!> \brief Compute Adams-Bashforth coefficients (order NAB, less or equal to 3).
!! NBD .EQ. 1
!! Standard Adams-Bashforth coefficients
!! NBD .GT. 1
!! Modified Adams-Bashforth coefficients to be used in con-
!! junction with Backward Differentiation schemes (order NBD)
!-----------------------------------------------------------------------
subroutine setabbd (ab,dtlag,nab,nbd)
  use kinds, only : DP
  implicit none

  REAL(DP), intent(out) :: AB(NAB)    !>!< Adams-Bashforth coefficients
  real(DP), intent(in)  :: DTLAG(nbd) !>!< Time-step history
  integer,  intent(in)  :: nab        !>!< Order of AB scheme
  integer,  intent(in)  :: nbd        !>!< Order of accompanying BDF scheme

  real(DP) :: dt0, dt1, dt2, dta, dts, dtb, dtc, dtd, dte

  IF ( NAB == 1 ) THEN
  
      AB(1) = 1.0
  
  ELSEIF ( NAB == 2 ) THEN
      DT0 = DTLAG(1)
      DT1 = DTLAG(2)
  
      DTA =  DT0/DT1
  
      IF ( NBD == 1 ) THEN
      
          AB(2) = -0.5*DTA
          AB(1) =  1.0 - AB(2)
      
      ELSEIF ( NBD == 2 ) THEN
      
          AB(2) = -DTA
          AB(1) =  1.0 - AB(2)
      
      ENDIF
  
  ELSEIF ( NAB == 3 ) THEN
      DT0 = DTLAG(1)
      DT1 = DTLAG(2)  
      DT2 = DTLAG(3)

      DTS =  DT1 + DT2
      DTA =  DT0 / DT1
      DTB =  DT1 / DT2
      DTC =  DT0 / DT2
      DTD =  DTS / DT1
      DTE =  DT0 / DTS
  
      IF ( NBD == 1 ) THEN
      
          AB(3) =  DTE*( 0.5*DTB + DTC/3. )
          AB(2) = -0.5*DTA - AB(3)*DTD
          AB(1) =  1.0 - AB(2) - AB(3)
      
      ELSEIF ( NBD == 2 ) THEN
      
          AB(3) =  2./3.*DTC*(1./DTD + DTE)
          AB(2) = -DTA - AB(3)*DTD
          AB(1) =  1.0 - AB(2) - AB(3)
      
      ELSEIF ( NBD == 3 ) THEN
      
          AB(3) =  DTE*(DTB + DTC)
          AB(2) = -DTA*(1.0 + DTB + DTC)
          AB(1) =  1.0 - AB(2) - AB(3)
      
      ENDIF
  
  ENDIF

  return
end subroutine setabbd

!-----------------------------------------------------------------------
!> \brief Compute backwards-difference (BDF) coefficients, order NBD
!-----------------------------------------------------------------------
subroutine setbd (bd,dtbd,nbd)
  use kinds, only : DP
  implicit none

  REAL(dp), intent(out) :: BD(*)   !>!< BDF coefficients
  real(dp), intent(in)  :: DTBD(*) !>!< Time-step history
  integer , intent(in)  :: nbd     !>!< Order of BDF scheme

  integer, PARAMETER :: NDIM = 10
  REAL(DP) :: BDMAT(NDIM,NDIM),BDRHS(NDIM), BDF
  INTEGER :: IR(NDIM),IC(NDIM)
  integer :: nsys, i, ibd

  BD(1:ndim) = 0._dp; bdf = -1
  ! BDF(1) is trivial
  IF (NBD == 1) THEN
      BD(1) = 1.
      BDF   = 1.
  ! BDF(>1) computed using a linear system
  ELSEIF (NBD >= 2) THEN
      NSYS = NBD+1
      CALL BDSYS (BDMAT,BDRHS,DTBD,NBD,NDIM)
      CALL LU    (BDMAT,NSYS,NDIM,IR,IC)
      CALL SOLVE (BDRHS,BDMAT,1,NSYS,NDIM,IR,IC)
      DO 30 I=1,NBD
          BD(I) = BDRHS(I)
      30 END DO
      BDF = BDRHS(NBD+1)
  ENDIF

  !   Normalize
  DO IBD=NBD,1,-1
      BD(IBD+1) = BD(IBD)
  END DO
  BD(1) = 1.
  DO IBD=1,NBD+1
      BD(IBD) = BD(IBD)/BDF
  END DO
!   write(6,1) (bd(k),k=1,nbd+1)
! 1 format('bd:',1p8e13.5)

  return
end subroutine setbd

!> Setup the linear system that defines BDF coefficients
subroutine bdsys (a,b,dt,nbd,ndim)
  use kinds, only : DP
  implicit none

  integer :: nbd, ndim
  REAL(DP) :: A(NDIM,9),B(9),DT(9)

  integer :: n, j, k, i
  real(DP) :: sumdt
  a = 0._dp
  N = NBD+1
  DO J=1,NBD
      A(1,J) = 1.
  END DO
  A(1,NBD+1) = 0.
  B(1) = 1.
  DO J=1,NBD
      SUMDT = 0.
      DO K=1,J
          SUMDT = SUMDT+DT(K)
      END DO
      A(2,J) = SUMDT
  END DO
  A(2,NBD+1) = -DT(1)
  B(2) = 0.
  DO I=3,NBD+1
      DO J=1,NBD
          SUMDT = 0.
          DO K=1,J
              SUMDT = SUMDT+DT(K)
          END DO
          A(I,J) = SUMDT**(I-1)
      END DO
      A(I,NBD+1) = 0.
      B(I) = 0.
  END DO
  return
end subroutine bdsys

!-------------------------------------------------------------------
!> \brief Set initial time for subintegration
!-------------------------------------------------------------------
subroutine tauinit (tau,ilag)
  use kinds, only : DP
  use tstep, only : nbd, dtlag
  implicit none

  real(DP) :: tau
  integer :: ilag
  integer :: i

  TAU   = 0.
  DO 10 I=NBD,ILAG+1,-1
      TAU = TAU+DTLAG(I)
  10 END DO
  return
end subroutine tauinit

!-----------------------------------------------------------------------
!> \brief Keep old velocity field(s)
!-----------------------------------------------------------------------
subroutine lagvel
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use soln, only : vxlag, vylag, vzlag, vx, vy, vz
  implicit none

  integer :: ntot1, ilag
  NTOT1 = NX1*NY1*NZ1*NELV

!    DO 100 ILAG=NBDINP-1,2,-1
  DO 100 ILAG=3-1,2,-1
      CALL COPY (VXLAG (1,1,1,1,ILAG),VXLAG (1,1,1,1,ILAG-1),NTOT1)
      CALL COPY (VYLAG (1,1,1,1,ILAG),VYLAG (1,1,1,1,ILAG-1),NTOT1)
      IF (NDIM == 3) &
      CALL COPY (VZLAG (1,1,1,1,ILAG),VZLAG (1,1,1,1,ILAG-1),NTOT1)
  100 END DO

  CALL OPCOPY (VXLAG,VYLAG,VZLAG,VX,VY,VZ)

  return
end subroutine lagvel

!----------------------------------------------------------------------
!> \brief Set up parameters for backward differentiation scheme.
!----------------------------------------------------------------------
subroutine setordbd
  use tstep, only : nbdinp, nbd, istep
  implicit none

!   IF (IFSPLIT .OR. NBDINP.EQ.0) THEN     undid hardwire, 3/6/92 pff
  IF ( NBDINP < 1) THEN
      NBD = 1
  ELSE
      IF ((ISTEP == 0) .OR. (ISTEP == 1))        NBD = 1
      IF ((ISTEP > 1) .AND. (ISTEP <= NBDINP))  NBD = ISTEP
      IF (ISTEP > NBDINP)                     NBD = NBDINP
  ENDIF

  return
end subroutine setordbd

!---------------------------------------------------------------
!> \brief Compute error norms of a (scalar) field variable X
!! defined on mesh 1 or mesh 2.
!! The error norms are normalized with respect to the volume
!! (except for Linf).
!!
!! \f$ l2 = |x|_2 \f$
!! \f$ semi = \sqrt{|<x|A|x>|/V} \f$ where A is the stiffness matrix
!---------------------------------------------------------------
subroutine normsc (h1,semi,l2,linf,x,imesh)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelv, nelt, ndim
  use geom, only : bm1, voltm1, volvm1
  use ctimer, only : tnmsc, nnmsc, dnekclock
  implicit none

  real(DP) :: h1, semi, l2, linf
  real(DP), intent(in) :: X  (LX1,LY1,LZ1,*)
  integer :: imesh

  real(DP), allocatable, dimension(:,:,:,:) :: y, ta1, ta2
  REAL(DP) :: LENGTH, vol
  integer :: nel, nxyz1, ntot1
  real(DP), external :: glamax, glsum
  real(DP) :: etime

  nnmsc = nnmsc + 1
  etime = dnekclock()

  allocate(Y(LX1,LY1,LZ1,LELT), TA1(LX1,LY1,LZ1,LELT), TA2(LX1,LY1,LZ1,LELT))

  IF (IMESH == 1) THEN
      NEL = NELV
      VOL = VOLVM1
  ELSEIF (IMESH == 2) THEN
      NEL = NELT
      VOL = VOLTM1
  else
    nel = -1
    vol = -1.
    write(*,*) "IMESH \notin {1,2}"
  ENDIF

  LENGTH = VOL**(1./(NDIM))
  NXYZ1  = NX1*NY1*NZ1
  NTOT1  = NXYZ1*NEL

  H1     = 0.
  SEMI   = 0.
  L2     = 0.
  LINF   = 0.

  LINF = GLAMAX (X,NTOT1)

  ta1 = x(:,:,:,1:nel)*x(:,:,:,1:nel) * bm1
  L2   = GLSUM  (TA1,NTOT1)
  IF (L2 < 0.0) L2 = 0.

  ta1 = 1._dp; ta2 = 0._dp
  etime = etime - dnekclock()
  CALL AXHELM (Y,X,TA1,TA2,IMESH,1)
  etime = etime + dnekclock()
  ta1 = y * x(:,:,:,1:nel)
  SEMI = GLSUM  (TA1,NTOT1)
  IF (SEMI < 0.0) SEMI = 0.

  H1   = SQRT((SEMI*LENGTH**2+L2)/VOL)
  SEMI = SQRT(SEMI/VOL)
  L2   = SQRT(L2/VOL)
  IF (H1 < 0.) H1 = 0.
  tnmsc = tnmsc + (dnekclock() - etime)

  return
end subroutine normsc

!---------------------------------------------------------------
!> \brief Compute error norms of a (vector) field variable (X1,X2,X3)
! defined on mesh 1 (velocity mesh only !)
! The error norms are normalized with respect to the volume
! (except for Linf).
!---------------------------------------------------------------
subroutine normvc (h1,semi,l2,linf,x1,x2,x3)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use size_m, only : lx1, ly1, lz1, lelt
  use geom, only : volvm1, bm1
  use ctimer, only : tnmvc, nnmvc, dnekclock
  implicit none

  REAL(DP) :: H1,SEMI,L2,LINF
  REAL(DP) :: X1 (LX1,LY1,LZ1,lelt)
  REAL(DP) :: X2 (LX1,LY1,LZ1,lelt)
  REAL(DP) :: X3 (LX1,LY1,LZ1,lelt)

  real(DP), allocatable, dimension(:,:,:,:) :: Y1, Y2, Y3, TA1, TA2
  REAL(DP) :: LENGTH
  integer :: imesh, nel, nxyz1, ntot1
  real(DP) :: vol
  real(DP), external :: glamax, glsum
  real(DP) :: etime

  nnmvc = nnmvc + 1
  etime = dnekclock()

  IMESH  = 1
  NEL    = NELV
  VOL    = VOLVM1
  LENGTH = VOL**(1./(NDIM))
  NXYZ1  = NX1*NY1*NZ1
  NTOT1  = NXYZ1*NEL

  H1     = 0.
  SEMI   = 0.
  L2     = 0.
  LINF   = 0.

  allocate(Y1 (LX1,LY1,LZ1,LELT) &
  ,Y2 (LX1,LY1,LZ1,LELT) &
  ,Y3 (LX1,LY1,LZ1,LELT) &
  ,TA1(LX1,LY1,LZ1,LELT) &
  ,TA2(LX1,LY1,LZ1,LELT) )


  IF (NDIM == 3) THEN
    ta1 = x1*x1 + x2*x2 + x3*x3
  else
    ta1 = x1*x1 + x2*x2
  ENDIF
  LINF = GLAMAX (TA1,NTOT1)
  LINF = SQRT( LINF )

  ta1 = ta1 * bm1
  L2 = GLSUM  (TA1,NTOT1)
  IF (L2 < 0.0) L2 = 0.

  ta1 = 1._dp; ta2 = 0._dp
  etime = etime - dnekclock()
  CALL OPHX  (Y1,Y2,Y3,X1,X2,X3,TA1,TA2)
  etime = etime + dnekclock()
  IF (NDIM == 3) THEN
    ta1 = y1*x1 + y2*x2 + y3*x3
  else
    ta1 = y1*x1 + y2*x2
  ENDIF

  SEMI = GLSUM (TA1,NTOT1)
  IF (SEMI < 0.0) SEMI = 0.

  H1   = SQRT((SEMI*LENGTH**2+L2)/VOL)
  SEMI = SQRT(SEMI/VOL)
  L2   = SQRT(L2  /VOL)
  IF (H1 < 0.) H1 = 0.

  tnmvc = tnmvc + (dnekclock() - etime)

  return
end subroutine normvc

subroutine opcolv (a1,a2,a3,c)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim
  use opctr, only : isclld, nrout, myrout, rname, dct, ncall, dcount
  implicit none
  REAL(DP) :: A1(1),A2(1),A3(1),C(1)
  integer :: ntot1, i, isbcnt

  NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
  if (isclld == 0) then
      isclld=1
      nrout=nrout+1
      myrout=nrout
      rname(myrout) = 'opcolv'
  endif

  isbcnt = ntot1*ndim
  dct(myrout) = dct(myrout) + (isbcnt)
  ncall(myrout) = ncall(myrout) + 1
  dcount      =      dcount + (isbcnt)
#endif

  IF (NDIM == 3) THEN
      DO 100 I=1,NTOT1
          A1(I)=A1(I)*C(I)
          A2(I)=A2(I)*C(I)
          A3(I)=A3(I)*C(I)
      100 END DO
  ELSE
      DO 200 I=1,NTOT1
          A1(I)=A1(I)*C(I)
          A2(I)=A2(I)*C(I)
      200 END DO
  ENDIF
  return
end subroutine opcolv

subroutine opcopy (a1,a2,a3,b1,b2,b3)
  use kinds, only : DP
  use size_m
  implicit none
  REAL(DP) :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
  integer :: ntot1
  NTOT1=NX1*NY1*NZ1*NELV
  CALL COPY(A1,B1,NTOT1)
  CALL COPY(A2,B2,NTOT1)
  IF(NDIM == 3)CALL COPY(A3,B3,NTOT1)
  return
end subroutine opcopy

!-----------------------------------------------------------------------
subroutine opdssum (a,b,c)! NOTE: opdssum works on FLUID/MHD arrays only!
  use kinds, only : DP
  use input, only : ifcyclic
  implicit none

  real(DP) :: a(1),b(1),c(1)

  if (ifcyclic) then
    write(*,*) "Oops: ifcyclic"
#if 0
      call rotate_cyc  (a,b,c,1)
      call vec_dssum   (a,b,c,nx1,ny1,nz1)
      call rotate_cyc  (a,b,c,0)
#endif
  else
      call vec_dssum   (a,b,c)
  endif

  return
end subroutine opdssum

!-----------------------------------------------------------------------
subroutine opdsop (a,b,c,op)! opdsop works on FLUID/MHD arrays only!
  use kinds, only : DP
  use input, only : ifcyclic
  implicit none

  real(DP) :: a(1),b(1),c(1)
  character(3) :: op

  if (ifcyclic) then
      if (op == '*  ' .OR. op == 'mul' .OR. op == 'MUL') then
          call vec_dsop    (a,b,c,op)
      else
          write(*,*) "Oops: op"
#if 0
          call rotate_cyc  (a,b,c,1)
          call vec_dsop    (a,b,c,nx1,ny1,nz1,op)
          call rotate_cyc  (a,b,c,0)
#endif
      endif
  else
      call vec_dsop    (a,b,c,op)
  endif

  return
end subroutine opdsop

!-----------------------------------------------------------------------
subroutine transpose(a,lda,b,ldb)
  use kinds, only : DP
  implicit none
  integer :: lda, ldb
  real(DP) :: a(lda,*),b(ldb,*)
  integer :: i, j

  do j=1,ldb
      do i=1,lda
          a(i,j) = b(j,i)
      enddo
  enddo
  return
end subroutine transpose

!-----------------------------------------------------------------------
!> \brief Compute the convective term CONV for a passive scalar field FI
! using the skew-symmetric formulation.
! The field variable FI is defined on mesh M1 (GLL) and
! the velocity field is assumed given.
subroutine convop(conv,fi)
  use kinds, only : DP
  use ctimer, only : icalld, tadvc, nadvc, etime1, dnekclock
  use size_m, only : lx1, ly1, lz1
  use size_m, only : nx1, ny1, nz1, nelv
  use dealias, only : vxd, vyd, vzd
  use input, only : param, ifpert
  use geom, only : bm1
  use soln, only : vx, vy, vz
  use tstep, only : nelfld, ifield
  implicit none

!     Arrays in parameter list
  REAL(DP) :: CONV (LX1,LY1,LZ1,*)
  REAL(DP) :: FI   (LX1,LY1,LZ1,*)

  integer :: nxyz1, ntot1, ntotz

#ifndef NOTIMER
  nadvc=nadvc + 1
  etime1=dnekclock()
#endif

  NXYZ1 = NX1*NY1*NZ1
  NTOT1 = NX1*NY1*NZ1*NELV
  NTOTZ = NX1*NY1*NZ1*nelfld(ifield)

  conv(:,:,:,1:nelfld(ifield)) = 0._dp

  if (param(86) /= 0.0) then  ! skew-symmetric form
!max        call convopo(conv,fi)
      goto 100
  endif

!     write(6,*) istep,param(99),' CONVOP',ifpert
!     ip99 = param(99)
!     if (istep.gt.5) call exitti(' CONVOP dbg: $',ip99)

  if (param(99) == 2 .OR. param(99) == 3) then
!max        call conv1d(conv,fi)  !    use dealiased form
  elseif (param(99) == 4) then
      if (ifpert) then
          etime1 = etime1 - dnekclock()
          call convect_new (conv,fi, .FALSE. ,vx,vy,vz, .FALSE. )
          etime1 = etime1 + dnekclock()
      else
          etime1 = etime1 - dnekclock()
          call convect_new (conv,fi, .FALSE. ,vxd,vyd,vzd, .TRUE. )
          etime1 = etime1 + dnekclock()
      endif
      conv(:,:,:,1:nelv) = conv(:,:,:,1:nelv) / bm1 ! local mass inverse
  elseif (param(99) == 5) then
!max        call convect_cons(conv,fi, .FALSE. ,vx,vy,vz, .FALSE. )
      conv(:,:,:,1:nelv) = conv(:,:,:,1:nelv) / bm1  ! local mass inverse
  else
!max        call conv1 (conv,fi)  !    use the convective form
  endif

  100 continue

#ifndef NOTIMER
  tadvc=tadvc+(dnekclock()-etime1)
#endif

  return
end subroutine convop

!-----------------------------------------------------------------------
!> \brief Compute OUTFLD = SUMi Di*INPi.
!! the divergence of the vector field (INPX,INPY,INPZ)
!---------------------------------------------------------------------
subroutine opdiv(outfld,inpx,inpy,inpz)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, lelv
  use size_m, only : nx2, ny2, nz2, nelv, ndim, nx1, ny1, nz1
  use geom, only : rxm2, rym2, rzm2, sxm2, sym2, szm2, txm2, tym2, tzm2
  use mesh, only : if_ortho
  use geom, only : bm1, jacmi
  use dxyz, only : dxm12, dytm12, dztm12
  use ixyz, only : ixm12, iytm12, iztm12

  implicit none

  real(DP) :: outfld (lx2,ly2,lz2,nelv)
  real(DP) :: inpx   (lx1,ly1,lz1,nelv)
  real(DP) :: inpy   (lx1,ly1,lz1,nelv)
  real(DP) :: inpz   (lx1,ly1,lz1,nelv)
  
  real(DP), allocatable :: work (:,:,:,:)

  real(DP) ::  ta1 (lx1*ly1*lz1) &
  ,             ta2 (lx1*ly1*lz1) &
  ,             ta3 (lx1*ly1*lz1)
  integer :: iflg, ntot2, i1, i2, e, iz, n1, n2, nxy2, nyz1

  allocate(work(lx2,ly2,lz2,lelv))

  iflg = 1

  ntot2 = nx2*ny2*nz2*nelv
  if (if_ortho) then
    nxy2  = nx2*ny2
    n1    = nx2*ny1
    n2    = nx2*ny2
    nyz1  = ny1*nz1
    do e=1,nelv
      call mxm (dxm12,nx2,inpx(:,:,:,e),nx1,ta1,nyz1)
      i1=1
      i2=1
      do iz=1,nz1
          call mxm (ta1(i1),nx2,iytm12,ny1,ta2(i2),ny2)
          i1=i1+n1
          i2=i2+n2
      enddo
      call mxm  (ta2,nxy2,iztm12,nz1,outfld(:,:,:,e),nz2)
      outfld(:,:,:,e) = outfld(:,:,:,e) * rxm2(:,:,:,e)
 
      call mxm  (ixm12,nx2,inpy(:,:,:,e),nx1,ta3,nyz1) ! reuse ta3 below
      i1=1
      i2=1
      do iz=1,nz1
          call mxm (ta3(i1),nx2,dytm12,ny1,ta2(i2),ny2)
          i1=i1+n1
          i2=i2+n2
      enddo
      call mxm     (ta2,nxy2,iztm12,nz1,ta1,nz2)
      outfld(:,:,:,e) = outfld(:,:,:,e) + reshape(ta1,(/lx2,ly2,lz2/)) * sym2(:,:,:,e)
 
      call mxm (ixm12,nx2,inpz(:,:,:,e),nx1,ta1,nyz1) 
      i1=1
      i2=1
      do iz=1,nz1
          call mxm (ta3(i1),nx2,iytm12,ny1,ta2(i2),ny2)
          i1=i1+n1
          i2=i2+n2
      enddo
      call mxm (ta2,nxy2,dztm12,nz1,ta3,nz2)
      outfld(:,:,:,e) = outfld(:,:,:,e) + reshape(ta3,(/lx2,ly2,lz2/)) * tzm2(:,:,:,e)
    enddo
    outfld = outfld * bm1 * jacmi

  else
    call multd (work,inpx,rxm2,sxm2,txm2,1,iflg)
    call copy  (outfld,work,ntot2)
    call multd (work,inpy,rym2,sym2,tym2,2,iflg)
    outfld = outfld + work
    if (ndim == 3) then
        call multd (work,inpz,rzm2,szm2,tzm2,3,iflg)
        outfld = outfld + work
    endif
  endif

  return
end subroutine opdiv

!-----------------------------------------------------------------------
!> \brief Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)
subroutine wgradm1(ux,uy,uz,u,nel) ! weak form of grad
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1
  use dxyz, only : dxm1, dxtm1
  use geom, only : rxm1, sxm1, txm1, rym1, sym1, tym1, rzm1, szm1, tzm1
  use wz_m, only : w3m1
  use mesh, only : if_ortho
  implicit none

  integer, parameter :: lxyz=lx1*ly1*lz1
  real(DP) :: ux(lx1,ly1,lz1,*),uy(lx1,ly1,lz1,*),uz(lx1,ly1,lz1,*),u(lxyz,*)

  real(DP) :: ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)

  integer :: e, n, nel

  N = nx1-1
  do e=1,nel
      call local_grad3(ur,us,ut,u(1,e),N,dxm1,dxtm1)
      if (if_ortho) then
        ux(:,:,:,e) = w3m1*(ur*rxm1(:,:,:,e))
        uy(:,:,:,e) = w3m1*(us*sym1(:,:,:,e))
        uz(:,:,:,e) = w3m1*(ut*tzm1(:,:,:,e))
      else
        ux(:,:,:,e) = w3m1*(ur*rxm1(:,:,:,e) + us*sxm1(:,:,:,e) + ut*txm1(:,:,:,e))
        uy(:,:,:,e) = w3m1*(ur*rym1(:,:,:,e) + us*sym1(:,:,:,e) + ut*tym1(:,:,:,e))
        uz(:,:,:,e) = w3m1*(ur*rzm1(:,:,:,e) + us*szm1(:,:,:,e) + ut*tzm1(:,:,:,e))
      endif
  enddo

  return
end subroutine wgradm1

!-----------------------------------------------------------------------
