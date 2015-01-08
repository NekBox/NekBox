!-----------------------------------------------------------------------
!> \brief Solve the convection-diffusion equation for passive scalar IPSCAL
subroutine cdscal (igeom)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, nfield, nid, mxprev
  use helmholtz, only : hsolve
  use input, only : ifmodel, ifkeps, ifaxis, ifaziv, iftran, iftmsh, ifprint
  use geom, only : bintm1, binvm1
  use soln, only : t, bq, tmask, tmult
  use tstep, only : nelfld, ifield, nmxnl, imesh, tolht, nmxh
  implicit none

  integer, intent(in) :: igeom

  LOGICAL ::          IFCONV

  real(DP), allocatable :: TA(:,:,:,:), TB(:,:,:,:)
  real(DP), allocatable :: H1(:,:,:,:), H2(:,:,:,:)

  integer, parameter :: laxt = mxprev

  integer, save :: napprox(2) = 0
  real(DP), allocatable, save :: approx(:,:,:)
  character(4) ::     name4

  integer :: nel, ntot, nfldt, if1, isd, iter, intype
!max    include 'ORTHOT'

!max  if (.not. allocated(approx)) allocate(approx(ktot,0:laxt))
  if (.not. allocated(approx)) allocate(approx(1,0:1,2))



  napprox(1) = laxt

  nel    = nelfld(ifield)
  ntot   = nx1*ny1*nz1*nel

  if (igeom == 1) then   ! geometry at t^{n-1}
      call makeq
      call lagscal
  else                   ! geometry at t^n

      IF (IFPRINT) THEN
          IF (IFMODEL .AND. IFKEPS) THEN
              NFLDT = NFIELD - 1
              IF (IFIELD == NFLDT .AND. NID == 0) THEN
                  WRITE (6,*) ' Turbulence Model - k/epsilon solution'
              ENDIF
          ELSE
              IF (IFIELD == 2 .AND. NID == 0) THEN
                  WRITE (6,*) ' Temperature/Passive scalar solution'
              ENDIF
          ENDIF
      ENDIF
      if1=ifield-1
      write(name4,1) if1-1
      1 format('PS',i2)
      if(ifield == 2) write(name4,'(A4)') 'TEMP'

  
  !        New geometry
  
      isd = 1
      if (ifaxis .AND. ifaziv .AND. ifield == 2) isd = 2
  !        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

      allocate(TA(LX1,LY1,LZ1,LELT), TB(LX1,LY1,LZ1,LELT))
      allocate(H1(LX1,LY1,LZ1,LELT), H2(LX1,LY1,LZ1,LELT))
      do iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

          INTYPE = 0
          IF (IFTRAN) INTYPE = -1
          CALL SETHLM  (H1,H2,INTYPE)
          CALL BCNEUSC (TA,-1)
          h2 = h2 + ta
          CALL BCDIRSC (T(1,1,1,1,IFIELD-1))
          CALL AXHELM  (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMESH,isd)
          tb = bq(:,:,:,:,ifield-1) - ta
          CALL BCNEUSC (TA,1)
          tb = tb + ta

      !        CALL HMHOLTZ (name4,TA,TB,H1,H2
      !    $                 ,TMASK(1,1,1,1,IFIELD-1)
      !    $                 ,TMULT(1,1,1,1,IFIELD-1)
      !    $                 ,IMESH,TOLHT(IFIELD),NMXH,isd)

          if(iftmsh(ifield)) then
              call hsolve  (name4,TA,TB,H1,H2 &
              ,tmask(1,1,1,1,ifield-1) &
              ,tmult(1,1,1,1,ifield-1) &
              ,imesh,tolht(ifield),nmxh,1 &
              ,approx,napprox,bintm1)
          else
              call hsolve  (name4,TA,TB,H1,H2 &
              ,tmask(1,1,1,1,ifield-1) &
              ,tmult(1,1,1,1,ifield-1) &
              ,imesh,tolht(ifield),nmxh,1 &
              ,approx,napprox,binvm1)
          endif

          t(:,:,:,:,ifield-1) = t(:,:,:,:,ifield-1) + ta

          call cvgnlps (ifconv) ! Check convergence for nonlinear problem
          if (ifconv) goto 2000

      !        Radiation case, smooth convergence, avoid flip-flop (ER).
          T(:,:,:,:,IFIELD-1) = T(:,:,:,:,IFIELD-1) - 0.5_dp*TA

      END DO
      2000 CONTINUE
      deallocate(h1,h2,tb)
      CALL BCNEUSC (TA,1)
      bq(:,:,:,:,ifield-1) = bq(:,:,:,:,ifield-1) + ta ! no idea why... pf
      deallocate(ta)

  endif

  return
end subroutine cdscal

!-----------------------------------------------------------------------
!> \brief Fill up user defined forcing function and collocate will the
!! mass matrix on the Gauss-Lobatto mesh.
subroutine makeuq
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use geom, only : bm1
  use soln, only : bq
  use tstep, only : time, dt, ifield, nelfld
  implicit none

  integer :: ntot

  ntot = nx1*ny1*nz1*nelfld(ifield)

  time = time-dt        ! Set time to t^n-1 for user function

  bq(:,:,:,:,ifield-1) = 0._dp
  call setqvol ( bq(1,1,1,1,ifield-1)          )
  bq(:,:,:,:,ifield-1) = bq(:,:,:,:,ifield-1) * bm1

  time = time+dt        ! Restore time

  return
end subroutine makeuq

!-----------------------------------------------------------------------
!> \brief  Set user specified volumetric forcing function (e.g. heat source).
subroutine setqvol(bql)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1
  use input, only : igroup, matype, cpgrp
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: bql(lx1*ly1*lz1,lelt)
  real(DP) :: cqvol
  integer :: nel, nxyz1, ntot1, iel, igrp

#ifndef MOAB
  nel   = nelfld(ifield)
  nxyz1 = nx1*ny1*nz1
  ntot1 = nxyz1*nel

  do iel=1,nel
      igrp = igroup(iel)
      if (matype(igrp,ifield) == 1) then ! constant source within a group
          cqvol = cpgrp(igrp,ifield,3)
          bql(:,iel) = cqvol
      else  !  pff 2/6/96 ............ default is to look at userq
          call nekuq (bql,iel)
      endif
  enddo

! 101 FORMAT(' Wrong material type (',I3,') for group',I3,', field',I2
!    $    ,/,' Aborting in SETQVOL.')
#else
! pulling in temperature right now, since we dont have anything else
  call userq2(bql)
#endif

  return
end subroutine setqvol

!------------------------------------------------------------------
!> \brief Generate user-specified volumetric source term (temp./p.s.)
!------------------------------------------------------------------
subroutine nekuq (bql,iel)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1
  use nekuse, only : qvol
  use parallel, only : lglel
  implicit none

  real(DP) :: bql(lx1,ly1,lz1,lelt)
  integer :: iel

  integer :: ielg, k, j, i

  ielg = lglel(iel)
  do k=1,nz1
      do j=1,ny1
          do i=1,nx1
              call nekasgn (i,j,k,iel)
              qvol = 0.0
              call userq   (i,j,k,ielg)
              bql(i,j,k,iel) = qvol
          enddo
      enddo
  END DO

  return
end subroutine nekuq

!-----------------------------------------------------------------------
!> \brief Eulerian scheme, add convection term to forcing function
!!  at current time step.
!---------------------------------------------------------------
subroutine convab()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1, lelt
  use geom, only : bm1
  use soln, only : t, vtrans, bq
  use tstep, only : ifield, nelfld
  implicit none

  real(DP), allocatable :: TA (:,:,:,:)
  integer :: nel, ntot1

  allocate(TA (LX1,LY1,LZ1,LELT))

  NEL = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL
  CALL CONVOP  (TA,T(1,1,1,1,IFIELD-1))
  bq(:,:,:,:,ifield-1) = bq(:,:,:,:,ifield-1) - bm1*ta *vtrans(:,:,:,:,ifield)

  return
end subroutine convab

!-----------------------------------------------------------------------
!> \brief Sum up contributions to 3rd order Adams-Bashforth scheme.
subroutine makeabq
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1, lelt
  use soln, only : vgradt1, vgradt2, bq
  use tstep, only : ab, ifield, nelfld
  implicit none

  real(DP), allocatable :: TA (:,:,:,:)
  real(DP) :: ab0, ab1, ab2
  integer :: nel, ntot1

  allocate(TA(LX1,LY1,LZ1,LELT))

  AB0   = AB(1)
  AB1   = AB(2)
  AB2   = AB(3)
  NEL   = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL

  ta = ab1*vgradt1(:,:,:,:,ifield-1) + ab2*vgradt2(:,:,:,:,ifield-1)
  vgradt2(:,:,:,:,ifield-1) = vgradt1(:,:,:,:,ifield-1)
  vgradt1(:,:,:,:,ifield-1) = bq(:,:,:,:,ifield-1)
  bq(:,:,:,:,ifield-1) = ab0*bq(:,:,:,:,ifield-1) + ta

  return
end subroutine makeabq

!-----------------------------------------------------------------------
!> \brief Add contributions to F from lagged BD terms.
!-----------------------------------------------------------------------
subroutine makebdq()
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1
  use geom, only : ifgeom
  use geom, only : bm1, bm1lag
  use soln, only : vtrans, t, tlag, bq
  use tstep, only : ifield, nelfld, dt, bd, nbd
  implicit none

  real(DP), allocatable :: TB(:,:,:,:)

  integer :: nel, ntot1, ilag

  allocate(TB(LX1,LY1,LZ1,LELT))

  NEL   = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL

  tb = bd(2) * bm1 * t(:,:,:,:,ifield-1)

  DO ILAG=2,NBD
      IF (IFGEOM) THEN
          tb = tb + bd(ilag+1) * bm1lag(:,:,:,:,ilag-1) * tlag(:,:,:,:,ilag-1, ifield-1)
      ELSE
          tb = tb + bd(ilag+1) * bm1 * tlag(:,:,:,:,ilag-1,ifield-1)
      ENDIF
  END DO

  tb = (1./DT) * tb * vtrans(:,:,:,:,ifield)
  bq(:,:,:,:,ifield-1) = bq(:,:,:,:,ifield-1) + tb

  return
end subroutine makebdq

!-----------------------------------------------------------------------
!> \brief Keep old passive scalar field(s)
!-----------------------------------------------------------------------
subroutine lagscal
  use size_m, only : nx1, ny1, nz1
  use soln, only : t, tlag
  use tstep, only : ifield, nelfld, nbdinp
  implicit none

  integer :: ntot1, ilag
  NTOT1 = NX1*NY1*NZ1*NELFLD(IFIELD)

  DO 100 ILAG=NBDINP-1,2,-1
      CALL COPY (TLAG(1,1,1,1,ILAG  ,IFIELD-1), &
      TLAG(1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
  100 END DO

  CALL COPY (TLAG(1,1,1,1,1,IFIELD-1),T(1,1,1,1,IFIELD-1),NTOT1)

  return
  end subroutine lagscal

!-----------------------------------------------------------------------
