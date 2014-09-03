!-----------------------------------------------------------------------
!> \brief Solve the convection-diffusion equation for passive scalar IPSCAL
subroutine cdscal (igeom)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, nfield, nid, mxprev
  use input, only : ifmodel, ifkeps, ifaxis, ifaziv, iftran, iftmsh, ifprint
  use mass, only : bintm1, binvm1
  use soln, only : t, bq, tmask, tmult
  use tstep, only : nelfld, ifield, nmxnl, imesh, tolht, nmxh
  implicit none

  integer, intent(in) :: igeom

  LOGICAL ::          IFCONV

  real(DP) :: TA(LX1,LY1,LZ1,LELT), TB(LX1,LY1,LZ1,LELT)
  real(DP) :: H1(LX1,LY1,LZ1,LELT), H2(LX1,LY1,LZ1,LELT)

  integer, parameter :: ktot = lx1*ly1*lz1*lelt
  integer, parameter :: laxt = mxprev

  integer, save :: napprox(2) = 0
  real(DP), allocatable, save :: approx(:,:)
  character(4) ::     name4

  integer :: nel, ntot, nfldt, if1, isd, iter, intype
!max    include 'ORTHOT'

  if (.not. allocated(approx)) allocate(approx(ktot,0:laxt))

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

      do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

          INTYPE = 0
          IF (IFTRAN) INTYPE = -1
          CALL SETHLM  (H1,H2,INTYPE)
          CALL BCNEUSC (TA,-1)
          CALL ADD2    (H2,TA,NTOT)
          CALL BCDIRSC (T(1,1,1,1,IFIELD-1))
          CALL AXHELM  (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMESH,isd)
          CALL SUB3    (TB,BQ(1,1,1,1,IFIELD-1),TA,NTOT)
          CALL BCNEUSC (TA,1)
          CALL ADD2    (TB,TA,NTOT)

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

          call add2    (t(1,1,1,1,ifield-1),ta,ntot)

          call cvgnlps (ifconv) ! Check convergence for nonlinear problem
          if (ifconv) goto 2000

      !        Radiation case, smooth convergence, avoid flip-flop (ER).
          CALL CMULT (TA,0.5,NTOT)
          CALL SUB2  (T(1,1,1,1,IFIELD-1),TA,NTOT)

      1000 END DO
      2000 CONTINUE
      CALL BCNEUSC (TA,1)
      CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TA,NTOT) ! no idea why... pf

  endif

  return
end subroutine cdscal

!-----------------------------------------------------------------------
!> \brief Fill up user defined forcing function and collocate will the
!! mass matrix on the Gauss-Lobatto mesh.
subroutine makeuq
  use size_m, only : nx1, ny1, nz1
  use mass, only : bm1
  use soln, only : bq
  use tstep, only : time, dt, ifield, nelfld
  implicit none

  integer :: ntot

  ntot = nx1*ny1*nz1*nelfld(ifield)

  time = time-dt        ! Set time to t^n-1 for user function

  call rzero   ( bq(1,1,1,1,ifield-1) ,    ntot)
  call setqvol ( bq(1,1,1,1,ifield-1)          )
  call col2    ( bq(1,1,1,1,ifield-1) ,bm1,ntot)

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
          call cfill (bql(1,iel),cqvol,nxyz1)
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
  do 10 k=1,nz1
      do 10 j=1,ny1
          do 10 i=1,nx1
              call nekasgn (i,j,k,iel)
              qvol = 0.0
              call userq   (i,j,k,ielg)
              bql(i,j,k,iel) = qvol
  10 END DO

  return
end subroutine nekuq

!-----------------------------------------------------------------------
!> \brief Eulerian scheme, add convection term to forcing function
!!  at current time step.
!---------------------------------------------------------------
subroutine convab()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1, lelt
  use mass, only : bm1
  use soln, only : t, vtrans, bq
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: TA (LX1,LY1,LZ1,LELT)
  integer :: nel, ntot1

  NEL = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL
  CALL CONVOP  (TA,T(1,1,1,1,IFIELD-1))
  CALL COL2    (TA,VTRANS(1,1,1,1,IFIELD),NTOT1)
  CALL SUBCOL3 (BQ(1,1,1,1,IFIELD-1),BM1,TA,NTOT1)

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

  real(DP) :: TA (LX1,LY1,LZ1,LELT)
  real(DP) :: ab0, ab1, ab2
  integer :: nel, ntot1

  AB0   = AB(1)
  AB1   = AB(2)
  AB2   = AB(3)
  NEL   = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL

  CALL ADD3S2 (TA,VGRADT1(1,1,1,1,IFIELD-1), &
  VGRADT2(1,1,1,1,IFIELD-1),AB1,AB2,NTOT1)
  CALL COPY   (   VGRADT2(1,1,1,1,IFIELD-1), &
  VGRADT1(1,1,1,1,IFIELD-1),NTOT1)
  CALL COPY   (   VGRADT1(1,1,1,1,IFIELD-1), &
  BQ(1,1,1,1,IFIELD-1),NTOT1)
  CALL ADD2S1 (BQ(1,1,1,1,IFIELD-1),TA,AB0,NTOT1)

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
  use mass, only : bm1, bm1lag
  use soln, only : vtrans, t, tlag, bq
  use tstep, only : ifield, nelfld, dt, bd, nbd
  implicit none

  real(DP) ::  TA (LX1,LY1,LZ1,LELT), TB(LX1,LY1,LZ1,LELT) &
  ,             H2 (LX1,LY1,LZ1,LELT)

  integer :: nel, ntot1, ilag
  real(DP) :: const

  NEL   = NELFLD(IFIELD)
  NTOT1 = NX1*NY1*NZ1*NEL
  CONST = 1./DT
  CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
  CALL CMULT (H2,CONST,NTOT1)

  CALL COL3  (TB,BM1,T(1,1,1,1,IFIELD-1),NTOT1)
  CALL CMULT (TB,BD(2),NTOT1)

  DO 100 ILAG=2,NBD
      IF (IFGEOM) THEN
          CALL COL3 (TA,BM1LAG(1,1,1,1,ILAG-1), &
          TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
      ELSE
          CALL COL3 (TA,BM1, &
          TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
      ENDIF
      CALL CMULT (TA,BD(ILAG+1),NTOT1)
      CALL ADD2  (TB,TA,NTOT1)
  100 END DO

  CALL COL2 (TB,H2,NTOT1)
  CALL ADD2 (BQ(1,1,1,1,IFIELD-1),TB,NTOT1)

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
