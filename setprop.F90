!> \file setprop.F90 \copybrief setprop

!------------------------------------------------------------------------
!> \brief Set variable property arrays
subroutine setprop
  use kinds, only : DP
  use size_m, only :nx1, ny1, nz1, nfield
  use ctimer, only : icalld, tspro, nspro, etime1, dnekclock
  use input, only : ifflow, ifmhd
  use mass, only : bm1
  use soln, only : vdiff, vtrans
  use tstep, only : ifield, nelfld, volfld, avdiff, avtran
  implicit none

  integer :: nxyz1, mfield, nfldt, ifld, nel, ntot1
  real(DP) :: vol
  real(DP), external :: glsc2

#ifndef NOTIMER
  if (icalld == 0) tspro=0.0
  icalld=icalld+1
  nspro=icalld
  etime1=dnekclock()
#endif

  NXYZ1 = NX1*NY1*NZ1
  MFIELD=2
  IF (IFFLOW) MFIELD=1
  nfldt = nfield
  if (ifmhd) nfldt = nfield+1

  ifld = ifield

  DO IFIELD=MFIELD,nfldt
!max      IF (IFSTRS .AND. IFIELD == 1) CALL STNRINV ! expensive !

      CALL VPROPS

      nel = nelfld(ifield)
      vol = volfld(ifield)
      ntot1 = nxyz1*nel

      avdiff(ifield) = glsc2 (bm1,vdiff (1,1,1,1,ifield),ntot1)/vol
      avtran(ifield) = glsc2 (bm1,vtrans(1,1,1,1,ifield),ntot1)/vol

  ENDDO

  ifield = ifld

#ifndef NOTIMER
  tspro=tspro+(dnekclock()-etime1)
#endif


  RETURN
end subroutine setprop
