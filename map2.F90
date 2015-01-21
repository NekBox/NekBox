!-----------------------------------------------------------------------
subroutine mapelpr()
  use size_m, only : nid, nelt
  use input, only : ifflow, ifmvbd, ifheat, iftmsh, npscal
  use parallel, only : np, nelgt, nelg, nelgv, node, lglel
  use tstep, only : ifield
  implicit none

  logical :: ifverbm
  integer :: mfield, nfldt, idum, n8, ie, inid, mtype, inelt, nn, nm
  integer, external :: iglmin, iglmax

  if(nid == 0) write(6,'(/,A)') ' mapping elements to processors'

  MFIELD=2
  IF (IFFLOW) MFIELD=1
  IF (IFMVBD) MFIELD=0

!   Set up TEMPORARY value for NFIELD - NFLDT
  NFLDT = 1
  IF (IFHEAT) NFLDT = 2 + NPSCAL

!   Distributed memory processor mapping
  IF (NP > NELGT) THEN
      IF(NID == 0) THEN
          WRITE(6,1000) NP,NELGT
          1000 FORMAT(2X,'ABORT: Too many processors (',I8 &
          ,') for to few elements (',I8,').' &
          ,/,2X,'ABORTING IN MAPELPR.')
      ENDIF
      call exitt
  ENDIF
  call set_proc_map()

  DO 1200 IFIELD=MFIELD,NFLDT
      IF (IFTMSH(IFIELD)) THEN
          NELG(IFIELD)      = NELGT
      ELSE
          NELG(IFIELD)      = NELGV
      ENDIF
  1200 END DO

!   Output the processor-element map:
  ifverbm= .TRUE. 
  if (np > 2050 .OR. nelgt > 40000) ifverbm= .FALSE. 

  if(ifverbm) then
      idum = 1
      if(nid == 0) then
          N8 = min(8,nelt)
          write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
          if (NELT > 8) write(6 ,1315) (lglel(ie),ie=9,NELT)
          DO inid=1,NP-1
              mtype = inid
              call csend(mtype,idum,4,inid,0)            ! handshake
              call crecv(mtype,inelt,4)               ! nelt of other cpus
              N8 = min(8,inelt)
          !             write(6 ,1310) inid+1,(lglel(ie,inid+1),ie=1,n8)
          !             IF (inelt.gt.8)
          !    &           write(6 ,1315) (lglel(ie,inid+1),ie=9,inelt)
          ENDDO
          1310 FORMAT(' RANK',I6,' IEG',8I8)
          1315 FORMAT('     ',6X,'    ',8I8)
      else
          mtype = nid
          call crecv(mtype,idum,4)                ! hand-shake
          call csend(mtype,nelt,4,0,0)            ! nelt
      endif
  endif

!   Check elemental distribution

!    IF (IPASS.EQ.2.AND.PARAM(156).eq.9) THEN
!       NXYZ=NX1*NY1*NZ1
!       DO 1400 IE=1,NELT
!          VTMP1=NODE
!          VTMP2=IE
!          CALL CFILL(VX(1,1,1,IE) ,VTMP1,NXYZ)
!          CALL CFILL(VY(1,1,1,IE) ,VTMP2,NXYZ)
!          CALL CFILL(T(1,1,1,IE,1),VTMP1,NXYZ)
! 1400    CONTINUE
!       call prepost(.true.,'   ')
!    ENDIF

  nn = iglmin(nelt,1)
  nm = iglmax(nelt,1)
  if(nid == 0) write(6,*) 'element load imbalance: ',nm-nn,nn,nm

  if(nid == 0) then
      write(6,*) 'done :: mapping elements to processors'
      write(6,*) ' '
  endif

  return
end subroutine mapelpr

!-----------------------------------------------------------------------
!> \brief Compute element to processor distribution according to (weighted)
!! physical distribution in an attempt to minimize exposed number of
!! element interfaces.
subroutine set_proc_map()
  implicit none

  call get_vert

  return
end subroutine set_proc_map

