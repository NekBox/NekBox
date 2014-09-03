!-----------------------------------------------------------------------
subroutine get_session_info(intracomm)
  use kinds, only : DP
  use size_m, only : nid
  use input, only : session, path, ifneknek
  use parallel, only : csize
  implicit none
  include 'mpif.h'

  real(DP) :: mpi_is_initialized
  integer :: ierr, intracomm

!   Find out the session name:
  call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
  if ( mpi_is_initialized == 0 ) call mpi_init (ierr)

  call mpi_comm_dup(mpi_comm_world,intracomm,ierr)
  call iniproc(intracomm)

  CALL BLANK(SESSION,132)
  CALL BLANK(PATH   ,132)

  ierr = 0
  IF(NID == 0) THEN
      OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD',ERR=24)
      READ(8,10) SESSION
      READ(8,10) PATH
      10 FORMAT(A132)
      CLOSE(UNIT=8)
      GOTO 23
      24 CONTINUE
      write(6,*) 'No file SESSION.NAME; using defaults of '
      write(6,*) 'PATH=. and SESSION=NEK.'
      PATH='.'
      SESSION='NEK'
  23 ENDIF

  call err_chk(ierr,' Cannot open SESSION.NAME!$')
     
  call bcast(SESSION,132*CSIZE)
  call bcast(PATH,132*CSIZE)

  IFNEKNEK  = .FALSE. 

  return
end subroutine get_session_info

!-----------------------------------------------------------------------
subroutine userchk_set_xfer
!   Dummy for singlmesh
  return
end subroutine userchk_set_xfer

!-----------------------------------------------------------------------
subroutine setintercomm(nekcommtrue,nptrue)
!   Dummy for singlmesh
  return
end subroutine setintercomm

!-----------------------------------------------------------------------
subroutine unsetintercomm(nekcommtrue,nptrue)
!   Dummy for singlmesh
  return
end subroutine unsetintercomm
!------------------------------------------------------------------------
subroutine uglmin(a,n)
!   Dummy for singlmesh
  return
end subroutine uglmin
!------------------------------------------------------------------------
subroutine happy_check(iflag)
!   Dummy for singlmesh
  return
end subroutine happy_check
!------------------------------------------------------------------------


