!-----------------------------------------------------------------------
subroutine get_session_info(intracomm)
  use mpif, only : mpi_comm_world
  use size_m, only : nid
  use input, only : session, path, series, ifneknek
  use parallel, only : csize
  implicit none

  integer :: mpi_is_initialized
  integer :: ierr, intracomm

!   Find out the session name:
  call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
  if ( mpi_is_initialized == 0 ) call mpi_init (ierr)

  call mpi_comm_dup(mpi_comm_world,intracomm,ierr)
  call iniproc(intracomm)

  CALL BLANK(SESSION,132)
  CALL BLANK(PATH   ,132)
  CALL BLANK(SERIES,132)

  ierr = 0
  IF(NID == 0) THEN
      OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD',ERR=24)
      READ(8,10) SESSION
      READ(8,10) PATH
      READ(8,10, iostat=ierr) SERIES
      if (ierr < 0) call chcopy(SERIES, SESSION, 132)
      ierr = 0
      10 FORMAT(A132)
      CLOSE(UNIT=8)
      GOTO 23
      24 CONTINUE
      write(6,*) 'No file SESSION.NAME; using defaults of '
      write(6,*) 'PATH=. and SESSION=NEK.'
      PATH='.'
      SESSION='NEK'
  ENDIF
  23 continue

  call err_chk(ierr,' Cannot open SESSION.NAME!$')
     
  call bcast(SESSION,132*CSIZE)
  call bcast(PATH,132*CSIZE)
  call bcast(SERIES,132*CSIZE)

  IFNEKNEK  = .FALSE. 

  return
end subroutine get_session_info

!-----------------------------------------------------------------------
!> \brief  Dummy for singlmesh
subroutine userchk_set_xfer
  return
end subroutine userchk_set_xfer

!-----------------------------------------------------------------------
!> \brief  Dummy for singlmesh
subroutine setintercomm(nekcommtrue,nptrue)
  implicit none
  integer :: nekcommtrue, nptrue
  return
end subroutine setintercomm

!-----------------------------------------------------------------------
!> \brief  Dummy for singlmesh
subroutine unsetintercomm(nekcommtrue,nptrue)
  implicit none
  integer :: nekcommtrue, nptrue
  return
end subroutine unsetintercomm
!------------------------------------------------------------------------
!> \brief  Dummy for singlmesh
subroutine uglmin(a,n)
  use kinds, only : DP
  implicit none
  real(DP) :: a
  integer :: n
  return
end subroutine uglmin
!------------------------------------------------------------------------
!> \brief  Dummy for singlmesh
subroutine happy_check(iflag)
  implicit none
  integer :: iflag
  return
end subroutine happy_check
!------------------------------------------------------------------------


