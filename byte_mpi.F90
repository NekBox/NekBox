subroutine byte_sync_mpi(mpi_fh)
  implicit none
  integer, intent(in) :: mpi_fh
#ifdef MPIIO
  include 'mpif.h'
  call MPI_file_sync(mpi_fh,ierr)
#endif

  return
end subroutine byte_sync_mpi

!--------------------------------------------------------------------------
subroutine byte_open_mpi(fname,mpi_fh,ierr)
  implicit none

  character(132) :: fname
  integer :: mpi_fh, ierr

#ifdef MPIIO
  include 'mpif.h'

  if(nid == pid0 .OR. nid == pid0r) then
  !        write(*,*) nid, 'call MPI_file_open',fname
      call MPI_file_open(nekcomm_io,fname, &
      MPI_MODE_RDWR+MPI_MODE_CREATE, &
      MPI_INFO_NULL,mpi_fh,ierr)
      if(ierr /= 0) then
          write(6,*) 'ABORT: Error in byte_open_mpi ', ierr
          return
      endif
  endif
#else
  write(6,*) 'byte_open_mpi: No MPI-IO support!'
  ierr=1
  return
#endif

  ierr=0
  return
end subroutine byte_open_mpi

!--------------------------------------------------------------------------
subroutine byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)
  use kinds, only : r4
  implicit none

  real(r4) :: buf(1)          ! buffer
  integer :: icount, iorank, mpi_fh, ierr

#ifdef MPIIO
  include 'mpif.h'

  if(nid == pid0 .OR. nid == pid0r) then
      iout = 4*icount ! icount is in 4-byte words
      if(iorank >= 0 .AND. nid /= iorank) iout = 0
  !        write(*,*) 'byte_read_mpi', nid, iout/4
#ifdef MPIIO_NOCOL
      call MPI_file_read(mpi_fh,buf,iout,MPI_BYTE, &
      MPI_STATUS_IGNORE,ierr)
#else
      call MPI_file_read_all(mpi_fh,buf,iout,MPI_BYTE, &
      MPI_STATUS_IGNORE,ierr)
#endif
      if(ierr /= 0) then
          write(6,*) 'ABORT: Error in byte_read_mpi ', ierr
          return
      endif
  endif
#else
  write(6,*) 'byte_read_mpi: No MPI-IO support!'
  ierr=1
  return
#endif
         
  ierr=0

  return
end subroutine byte_read_mpi

!--------------------------------------------------------------------------
subroutine byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
  use kinds, only : r4
  implicit none

  real(r4) :: buf(1)          ! buffer
  integer :: icount, iorank, mpi_fh, ierr

#ifdef MPIIO
  include 'mpif.h'

  if(nid == pid0 .OR. nid == pid0r) then
      iout = 4*icount ! icount is in 4-byte words
      if(iorank >= 0 .AND. nid /= iorank) iout = 0
  !        write(*,*) 'byte_write', nid, iout/4
#ifdef MPIIO_NOCOL
      call MPI_file_write(mpi_fh,buf,iout,MPI_BYTE, &
      MPI_STATUS_IGNORE,ierr)
#else
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_BYTE, &
      MPI_STATUS_IGNORE,ierr)
#endif
      if(ierr /= 0) then
          write(6,*) 'ABORT: Error in byte_write_mpi ', ierr
          return
      endif
  endif
#else
  write(6,*) 'byte_write_mpi: No MPI-IO support!'
  ierr=1
  return
#endif
  ierr=0
  return
end subroutine byte_write_mpi

!--------------------------------------------------------------------------
subroutine byte_close_mpi(mpi_fh,ierr)
  use size_m, only : nid
  implicit none

  integer :: mpi_fh, ierr

#ifdef MPIIO
  include 'mpif.h'
  if(nid == pid0 .OR. nid == pid0r) then
      call MPI_file_close(mpi_fh,ierr)
  endif
  if(ierr /= 0) then
      write(6,*) 'ABORT: Error in byte_close_mpi ', ierr
      return
  endif
#else
  if(nid == 0) write(6,*) 'byte_close_mpi: No MPI-IO support!'
  ierr=1
  return
#endif

  return
end subroutine byte_close_mpi

!--------------------------------------------------------------------------
subroutine byte_set_view(ioff_in,mpi_fh)
  use kinds, only : i8
  implicit none

  integer(i8) :: ioff_in
  integer :: mpi_fh

#ifdef MPIIO
  include 'mpif.h'
      
  if(nid == pid0 .OR. nid == pid0r) then
      if(ioff_in < 0) then
          write(6,*) 'byte_set_view: offset<0!'
          call exitt
      endif
  !         write(*,*) 'dataoffset', nid, ioff_in
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE, &
      'native',MPI_INFO_NULL,ierr)
      if(ierr /= 0) then
          write(6,*) 'ABORT: Error in byte_set_view ', ierr
          call exitt
      endif
  endif
#endif

  return
end subroutine byte_set_view

!--------------------------------------------------------------------------
subroutine nek_comm_io(nn)
  implicit none

  integer :: nn

#ifdef MPIIO
  include 'mpif.h'
!    c!ommon /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
!    c!ommon /scrns/  irank_io(0:lp-1)

#ifdef MPIIO_NOCOL
  if(nid == 0) then
      j = 0
      if(nid == pid0 .OR. nid == pid0r) then
          irank_io(j) = nid
          j = j + 1
      endif
      do ir = 1,np-1
          call csend(ir,idum,4,ir,0)           ! handshake
          call crecv(ir,ibuf,4)
          if(ibuf > 0) then
              irank_io(j) = ibuf
              j = j + 1
          endif
      enddo
  else
      mtype = nid
      ibuf = -1
      if(nid == pid0) then
          ibuf = nid
      endif
      call crecv(mtype,idum,4)                ! hand-shake
      call csend(mtype,ibuf,4,0,0)            ! u4 :=: u8
  endif

  call bcast(irank_io,isize*nn)

!      write(6,*) 'nid', nid, (irank_io(i),i=0,nn-1)

  call mpi_comm_group (nekcomm,nekgroup,ierr)
  if(ierr > 0) call exitt
  call mpi_group_incl (nekgroup,nn,irank_io,nekgroup_io,ierr)
  if(ierr > 0) call exitt
  call mpi_comm_create(nekcomm,nekgroup_io,nekcomm_io,ierr)
  if(ierr > 0) call exitt
  call mpi_group_free (nekgroup_io,ierr)
  if(ierr > 0) call exitt
  call mpi_group_free (nekgroup,ierr)
  if(ierr > 0) call exitt
#else
  nekcomm_io = nekcomm
  return
#endif

#endif

  return
end subroutine nek_comm_io
