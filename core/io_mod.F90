!==============================================================================
!> \file io_mod.F90
!! \brief Input/Output module with simplified interface
!! \date November 2014
!! \author Max Hutchinson
!!
!! This module should hold input/output routines, giving them implicit
!! interfaces.
module io
  implicit none
 
  character(132) :: load_name = 'NONE' !>!< name of output to load

  public load_ic, load_name
  private

contains

!--------------------------------------------------------------------
!> \brief Load initial condition from a previous multi-file output
subroutine load_ic()
  use kinds, only : DP
  use parallel, only : nid
  use soln, only : vx, vy, vz, t
  use restart, only : pid0, fid0
  use size_m, only : nx1, ny1, nz1
  use tstep, only : time

  integer :: nelo !>!< number of i/o elements per io-node
  integer :: word_size_load !>!< number of bytes per word

  integer :: ierr, i
  integer, parameter :: pad_size = 1
  real(DP), allocatable :: padding(:,:,:,:)

  if (load_name == 'NONE') then
    call get_restart_name(load_name)
  endif

  if (nid == pid0) then
    call mbyte_open(load_name, fid0, ierr)
  endif

  ! read and seek past header 
  call mfo_read_header(nelo, word_size_load, time)

  ! seek past positions
  if (nid == pid0) then
    allocate(padding(nx1, ny1, nz1, pad_size))
    do i = 1, nelo, pad_size
      call byte_read(padding, word_size_load * size(padding) / 4, ierr)
      call byte_read(padding, word_size_load * size(padding) / 4, ierr)
      call byte_read(padding, word_size_load * size(padding) / 4, ierr)
    enddo
  endif

  ! read velocities
  call mfo_read_vector(vx, vy, vz, size(vx,4), size(vx,1), size(vx,2), size(vx,3), word_size_load) 

  ! seek past pressure
  if (nid == pid0) then
    do i = 1, nelo, pad_size
      call byte_read(padding, word_size_load * size(padding) / 4, ierr)
    enddo
  endif
  call mfo_read_scalar(t(:,:,:,:,1), size(t, 4), size(t,1), size(t,2), size(t,3), word_size_load)

  if (nid == pid0) then
    call byte_close(ierr)
  endif

end subroutine load_ic

subroutine get_restart_name(fname)
  use kinds, only : DP
  use input, only : ifreguo, session, param
  use restart, only : nfileo, ifdiro
  use string, only : ltrunc
  implicit none

  character(132) ::  fname 
  character(1) ::   fnam1(132)

  character(6), save ::  six = "??????"
  character(6) :: str

  character(1), save :: slash = '/', dot = '.'

  integer :: k, len, ndigit
  integer, external :: i_find_prefix, mod1
  real(DP) :: rfileo

  fname = ''

#ifdef MPIIO
  rfileo = 1
#else
  rfileo = nfileo
#endif
  ndigit = int(log10(rfileo) + 1)

  k = 1
  if (ifdiro) then                                  !  Add directory
      call chcopy(fnam1(1),'A',1)
      call chcopy(fnam1(2),six,ndigit)  ! put ???? in string
      k = 2 + ndigit
      call chcopy(fnam1(k),slash,1)
      k = k+1
  endif

  len=ltrunc(session,132)                           !  Add SESSION
  call chcopy(fnam1(k),session,len)
  k = k+len

  if (ifreguo) then
      len=4
      call chcopy(fnam1(k),'_reg',len)
      k = k+len
  endif

  call chcopy(fnam1(k),six,ndigit)                  !  Add file-id holder
  k = k + ndigit

  call chcopy(fnam1(k  ),dot,1)                     !  Add .f appendix
  call chcopy(fnam1(k+1),'f',1)
  k = k + 2

  write(str,4) int(param(69))                                 !  Add nfld number
  4 format(i5.5)
  call chcopy(fnam1(k),str,5)
  k = k + 5

  call chcopy(fname(1:132),fnam1(1),k-1)

end subroutine get_restart_name

!-----------------------------------------------------------------------
!> \brief Read header and return number of elements and word size
subroutine mfo_read_header(nelo, word_size_file, time)
  use kinds, only : r4, DP
  use size_m, only : nid, nelt, lelt
  use parallel, only : lglel, isize
  use restart, only : nfileo, pid0, pid1
  use restart, only : iHeaderSize
  implicit none

  integer, intent(out) :: nelo
  integer, intent(out) :: word_size_file ! (intent out)
  real(DP), intent(out) :: time

  integer :: nelo_file !>!< number of i/o elements per io-node (from file)
  real(r4) :: test_pattern
  real(r4), allocatable :: padding(:)
  integer :: lglist(0:lelt)
  character(1) :: rdcode1(10)
  integer :: fid0, istep, nelgt
  character(132) :: hdr

  integer :: idum, nfileoo, j, mtype, inelp, ierr, i
  integer :: len
  integer :: pad_size
  integer :: nxo, nyo, nzo

  call nekgsync()
  idum = 1

  nfileoo = nfileo
  if(nid == pid0) then                ! how many elements to dump
      nelo = nelt
      do j = pid0+1,pid1
          mtype = j
          call csend(mtype,idum,4,j,0)   ! handshake
          call crecv(mtype,inelp,4)
          nelo = nelo + inelp
      enddo
  else
      mtype = nid
      call crecv(mtype,idum,4)          ! hand-shake
      call csend(mtype,nelt,4,pid0,0)   ! u4 :=: u8
  endif

  ierr = 0
  if(nid == pid0) then
      pad_size = (8 * (2**20) - (iHeaderSize + 4) ) / 4
      allocate(padding(pad_size)); padding = 0.
      call byte_read(hdr,iHeaderSize/4,ierr)
      read(hdr, 1) word_size_file,nxo,nyo,nzo,nelo_file,nelgt,time,istep,fid0,nfileoo &
      ,   (rdcode1(i),i=1,10)
      1 format(5x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13, &
      &        1x,i9,1x,i6,1x,i6,1x,10a)
      call byte_read(test_pattern,1,ierr)
      ! pad up to 8MB
      call byte_read(padding, pad_size, ierr)
      deallocate(padding)
  endif
  call bcast(word_size_file, isize)
  call bcast(time, 8)

  call err_chk(ierr,'Error writing header in mfo_write_hdr. $')

  if(nid == pid0) then
      call byte_read(lglist,nelt,ierr)
      pad_size = -nelt
      do j = pid0+1,pid1
          mtype = j
          call csend(mtype,idum,4,j,0)   ! handshake
          len = 4*(lelt+1)
          call crecv(mtype,lglist,len)
          if(ierr == 0) then
              call byte_read(lglist(1),lglist(0),ierr)
          pad_size = pad_size - lglist(0)
          endif
      enddo

    ! pad up to 8MB
    do while (pad_size < 0) 
      pad_size = pad_size + (8 * (2**20)) / 4
    enddo
    allocate(padding(pad_size)); padding = 0.
    call byte_read(padding, pad_size, ierr)
    deallocate(padding)

  else
      mtype = nid
      call crecv(mtype,idum,4)          ! hand-shake
              
      lglist(0) = nelt
      do j = 1, nelt
        lglist(j) = lglel(j)
      enddo

      len = 4*(nelt+1)
      call csend(mtype,lglist,len,pid0,0)
  endif

  call err_chk(ierr,'Error reading global nums in mfo_write_hdr$')
  return
end subroutine mfo_read_header

!-----------------------------------------------------------------------
!> \brief Read a scalar field
subroutine mfo_read_scalar(u,nel,mx,my,mz, wdsizo)
  use kinds, only : DP, r4
  use size_m, only : nid, lelt, lxo
  use restart, only : pid0, pid1
  implicit none

  integer, intent(in) :: nel, mx, my, mz, wdsizo
  real(DP), intent(in) :: u(mx,my,mz,1)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, ntot, idum, ierr, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  ntot = nxyz*nel

  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*2*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*1*lelt))
  endif

  if (nid == pid0) then
      idum = nel
      nout = wdsizo/4 * nxyz * idum
      if(wdsizo == 4 .and. ierr == 0) then
        call byte_read(u4,nout,ierr)          ! u4 :=: u8
      elseif(ierr == 0) then
        call byte_read(u8,nout,ierr)          ! u4 :=: u8
      endif

      if (wdsizo == 4) then             ! 32-bit output
          call copy4r (u,u4,nxyz * idum)
      else
          call copy   (u,u8,nxyz * idum)
      endif


  ! read in the data of my childs
      idum  = 1
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)       ! handshake
          call crecv(mtype,idum,4)       ! handshake

          nout  = wdsizo/4 * nxyz * idum
          if (wdsizo == 4 .AND. ierr == 0) then
            call byte_read(u4,nout,ierr)
            call csend(mtype, u4, nout*4, k, 0)
          elseif(ierr == 0) then
            call byte_read(u8,nout,ierr)
            call csend(mtype, u8, nout*4, k, 0)
          endif
      enddo

  else
      mtype = nid
      call crecv(mtype,idum,4)            ! hand-shake
      call csend(mtype, nel, 4, pid0, 0)

      if (wdsizo == 4) then             ! 32-bit output
        call crecv(mtype, u4, nxyz * nel *wdsizo)
        call copy4r (u, u4, nxyz * nel)
      else
        call crecv(mtype, u8, nxyz * nel *wdsizo)
        call copy   (u, u8, nxyz * nel)
      endif

  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_outs. $')

  return
end subroutine mfo_read_scalar

!-----------------------------------------------------------------------
!> \brief Read a vector field
subroutine mfo_read_vector(u,v,w,nel,mx,my,mz, wdsizo) 
  use kinds, only : DP, r4
  use size_m, only : nid, ndim, lxo, lelt
  use input, only : if3d
  use restart, only : pid0, pid1
  implicit none
 
  integer, intent(in) :: mx, my, mz, wdsizo
  real(DP), intent(in) :: u(mx*my*mz,*),v(mx*my*mz,*),w(mx*my*mz,*)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, nel, idum, ierr
  integer :: j, iel, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*6*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*3*lelt))
  endif
  
  if (nid == pid0) then
      nout = wdsizo/4 * ndim * nel * nxyz
      if (wdsizo == 4 .and. ierr == 0) then
        call byte_read(u4,nout,ierr)          ! u4 :=: u8
      elseif (ierr == 0) then
        call byte_read(u8,nout,ierr)          ! u4 :=: u8
      endif

      j = 0
      if (wdsizo == 4) then             ! 32-bit output
          do iel = 1,nel
              call copy4r   (u(1,iel), u4(j+1),nxyz)
              j = j + nxyz
              call copy4r   (v(1,iel), u4(j+1),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copy4r (w(1, iel), u4(j+1),nxyz)
                  j = j + nxyz
              endif
          enddo
      else
          do iel = 1,nel
              call copy     (u(1,iel), u8(j+1),nxyz)
              j = j + nxyz
              call copy     (v(1,iel), u8(j+1),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copy   (w(1,iel), u8(j+1),nxyz)
                  j = j + nxyz
              endif
          enddo
      endif

  ! read in the data of my childs
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)  ! handshake
          call crecv(mtype,idum,4)      ! hand-shake
              
          nout  = wdsizo/4 * ndim*nxyz * idum
          if (wdsizo == 4 .AND. ierr == 0) then
              call byte_read(u4,nout,ierr)
              call csend(mtype,u4,nout*4, k, 0)
          elseif(ierr == 0) then
              call byte_read(u8,nout,ierr)
              call csend(mtype,u8,nout*4, k, 0)
          endif
      enddo
  else
      mtype = nid
      call crecv(mtype,idum,4)            ! hand-shake
      call csend(mtype,nel,4,pid0,0)     ! u4 :=: u8

      if (wdsizo == 4) then             ! 32-bit output
        call crecv(mtype,u4,wdsizo*(nel*nxyz*ndim)) ! u4 :=: u8

        j = 0
        do iel = 1,nel
            call copy4r   (u(1,iel), u4(j+1),nxyz)
            j = j + nxyz
            call copy4r   (v(1,iel), u4(j+1),nxyz)
            j = j + nxyz
            if(if3d) then
                call copy4r (w(1,iel), u4(j+1),nxyz)
                j = j + nxyz
            endif
        enddo

      else
        call crecv(mtype,u8,wdsizo*(nel*nxyz*ndim))     ! u4 :=: u8
          j = 0
          do iel = 1,nel
              call copy     (u(1,iel), u8(j+1),nxyz)
              j = j + nxyz
              call copy     (v(1,iel), u8(j+1),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copy   (w(1,iel), u8(j+1),nxyz)
                  j = j + nxyz
              endif
          enddo
      endif
  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_outv. $')
  return
end subroutine mfo_read_vector

end module io

