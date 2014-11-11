module io
  implicit none


contains

subroutine load_ic(word_size_load)
  use kinds, only : r4, DP
  use input, only : session
  use parallel, only : nid
  use soln, only : vx, vy, vz, t
  use restart, only : pid0, fid0, wdsizo
  use size_m, only : nx1, ny1, nz1
  integer, intent(in) :: word_size_load
  character(132) :: session_pop
  character(132) :: hname
  integer :: nelo !>!< number of i/o elements per io-node
  integer :: wdsiztmp, ierr
  real(DP), allocatable :: padding(:,:,:,:)

#if 0
  session_pop = session
  session = prefix
  if (nid == pid0) then
    call mfo_open_files('   ', ierr)
  endi
  session = session_pop
#else
  if (nid == pid0) then
    hname = "io_test?.f00001"
    call mbyte_open(hname, fid0, ierr)
  endif
#endif

  wdsiztmp = wdsizo
  wdsizo = word_size_load

  ! read and seek past header 
  call mfo_read_hdr(nelo)
  ! seek past positions
  if (nid == pid0) then
    allocate(padding(nx1, ny1, nz1, nelo))
    call byte_read(padding, wdsizo * size(padding) / 4, ierr)
    call byte_read(padding, wdsizo * size(padding) / 4, ierr)
    call byte_read(padding, wdsizo * size(padding) / 4, ierr)
  endif
  ! read velocities
  call mfo_inv(vx, vy, vz, size(vx,4), size(vx,1), size(vx,2), size(vx,3))    
  ! seek past pressure
  if (nid == pid0) then
    call byte_read(padding, wdsizo * size(padding) / 4, ierr)
  endif
  call mfo_ins(t(:,:,:,:,1), size(t, 4), size(t,1), size(t,2), size(t,3))

  if (nid == pid0) then
    call byte_close(ierr)
  endif

  wdsizo = wdsiztmp

end subroutine load_ic

!-----------------------------------------------------------------------
!> \brief read a scalar field
subroutine mfo_ins(u,nel,mx,my,mz)
  use kinds, only : DP, r4
  use size_m, only : nid, lelt, lxo
  use restart, only : wdsizo, pid0, pid1
  implicit none

  integer, intent(in) :: nel, mx, my, mz
  real(DP), intent(in) :: u(mx,my,mz,1)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, len, leo, ntot, idum, ierr, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  len  = 8 + 8*(lelt*nxyz)  ! recv buffer size
  leo  = 8 + wdsizo*(nel*nxyz)
  ntot = nxyz*nel

  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*2*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*1*lelt))
  endif

  if (nid == pid0) then

      if(wdsizo == 4 .and. ierr == 0) then
        nout = wdsizo/4 * ntot
        call byte_read(u4,nout,ierr)          ! u4 :=: u8
      elseif(ierr == 0) then
        nout = wdsizo/4 * ntot
        call byte_read(u8,nout,ierr)          ! u4 :=: u8
      endif

      if (wdsizo == 4) then             ! 32-bit output
          call copy4r (u,u4,ntot)
      else
          call copy   (u,u8,ntot)
      endif


  ! read in the data of my childs
      idum  = 1
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)       ! handshake
          call crecv(mtype,idum,4)       ! handshake

          if (wdsizo == 4 .AND. ierr == 0) then
            nout  = wdsizo/4 * nxyz * idum
            call byte_read(u4,nout,ierr)
            call csend(mtype, u4, nout*4, k, 0)
          elseif(ierr == 0) then
            nout  = wdsizo/4 * nxyz * idum 
            call byte_read(u8,nout,ierr)
            call csend(mtype, u8, nout*4, k, 0)
          endif
      enddo

  else
      mtype = nid
      call crecv(mtype,idum,4)            ! hand-shake
      call csend(mtype, nel, 4, pid0, 0)

      if (wdsizo == 4) then             ! 32-bit output
        call crecv(mtype, u4, ntot*wdsizo)
        call copy4r (u, u4,ntot)
      else
        call crecv(mtype, u8, ntot*wdsizo)
        call copy   (u,u8,ntot)
      endif

  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_outs. $')

  return
end subroutine mfo_ins

!-----------------------------------------------------------------------
!> \brief output a vector field
subroutine mfo_inv(u,v,w,nel,mx,my,mz) 
  use kinds, only : DP, r4
  use size_m, only : nid, ndim, lxo, lelt
  use input, only : if3d
  use restart, only : wdsizo, pid0, pid1
  implicit none
 
  integer, intent(in) :: mx, my, mz
  real(DP), intent(in) :: u(mx*my*mz,*),v(mx*my*mz,*),w(mx*my*mz,*)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, len, leo, nel, idum, ierr
  integer :: j, iel, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  len  = 8 + 8*(lelt*nxyz*ndim)   ! recv buffer size (u4)
  leo  = 8 + wdsizo*(nel*nxyz*ndim)
  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*6*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*3*lelt))
  endif
  
  if (nid == pid0) then
      nout = wdsizo/4 * ndim*nel * nxyz
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
end subroutine mfo_inv

!-----------------------------------------------------------------------
!> \brief write hdr, byte key, els.
subroutine mfo_read_hdr(nelo)
  use kinds, only : r4
  use size_m, only : nid, nelt, lelt, ldimt
  use input, only : ifxyo, ifvo, ifpo, ifto, ifpsco
  use parallel, only : nelgt, lglel
  use restart, only : nfileo, pid0, pid1, rdcode1, wdsizo, nxo, nyo, nzo
  use restart, only : fid0, iHeaderSize
  use tstep, only : istep, time
  implicit none

  integer, intent(out) :: nelo

  real(r4) :: test_pattern
  real(r4), allocatable :: padding(:)
  integer :: lglist(0:lelt)

  character(132) :: hdr

  integer :: idum, nfileoo, j, mtype, inelp, ierr, i, npscalo, k
  integer :: ibsw_out, len
  integer :: pad_size

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
      call byte_read(test_pattern,1,ierr)
      ! pad up to 8MB
      call byte_read(padding, pad_size, ierr)
      deallocate(padding)
  endif

  call err_chk(ierr,'Error writing header in mfo_write_hdr. $')

  if(nid == pid0) then
      call byte_read(lglel,nelt,ierr)
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
      lglist(1:nelt) = lglel(1:nelt)

      len = 4*(nelt+1)
      call csend(mtype,lglist,len,pid0,0)
  endif

  call err_chk(ierr,'Error reading global nums in mfo_write_hdr$')
  return
end subroutine mfo_read_hdr


end module io

