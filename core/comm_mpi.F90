!---------------------------------------------------------------------
  subroutine iniproc(intracomm)
    use kinds, only : DP
    use mpif, only : mpi_comm_world, mpi_double_precision, mpi_real
    use mpif, only : mpi_tag_ub
    use size_m, only : nid, lp, lelg
    use parallel, only : np, wdsize, ifdblas, isize, lsize, csize, pid
    use parallel, only : nullpid, node0, node, cr_h
    use parallel, only : np_=>np,nekcomm,nekreal
    implicit none

    integer :: intracomm
    logical :: flag

    integer :: nval, ierr
    real(DP) :: eps, oneeps

!      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
!      if ( mpi_is_initialized .eq. 0 ) then
!         call mpi_init (ierr)
!      endif

! set nek communicator
    call init_nek_comm(intracomm)

    if(nid == 0) call printHeader

! check upper tag size limit
    call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)
    if (nval < (10000+max(lp,lelg))) then
        if(nid == 0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
        call exitt
    endif

    IF (NP > LP) THEN
        WRITE(6,*) &
        'ERROR: Code compiled for a max of',LP,' processors.'
        WRITE(6,*) &
        'Recompile with LP =',NP,' or run with fewer processors.'
        WRITE(6,*) &
        'Aborting in routine INIPROC.'
        call exitt
    endif

! set word size for REAL
    wdsize=4
    eps=1.0e-12
    oneeps = 1.0+eps
    if (oneeps /= 1.0) then
        wdsize=8
    else
        if(nid == 0) &
        write(6,*) 'ABORT: single precision mode not supported!'
        call exitt
    endif
    nekreal = mpi_real
    if (wdsize == 8) nekreal = mpi_double_precision

    ifdblas = .FALSE. 
    if (wdsize == 8) ifdblas = .TRUE. 

! set word size for INTEGER
! HARDCODED since there is no secure way to detect an int overflow
    isize = 4

! set word size for LOGICAL
    lsize = 4

! set word size for CHARACTER
    csize = 1

    PID = 0
    NULLPID=0
    NODE0=0
    NODE= NID+1

    if (nid == 0) then
        write(6,*) 'Number of processors:',np
        WRITE(6,*) 'REAL    wdsize      :',WDSIZE
        WRITE(6,*) 'INTEGER wdsize      :',ISIZE
    endif

    call crystal_setup(cr_h,nekcomm,np)  ! set cr handle to new instance

    return
  end subroutine iniproc

!-----------------------------------------------------------------------
  subroutine init_nek_comm(intracomm)
    use parallel, only : nid,np,nekcomm
    implicit none

    integer :: intracomm
    integer, external :: mynode, numnodes

    nekcomm = intracomm
    nid     = mynode()
    np      = numnodes()

    return
  end subroutine init_nek_comm
!-----------------------------------------------------------------------
!> \brief Global vector commutative operation
  subroutine gop( x, w, op, n)
    use kinds, only : DP
    use mpif, only : mpi_max, mpi_min, mpi_prod, mpi_sum
    use ctimer, only : ifsync, icalld, tgop, ngop, etime1, dnekclock
    use parallel, only :nid,nekcomm,nekreal 
    implicit none

    integer :: n
    real(DP) :: x(n), w(n)
    character(3) :: op

    integer :: ierr

    if (ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) then
        tgop =0.0d0
        ngop =0
        icalld=1
    endif
    ngop = ngop + 1
    etime1=dnekclock()
#endif

    if (op == '+  ') then
        call mpi_allreduce (x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
    elseif (op == 'M  ') then
        call mpi_allreduce (x,w,n,nekreal,mpi_max ,nekcomm,ierr)
    elseif (op == 'm  ') then
        call mpi_allreduce (x,w,n,nekreal,mpi_min ,nekcomm,ierr)
    elseif (op == '*  ') then
        call mpi_allreduce (x,w,n,nekreal,mpi_prod,nekcomm,ierr)
    else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
        call exitt
    endif

    call copy(x,w,n)

#ifndef NOTIMER
    tgop =tgop +(dnekclock()-etime1)
#endif

    return
  end subroutine gop
!-----------------------------------------------------------------------
!> \brief Global vector commutative operation
  subroutine igop( x, w, op, n)
    use mpif, only : mpi_integer, mpi_max, mpi_min, mpi_prod, mpi_sum
    use parallel, only : nid,nekcomm
    implicit none

    integer, intent(in) :: n
    integer, intent(inout) :: x(n), w(n)
    character(3), intent(in) :: op
    integer :: ierr

    if     (op == '+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
    elseif (op == 'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
    elseif (op == 'm  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
    elseif (op == '*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
    else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
    endif

    x = w

    return
    end subroutine igop
!-----------------------------------------------------------------------
!> \brief Global vector commutative operation
  subroutine i8gop( x, w, op, n)
    use kinds, only : i8
    use mpif, only : mpi_integer8, mpi_max, mpi_min, mpi_prod, mpi_sum
    use parallel, only : nid,nekcomm
    implicit none
    integer :: n
    integer(i8) :: x(n), w(n)
    character(3) :: op
    integer :: ierr

    if     (op == '+  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_sum ,nekcomm,ierr)
    elseif (op == 'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_max ,nekcomm,ierr)
    elseif (op == 'm  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_min ,nekcomm,ierr)
    elseif (op == '*  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_prod,nekcomm,ierr)
    else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
    endif

    x = w

    return
    end subroutine i8gop
!-----------------------------------------------------------------------
  subroutine csend(mtype,buf,len,jnid,jpid)
    use kinds, only : r4
    use mpif, only : mpi_byte
    use parallel, only : nekcomm
    implicit none
    real(r4) :: buf(1)
    integer :: mtype, len, jnid, jpid
    integer :: ierr

    call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)

    return
    end subroutine csend
!-----------------------------------------------------------------------
  subroutine crecv(mtype,buf,lenm)
    use kinds, only : r4
    use mpif, only : mpi_any_source, mpi_byte, mpi_status_size
    use parallel, only : nekcomm, nid
    implicit none
    integer :: mtype,  lenm
    real(r4) :: buf(1)

    integer :: status(mpi_status_size)
    integer :: len, jnid, ierr

    len = lenm
    jnid = mpi_any_source

    call mpi_recv (buf,len,mpi_byte &
    ,jnid,mtype,nekcomm,status,ierr)

    if (len > lenm) then
        write(6,*) nid,'long message in mpi_crecv:',len,lenm
        call exitt
    endif

    return
  end subroutine crecv
!-----------------------------------------------------------------------
  integer function numnodes()
    use parallel, only : nekcomm
    implicit none
    integer :: ierr

    call mpi_comm_size (nekcomm, numnodes , ierr)

    return
  end function numnodes
!-----------------------------------------------------------------------
  integer function mynode()
    use parallel, only : nekcomm
    implicit none
    integer :: myid, ierr

    call mpi_comm_rank (nekcomm, myid, ierr)
    mynode = myid

    return
  end function mynode
!-----------------------------------------------------------------------
!> \brief Broadcast logical variable to all processors.
  subroutine lbcast(ifif)
    use kinds, only : r4
    use parallel, only : np, isize
    implicit none

    logical :: ifif
    integer :: item

    if (np == 1) return

    item=0
    if (ifif) item=1
    call bcast(real(item, kind=r4),isize)
    ifif= .FALSE. 
    if (item == 1) ifif= .TRUE. 

    return
  end subroutine lbcast
!-----------------------------------------------------------------------
  subroutine bcast(buf,len)
    use kinds, only : r4
    use mpif, only : mpi_byte
    use parallel, only : nekcomm

    implicit none
    real(r4) :: buf
    integer :: len, ierr

    call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)

    return
  end subroutine bcast
!-----------------------------------------------------------------------
!     Note: len in bytes
  integer function irecv(msgtag,x,len)
    use mpif, only : mpi_any_source, mpi_byte
    use parallel, only : nekcomm
    implicit none
    integer :: msgtag, x(1), len
    integer :: ierr, imsg

    call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag &
    ,nekcomm,imsg,ierr)
    irecv = imsg

    return
  end function irecv
!-----------------------------------------------------------------------
  subroutine nekgsync()
    use parallel, only : nekcomm
    implicit none
    integer :: ierr

    call mpi_barrier(nekcomm,ierr)

    return
  end subroutine nekgsync
!-----------------------------------------------------------------------
subroutine exitti(stringi,idata)
  use size_m, only : nid
  use string, only : indx1
  implicit none

  integer :: idata

  character(132) :: stringi
  character(132) :: stringo
  character(11) :: s11
  integer :: len, k

  call blank(stringo,132)
  call chcopy(stringo,stringi,132)
  len = indx1(stringo,'$',1)
  write(s11,11) idata
  11 format(1x,i10)
  call chcopy(stringo(len:len),s11,11)

  if (nid == 0) write(6,1) (stringo(k:k),k=1,len+10)
  1 format('EXIT: ',132a1)

  call exitt

  return
end subroutine exitti

!-----------------------------------------------------------------------
subroutine err_chk(ierr,istring)
  use size_m, only : nid
  use string, only : indx1
  implicit none

  integer :: ierr

  character(*) :: istring
  character(132) :: ostring
  character(10) :: s10

  integer :: len, k
  integer, external :: iglsum

  ierr = iglsum(ierr,1)
  if(ierr == 0) return

  len = indx1(istring,'$',1)
  call blank(ostring,132)
  write(s10,11) ierr
  11 format(1x,' ierr=',i3)

  call chcopy(ostring,istring,len-1)
  call chcopy(ostring(len:len),s10,10)

  if (nid == 0) write(6,1) (ostring(k:k),k=1,len+10)
  1 format('ERROR: ',132a1)

  call exitt

  return
end subroutine err_chk

!-----------------------------------------------------------------------
subroutine exitt0
  implicit none

  integer :: ierr

  write(6,*) 'Emergency exit'

  call print_stack()
  call flush_io

  call mpi_finalize (ierr)
#ifdef EXTBAR
  call exit_(0)
#else
  call exit(0)
#endif
     

  return
end subroutine exitt0
!-----------------------------------------------------------------------
subroutine exitt
  use kinds, only : DP, i8
  use size_m, only : nid, nx1, ny1, nz1
  use ctimer, only : dnekclock, ttotal, etimes, ttime
  use input, only : ifneknek
  use parallel, only : nvtot, np
  use tstep, only : istep
  implicit none

#ifdef PAPI
  real(r4) :: papi_mflops
#endif
  integer(i8) :: papi_flops

  logical :: ifopen              !check for opened files

  real(DP) :: tstop, dtmp1, dtmp2, dtmp3, dgp
  integer :: nxyz, ierr


!   Communicate unhappiness to the other session
!  if (ifneknek .AND. icall == 0) call happy_check(0)
  if (ifneknek) call happy_check(0)

  call nekgsync()

  papi_flops = 0
#ifdef PAPI
  call nek_flops(papi_flops,papi_mflops)
#endif

  tstop  = dnekclock()
  ttotal = tstop-etimes
  nxyz   = nx1*ny1*nz1

  if (nid == 0) then
      inquire(unit=50,opened=ifopen)
      if(ifopen) close(50)
      dtmp1 = 0
      dtmp2 = 0
      dtmp3 = 0
      if(istep > 0) then
          dgp   = nvtot
          dgp   = max(dgp,1._dp)
          dtmp1 = np*ttime/(dgp*max(istep,1))
          dtmp2 = ttime/max(istep,1)
          dtmp3 = 1.*papi_flops/1e6
      endif
      write(6,*) ' '
      write(6,'(A)') 'call exitt: dying ...'
      write(6,*) ' '
      call print_stack()
      write(6,*) ' '
      write(6,'(4(A,1p1e13.5,A,/))') &
      'total elapsed time             : ',ttotal, ' sec' &
      ,'total solver time incl. I/O    : ',ttime , ' sec' &
      ,'time/timestep                  : ',dtmp2 , ' sec' &
      ,'CPU seconds/timestep/gridpt    : ',dtmp1 , ' sec'
#ifdef PAPI
      write(6,'(2(A,1g13.5,/))') &
      'Gflops                         : ',dtmp3/1000. &
      ,'Gflops/s                       : ',papi_mflops/1000.
#endif
  endif
  call flush_io

  call mpi_finalize (ierr)
#ifdef EXTBAR 
  all exit_(0)
#else
  call exit(0)
#endif
  return
end subroutine exitt
!-----------------------------------------------------------------------
subroutine printHeader
  implicit none
  INCLUDE 'HEADER'

 return
end subroutine printHeader
!-----------------------------------------------------------------------
integer function igl_running_sum(in)
  use mpif, only : mpi_integer, mpi_sum
  use parallel, only : nekcomm
  implicit none
  integer :: in
  integer :: x,w,r, ierr

  x = in  ! running sum
  w = in  ! working buff
  r = 0   ! recv buff

  call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
  igl_running_sum = r

  return
end function igl_running_sum
!-----------------------------------------------------------------------
subroutine msgwait(imsg)
  use mpif, only : mpi_status_size
  implicit none
  integer :: status(mpi_status_size)
  integer :: imsg, ierr

  call mpi_wait (imsg,status,ierr)

  return
end subroutine msgwait
