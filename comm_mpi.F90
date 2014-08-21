!---------------------------------------------------------------------
    subroutine iniproc(intracomm)
    use size_m
    use parallel
    include 'mpif.h'

    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

    logical :: flag

!      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
!      if ( mpi_is_initialized .eq. 0 ) then
!         call mpi_init (ierr)
!      endif

! set nek communicator
    call init_nek_comm(intracomm)
    nid  = nid_
    np   = np_

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
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    nekcomm = intracomm
    nid     = mynode()
    np      = numnodes()

    return
    end subroutine init_nek_comm
!-----------------------------------------------------------------------
    subroutine gop( x, w, op, n)

!     Global vector commutative operation

    use ctimer

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    real :: x(n), w(n)
    character(3) :: op

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
    subroutine igop( x, w, op, n)

!     Global vector commutative operation

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    integer :: x(n), w(n)
    character(3) :: op

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

    call icopy(x,w,n)

    return
    end subroutine igop
!-----------------------------------------------------------------------
    subroutine i8gop( x, w, op, n)

!     Global vector commutative operation

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    integer*8 :: x(n), w(n)
    character(3) :: op

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

    call i8copy(x,w,n)

    return
    end subroutine i8gop
!-----------------------------------------------------------------------
    subroutine csend(mtype,buf,len,jnid,jpid)
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    real*4 :: buf(1)

    call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)

    return
    end subroutine csend
!-----------------------------------------------------------------------
    subroutine crecv(mtype,buf,lenm)
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)

    real*4 :: buf(1)
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
    subroutine crecv3(mtype,buf,len,lenm)
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)

    real*4 :: buf(1)
    len = lenm
    jnid = mpi_any_source

    call mpi_recv (buf,len,mpi_byte &
    ,jnid,mtype,nekcomm,status,ierr)
    call mpi_get_count (status,mpi_byte,len,ierr)

    if (len > lenm) then
        write(6,*) nid,'long message in mpi_crecv:',len,lenm
        call exitt
    endif

    return
    end subroutine crecv3
!-----------------------------------------------------------------------
    integer function numnodes()
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    call mpi_comm_size (nekcomm, numnodes , ierr)

    return
    end function numnodes
!-----------------------------------------------------------------------
    integer function mynode()
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: myid

    call mpi_comm_rank (nekcomm, myid, ierr)
    mynode = myid

    return
    end function mynode
!-----------------------------------------------------------------------
    real*8 function dnekclock()
    include 'mpif.h'

    dnekclock=mpi_wtime()

    return
    END function
!-----------------------------------------------------------------------
    real*8 function dnekclock_sync()
    include 'mpif.h'

    call nekgsync()
    dnekclock_sync=mpi_wtime()

    return
    END function
!-----------------------------------------------------------------------
    subroutine lbcast(ifif)

!     Broadcast logical variable to all processors.

    use size_m
    use parallel
    include 'mpif.h'

    logical :: ifif

    if (np == 1) return

    item=0
    if (ifif) item=1
    call bcast(item,isize)
    ifif= .FALSE. 
    if (item == 1) ifif= .TRUE. 

    return
    end subroutine lbcast
!-----------------------------------------------------------------------
    subroutine bcast(buf,len)
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    real*4 :: buf(1)

    call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)

    return
    end subroutine bcast
!-----------------------------------------------------------------------
    subroutine create_comm(inewcomm)
    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

!      call mpi_comm_group (mpi_comm_world,itmp,ierr)
!      call mpi_comm_create (mpi_comm_world,itmp,icomm,ierr)
!      call mpi_group_free (itmp,ierr)

    call mpi_comm_dup(nekcomm,inewcomm,ierr)

    return
    end subroutine create_comm
!-----------------------------------------------------------------------
    function isend(msgtag,x,len,jnid,jpid)

!     Note: len in bytes

    integer :: x(1)

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    call mpi_isend (x,len,mpi_byte,jnid,msgtag &
    ,nekcomm,imsg,ierr)
    isend = imsg
!     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)

    return
    end function isend
!-----------------------------------------------------------------------
    function irecv(msgtag,x,len)

!     Note: len in bytes

    integer :: x(1)

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag &
    ,nekcomm,imsg,ierr)
    irecv = imsg
!     write(6,*) nid,' irecv:',imsg,msgtag,len


    return
    end function irecv
!-----------------------------------------------------------------------
    subroutine msgwait(imsg)

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)

!     write(6,*) nid,' msgwait:',imsg

    call mpi_wait (imsg,status,ierr)

    return
    end subroutine msgwait
!-----------------------------------------------------------------------
    subroutine nekgsync()

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

    call mpi_barrier(nekcomm,ierr)

    return
    end subroutine nekgsync
!-----------------------------------------------------------------------
    subroutine exittr(stringi,rdata,idata)
    use ctimer
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    character(1) :: stringi(132)
    character(1) :: stringo(132)
    character(25) :: s25

    call blank(stringo,132)
    call chcopy(stringo,stringi,132)
    len = indx1(stringo,'$',1)
    write(s25,25) rdata,idata
    25 format(1x,1p1e14.6,i10)
    call chcopy(stringo(len),s25,25)

    if (nid == 0) write(6,1) (stringo(k),k=1,len+24)
    1 format('EXIT: ',132a1)

    call exitt

    return
    end subroutine exittr
!-----------------------------------------------------------------------
    subroutine exitti(stringi,idata)
    use ctimer
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    character(1) :: stringi(132)
    character(1) :: stringo(132)
    character(11) :: s11

    call blank(stringo,132)
    call chcopy(stringo,stringi,132)
    len = indx1(stringo,'$',1)
    write(s11,11) idata
    11 format(1x,i10)
    call chcopy(stringo(len),s11,11)

    if (nid == 0) write(6,1) (stringo(k),k=1,len+10)
    1 format('EXIT: ',132a1)

    call exitt

    return
    end subroutine exitti
!-----------------------------------------------------------------------
    subroutine err_chk(ierr,string)
    use size_m
!     use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf
!     use ctimer
    character(1) :: string(132)
    character(1) :: ostring(132)
    character(10) :: s10

    ierr = iglsum(ierr,1)
    if(ierr == 0) return

    len = indx1(string,'$',1)
    call blank(ostring,132)
    write(s10,11) ierr
    11 format(1x,' ierr=',i3)

    call chcopy(ostring,string,len-1)
    call chcopy(ostring(len),s10,10)

    if (nid == 0) write(6,1) (ostring(k),k=1,len+10)
    1 format('ERROR: ',132a1)

    call exitt

    return
    end subroutine err_chk

!-----------------------------------------------------------------------
    subroutine exitt0
    use ctimer
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf
    include 'mpif.h'

    real*4 :: papi_mflops
    integer*8 :: papi_flops

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
    use ctimer
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf
    include 'mpif.h'
    common /happycallflag/ icall

    real*4 :: papi_mflops
    integer*8 :: papi_flops
    logical :: ifopen              !check for opened files


!     Communicate unhappiness to the other session
    if (ifneknek .AND. icall == 0) call happy_check(0)

    call nekgsync()


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
            dgp   = max(dgp,1.)
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

    INCLUDE 'HEADER'

    return
    end subroutine printHeader
!-----------------------------------------------------------------------
    function igl_running_sum(in)

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)
    integer :: x,w,r

    x = in  ! running sum
    w = in  ! working buff
    r = 0   ! recv buff

    call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
    igl_running_sum = r

    return
    end function igl_running_sum
!-----------------------------------------------------------------------
    subroutine platform_timer(ivb) ! mxm, ping-pong, and all_reduce timer

    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf


    call mxm_test_all(nid,ivb)  ! measure mxm times
!     call exitti('done mxm_test_all$',ivb)

    call comm_test(ivb)         ! measure message-passing and all-reduce times

    return
    end subroutine platform_timer
!-----------------------------------------------------------------------
    subroutine comm_test(ivb) ! measure message-passing and all-reduce times
! ivb = 0 --> minimal verbosity
! ivb = 1 --> fully verbose
! ivb = 2 --> smaller sample set(shorter)

    use size_m
    use parallel

    call gop_test(ivb)   ! added, Jan. 8, 2008

    log_np=log2(np)
    np2 = 2**log_np
    if (np2 == np) call gp2_test(ivb)   ! added, Jan. 8, 2008

    io = 6
    n512 = min(512,np-1)

    do nodeb=1,n512
        call pingpong(alphas,betas,0,nodeb,.0005,io,ivb)
        if (nid == 0) write(6,1) nodeb,np,alphas,betas
        1 format(2i10,1p2e15.7,' alpha beta')
    enddo

    return
    end subroutine comm_test
!-----------------------------------------------------------------------
    subroutine pingpong(alphas,betas,nodea,nodeb,dt,io,ivb)

    use size_m
    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

    parameter  (lt=lx1*ly1*lz1*lelt)
    parameter (mwd = 3*lt)
    common /scrns/ x(mwd),y(mwd)

    include 'mpif.h'
    integer :: status(mpi_status_size)

    character(10) :: fname

    if (nid == nodea) then
        write(fname,3) np,nodeb
        3 format('t',i4.4,'.',i4.4)
        if (io /= 6) open (unit=io,file=fname)
    endif

    call nekgsync
    call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

    nwds = 0
    if (nid == nodea .AND. ivb > 0) write(io,*)

    betas = 0  ! Reported inverse bandwidth
    count = 0

    do itest = 1,500

        nloop = msg_vol/(nwds+2)
        nloop = min(nloop,1000)
        nloop = max(nloop,1)

        len   = 8*nwds
             
        call ping_loop(t1,t0,len,nloop,nodea,nodeb,nid,x,y,x,y)

        if (nid == nodea) then
            tmsg = (t1-t0)/(2*nloop)   ! 2*nloop--> Double Buffer
            tmsg = tmsg / 2.           ! one-way cost = 1/2 round-trip
            tpwd = tmsg                ! time-per-word
            if (nwds > 0) tpwd = tmsg/nwds
            if (ivb > 0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd
            1 format(3i6,i12,1p2e16.8,' pg')

            if (nwds == 1) then
                alphas = tmsg
            elseif (nwds > 10000) then   ! "average" beta
                betas = (betas*count + tpwd)/(count+1)
                count = count + 1
            endif
        endif

        if (ivb == 2) then
            nwds = (nwds+1)*1.25
        else
            nwds = (nwds+1)*1.016
        endif
        if (nwds > mwd) then
        !        if (nwds.gt.1024) then
            if (nid == nodea .AND. io /= 6) close(unit=io)
            call nekgsync
            return
        endif

    enddo

    if (nid == nodea .AND. io /= 6) close(unit=io)
    call nekgsync

    return
    end subroutine pingpong
!-----------------------------------------------------------------------
    subroutine pingpongo(alphas,betas,nodea,nodeb,dt,io,ivb)

    use size_m
    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

    parameter  (lt=lx1*ly1*lz1*lelt)
    parameter (mwd = 3*lt)
    common /scrns/ x(mwd),y(mwd)

    include 'mpif.h'
    integer :: status(mpi_status_size)

    character(10) :: fname

    if (nid == nodea) then
        write(fname,3) np,nodeb
        3 format('t',i4.4,'.',i4.4)
        if (io /= 6) open (unit=io,file=fname)
    endif

    call nekgsync
    call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

    nwds = 0
    if (nid == nodea .AND. ivb > 0) write(io,*)

    betas = 0  ! Reported inverse bandwidth
    count = 0

    do itest = 1,500
        call nekgsync
        nloop = msg_vol/(nwds+2)
        nloop = min(nloop,1000)
        nloop = max(nloop,1)

        len   = 8*nwds
        jnid = mpi_any_source

        if (nid == nodea) then

            msg  = irecv(itest,y,1)
            call csend(itest,x,1,nodeb,0)   ! Initiate send, to synch.
            call msgwait(msg)

            t0 = mpi_wtime ()
            do i=1,nloop
                call mpi_irecv(y,len,mpi_byte,mpi_any_source,i &
                ,nekcomm,msg,ierr)
                call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)
                call mpi_wait (msg,status,ierr)
            enddo
            t1 = mpi_wtime ()
            tmsg = (t1-t0)/nloop
            tmsg = tmsg / 2.       ! Round-trip message time = twice one-way
            tpwd = tmsg
            if (nwds > 0) tpwd = tmsg/nwds
            if (ivb > 0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd
            1 format(3i6,i12,1p2e16.8,' pg')

            if (nwds == 1) then
                alphas = tmsg
            elseif (nwds > 10000) then
                betas = (betas*count + tpwd)/(count+1)
                count = count + 1
            endif

        elseif (nid == nodeb) then

            call crecv(itest,y,1)           ! Initiate send, to synch.
            call csend(itest,x,1,nodea,0)

            t0 = dnekclock()
            do i=1,nloop
                call mpi_recv (y,len,mpi_byte &
                ,jnid,i,nekcomm,status,ierr)
                call mpi_send (x,len,mpi_byte,nodea,i,nekcomm,ierr)
            enddo
            t1 = dnekclock()
            tmsg = (t1-t0)/nloop

        endif

        nwds = (nwds+1)*1.016
        if (nwds > mwd) then
            if (nid == nodea .AND. io /= 6) close(unit=io)
            call nekgsync
            return
        endif

    enddo

    if (nid == nodea .AND. io /= 6) close(unit=io)
    call nekgsync

    return
    end subroutine pingpongo
!-----------------------------------------------------------------------
    subroutine get_msg_vol(msg_vol,dt,nodea,nodeb)
    use size_m
    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
    parameter (lt=lx1*ly1*lz1*lelt)
    common /scrns/ x(3*lt),y(3*lt)

!     Est. msg vol for dt s

    msg_vol = 1000

    nwds  = min(1000,lt)
    nloop = 50
     
    tmsg = 0.
    call gop(tmsg,t1,'+  ',1)

    len = 8*nwds
    if (nid == nodea) then

        msg  = irecv(1,y,1)
        call csend(1,x,1,nodeb,0)   ! Initiate send, to synch.
        call msgwait(msg)

        t0 = dnekclock()
        do i=1,nloop
            msg  = irecv(i,y,len)
            call csend(i,x,len,nodeb,0)
            call msgwait(msg)
        enddo
        t1   = dnekclock()
        tmsg = (t1-t0)/nloop
        tpwd = tmsg/nwds

    elseif (nid == nodeb) then

        call crecv(1,y,1)           ! Initiate send, to synch.
        call csend(1,x,1,nodea,0)

        t0 = dnekclock()
        do i=1,nloop
            call crecv(i,y,len)
            call csend(i,x,len,nodea,0)
        enddo
        t1   = dnekclock()
        tmsg = (t1-t0)/nloop
        tmsg = 0.

    endif

    call gop(tmsg,t1,'+  ',1)
    msg_vol = nwds*(dt/tmsg)
!     if (nid.eq.nodea) write(6,*) nid,msg_vol,nwds,dt,tmsg,' msgvol'

    return
    end subroutine get_msg_vol
!-----------------------------------------------------------------------
    subroutine gop_test(ivb)
    use size_m
    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
    include 'mpif.h'
    integer :: status(mpi_status_size)

    parameter  (lt=lx1*ly1*lz1*lelt)
    parameter (mwd = 3*lt)
    common /scrns/ x(mwd),y(mwd)
    common /scruz/ times(2,500)
    common /scrcg/ nwd(500)

    nwds  = 1
    mtest = 0
    do itest = 1,500
        nwds = (nwds+1)*1.016
        if (nwds > mwd) goto 100
        mtest = mtest+1
        nwd(mtest) = nwds
    enddo
    100 continue

    nwds = 1
    do itest = mtest,1,-1

        tiny = 1.e-27
        call cfill(x,tiny,mwd)
        nwds = nwd(itest)
        call nekgsync

        t0 = mpi_wtime ()
        call gop(x,y,'+  ',nwds)
        call gop(x,y,'+  ',nwds)
        call gop(x,y,'+  ',nwds)
        call gop(x,y,'+  ',nwds)
        call gop(x,y,'+  ',nwds)
        call gop(x,y,'+  ',nwds)
        t1 = mpi_wtime ()

        tmsg = (t1-t0)/6 ! six calls
        tpwd = tmsg
        if (nwds > 0) tpwd = tmsg/nwds
        times(1,itest) = tmsg
        times(2,itest) = tpwd

    enddo
    101 continue


    if (nid == 0) then
        nwds = 1
        do itest=1,500
            if (ivb > 0 .OR. itest == 1) &
            write(6,1) np,nwds,(times(k,itest),k=1,2)
            1 format(i12,i12,1p2e16.8,' gop')
            nwds = (nwds+1)*1.016
            if (nwds > mwd) goto 102
        enddo
        102 continue
    endif

    return
    end subroutine gop_test
!-----------------------------------------------------------------------
    subroutine gp2_test(ivb)

    use size_m
    include 'mpif.h'

    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)

    parameter  (lt=lx1*ly1*lz1*lelt)
    parameter (mwd = 3*lt)
    common /scrns/ x(mwd),y(mwd)
    common /scruz/ times(2,500)

    call rzero(x,mwd)

    nwds = 1
    do itest = 1,500
        call gp2(x,y,'+  ',1,nid,np)

        t0 = mpi_wtime ()
        call gp2(x,y,'+  ',nwds,nid,np)
        call gp2(x,y,'+  ',nwds,nid,np)
        call gp2(x,y,'+  ',nwds,nid,np)
        call gp2(x,y,'+  ',nwds,nid,np)
        t1 = mpi_wtime ()

        tmsg = (t1-t0)/4 ! four calls
        tpwd = tmsg
        if (nwds > 0) tpwd = tmsg/nwds
        times(1,itest) = tmsg
        times(2,itest) = tpwd

        nwds = (nwds+1)*1.016
        if (nwds > mwd) goto 101
    enddo
    101 continue


    if (nid == 0) then
        nwds = 1
        do itest=1,500
            if (ivb > 0 .OR. itest == 1) &
            write(6,1) np,nwds,(times(k,itest),k=1,2)
            1 format(i12,i12,1p2e16.8,' gp2')
            nwds = (nwds+1)*1.016
            if (nwds > mwd) goto 102
        enddo
        102 continue
    endif

    return
    end subroutine gp2_test
!-----------------------------------------------------------------------
    integer function xor(m,n)

!  If NOT running on a parallel processor, it is sufficient to
!  have this routine return a value of XOR=1.

!  Pick one of the following:

!  UNIX 4.2, f77:
    XOR = OR(M,N)-AND(M,N)

!  Intel FTN286:
!     XOR = M.NEQV.N

!  Ryan-McFarland Fortran
!      XOR = IEOR(M,N)

!     XOR = 0
!     IF(M.EQ.1 .OR.  N.EQ.1) XOR=1
!     IF(M.EQ.0 .AND. N.EQ.0) XOR=0
!     IF(M.EQ.1 .AND. N.EQ.1) XOR=0
!     IF(M.GT.1 .OR.N.GT.1 .OR.M.LT.0.OR.N.LT.0) THEN
!        PRINT*,'ERROR IN XOR'
!        STOP
!     ENDIF

    return
    end function xor
!-----------------------------------------------------------------------
    subroutine gp2( x, w, op, n, nid, np)

!     Global vector commutative operation using spanning tree.

!     Std. fan-in/fan-out

    real :: x(n), w(n)
    character(3) :: op

    integer :: bit, bytes, cnt, diff, spsize, i, &
    parent, troot, xor, root, lnp, log2
    logical :: ifgot

    integer :: type
    save    type
    data    type  /998/

    type  = type+100
    if (type > 9992) type=type-998
    typer = type-1
    bytes = 8*n

    root    = 0
    troot   = max0((nid/np)*np, root)
    diff    = xor(nid,troot)
    nullpid = 0

!     Accumulate contributions from children, if any
    level2=1
    5 continue
    level=level2
    level2=level+level
    if (mod(nid,level2) /= 0) goto 20
    call crecv(type,w,bytes)
    if (op == '+  ') then
        do i=1,n
            x(i) = x(i) + w(i)
        enddo
    elseif (op == '*  ') then
        do i=1,n
            x(i) = x(i) * w(i)
        enddo
    elseif (op == 'M  ') then
        do i=1,n
            x(i) = max(x(i),w(i))
        enddo
    elseif (op == 'm  ') then
        do i=1,n
            x(i) = min(x(i),w(i))
        enddo
    endif
    if (level2 < np) goto 5

!     Pass result back to parent
    20 parent = nid-level
    if (nid /= 0) call csend(type,x,bytes,parent,nullpid)

!     Await final answer from node 0 via log_2 fan out
    level=np/2
    ifgot= .FALSE. 
    if (nid == root) ifgot= .TRUE. 

    lnp = log2(np)
    do i=1,lnp
        if (ifgot) then
            jnid=nid+level
            call csend(typer,x,bytes,jnid,nullpid)
        elseif (mod(nid,level) == 0) then
            call crecv(typer,x,bytes)
            ifgot= .TRUE. 
        endif
        level=level/2
    enddo

    return
    end subroutine gp2
!-----------------------------------------------------------------------
    subroutine ping_loop1(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

    real :: x(1),y(1)

    include 'mpif.h'
    integer :: status(mpi_status_size)

    i=0
    if (nid == nodea) then
        call nekgsync
        call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
        call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
    !        call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
        call msgwait(msg)                                          ! 1b

        t0 = mpi_wtime ()
        do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
        !           call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2b
        enddo
        t1 = mpi_wtime ()

    elseif (nid == nodeb) then

        call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
        call nekgsync
        call mpi_wait (msg,status,ierr)                            ! 1a

        j=i
        do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
        !           call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
        enddo
    !        call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb
        call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

    else
        call nekgsync
    endif

    return
    end subroutine ping_loop1
!-----------------------------------------------------------------------
    subroutine ping_loop2(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

    real :: x(1),y(1)

    include 'mpif.h'
    integer :: status(mpi_status_size)

    i=0
    if (nid == nodea) then
        call nekgsync
        call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
        call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
        call msgwait(msg)                                          ! 1b

        t0 = mpi_wtime ()
        do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_wait (msg,status,ierr)                         ! 2b
        enddo
        t1 = mpi_wtime ()

    elseif (nid == nodeb) then

        call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
        call nekgsync
        call mpi_wait (msg,status,ierr)                            ! 1a

        j=i
        do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
        enddo
        call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

    else
        call nekgsync
    endif

    return
    end subroutine ping_loop2
!-----------------------------------------------------------------------
    subroutine ping_loop(t1,t0,len,nloop,nodea,nodeb,nid,x1,y1,x2,y2)
!     Double Buffer : does 2*nloop timings

    common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

    real :: x1(1),y1(1),x2(1),y2(1)

    include 'mpif.h'
    integer :: status(mpi_status_size)

    itag=1
    if (nid == nodea) then
        call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)   ! 1b
        call nekgsync


        t0 = mpi_wtime ()
        do i=1,nloop
            call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 1a
            call mpi_irecv(y2,len,mpi_byte,nodeb,itag,nekcomm,msg2,ierr)! 2b
            call mpi_wait (msg1,status,ierr)                            ! 1b
            call mpi_send (x2,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 2a
            call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)! 3b
            call mpi_wait (msg2,status,ierr)                            ! 2b
        enddo
        t1 = mpi_wtime ()
        call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)        ! nb
        call mpi_wait (msg1,status,ierr)                              ! nb

    elseif (nid == nodeb) then

        call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)   ! nb
        call nekgsync


        do i=1,nloop
            call mpi_wait (msg1,status,ierr)                            ! 1a
            call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 1b
            call mpi_irecv(y2,len,mpi_byte,nodea,itag,nekcomm,msg2,ierr)! 2a
            call mpi_wait (msg2,status,ierr)                            ! 2a
            call mpi_send (x2,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 2b
            call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)! 3a
        enddo
        call mpi_wait (msg1,status,ierr)                            ! 2a
        call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)        ! nb

    else
        call nekgsync
    endif

    return
    end subroutine ping_loop
!-----------------------------------------------------------------------
    integer*8 function i8gl_running_sum(in)

    include 'mpif.h'
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)
    integer*8 :: x,r

    x = in  ! running sum
    r = 0   ! recv buff

    call mpi_scan(x,r,1,mpi_integer8,mpi_sum,nekcomm,ierr)
    i8gl_running_sum = r

    return
    END function
!-----------------------------------------------------------------------
