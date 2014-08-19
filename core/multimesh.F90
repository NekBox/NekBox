!-----------------------------------------------------------------------
    subroutine get_session_info(intracomm)
    include 'mpif.h'
    include 'SIZE'
    include 'GLOBALCOM'
    include 'TSTEP'
    include 'INPUT'
          
    common /happycallflag/ icall
    common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
    integer :: nid_global_root(0:nsessmax-1)
    character(132) :: session_mult(0:nsessmax-1), path_mult(0:nsessmax-1)

    logical :: ifhigh

!     nsessmax = upper limit for number of sessions
!     npmax = upper limit for number of processors in each session


!     Read from a file: number of sessions (nsessions),
!     session name (SESSION_MULT(n-1)),
!     number of pocessors in each session (npsess(n-1))
!     and path (PATH_MULT(n-1))


    call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
    if ( mpi_is_initialized == 0 ) call mpi_init (ierr)

    call mpi_comm_size(mpi_comm_world,np_global,ierr)
    call mpi_comm_rank(mpi_comm_world,nid_global,ierr)

    nid=nid_global
    nekcomm=mpi_comm_world

    ierr = 0
    if (nid_global == 0) then
        open (unit=8,file='SESSION.NAME',status='old',err=24)
        read(8,*) nsessions
        do n=0,nsessions-1
            call blank(session_mult(n),132)
            call blank(path_mult(n)   ,132)
            read(8,10) session_mult(n)
            read(8,10) path_mult(n)
            read(8,*)  npsess(n)
        enddo
        10 format(a132)
        close(unit=8)
        goto 23
        24 ierr = 1
    endif
    23 continue
          
    call err_chk(ierr,' Cannot open SESSION.NAME!$')

    call bcast(nsessions,4)

    if (nsessions > 2) &
    call exitti('More than 2 sessions are currently &
    not supported!$',1)

    do n=0,nsessions-1
        call bcast(npsess(n),4)
        call bcast(session_mult(n),132)
        call bcast(path_mult(n),132)
    enddo

    npall=0
    do n=0,nsessions-1
        npall=npall+npsess(n)
    enddo
         
!     Check if number of processors in each session is consistent
!     with the total number of processors

    if (npall /= np_global) &
    call exitti('Wrong number of processors!$',1)

!     Assign key for splitting into multiple groups

    nid_global_root_next=0
    do n=0,nsessions-1
        nid_global_root(n)=nid_global_root_next
        nid_global_root_next=nid_global_root(n)+npsess(n)
        if (nid_global >= nid_global_root(n) .AND. &
        nid_global < nid_global_root_next) &
        idsess=n
    enddo
             
    call mpi_comm_split(mpi_comm_world,idsess,nid,intracomm,ierr)
     
    session = session_mult(idsess)
    path    = path_mult   (idsess)

!     Intercommunications set up only for 2 sessions

    if (nsessions > 1) then

        if (idsess == 0) idsess_neighbor=1
        if (idsess == 1) idsess_neighbor=0
         
        call mpi_intercomm_create(intracomm,0,mpi_comm_world, &
        nid_global_root(idsess_neighbor), 10,intercomm,ierr)

        np_neighbor=npsess(idsess_neighbor)
              
        call iniproc(intracomm)

        ifhigh= .TRUE. 
        call mpi_intercomm_merge(intercomm, ifhigh, iglobalcomm, ierr)
              
        ifneknek  = .TRUE. 

        ninter = 1 ! Initialize NEKNEK interface extrapolation order to 1.

        icall = 0  ! Emergency exit call flag

    endif

    return
    end subroutine get_session_info
!-----------------------------------------------------------------------
    subroutine multimesh_create

    include 'SIZE'
    include 'TOTAL'
    integer :: iList_all(ldim+2,nmaxcom)
    real ::    pts(ldim,nmaxcom)

    integer :: icalld
    save    icalld
    data    icalld  /0/

!     Set interpolation flag: points with bc = 'int' get intflag=1.
!     Boundary conditions are changed back to 'v' or 't'.

    if (icalld == 0) then
        call set_intflag
        icalld = icalld + 1
    endif

!     Every processor sends all the points with intflag=1 to the processors
!     of remote session
!     (send_points)

!    Root processor receives all the boundary points from remote session
!    and redistributes it among its processors for efficient
!    localization of points (recv_points)

    call exchange_points(pts,iList_all,npoints_all)
          
!     Find which boundary points of remote session are
!     within the local mesh (intpts_locate).
!     Communicate the point info(iList) to the remote processors
!     Points which are inside the mesh are masked and
!     their identities are stored in array

    call intpts_locate (pts,iList_all,npoints_all)

    return

    end subroutine multimesh_create
!-------------------------------------------------------------
    subroutine set_intflag
    include 'SIZE'
    include 'TOTAL'
    include 'NEKNEK'
    CHARACTER CB*3

!     Set interpolation flag: points with boundary condition = 'int'
!     get intflag=1.

!     Boundary conditions are changed back to 'v' or 't'.

    if (ifflow) then
        ifield = 1
    elseif (ifheat) then
        ifield   = 2
    endif

    nfaces = 2*ndim
    nel    = nelfld(ifield)
          
    nflag=nel*nfaces
    call izero(intflag,nflag)

    do 2010 iel=1,nel
        do 2010 iface=1,nfaces
            cb=cbc(iface,iel,ifield)
            if (cb == 'int') then
                intflag(iface,iel)=1
                if (ifield == 1) cbc(iface,iel,ifield)='v'
                if (ifield == 2) cbc(iface,iel,ifield)='t'
            endif
    2010 END DO

    return
    end subroutine set_intflag

!-------------------------------------------------------------
    subroutine exchange_points(pts,iList_all,npoints_all)
    include 'SIZE'
    include 'TOTAL'
    include 'NEKUSE'
    include 'NEKNEK'
    include 'mpif.h'
    common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
    integer :: status(mpi_status_size)
    real ::    pts(ldim,nmaxcom)
    integer :: jsend((ldim+2)*nmaxl),jrecv((ldim+2)*nmaxl)
    common /exchr/ rsend(ldim*nmaxl), rrecv(ldim*nmaxl)
    integer :: iList_all(ldim+2,nmaxcom)
        
!     Look for boundary points with Diriclet b.c. (candidates for
!     interpolation)

          
    if (ifflow) then
        ifield = 1
    elseif (ifheat) then
        ifield   = 2
    endif

    nfaces = 2*ndim
    nel    = nelfld(ifield)

    ip = 0
    do 2010 iel=1,nel
        do 2010 iface=1,nfaces
            if (intflag(iface,iel) == 1) then
                call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                do 100 iz=kz1,kz2
                    do 100 iy=ky1,ky2
                        do 100 ix=kx1,kx2
                            call nekasgn (ix,iy,iz,iel)
                            ip=ip+1

                            if (if3d) then
                                jsend((ldim+2)*(ip-1)+1)=ix
                                jsend((ldim+2)*(ip-1)+2)=iy
                                jsend((ldim+2)*(ip-1)+3)=iz
                                jsend((ldim+2)*(ip-1)+4)=iel
                                jsend((ldim+2)*(ip-1)+5)=nid

                                rsend(ldim*(ip-1)+1)=x
                                rsend(ldim*(ip-1)+2)=y
                                rsend(ldim*(ip-1)+3)=z
                            else
                                jsend((ldim+2)*(ip-1)+1)=ix
                                jsend((ldim+2)*(ip-1)+2)=iy
                                jsend((ldim+2)*(ip-1)+3)=iel
                                jsend((ldim+2)*(ip-1)+4)=nid

                                rsend(ldim*(ip-1)+1)=x
                                rsend(ldim*(ip-1)+2)=y
                            endif


                            if (ip > nmaxl) then
                                write(6,*) nid, &
                                ' ABORT: nbp (current ip) too large',ip,nmaxl
                                call exitt
                            endif

                100 END DO
            endif
    2010 END DO

    nbp = ip

    call izero(iList_all,(ldim+2)*nmaxcom)
    call rzero(pts,  ldim*nmaxcom)

!     Send amount of points boundary points, its coordinates and
!     its local identities
!     to the the processors of remote session (for load balancing)

!     Determine amount of points to send to each processor

    ndistrib = nbp/np_neighbor
    iremainder = nbp - ndistrib*np_neighbor

    ip=0
    iaddress=1
    iraddress=1

    do id=0,np_neighbor-1
        if(id <= iremainder-1) then
            nsend = ndistrib + 1
        else
            nsend = ndistrib
        endif

        len=isize
        call mpi_irecv (nrecv,len,mpi_byte, id, id, &
        intercomm, msg,ierr)
        call mpi_send  (nsend,len,mpi_byte, id, nid, &
        intercomm, ierr)

        call mpi_wait (msg,status,ierr)

        len=(ldim+2)*nrecv*isize
        call mpi_irecv (jrecv, len ,mpi_byte, id, 100+id, &
        intercomm, msg, ierr)
        len=(ldim+2)*nsend*isize
        call mpi_send (jsend(iaddress),len,mpi_byte, id, 100+nid, &
        intercomm, ierr)
        iaddress=iaddress+(ldim+2)*nsend

        call mpi_wait (msg, status, ierr)

        len=ldim*nrecv*wdsize
        call mpi_irecv (rrecv, len ,mpi_byte, id, 200+id, &
        intercomm, msg, ierr)
        len=ldim*nsend*wdsize
        call mpi_send (rsend(iraddress),len,mpi_byte, id, 200+nid, &
        intercomm, ierr)
        iraddress=iraddress+ldim*nsend

        call mpi_wait (msg, status, ierr)

        do n=1,nrecv
            ip = ip + 1

            if (ip > nmaxcom) then
                write(6,*) nid,'ABORT: increase nmaxcom',ip,nmaxcom
                call exitt
            endif

                        
            do j=1,ldim+2
                iList_all(j,ip)=jrecv((ldim+2)*(n-1)+j)
            enddo
            do j=1,ldim
                pts(j,ip)=rrecv(ldim*(n-1)+j)
            enddo

        enddo ! n

    enddo    ! id

    call neknekgsync()

    npoints_all = ip

    return
    end subroutine exchange_points
!-------------------------------------------------------------------------
    subroutine intpts_locate (pts,iList_all,npoints_all)

    include 'SIZE'
    include 'TOTAL'
    include 'NEKNEK'
    include 'mpif.h'
    integer :: jsend((ldim+1)*nmaxl),jrecv((ldim+1)*nmaxl)
    integer :: rcode_all(nmaxcom),elid_all(nmaxcom),proc_all(nmaxcom)
    real ::    dist(nmaxcom)
    real ::    pts(ldim,nmaxcom)
    real ::    dist_all(nmaxcom)
    real ::    rst_all(nmaxcom*ldim)
    integer :: iList_all(ldim+2,nmaxcom)
    character(3) :: CB
    integer :: status(mpi_status_size)
    logical :: ifcomm

    call intpts_setup(-1.0,inth_multi) ! use default tolerance
            

    call findpts(inth_multi,rcode_all,1, &
    proc_all,1, &
    elid_all,1, &
    rst_all,ndim, &
    dist_all,1, &
    pts(1,1),ndim, &
    pts(2,1),ndim, &
    pts(3,1),ndim,npoints_all)


!     Rearrange arrays so that only the points which are found within
!     the mesh are marked for interpolation: those are truly internal
!     points (rcode_all=0) plus some of the "boundary points"
!     (rcode_all=1) which fall on the boundaries of the elements
!     of another domain (essentially also internal points) or
!     wall or periodic boundaries

!     Points with rcode_all(i)=1 which fall on the "interpolation"
!     boundaries of another domain are treated with boundary
!     conditions and not with interpolation routines, and are
!     therefore excluded.

!     rcode_all, proc_all ... npoints_all represent all boundary points,
!     rcode, proc ... npoints - only interpolation points
!     (used by interpolation routines)
!     Contained in a common block /multipts_r/, /multipts_i/ (in NEKNEK)



!     Exclude "true" boundary points

    if (ifflow) then
        ifield = 1
    elseif (ifheat) then
        ifield   = 2
    endif

    NEL    = NELFLD(IFIELD)
    NXYZ  =NX1*NY1*NZ1
    NTOT  =NXYZ*NEL

    call izero(imask,ntot)

    ip=0
    icount=0

    do 100 i=1,npoints_all
                 
        if (rcode_all(i) < 2) then

            iel  = elid_all(i) + 1

            if (rcode_all(i) == 1 .AND. dist_all(i) > 1e-12) then
                if (ndim == 2) &
                write(6,'(A,3E15.7)') &
                ' WARNING: point on boundary or outside the mesh xy[z]d^2: ', &
                (pts(k,i),k=1,ndim),dist_all(i)
                if (ndim == 3) &
                write(6,'(A,4E15.7)') &
                ' WARNING: point on boundary or outside the mesh xy[z]d^2: ', &
                (pts(k,i),k=1,ndim),dist_all(i)
                GOTO 100
            endif
            ip=ip+1
            rcode(ip) = rcode_all(i)
            elid(ip)  = elid_all(i)
            proc(ip)  = proc_all(i)
            do j=1,ndim
                rst(ndim*(ip-1)+j)   = rst_all(ndim*(i-1)+j)
            enddo
            do n=1,ndim+1
                iList(n,ip) = iList_all(n,i)
            end do

        !    Check receiving remote processor id information
            n=ndim+2
        !     First point
            if (ip == 1) then
                icount=1
                infosend(icount,1)=iList_all(n,i)
                infosend(icount,2)=ip
                iprocp=iList_all(n,i)
            else
                iproc  = iList_all(n,i)
                if (iproc /= iprocp) then
                    if (iproc > iprocp) then
                        icount=icount+1
                        infosend(icount,1)=iproc
                        infosend(icount,2)=ip
                        iprocp=iproc
                    else
                        write(6,*) 'Attention: intrinsic sorting is not achieved'
                        write(6,'(a7,3i7,1x,a10)') 'switch', ip, &
                        iproc,iprocp, session
                        ierror = 1
                        goto 100
                    end if
                end if
            end if
                  
        end if  !  rcode_all

    100 END DO

    call iglmax(ierror,1)
    if (ierror == 1) call exitt

    npoints=ip
    npsend=icount

!     Store point counts (infosend(:,2)) instead of pointers

    do n=2,npsend
        ip=infosend(n-1,2)
        ic=infosend(n,2)
        infosend(n-1,2)=ic-ip
    end do
    infosend(npsend,2)=npoints-infosend(npsend,2)+1

    write(6,'(a7,i12,1x,a10)') 'found', npoints, session

!     Sending all processors information about the number
!     of receiving points, including zero points
!     ifcomm=.false. means there are no more processors
!     with non-zero receiving points (either if npsend=0
!     or if we already notified all the relevant processors)

    if (npsend == 0) then
        ifcomm= .FALSE. 
    else
        ifcomm= .TRUE. 
        icount=1
    end if

    il=0
    do id=0,np_neighbor-1

        if (ifcomm) then

            idsend=infosend(icount,1)

            if (id < idsend) then
                nsend=0
            else if (id == idsend) then
                nsend=infosend(icount,2)
                icount =icount+1
            end if

        else   !  not(ifcomm)
            nsend=0
        end if

        if (ifcomm .AND. (icount > npsend)) ifcomm= .FALSE. 

        len=isize
        call mpi_send (nsend,len,mpi_byte, id, nid, intercomm, ierr)

        if (nsend /= 0) then
        !     Send the points identity to communicating processors
            len=(ldim+1)*nsend*isize
            do i=1,nsend
                il=il+1
                do j=1,ldim+1
                    jsend((ldim+1)*(i-1)+j)=iList(j,il)
                enddo
            end do

            call mpi_send(jsend,len,mpi_byte,id,100+nid, intercomm,ierr)
        end if

    end do


!     Receive information about sending processors
    6 irecvcnt=0
    il=0
    do id=0,np_neighbor-1
        len=isize
        call mpi_recv (nrecv,len,mpi_byte,id,id, intercomm, status, ierr)
        if (nrecv /= 0) then
            irecvcnt=irecvcnt+1
            inforecv(irecvcnt,1)=id
            inforecv(irecvcnt,2)=nrecv

        !     Receive information about point identities
            len=(ldim+1)*nrecv*isize
            call mpi_recv(jrecv,len,mpi_byte,id,100+id,intercomm,status,ierr)
                     
        !     Receiving processor nid masks which points he receives from
        !     remote processor id as interpolation points using the identity
        !     information contained in array 'jrecv'. Those points are masked
        !     with imask=1 (imask=0 for all other points).

            do ii=1,nrecv
                ix=jrecv((ldim+1)*(ii-1)+1)
                iy=jrecv((ldim+1)*(ii-1)+2)
                iz = 1
                if (if3d) iz=jrecv((ldim+1)*(ii-1)+3)
                ie=jrecv((ldim+1)*(ii-1)+ldim+1)
                il=il+1
                iden(1,il)=ix
                iden(2,il)=iy
                if (if3d) iden(3,il)=iz
                iden(ldim+1,il)=ie
                imask(ix,iy,iz,ie)=1
            enddo

        end if  !  jrecv /= 0
              
    end do  !  all remote processors

    nprecv=irecvcnt
    1 format (a7,3i7,a7)

    return
          
    call neknekgsync()

    return
    end subroutine intpts_locate
      
!----------------------------------------------------------------
    subroutine get_values(which_field)
    include 'SIZE'
    include 'TOTAL'
    include 'NEKNEK'
    include 'mpif.h'

    character(3) :: which_field(nfld)
    real :: field(lx1*ly1*lz1*lelt,nfldmax)
    real :: rsend(nfldmax*nmaxl), rrecv(nfldmax*nmaxl)
    real :: fieldout(nfldmax,nmaxl)

    integer :: status(mpi_status_size)

!     Information about communication of points is contained in common
!     block /proclist/.
     
!     iden(1,:) = ix
!     iden(2,:) = iy
!     iden(3,:) = iz
!     iden(4,:) = iel


!     Put field values for the field ifld=1,nfld to the working array
!     field (:,:,:,:,nfld) used byt the findpts_value routine according
!     to the field identificator which_field(ifld)


    nv = nx1*ny1*nz1*nelv
    nt = nx1*ny1*nz1*nelt
    do ifld=1,nfld
        if (which_field(ifld) == 't' ) call copy(field(1,ifld),t ,nt)
        if (which_field(ifld) == 'vx') call copy(field(1,ifld),vx,nt)
        if (which_field(ifld) == 'vy') call copy(field(1,ifld),vy,nt)
        if (which_field(ifld) == 'vz') call copy(field(1,ifld),vz,nt)

    !     Find interpolation values

        call findpts_eval(inth_multi,fieldout(ifld,1),nfldmax, &
        rcode,1, &
        proc,1, &
        elid,1, &
        rst,ndim,npoints, &
        field(1,ifld))

    enddo

!     Send interpolation values to the corresponding processors
!     of remote session

    call neknekgsync()

    il=0
    do n=1,npsend
        id       = infosend(n,1)
        nsend    = infosend(n,2)
        do i=1,nsend
            il=il+1
            do ifld=1,nfld
                rsend(nfld*(i-1)+ifld)=fieldout(ifld,il)
            enddo
        enddo
        len=nfld*nsend*wdsize
        call mpi_send (rsend,len,mpi_byte, id, nid, intercomm, ierr)
    end do ! n=1,npsend

    il=0
    do n=1,nprecv
        id    = inforecv(n,1)
        nrecv = inforecv(n,2)
        len=nfld*nrecv*wdsize
        call mpi_recv (rrecv,len,mpi_byte,id,id,intercomm,status,ierr)

    !           Extract point identity
            
        do i=1,nrecv
            il=il+1
            ix = iden(1,il)
            iy = iden(2,il)
            iz = 1
            if (if3d) iz=iden(3,il)
            ie=iden(ldim+1,il)
            do ifld=1,nfld
                valint(ix,iy,iz,ie,ifld)=rrecv(nfld*(i-1)+ifld)
            enddo
        enddo      ! i=1,nrecv

    enddo ! n=1,nprecv

    call neknekgsync()

    return
    end subroutine get_values
!--------------------------------------------------------------------------
    subroutine userchk_set_xfer
    include 'SIZE'
    include 'TOTAL'
    include 'NEKNEK'
    real :: l2,linf
    character(3) :: which_field(nfldmax)

!     nfld is the number of fields to interpolate.
!     nfld = 3 for just veliocities, nfld = 4 for velocities + temperature

    which_field(1)='vx'
    which_field(2)='vy'
    which_field(3)='vz'
    if (nfld > 3) which_field(4)='t'

    if (nsessions > 1) call get_values(which_field)

    return
    end subroutine userchk_set_xfer
!------------------------------------------------------------------------
    subroutine bcopy
    include 'SIZE'
    include 'TOTAL'
    include 'NEKNEK'

    n    = nx1*ny1*nz1*nelt

    do k=1,nfld
        call copy(bdrylg(1,k,2),bdrylg(1,k,1),n)
        call copy(bdrylg(1,k,1),bdrylg(1,k,0),n)
        call copy(bdrylg(1,k,0),valint(1,1,1,1,k),n)
    enddo

!     Order of extrpolation is contolled by the parameter NINTER contained
!     in NEKNEK. First order interface extrapolation, NINTER=1 (time lagging)
!     is activated. It is unconditionally stable.  If you want to use
!     higher-order interface extrapolation schemes, you need to increase
!     ngeom to ngeom=3-5 for scheme to be stable.

    if (NINTER == 1 .OR. istep == 0) then
        c0=1
        c1=0
        c2=0
    else if (NINTER == 2 .OR. istep == 1) then
        c0=2
        c1=-1
        c2=0
    else
        c0=3
        c1=-3
        c2=1
    endif
         
    do k=1,nfld
        do i=1,n
            ubc(i,1,1,1,k) = &
            c0*bdrylg(i,k,0)+c1*bdrylg(i,k,1)+c2*bdrylg(i,k,2)
        enddo
    enddo

    return
    end subroutine bcopy
!---------------------------------------------------------------------
    subroutine setintercomm(nekcommtrue,nptrue)
    include 'SIZE'
    include 'GLOBALCOM'
    common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
          
    nekcommtrue=nekcomm
    nekcomm=iglobalcomm

    nptrue=np
    np=np_global

    return
    end subroutine setintercomm

!-----------------------------------------------------------------------
    subroutine unsetintercomm(nekcommtrue,nptrue)
    include 'SIZE'
    include 'GLOBALCOM'
    common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal

    nekcomm=nekcommtrue
    np=nptrue

    return
    end subroutine unsetintercomm

!-----------------------------------------------------------------------
    function uglmin(a,n)
    REAL :: A(1)
    DIMENSION TMP(1),WORK(1)

    call happy_check(1)
    call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
    uglmin=glmin(a,n)
    call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

    return
    end function uglmin

!-----------------------------------------------------------------------
    function uglamax(a,n)
    REAL :: A(1)
    DIMENSION TMP(1),WORK(1)

    call happy_check(1)
    call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
    uglamax=glamax(a,n)
    call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

    return
    end function uglamax
!------------------------------------------------------------------------
    subroutine neknekgsync()
    include 'SIZE'
    include 'GLOBALCOM'

    call happy_check(1)
    call mpi_barrier(intercomm,ierr)
    return
    end subroutine neknekgsync
!------------------------------------------------------------------------
    subroutine happy_check(ihappy)
    include 'SIZE'
    include 'TOTAL'
    common /happycallflag/ icall

!     Happy check
    call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
    iglhappy=iglmin(ihappy,1)
    call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original
    if (ihappy == 1 .AND. iglhappy == 0) then
        if (nid == 0) then
            WRITE (6,*) '       '
            WRITE (6,'(A,1i7,A,1e13.5)') &
            ' Emergency exit due to the other session:', &
            ISTEP,'   time =',TIME
            WRITE (6,*)
        end if
        icall=1
        call exitt
    end if

    return
    end subroutine happy_check
          
