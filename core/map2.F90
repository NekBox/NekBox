!-----------------------------------------------------------------------
    subroutine mapelpr()
    use size_m
    include 'INPUT'
    include 'PARALLEL'
    include 'SCRCT'
    include 'SOLN'
    include 'TSTEP'

    logical :: ifverbm

    if(nid == 0) write(6,'(/,A)') ' mapping elements to processors'

    MFIELD=2
    IF (IFFLOW) MFIELD=1
    IF (IFMVBD) MFIELD=0

!     Set up TEMPORARY value for NFIELD - NFLDT
    NFLDT = 1
    IF (IFHEAT) NFLDT = 2 + NPSCAL

!     Distributed memory processor mapping
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

!     Output the processor-element map:
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

!     Check elemental distribution

!      IF (IPASS.EQ.2.AND.PARAM(156).eq.9) THEN
!         NXYZ=NX1*NY1*NZ1
!         DO 1400 IE=1,NELT
!            VTMP1=NODE
!            VTMP2=IE
!            CALL CFILL(VX(1,1,1,IE) ,VTMP1,NXYZ)
!            CALL CFILL(VY(1,1,1,IE) ,VTMP2,NXYZ)
!            CALL CFILL(T(1,1,1,IE,1),VTMP1,NXYZ)
! 1400    CONTINUE
!         call prepost(.true.,'   ')
!      ENDIF

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
    subroutine set_proc_map()

!     Compute element to processor distribution according to (weighted)
!     physical distribution in an attempt to minimize exposed number of
!     element interfaces.

    use size_m
    include 'INPUT'
    include 'PARALLEL'
    include 'SOLN'
    include 'SCRCT'
    include 'TSTEP'
    include 'ZPER'
    common /ctmp0/ iwork(lelt)

    REAL*8 :: dnekclock,t0

    t0 = dnekclock()
!     if (.not.(ifgtp.or.ifgfdm)) then
    if ( .NOT. ifgtp) then
    
    !        rsb element to processor mapping
    
        if (ifgfdm)       call gfdm_elm_to_proc(gllnid,np) ! gfdm w/ .map

        call get_map

    endif

    if(ifzper .OR. ifgtp) call gfdm_elm_to_proc(gllnid,np) ! special processor map

!     compute global to local map (no processor info)

    if ( .NOT. ifmoab) then
        IEL=0
        CALL IZERO(GLLEL,NELGT)
        DO IEG=1,NELGT
            IF (GLLNID(IEG) == NID) THEN
                IEL = IEL + 1
                GLLEL(IEG)=IEL
                NELT = IEL
                if (ieg <= nelgv) NELV = IEL
            ENDIF
        !        write(6,*) 'map2 ieg:',ieg,nelv,nelt,nelgv,nelgt
        ENDDO
    
    !     dist. global to local map to all processors
    
        npass = 1 + nelgt/lelt
        k=1
        do ipass = 1,npass
            m = nelgt - k + 1
            m = min(m,lelt)
            if (m > 0) call igop(gllel(k),iwork,'+  ',m)
            k = k+m
        enddo
    endif

!     compute local to global map
!     (i.e. returns global element number given local index and proc id)

    do ieg=1,nelgt
        mid  =gllnid(ieg)
        ie   =gllel (ieg)
        if (mid == nid) lglel(ie)=ieg
    enddo

!     All Done.

    return
    end subroutine set_proc_map
!-----------------------------------------------------------------------
    subroutine gfdm_elm_to_proc(gllnid,np)


    use size_m
    include 'ZPER'

    integer :: gllnid(1)

    common /ctmp1/  map_st(lelg_sm)
    common /vptsol/ iwork(0:lp)
    integer :: nelbox(3),nstride_box(3)

    call gfdm_set_pst(ip,is,it,nelbox,nstride_box,nx2,ny2,nz2)

    nep = nelbox(ip)
    nes = nelbox(is)
         
    if(nelbox(it) == 0) nelbox(it)=1
    net = nelbox(it)

    nst = nes*net
    if (nst < np) then
        if (nid == 0) &
        write(6,*) 'ERROR, number of elements in plane must be > np' &
        ,nst,np,nep,nes,net
        call exitt
    endif


    call gfdm_map_2d(map_st,nes,net,iwork,np)
    call gfdm_build_global_el_map (gllnid,map_st,nes,net &
    ,nelbox,nstride_box,ip,is,it)

    return
    end subroutine gfdm_elm_to_proc
!-----------------------------------------------------------------------
    subroutine gfdm_map_2d(map_st,nes,net,num_el,np)

!     Set up a 2D processor decomposition of an NES x NET array.

    integer :: map_st(nes,net),num_el(0:np)

!     First a stupid dealing algorithm to determine how many
!     elements on each processor

    do i=0,np
        num_el(i) = 0
    enddo

    k = np-1
    do i=1,nes*net
        num_el(k) = num_el(k)+1
        k=k-1
        if (k < 0) k = np-1
    enddo

    jnid = 0
    nel_cnt = 0
    nel_cur = num_el(jnid)
    do j=1,net,2
        do i=1,nes                 ! Count down
            nel_cnt = nel_cnt + 1
            if (nel_cnt > nel_cur) then
                jnid=jnid+1
                nel_cur = num_el(jnid)
                nel_cnt = 1
            endif
            map_st(i,j) = jnid
        enddo
    
        j1 = j+1
        if (j1 <= net) then
            do i=nes,1,-1                ! Count up
                nel_cnt = nel_cnt + 1
                if (nel_cnt > nel_cur) then
                    jnid=jnid+1
                    nel_cur = num_el(jnid)
                    nel_cnt = 1
                endif
                map_st(i,j1) = jnid
            enddo
        endif
    enddo

    return
    end subroutine gfdm_map_2d
!-----------------------------------------------------------------------
    subroutine gfdm_set_pst(ip,is,it,nelbox,nstride_box,nxp,nyp,nzp)

    use size_m
    include 'INPUT'
    include 'ZPER'

    integer :: nelbox(3),nstride_box(3)

    if (if3d) then
        if (param(118) < 0) then
            ip = 3
            is = 1
            it = 2
        elseif (param(117) < 0) then
            ip = 2
            is = 3
            it = 1
        else
            ip = 1
            is = 2
            it = 3
        endif
    else
        if (param(117) < 0) then
            ip = 2
            is = 1
            it = 3
        else
            ip = 1
            is = 2
            it = 3
        endif
    endif

    pst2lex(1)=ip       ! identify x-, y- or z with primary direction
    pst2lex(2)=is
    pst2lex(3)=it

    lex2pst(ip)=1
    lex2pst(is)=2
    lex2pst(it)=3

    nelbox(1) = nelx
    nelbox(2) = nely
    nelbox(3) = nelz

    nstride_box(1) = 1
    nstride_box(2) = nelx
    nstride_box(3) = nelx*nely

    ngfdm_p(1) = nelx*nxp
    ngfdm_p(2) = nely*nyp
    ngfdm_p(3) = nelz*nzp
    write(6,*) 'ngfdm:',(ngfdm_p(k),k=1,3)

    return
    end subroutine gfdm_set_pst
!-----------------------------------------------------------------------
    subroutine gfdm_build_global_el_map (gllnid,map_st,nes,net &
    ,nelbox,nstride_box,ip,is,it)

    use size_m
    integer :: gllnid(1)
    integer :: map_st(nes,net)
    integer :: nelbox(3),nstride_box(3)

    integer :: proc

    do jt=1,nelbox(it)
        do js=1,nelbox(is)
            proc = map_st(js,jt)
            do jp=1,nelbox(ip)
                ieg = 1 + nstride_box(ip)*(jp-1)   & ! nstride_p=nes*net
                + nstride_box(is)*(js-1) &
                + nstride_box(it)*(jt-1)
                gllnid(ieg) = proc
            enddo
        enddo
    enddo

    return
    end subroutine gfdm_build_global_el_map
!-----------------------------------------------------------------------
    subroutine outmati(u,m,n,name6)
    integer :: u(m,n)
    character(6) :: name6
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

!     Print out copies of a global matrix

    do mid=0,np-1
        call nekgsync
        if (mid == nid) then
            n20 = min(n,20)
            write(6,1) nid,m,n,name6
            1 format(//,3i6,'  Matrix:',2x,a6,/)
            do i=1,m
                write(6,2) nid,name6,(u(i,j),j=1,n20)
            enddo
            2 format(i3,1x,a6,20i6)
        endif
        call nekgsync
    enddo
    return
    end subroutine outmati
!-----------------------------------------------------------------------
    subroutine outmati8(u,m,n,name6)
    integer*8 :: u(m,n)
    character(6) :: name6
    common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

!     Print out copies of a global matrix

    do mid=0,np-1
        call nekgsync
        if (mid == nid) then
            n20 = min(n,20)
            write(6,1) nid,m,n,name6
            1 format(//,3i6,'  Matrix:',2x,a6,/)
            do i=1,m
                write(6,2) nid,name6,(u(i,j),j=1,n20)
            enddo
            2 format(i3,1x,a6,20i6)
        endif
        call nekgsync
    enddo
    return
    end subroutine outmati8
!-----------------------------------------------------------------------
    subroutine get_map

    call get_vert

    return
    end subroutine get_map
!-----------------------------------------------------------------------
