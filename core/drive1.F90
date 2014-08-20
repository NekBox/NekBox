

!-----------------------------------------------------------------------
!> Initialize nek
!!
!! This includes reading the input files
    subroutine nek_init(intracomm)

    use ctimer
    use size_m
    use dealias, only : init_dealias

    include 'TOTAL'
    include 'DOMAIN'
    include 'ZPER'

    include 'OPCTR'

    real :: kwave2
    real*8 :: t0, tpp

    logical :: ifemati,ifsync_

    call get_session_info(intracomm)
          
!     Initalize Nek (MPI stuff, word sizes, ...)
!     call iniproc (initalized in get_session_info)

    etimes = dnekclock()
    istep  = 0
    tpp    = 0.0

    call opcount(1)

!     Initialize and set default values.
    call init_dealias

    call initdim
    call initdat
    call files

!     Read .rea +map file
    etime1 = dnekclock()
    call readat

    ifsync_ = ifsync
    ifsync = .TRUE. 

!     Initialize some variables
    call setvar

!     Check for zero steps
    instep=1
    if (nsteps == 0 .AND. fintim == 0.) instep=0

!     Setup domain topology
    igeom = 2
    call setup_topo

!     Compute GLL stuff (points, weights, derivate op, ...)
    call genwz

!     Initalize io unit
    call io_init

!     Set size for CVODE solver
#if 0
    if(ifcvode .AND. nsteps > 0) call cv_setsize(0,nfield)
#endif

!     USRDAT
    if(nid == 0) write(6,*) 'call usrdat'
    call usrdat
    if(nid == 0) write(6,'(A,/)') ' done :: usrdat'

!     generate geometry (called after usrdat in case something changed)
    call gengeom (igeom)

!max    if (ifmvbd) call setup_mesh_dssum ! Set up dssum for mesh (needs geom)

!     USRDAT2
    if(nid == 0) write(6,*) 'call usrdat2'
    call usrdat2
    if(nid == 0) write(6,'(A,/)') ' done :: usrdat2'

    call geom_reset(1)    ! recompute Jacobians, etc.
    call vrdsmsh          ! verify mesh topology

    call echopar ! echo back the parameter stack
    call setlog  ! Initalize logical flags

!     Zero out masks corresponding to Dirichlet boundary points.
    call bcmask

!     Need eigenvalues to set tolerances in presolve (SETICS)
    if (fintim /= 0.0 .OR. nsteps /= 0) call geneig (igeom)

!     Verify mesh topology
    call vrdsmsh

!     Pressure solver initialization (uses "SOLN" space as scratch)
    if (ifflow .AND. (fintime /= 0 .OR. nsteps /= 0)) then
        call estrat
        if (iftran .AND. solver_type == 'itr') then
            call set_overlap
        elseif (solver_type == 'fdm' .OR. solver_type == 'pdm')then
            write(*,*) "Oops: gfdm"
#if 0
            ifemati = .TRUE. 
            kwave2  = 0.0
            if (ifsplit) ifemati = .FALSE. 
            call gfdm_init(nx2,ny2,nz2,ifemati,kwave2)
#endif
        elseif (solver_type == '25D') then
            write(*,*) "Oops"
#if 0
            call g25d_init
#endif
        endif
    endif

!     USRDAT3
    if(nid == 0) write(6,*) 'call usrdat3'
    call usrdat3
    if(nid == 0) write(6,'(A,/)') ' done :: usrdat3'

!     Set initial conditions + compute field properties
    call setics
    call setprop

!     USRCHK
    if(instep /= 0) then
        if(nid == 0) write(6,*) 'call userchk'
        call userchk
        if(nid == 0) write(6,'(A,/)') ' done :: userchk'
    endif

!     Initialize CVODE
#if 0
    if(ifcvode .AND. nsteps > 0) call cv_init
#endif

    call comment
!    call sstest (isss)

    call dofcnt

    jp = 0  ! Set perturbation field count to 0 for baseline flow

!    call in_situ_init()

!     Initalize timers to ZERO
    call time00
    call opcount(2)

    etims0 = dnekclock_sync()
    IF (NID == 0) THEN
        WRITE (6,*) ' '
        IF (TIME /= 0.0) WRITE (6,'(A,E14.7)') ' Initial time:',TIME
        WRITE (6,'(A,g13.5,A)') &
        ' Initialization successfully completed ', &
        etims0-etimes, ' sec'
    ENDIF

    ifsync = ifsync_ ! restore initial value

    return
    end subroutine nek_init
!-----------------------------------------------------------------------
    subroutine nek_solve

    use ctimer
    use size_m
    include 'TSTEP'
    include 'INPUT'

    real*4 :: papi_mflops
    integer*8 :: papi_flops

    call nekgsync()

    if (instep == 0) then
        if(nid == 0) write(6,'(/,A,/,A,/)') &
        ' nsteps=0 -> skip time loop', &
        ' running solver in post processing mode'
    else
        if(nid == 0) write(6,'(/,A,/)') 'Starting time loop ...'
    endif

    isyc  = 0
    itime = 0
    if(ifsync) isyc=1
    itime = 1
    call nek_comm_settings(isyc,itime)

    call nek_comm_startstat()

    istep  = 0
    msteps = 1

    do kstep=1,nsteps,msteps
        call nek__multi_advance(kstep,msteps)
        call userchk
        call prepost ( .FALSE. ,'his')
!        call in_situ_check()
        if (lastep == 1) goto 1001
    enddo
    1001 lastep=1


    call nek_comm_settings(isyc,0)

    call comment

!     check for post-processing mode
    if (instep == 0) then
        nsteps=0
        istep=0
        if(nid == 0) write(6,*) 'call userchk'
        call userchk
        if(nid == 0) write(6,*) 'done :: userchk'
        call prepost ( .TRUE. ,'his')
    else
        if (nid == 0) write(6,'(/,A,/)') &
        'end of time-step loop'
    endif


    RETURN
    end subroutine nek_solve
!-----------------------------------------------------------------------
    subroutine nek_advance

    use ctimer
    use size_m
    include 'TOTAL'

    common /cgeom/ igeom

    call nekgsync
    IF (IFTRAN) CALL SETTIME
    if (ifmhd ) call cfl_check
    CALL SETSOLV
    CALL COMMENT

    if (ifsplit) then   ! PN/PN formulation

        igeom = 1
        if (ifheat)          call heat     (igeom)
        call setprop
        call qthermal
        igeom = 1
        if (ifflow)          call fluid    (igeom)
        if (param(103) > 0) call q_filter(param(103))
        call setup_convect (2) ! Save convective velocity _after_ filter

    else                ! PN-2/PN-2 formulation
      write(*,*) "Oops! Pn-2/Pn-2"
#if 0
        call setprop
        do igeom=1,ngeom

            if (igeom > 2) call userchk_set_xfer

            if (ifgeom) then
                call gengeom (igeom)
                call geneig  (igeom)
            endif

            if (ifmhd) then
                if (ifheat)      call heat     (igeom)
                call induct   (igeom)
            elseif (ifpert) then
                write(*,*) "Oops! ifpert"
#if 0
                if (ifbase .AND. ifheat)  call heat          (igeom)
                if (ifbase .AND. ifflow)  call fluid         (igeom)
                if (ifflow)             call fluidp        (igeom)
                if (ifheat)             call heatp         (igeom)
#endif
            else  ! std. nek case
                if (ifheat)             call heat          (igeom)
                if (ifflow)             call fluid         (igeom)
                if (ifmvbd)             call meshv         (igeom)
            endif

            if (igeom == ngeom .AND. param(103) > 0) &
            call q_filter(param(103))

            call setup_convect (igeom) ! Save convective velocity _after_ filter

        enddo
#endif
    endif

    return
    end subroutine nek_advance
!-----------------------------------------------------------------------
    subroutine nek_end

    use size_m
    include 'TSTEP'
    include 'PARALLEL'
    include 'OPCTR'

    if(instep /= 0)  call runstat
    if(xxth(1) > 0) call crs_stats(xxth(1))

       
!    call in_situ_end()
    return
    end subroutine nek_end
!-----------------------------------------------------------------------
    subroutine nek__multi_advance(kstep,msteps)

    use size_m
    include 'TOTAL'

    do i=1,msteps
        istep = istep+i
        call nek_advance
    enddo

    return
    end subroutine nek__multi_advance
!-----------------------------------------------------------------------
