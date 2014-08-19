!---------------------------------------------------------------------
! VisIt Simulation Code for Nek5000

!   Provide a connection to VisIt simulation code. You will be able
! to connect to VisIt and drive the simulation though VisIt or
! have sim code send data to VisIt.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! visit_init

!   Initialize VisIt. This will create the sim file needed by Visit
! and create the initial VisIt setup.
!---------------------------------------------------------------------
    subroutine visit_init
    implicit none
    include "visitfortransimV2interface.inc"
    include "mpif.h"
!     // local variables
    integer :: err
!     // SIMSTATE common block
    integer :: runflag, endflag
    common /SIMSTATE/ runflag, endflag
    save /SIMSTATE/
!     // PARALLEL state common block
    integer :: par_rank, par_size
    common /PARALLEL/ par_rank, par_size
    save /PARALLEL/

    #ifdef VISIT_STOP
!     // The sim will wait for VisIt to connect after first step.
    runflag = 0
    #else
!     // Default is to run sim and VisIt can connect any time.
    runflag = 1
    #endif
    endflag = 0

!     // Determine the rank and size of this MPI task so we can tell
!     // VisIt's libsim about it.
    call MPI_COMM_RANK(MPI_COMM_WORLD, par_rank, err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, par_size, err)
    if(par_size > 1) then
        err = visitsetparallel(1)
    endif
    err = visitsetparallelrank(par_rank)

    call simulationarguments()
!     // TODO: look for the visitsetupenv2 function.
!     // Has better scaling, but has not been release for fortran.
    err = visitsetupenv()
!     // Have the master process create the sim file.
    if(par_rank == 0) then
        err = visitinitializesim("nek5000", 7, &
        "Nek5000 Simulation", 18, &
        "/no/useful/path", 15, &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN)
    endif

    end subroutine visit_init

!---------------------------------------------------------------------
! simulationarguments
!     This routine handles command line arguments
!     -dir <VisIt directory>
!     -options <VisIt Options>
!     -trace <VisIt trace file>
!---------------------------------------------------------------------
    subroutine simulationarguments()
    implicit none
    character (len=80) str
    integer :: err, i, N, len
    integer :: visitsetoptions, visitsetdirectory, visitopentracefile

    N = iargc()
    i = 1
    len = 80
    5 if (i <= N) then
        call getarg(i, str)
        if(str == "-dir") then
            call getarg(i+1, str)
            err = visitsetdirectory(str, len)
            i = i + 1
        elseif(str == "-options") then
            call getarg(i+1, str)
            err = visitsetoptions(str, len)
            i = i + 1
        elseif(str == "-trace") then
            call getarg(i+1, str)
            err = visitopentracefile(str, len)
            i = i + 1
        endif
        i = i + 1
        goto 5
    endif
    end subroutine simulationarguments

!---------------------------------------------------------------------
! processvisitcommand
!---------------------------------------------------------------------
    integer function processvisitcommand()
    implicit none
    include "mpif.h"
    include "visitfortransimV2interface.inc"
!     // PARALLEL state common block
    integer :: par_rank, par_size
    common /PARALLEL/ par_rank, par_size
    integer :: command, e, doloop, success, ret
    integer :: VISIT_COMMAND_PROCESS
    integer :: VISIT_COMMAND_SUCCESS
    integer :: VISIT_COMMAND_FAILURE
    parameter (VISIT_COMMAND_PROCESS = 0)
    parameter (VISIT_COMMAND_SUCCESS = 1)
    parameter (VISIT_COMMAND_FAILURE = 2)

    if(par_rank == 0) then
        success = visitprocessenginecommand()

        if(success > 0) then
            command = VISIT_COMMAND_SUCCESS
            ret = 1
        else
            command = VISIT_COMMAND_FAILURE
            ret = 0
        endif

        call MPI_BCAST(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,e)
    else
        doloop = 1
        2345 call MPI_BCAST(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,e)
        if(command == VISIT_COMMAND_PROCESS) then
            success = visitprocessenginecommand()
        elseif(command == VISIT_COMMAND_SUCCESS) then
            ret = 1
            doloop = 0
        else
            ret = 0
            doloop = 0
        endif

        if(doloop /= 0) then
            goto 2345
        endif
    endif
    processvisitcommand = ret
    end function processvisitcommand

!---------------------------------------------------------------------
! visit_check
!---------------------------------------------------------------------
    subroutine visit_check()
    implicit none
    include "mpif.h"
    include "visitfortransimV2interface.inc"
!     // functions
    integer :: processvisitcommand
!     // local variables
    integer :: visitstate, result, blocking, ierr
!     // SIMSTATE common block
    integer :: runflag, endflag
    common /SIMSTATE/ runflag, endflag
!     // PARALLEL state common block
    integer :: par_rank, par_size
    common /PARALLEL/ par_rank, par_size

!     // If at the end of sim, lets not force an update.
    if(endflag == 0) then
    !        // Check if we are connected to VisIt
        result = visitisconnected()
        if(result == 1) then
        !            // Tell VisIt that the timestep changed
            result = visittimestepchanged()
        !            // Tell VisIt to update its plots
            result = visitupdateplots()
        endif
    endif

!     write (6, 2000)
    2000 format('VisIt Check!')
    do 10
    !         // If we are running don't block
        if(runflag == 1) then
            blocking = 0
        else
            blocking = 1
        endif

    !         // Detect input from VisIt on processor 0 and then broadcast
    !         // the results of that input to all processors.
        if(par_rank == 0) then
            visitstate = visitdetectinput(blocking, -1)
        endif
        call MPI_BCAST(visitstate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (visitstate == 0) then
        !             // Okay - Process time step.
            goto 1234
        elseif (visitstate == 1) then
        !             // Attempt to Connect VisIt
            ierr = runflag
            runflag = 0
            result = visitattemptconnection()
            if (result == 1) then
                write (6, 2001)
                2001 format('VisIt connected!')
            else
                write (6, 2002)
                2002 format('VisIt did not connected!')
            endif
            flush( 6 )
            runflag = ierr
        elseif (visitstate == 2) then
        !             // Engine socket input
        !             ierr = runflag
        !             runflag = 0
            if (processvisitcommand() == 0) then
                result = visitdisconnect()
            !                 // If VisIt is disconnect lets run.
                runflag = 1
            endif
        !             // Check if user wants to exit sim.
            if(runflag == 2) then
                goto 1234
            endif
        elseif (visitstate == 3) then
        !             // Console socket input
        elseif (visitstate < 0) then
        !             // Error
            goto 1234
        endif
    10 END DO
    1234 end

!---------------------------------------------------------------------
! visit_end
!---------------------------------------------------------------------
    subroutine visit_end()
    implicit none
    include "visitfortransimV2interface.inc"
!     // local variables
    integer :: result
!     // SIMSTATE common block
    integer :: runflag, endflag
    common /SIMSTATE/ runflag, endflag

!     // Check if we are connected to VisIt
    result = visitisconnected()
    if(result == 1) then
    !        // Let VisIt exit the sim.

    !        // This will tell the visit_check function we are at the end.
        endflag = 1
        runflag = 0

        do 10
            call visit_check()

        !           // User asked to finish.
            if (endflag == 2) then
                EXIT
            endif
        10 END DO
    endif

    end subroutine visit_end

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! These functions must be defined to satisfy the
! visitfortransimV2interface lib.

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!---------------------------------------------------------------------
! visitcommandcallback
!   Handle User defined functions from VisIt user interface.
!---------------------------------------------------------------------
    subroutine visitcommandcallback (cmd, lcmd, args, largs)
    implicit none
    character(8) :: cmd, args
    integer ::     lcmd, largs
    include "visitfortransimV2interface.inc"
!     // SIMSTATE common block
    integer :: runflag, endflag
    common /SIMSTATE/ runflag, endflag

!     // Handle the commands that we define in visitgetmetadata.
    if(visitstrcmp(cmd, lcmd, "stop", 4) == 0) then
        runflag = 0
    elseif(visitstrcmp(cmd, lcmd, "step", 4) == 0) then
        runflag = 2
    elseif(visitstrcmp(cmd, lcmd, "run", 3) == 0) then
        runflag = 1
    elseif(visitstrcmp(cmd, lcmd, "exit", 4) == 0) then
        call exitt()
    elseif(visitstrcmp(cmd, lcmd, "finish", 6) == 0) then
        if(endflag == 1) then
            endflag = 2
            runflag = 2
        endif
    endif
    end subroutine visitcommandcallback

!---------------------------------------------------------------------
! visitbroadcastintfunction
!---------------------------------------------------------------------
    integer function visitbroadcastintfunction(value, sender)
    implicit none
    include "mpif.h"
    integer :: value, sender, ierr
    call MPI_BCAST(value,1,MPI_INTEGER,sender,MPI_COMM_WORLD,ierr)
    visitbroadcastintfunction = 0
    end function 

!---------------------------------------------------------------------
! visitbroadcaststringfunction
!---------------------------------------------------------------------
    integer function visitbroadcaststringfunction(str, lstr, sender)
    implicit none
    include "mpif.h"
    character(8) :: str
    integer ::     lstr, sender, ierr
    call MPI_BCAST(str,lstr,MPI_CHARACTER,sender,MPI_COMM_WORLD,ierr)
    visitbroadcaststringfunction = 0
    end function 

!---------------------------------------------------------------------
! visitslaveprocesscallback
!---------------------------------------------------------------------
    subroutine visitslaveprocesscallback ()
    implicit none
    include "mpif.h"
    integer :: c, ierr, VISIT_COMMAND_PROCESS
    parameter (VISIT_COMMAND_PROCESS = 0)
    c = VISIT_COMMAND_PROCESS
    call MPI_BCAST(c, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end subroutine visitslaveprocesscallback

!---------------------------------------------------------------------
! visitactivatetimestep
!---------------------------------------------------------------------
    integer function visitactivatetimestep()
    implicit none
    include "visitfortransimV2interface.inc"
    visitactivatetimestep = VISIT_OKAY
    end function visitactivatetimestep

!---------------------------------------------------------------------
! visitgetmetadata
!   This function tells VisIt about all of the data that this
!   simulation will output. Mesh, variables, expressions and commands.
!---------------------------------------------------------------------
    integer function visitgetmetadata()
    include 'SIZE'
    include 'TSTEP'
    include 'INPUT'
    INCLUDE 'PARALLEL'
    include "visitfortransimV2interface.inc"
!     // SIMSTATE common block
    integer :: runflag, endflag
    common /SIMSTATE/ runflag, endflag
!     // PARALLEL state common block
    integer :: par_rank, par_size
    common /PARALLEL/ par_rank, par_size
!     // Local variables
    integer :: md, mmd, vmd, cmd, emd, err, iDim, k
    character(3) :: seqstring

    if(visitmdsimalloc(md) == VISIT_OKAY) then
        err = visitmdsimsetcycletime(md, ISTEP, TIME)
        if(runflag == 1) then
            err = visitmdsimsetmode(md, VISIT_SIMMODE_RUNNING)
        else
            err = visitmdsimsetmode(md, VISIT_SIMMODE_STOPPED)
        endif

    !       // TODO: ask why this changes after first save?
    !       write(6,*) 'VISIT: IFXYO',IFXYO
    !       if(IFXYO) then
    !         // Add a mesh
        if(IF3D) then
            iDim = 3
        else
            iDim = 2
        endif

        if(visitmdmeshalloc(mmd) == VISIT_OKAY) then
            err = visitmdmeshsetname(mmd, "mesh", 4)
            err = visitmdmeshsetmeshtype(mmd, &
            VISIT_MESHTYPE_CURVILINEAR)
            err = visitmdmeshsettopologicaldim(mmd, iDim)
            err = visitmdmeshsetspatialdim(mmd, iDim)
            err = visitmdmeshsetnumdomains(mmd, NELGT)
            err = visitmdmeshsetdomaintitle(mmd, "Domains", 7)
            err = visitmdmeshsetdomainpiecename(mmd, "domain", 6)
        !           err = visitmdmeshsetxunits(mmd, "cm", 2)
        !           err = visitmdmeshsetyunits(mmd, "cm", 2)
            err = visitmdmeshsetxlabel(mmd, "X-Axis", 6)
            err = visitmdmeshsetylabel(mmd, "Y-Axis", 6)
            if(IF3D) then
            !             err = visitmdmeshsetzunits(mmd, "cm", 2)
                err = visitmdmeshsetzlabel(mmd, "Z-Axis", 6)
            endif
            err = visitmdsimaddmesh(md, mmd)
        endif
    !       endif

        if(IFVO) then
        !         // Add a X velocity variable on the mesh.
            if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                err = visitmdvarsetname(vmd, "x_velocity", 10)
                err = visitmdvarsetmeshname(vmd, "mesh", 4)
            !           err = visitmdvarsetunits(vmd, "cm", 2)
                err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
            !           err = visitmdvarsettreatasascii(vmd, 0)
                err = visitmdsimaddvariable(md, vmd)
            endif

        !         // Add a Y velocity variable on the mesh.
            if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                err = visitmdvarsetname(vmd, "y_velocity", 10)
                err = visitmdvarsetmeshname(vmd, "mesh", 4)
            !           err = visitmdvarsetunits(vmd, "cm", 2)
                err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
            !           err = visitmdvarsettreatasascii(vmd, 0)
                err = visitmdsimaddvariable(md, vmd)
            endif

            if(IF3D) then
            !           // Add a Z velocity variable on the mesh.
                if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                    err = visitmdvarsetname(vmd, "z_velocity", 10)
                    err = visitmdvarsetmeshname(vmd, "mesh", 4)
                !             err = visitmdvarsetunits(vmd, "cm", 2)
                    err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                    err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
                !             err = visitmdvarsettreatasascii(vmd, 0)
                    err = visitmdsimaddvariable(md, vmd)
                endif
            endif

        !         // Add a velocity expression.
            if(visitmdexpralloc(emd) == VISIT_OKAY) then
                err = visitmdexprsetname(emd, "velocity", 8)
                if(IF3D) then
                    err = visitmdexprsetdefinition(emd, &
                    "{x_velocity, y_velocity, z_velocity}", 36)
                else
                    err = visitmdexprsetdefinition(emd, &
                    "{x_velocity, y_velocity}", 24)
                endif
                err = visitmdexprsettype(emd, VISIT_VARTYPE_VECTOR)

                err = visitmdsimaddexpression(md, emd)
            endif

        !         // Add a velocity magnitude expression.
            if(visitmdexpralloc(emd) == VISIT_OKAY) then
                err = visitmdexprsetname(emd, "velocity_mag", 12)
                err = visitmdexprsetdefinition(emd, &
                "magnitude(velocity)", 19)
                err = visitmdexprsettype(emd, VISIT_VARTYPE_SCALAR)

                err = visitmdsimaddexpression(md, emd)
            endif
        endif

        if(IFPO) then
        !         // Add a pressure variable on the mesh.
            if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                err = visitmdvarsetname(vmd, "pressure", 8)
                err = visitmdvarsetmeshname(vmd, "mesh", 4)
            !           err = visitmdvarsetunits(vmd, "cm", 2)
                err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
            !           err = visitmdvarsettreatasascii(vmd, 0)
                err = visitmdsimaddvariable(md, vmd)
            endif
        endif

        if(IFTO) then
        !         // Add a temperature variable on the mesh.
            if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                err = visitmdvarsetname(vmd, "temperature", 11)
                err = visitmdvarsetmeshname(vmd, "mesh", 4)
            !           err = visitmdvarsetunits(vmd, "cm", 2)
                err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
            !           err = visitmdvarsettreatasascii(vmd, 0)
                err = visitmdsimaddvariable(md, vmd)
            endif
        endif

        do k=1,LDIMT-1
        !         // Add a user defined variable on the mesh.
            if(IFPSCO(k)) then
                if(visitmdvaralloc(vmd) == VISIT_OKAY) then
                    write (seqstring,'(I0)') k
                    err = visitmdvarsetname(vmd, "s"//trim(seqstring), &
                    &                                 1+len_trim(seqstring))
                    err = visitmdvarsetmeshname(vmd, "mesh", 4)
                !             err = visitmdvarsetunits(vmd, "cm", 2)
                    err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
                    err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
                !             err = visitmdvarsettreatasascii(vmd, 0)
                    err = visitmdsimaddvariable(md, vmd)
                endif
            endif
        enddo


    !       // Add simulation commands
        err = visitmdcmdalloc(cmd)
        if(err == VISIT_OKAY) then
            err = visitmdcmdsetname(cmd, "stop", 4)
            err = visitmdsimaddgenericcommand(md, cmd)
        endif
        err = visitmdcmdalloc(cmd)
        if(err == VISIT_OKAY) then
            err = visitmdcmdsetname(cmd, "step", 4)
            err = visitmdsimaddgenericcommand(md, cmd)
        endif
        err = visitmdcmdalloc(cmd)
        if(err == VISIT_OKAY) then
            err = visitmdcmdsetname(cmd, "run", 3)
            err = visitmdsimaddgenericcommand(md, cmd)
        endif
        err = visitmdcmdalloc(cmd)
        if(err == VISIT_OKAY) then
            err = visitmdcmdsetname(cmd, "exit", 4)
            err = visitmdsimaddgenericcommand(md, cmd)
        endif
        err = visitmdcmdalloc(cmd)
        if(err == VISIT_OKAY) then
            err = visitmdcmdsetname(cmd, "finish", 6)
            err = visitmdsimaddgenericcommand(md, cmd)
        endif
    endif
    visitgetmetadata = md
    end function visitgetmetadata

!---------------------------------------------------------------------
! visitgetmesh
!    Use this function to return mesh data to VisIt.
!---------------------------------------------------------------------
    integer function visitgetmesh(domain, name, lname)
    include 'SIZE'
    include 'TOTAL'
    character(8) :: name
    integer ::     domain, lname
    include "visitfortransimV2interface.inc"

!     // local variables
    integer :: vh, x, y, z, err, dl, dmdims(3)

    vh = VISIT_INVALID_HANDLE

!     // Fortran starting index 1, but VisIt is 0
!     // Also need to get local domain index
    domain = GLLEL(domain + 1)

    if(visitstrcmp(name, lname, "mesh", 4) == 0) then
        if(visitcurvmeshalloc(vh) == VISIT_OKAY) then
            err = visitvardataalloc(x)
            err = visitvardataalloc(y)
            if(IF3D) err = visitvardataalloc(z)

            dl = NX1 * NY1 * NZ1
            err = visitvardatasetd(x, VISIT_OWNER_SIM, 1, dl, &
            XM1(1,1,1,domain))
            err = visitvardatasetd(y, VISIT_OWNER_SIM, 1, dl, &
            YM1(1,1,1,domain))
            if(IF3D) then
                err = visitvardatasetd(z, VISIT_OWNER_SIM, 1, dl, &
                ZM1(1,1,1,domain))
            endif

            dmdims(1) = NX1
            dmdims(2) = NY1
            dmdims(3) = NZ1
            if(IF3D) then
                err = visitcurvmeshsetcoordsxyz(vh, dmdims, x, y, z)
            else
                err = visitcurvmeshsetcoordsxy(vh, dmdims, x, y)
            endif
        endif
    endif

    visitgetmesh = vh
    end function visitgetmesh

!---------------------------------------------------------------------
! visitgetmaterial
!    Use this function to return material data to VisIt.
!---------------------------------------------------------------------
    integer function visitgetmaterial(domain, name, lname)
    implicit none
    character(8) :: name
    integer ::     domain, lname
    include "visitfortransimV2interface.inc"
    visitgetmaterial = VISIT_ERROR
    end function visitgetmaterial

!---------------------------------------------------------------------
! visitgetvariable
!    Use this function to return variable data to VisIt.
!---------------------------------------------------------------------
    integer function visitgetvariable(gdomain, name, lname)
!     implicit none
    include 'SIZE'
    include 'TOTAL'
    character(8) :: name
    integer ::     gdomain, lname
    include "visitfortransimV2interface.inc"
!     // local vars
    integer :: h, nvals, err, domain, k

    nvals = NX1 * NY1 * NZ1

    h = VISIT_INVALID_HANDLE

!     // Fortran starting index 1, but VisIt is 0
!     // Also need to get local domain index
    gdomain = gdomain + 1
    domain = GLLEL(gdomain)

    if(visitvardataalloc(h) == VISIT_OKAY) then
        if(visitstrcmp(name, lname, "temperature", 11) == 0) then
            err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
            T(1,1,1,domain,1))
        elseif(visitstrcmp(name, lname, "pressure", 8) == 0) then
            if(gdomain <= NELGV) then
                err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
                PR(1,1,1,domain))
            endif
        elseif(visitstrcmp(name, lname, "x_velocity", 10) == 0) then
            if(gdomain <= NELGV) then
                err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
                VX(1,1,1,domain))
            endif
        elseif(visitstrcmp(name, lname, "y_velocity", 10) == 0) then
            if(gdomain <= NELGV) then
                err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
                VY(1,1,1,domain))
            endif
        elseif(visitstrcmp(name, lname, "z_velocity", 10) == 0) then
            if(gdomain <= NELGV) then
                err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
                VZ(1,1,1,domain))
            endif
        elseif(visitstrcmp("s", 1, name, 1) == 0) then
        !             // Handle the user define variables.
            read( name(2:lname), '(i10)' ) k
            k = k + 1
            if(IFPSCO(k)) then
                err = visitvardatasetd(h, VISIT_OWNER_SIM, 1, nvals, &
                T(1,1,1,domain,k))
            endif
        endif
    endif

    visitgetvariable = h
    end function visitgetvariable

!---------------------------------------------------------------------
! visitgetcurve
!    Use this function to return curve data to VisIt.
!---------------------------------------------------------------------
    integer function visitgetcurve(handle, name, lname)
    implicit none
    character(8) :: name
    integer ::     handle, lname
    include "visitfortransimV2interface.inc"
    visitgetcurve = VISIT_ERROR
    end function visitgetcurve

!---------------------------------------------------------------------
! visitgetdomainlist
!    This function returns a list of domains owned by this process.
!---------------------------------------------------------------------
    integer function visitgetdomainlist()
    INCLUDE 'SIZE'
    INCLUDE 'PARALLEL'
    include "visitfortransimV2interface.inc"
!     // local vars
    integer :: h, dl, err

    h = VISIT_INVALID_HANDLE

    if(visitdomainlistalloc(h) == VISIT_OKAY) then
        if(visitvardataalloc(dl) == VISIT_OKAY) then
        !             // Hack to work around the index difference between
        !             // fortran and VisIt C. Temp change and make copy.
            do i=1,NELT
                LGLEL(i) = LGLEL(i) - 1
            enddo
            err = visitvardataseti(dl, VISIT_OWNER_COPY,1,NELT,LGLEL)
            err = visitdomainlistsetdomains(h, NELGT, dl)
        !             // restore correct fortran values.
            do i=1,NELT
                LGLEL(i) = LGLEL(i) + 1
            enddo
        endif
    endif

    visitgetdomainlist = h
    end function visitgetdomainlist

!---------------------------------------------------------------------
! visitgetdomainbounds
!    This function allows VisIt to create ghost zones between domains.
!---------------------------------------------------------------------
    integer function visitgetdomainbounds(name, lname)
    implicit none
    character(8) :: name
    integer ::     lname
    include "visitfortransimV2interface.inc"
    visitgetdomainbounds = VISIT_INVALID_HANDLE
    end function visitgetdomainbounds

!---------------------------------------------------------------------
! visitgetdomainnesting
!    This is used to tell VisIt how AMR patches are nested.
!---------------------------------------------------------------------
    integer function visitgetdomainnesting(name, lname)
    implicit none
    character(8) :: name
    integer ::     lname
    include "visitfortransimV2interface.inc"
    visitgetdomainnesting = VISIT_INVALID_HANDLE
    end function visitgetdomainnesting

