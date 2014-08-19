!     Nek-MOAB interface notes
!     ------------------------
! Files, repo:
! Most of the Nek-MOAB interface is in the trunk/nek/3rdparty/moab.f source file, which is
! called into from nek's rdparam, readat, and usrchk subroutines.  From 37k ft, nek's datastructure
! stores fields on 3d (nxnxn) arrays on elements, corresponding to a (n-1)th order Guass-Lobatto
! quadrature.  Data is duplicated in those arrays between neighboring elements.

! Global-local ids, etc:
! - in Nek, vertices not really identified as entities, just by coordinates; fields and
!   vertex coordinates are stored in 3d arrays, in lexicographical order, on elements
! - corner vertex positions also held, in separate arrays, in fe order
! - each proc stores global to local processor map (gllnid) and element number (gllel) arrays
!   for element to processor/local id location
! - 2 possible materials (fluid, solid), with fluid elements before solid elements locally
! - side boundary conditions stored in 3d array (elem, side#=6, field), where field is temp or vel

! Flow:
! - rdparam (connect2.f): read ifmoab, ifcoup, ifvcoup
! - readat (connect2.f):
!   . read h5mfle, block/sideset info
!   . bcast h5mfle, block/sideset info
!   . call moab_dat
!     - call nekMOAB_load:
!       . instantiate iMesh, iMeshP, get root set
!       . load file
!       . get tag handles for material, neumann sets
!     - call nekMOAB_get_elems:
!       . foreach block (fluid + solid), init elem iter, elem counts, starts
!       . check against hard-set sizes from SIZE
!       . foreach block (fluid + solid), call nekMOAB_gllnid:
!         - set gllnid(e), global id to proc index
!     - call chk_nel(connect2.f): check local sizes against global hard-set sizes from SIZE
!     - call nekMOAB_create_tags: create results fields tags, set SEM_DIMS on root set
!     - call mapelpr(map2.f):
!       . check #procs against nelgt, exit if too few elems
!       . call set_proc_map(map2.f)-get_map(map2.f)-get_vert(navier8.f)-nekMOAB_loadConn:
!         - get connectivity in vertex(8,elem) in lex order in terms of global ids
!       . compute, distribute gllel (global to local elem map)
!     - call moab_geometry(xm1,ym1,zm1):
!       . call nekMOAB_loadCoord: for each block (fluid+solid):
!         - get connectivity pointer
!         - get vertex coordinates in blocked-storage array (xxx... yyy... zzz...)
!       . call permute_and_map for each coord array: permute coords, then compute spectral pts
!         using high-order (hex27) vertices
!     - call xml2xc: set corner vertex positions (xc, yc, zc), flip corners from lex to fe ordering

!     - call ifill: fill bc array moabbc with -1s
!     - call nekMOAB_BC(moabbc):
!       . for each sideset:
!         - get sideset id, field #
!         - call nekMOAB_intBC(moabbc->bcdata):
!           . get all faces in set recursively; foreach face:
!             - get connectivity of face, hexes adj to face; foreach hex:
!               . get hex connectivity, side # of face wrt hex
!               . call nekMOAB_getElNo: get local hex #, through global ids, gllel
!               . set bcdata(side_no, elem_no, field) to sideset id
!         - end (nekMOAB_intBC)
!     - end (nekMOAB_BC)
!   . end (moab_dat)
! - setlog (bdry.f): print values of ifmoab, ifcoup, ifvcoup
! (solve)
! - usrchk (zero.usr, xxx.usr):
!   . call usr_moab_output (zero.usr, xxx.usr):
!     - call nekMOAB_export_vars (moab.f): write fields to MOAB tags
!     - if an io step, call iMesh_save


    #define IMESH_NULL 0
    #define IMESH_ASSERT \
    if (ierr /= 0) call imesh_err(ierr,imeshh,'moab.f ',__LINE__)

!-----------------------------------------------------------------------
    subroutine nekMOAB_init(comm_f, imesh_instance, partn_handle, &
    fileset_handle, ierr)

    implicit none
    #include "NEKMOAB"
    #include "mpif.h"
    IBASE_HANDLE_T imesh_instance, comm_c
    iBase_EntitySetHandle fileset_handle, partn_handle
    integer :: comm_f, ierr, comm_sz

    if (imesh_instance == 0) then
    !     !Initialize imesh and load file
        imeshh = IMESH_NULL
        call iMesh_newMesh(" ", imeshh, ierr)
        IMESH_ASSERT
        iCreatedImesh = 1
    else
        imeshh = imesh_instance
        iCreatedImesh = 0
    endif

    #ifdef MPI
    call MPI_Comm_size(comm_f, comm_sz, ierr)
    if (ierr /= MPI_SUCCESS) return

    if (comm_sz > 1 .AND. partn_handle == 0) then
        call moab_comm_f2c(comm_f, comm_c)
        call iMeshP_createPartitionAll(%VAL(imeshh), &
        %VAL(comm_c), hPartn, ierr)
        IMESH_ASSERT
        iCreatedPartn = 1
    else
        hPartn = partn_handle
        iCreatedPartn = 0
    endif
    #endif

    if (fileset_handle /= 0) fileset = fileset_handle

    return
    end subroutine nekMOAB_init

!-----------------------------------------------------------------------

    subroutine print_tag_values(tagh, nvals, is_v)
    implicit none
    #include "NEKMOAB"

    iBase_TagHandle tagh
    iBase_EntityHandle elems, verts, adjs
    integer :: elem_size, vert_size, adj_size, ierr, tagv_size, vind
    integer :: offset, offset_size, nvals, maxprint, is_v
    real :: tag_vals
    IBASE_HANDLE_T tagv_ptr, verts_ptr, elems_ptr, offset_ptr, adj_ptr
    iBase_EntitySetHandle  tmp_set
    character*(80) tag_name

    pointer (tagv_ptr, tag_vals(1))
    pointer (adj_ptr, adjs(1))
    pointer (verts_ptr, verts(1))
    pointer (elems_ptr, elems(1))
    pointer (offset_ptr, offset(1))

    call iMesh_getTagName(%VAL(imeshh), %VAL(tagh), tag_name, ierr)
    IMESH_ASSERT

    elem_size = 0
    call iMesh_getEntitiesRec(%VAL(imeshh), %VAL(fileset), &
    %VAL(iBase_REGION), %VAL(iMesh_ALL_TOPOLOGIES), %VAL(1), &
    elems_ptr, elem_size, elem_size, ierr)
    IMESH_ASSERT

    if (is_v == 1) then

    ! Put all vertices bounding elements into a set, as a means to
    ! get unique list
        call iMesh_createEntSet(%VAL(imeshh), %VAL(0), &
        tmp_set, ierr)
        IMESH_ASSERT

    !     moab->get_adjacencies(srcelms,0,false,src_verts,Interface::UNION)
        offset_size = 0
        adj_size = 0
        call iMesh_getEntArrAdj(%VAL(imeshh), elems(1), &
        %VAL(elem_size), %VAL(iBase_VERTEX), adj_ptr, &
        adj_size, adj_size, offset_ptr, &
        offset_size, offset_size, ierr)
        IMESH_ASSERT

        call iMesh_addEntArrToSet(%VAL(imeshh), &
        adjs(1), %VAL(adj_size), %VAL(tmp_set), ierr)

        vert_size = 0
        call iMesh_getEntitiesRec(%VAL(imeshh), %VAL(tmp_set), &
        %VAL(iBase_VERTEX), %VAL(iMesh_ALL_TOPOLOGIES), &
        %VAL(0), verts_ptr, vert_size, vert_size, ierr)
        IMESH_ASSERT

        tagv_size = 0
        call iMesh_getDblArrData(%VAL(imeshh), &
        verts(1), %VAL(vert_size), %VAL(tagh), tagv_ptr, &
        tagv_size, tagv_size, ierr)
        IMESH_ASSERT
    else

        tagv_size = 0
        call iMesh_getDblArrData(%VAL(imeshh), &
        elems(1), %VAL(elem_size), %VAL(tagh), tagv_ptr, &
        tagv_size, tagv_size, ierr)
        IMESH_ASSERT

    endif

!      write (*,22) , elem_size, vert_size
!   22   format("Number of elements = ", I8, " Number of vertices=", I8)

    maxprint = nvals
    if (is_v == 1) then
        if (nvals > vert_size)  maxprint = vert_size
        write (*,*) tag_name, " -- Printing values", maxprint, &
        ' out of ', vert_size
    else
        if (nvals > elem_size)  maxprint = elem_size
        write (*,*) tag_name, " -- Printing values", maxprint, &
        ' out of ', elem_size
    endif

    do 10, vind = 1, maxprint
        if (is_v == 1) then
            write (*,50) verts(vind), tag_vals(vind)
        else
            write (*,50) elems(vind), tag_vals(vind)
        endif
        50 format(I12, "  ", F15.10)
    10 END DO

    return
    end subroutine print_tag_values

!-----------------------------------------------------------------------

    subroutine nekMOAB_import
    implicit none
    #include "NEKMOAB"
    include 'PARALLEL'
    include 'GEOM'

    integer :: nekcomm, nekgroup, nekreal, nid_, np_
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

    integer :: ierr

    if (imeshh == 0 .OR. hPartn == 0 .OR. &
    fileset == 0) then
        call nekMOAB_init(nekcomm, imeshh, hPartn, fileset, &
        ierr)
        IMESH_ASSERT
    endif
             
    if (fileset == 0) then
        call nekMOAB_load   ! read mesh using MOAB
    endif

    #ifdef MPI
!      partsSize = 0
!      rpParts = IMESH_NULL
!      call iMeshP_getLocalParts(%VAL(imeshh), %VAL(hPartn),
!     $     rpParts, partsSize, partsSize, ierr)
    #endif

    call nekMOAB_create_tags             ! allocate MOAB tags to reflect Nek variables

    call nekMOAB_get_elems              ! read material sets and establish mapping
    call chk_nel

    call mapelpr                        ! create gllel mapping
    call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
    call xml2xc                         ! fill xc,yc,zc

    call nekMOAB_BC             ! read MOAB BCs

!      call nekMOAB_compute_diagnostics

    return
    end subroutine nekMOAB_import
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_instance(imesh_instance, partn, &
    fileset_handle)

    implicit none
    #include "NEKMOAB"
    #include "mpif.h"
    IBASE_HANDLE_T imesh_instance, partn
    iBase_EntitySetHandle fileset_handle

    imesh_instance = imeshh
    fileset_handle = fileset
    partn = hPartn

    return
    end subroutine nekMOAB_get_instance

!-----------------------------------------------------------------------
    subroutine nekMOAB_create_tags

    implicit none
    #include "NEKMOAB"

    integer :: ierr, ntot
    integer :: semdim(3)
    iBase_TagHandle tagh

    ntot = nx1*ny1*nz1

! tags used for initializing model, should already be there
    globalIdTag = 0
    call iMesh_getTagHandle(%VAL(imeshh), &
    "GLOBAL_ID",        & !/*in*/ const char* tag_name,
    globalIdTag,        & !/*out*/ iBase_TagHandle *tag_handle,
    ierr)
    IMESH_ASSERT
          
    matsetTag = 0
    call iMesh_getTagHandle(%VAL(imeshh), &
    "MATERIAL_SET",  & !/*in*/ const char* tag_name,
    matSetTag,  & !/*out*/ iBase_TagHandle *tag_handle,
    ierr)
    IMESH_ASSERT

    neusetTag = 0
    call iMesh_getTagHandle(%VAL(imeshh), &
    "NEUMANN_SET",  & !/*in*/ const char* tag_name,
    neuSetTag,  & !/*out*/ iBase_TagHandle *tag_handle,
    ierr)
    IMESH_ASSERT

! create a tag to store SEM dimensions, and set it on the file set
    semdim(1) = nx1
    semdim(2) = ny1
    semdim(3) = nz1
    call nekMOAB_create_find_tag("SEM_DIMS", &
    " moab:TAG_STORAGE_TYPE=SPARSE", 3, iBase_INTEGER, tagh, &
     .TRUE. , ierr)
    IMESH_ASSERT
    call iMesh_setEntSetData(%VAL(imeshh), %VAL(fileset), %VAL(tagh), &
    semdim, 12, ierr)
    IMESH_ASSERT

! tags for results variables
    call nekMOAB_create_find_tag("SEM_X", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, xm1Tag, .TRUE. , ierr)
    IMESH_ASSERT
    call nekMOAB_create_find_tag("SEM_Y", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, ym1Tag, .TRUE. , ierr)
    IMESH_ASSERT
    call nekMOAB_create_find_tag("SEM_Z", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, zm1Tag, .TRUE. , ierr)
    IMESH_ASSERT

    call nekMOAB_create_find_tag("VX", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, vxTag, .TRUE. , ierr)
    IMESH_ASSERT
    call nekMOAB_create_find_tag("VY", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, vyTag, .TRUE. , ierr)
    IMESH_ASSERT
    call nekMOAB_create_find_tag("VZ", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    ntot, iBase_DOUBLE, vzTag, .TRUE. , ierr)
    IMESH_ASSERT

    tTag = 0
    call nekMOAB_create_find_tag("TEMP", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=600.0", &
    &      1, iBase_DOUBLE, tTag, .TRUE. , ierr)
    IMESH_ASSERT

    pTag = 0
    if (nx2 == nx1 .AND. ny2 == ny1 .AND. nz2 == nz1) then
        call nekMOAB_create_find_tag("PRESS", &
        " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
        &      1, iBase_DOUBLE, pTag, .TRUE. , ierr)
        IMESH_ASSERT
    endif

! may or may not have these tags, depending on coupler state
    dTag = 0
    call nekMOAB_create_find_tag("DENSITY", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", &
    &      1, iBase_DOUBLE, dTag, .TRUE. , ierr)
    IMESH_ASSERT

    powTag = 0
    call nekMOAB_create_find_tag("POWER", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", &
    &      1, iBase_DOUBLE, powTag, .TRUE. , ierr)
    IMESH_ASSERT

!     initialize these tags to zero since their use depends on user input
    vtTag = 0
    vpTag = 0
    vdTag = 0
    vpowTag = 0

    call nekMOAB_create_find_tag("VTEMP", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", &
    &      1, iBase_DOUBLE, vtTag, .TRUE. , ierr)
    IMESH_ASSERT

    if (nx2 == nx1 .AND. ny2 == ny1 .AND. nz2 == nz1) then
        call nekMOAB_create_find_tag("VPRESS", &
        " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", &
        &         1, iBase_DOUBLE, vpTag, .TRUE. , ierr)
        IMESH_ASSERT
    endif

! may or may not have these tags, depending on coupler state
    call nekMOAB_create_find_tag("VDENSITY", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", &
    &      1, iBase_DOUBLE, vdTag, .TRUE. , ierr)
    IMESH_ASSERT

    call nekMOAB_create_find_tag("VPOWER", &
    " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", &
    &      1, iBase_DOUBLE, vpowTag, .TRUE. , ierr)
    IMESH_ASSERT

    return
    end subroutine nekMOAB_create_tags
!-----------------------------------------------------------------------
    subroutine nekMOAB_create_find_tag(tag_name, tag_opts, tag_size, &
    tag_datatype, tag_handle, create_if_missing, ierr)
! attempt to get a tag handle of a specified name; if not found, create one; if found
! verify tag characteristics
    implicit none
    #include "NEKMOAB"
    include 'mpif.h'

    character*(*) tag_name, tag_opts
    integer :: tag_size, tag_datatype, ierr, dum_int
    logical :: create_if_missing
    iBase_TagHandle tag_handle

    call iMesh_getTagHandle(%VAL(imeshh), tag_name, tag_handle, &
    ierr)
    if (iBase_SUCCESS /= ierr .AND. create_if_missing) then
    ! need to create it
        call iMesh_createTagWithOptions(%VAL(imeshh), tag_name, &
        tag_opts, %VAL(tag_size), %VAL(tag_datatype), &
        tag_handle, ierr)
        IMESH_ASSERT
    else if (create_if_missing) then
    ! verify characteristics: size, datatype
        call iMesh_getTagSizeValues(%VAL(imeshh), %VAL(tag_handle), &
        dum_int, ierr)
        if (ierr /= iBase_SUCCESS .OR. dum_int /= tag_size) then
            ierr = iBase_INVALID_TAG_HANDLE
            return
        endif
        call iMesh_getTagType(%VAL(imeshh), %VAL(tag_handle), &
        dum_int, ierr)
        if (ierr /= iBase_SUCCESS .OR. dum_int /= tag_datatype) then
            ierr = iBase_INVALID_TAG_HANDLE
            return
        endif
    endif

    return
    end subroutine nekMOAB_create_find_tag

!-----------------------------------------------------------------------
    subroutine nekMOAB_load

!     Load "filename" into imesh/moab, store imesh handle
!     in /nekmoab/ common block

    implicit none
    #include "NEKMOAB"
    include 'mpif.h'
! two forms of load options, depending whether we\'re running serial or parallel
    character*(*) parLoadOpt, serLoadOpt
    parameter(parLoadOpt=" moab:PARALLEL=READ_PART   moab:PARTITION=PA &
    RALLEL_PARTITION moab:PARALLEL_RESOLVE_SHARED_ENTS moab:PARTITION_ &
    DISTRIBUTE moab:CPUTIME")
    parameter(serLoadOpt = " ")
    integer :: nekcomm, nekgroup, nekreal, nid_, np_
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

    integer :: ierr
    IBASE_HANDLE_T ccomm

!     create a file set to load into
    call iMesh_createEntSet(%VAL(imeshh), %VAL(0), fileset, ierr)
    IMESH_ASSERT

    #ifdef MPI
    if (np_ > 1) then
        call iMeshP_loadAll(%VAL(imeshh), %VAL(hPartn),%VAL(fileset), &
        H5MFLE, parLoadOpt, ierr)
        IMESH_ASSERT

    else
        #endif
        call iMesh_load(%VAL(imeshh), %VAL(fileset), &
        H5MFLE, serLoadOpt, ierr)
        IMESH_ASSERT
        #ifdef MPI
    endif
    #endif

    return
    end subroutine nekMOAB_load
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_elems()

! get fluid/solid elements and establish global id mapping

    implicit none
    #include "NEKMOAB"
    include 'PARALLEL'
    include 'SCRCT'
    include 'ZPER'

    include 'mpif.h'

    integer :: iglsum, i, p, tmp, ierr
    iBase_EntitySetHandle dumsets(numsts)
    IBASE_HANDLE_T valptr, setsptr
    integer :: dumval, dumalloc, dumsize, dumnum, ns, ilast, atend, &
    count, tmpcount, ietmp1(numsts), ietmp2(numsts), npass, &
    ipass, m, k, iwork, dumval2
    common /ctmp0/ iwork(lelt)
          
    integer :: nekcomm, nekgroup, nekreal, nid_, np_, j
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

    call nekMOAB_get_set_infos(imeshh, fileset, nekcomm, &
    matsetTag, numflu, numoth, numsts, iestart, iecount, ieiter, &
    matsets, matids, nelgv, nelgt)

! set material types
    tmp = 1
    do i = 1, numflu+numoth
        do j = iestart(i), iestart(i)+iecount(i)-1
            IMATIE(j) = matindx(i)
        enddo
    enddo

    call nekMOAB_get_gids(imeshh, nekcomm, nid_, globalIdTag, numflu, &
    numoth, nelgv, nelgt, iestart, iecount, ieiter, lglel, gllel, &
    gllnid)

    return
    end subroutine nekMOAB_get_elems
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_gids(imeshh, nekcomm, nid_, globalIdTag, &
    numflu, numoth, nelgv, nelgt, iestart, iecount, ieiter, &
    lglel, gllel, gllnid)

    implicit none

    #ifdef PTRSIZE8
    #define POINTER_SIZE 8
    #else
    #define POINTER_SIZE 4
    #endif
    #ifdef MPI
    #include "iMeshP_f.h"
    #include "mpif.h"
    #else
    #include "iMesh_f.h"
    #endif

    iMesh_Instance imeshh
    integer :: numflu, numoth, nid_, iestart(*), iecount(*), nekcomm, &
    lglel(*), gllel(*), gllnid(*), numsts, nelgv, nelgt
    iBase_EntityArrIterator ieiter(*)
    iBase_TagHandle globalIdTag

    include 'SIZE'

    common /ctmp0/ iwork(lelt)
    integer :: ietmp1(2), ietmp2(2), ierr, i, ipass, k, m, iwork, atend, &
    count, mcount, tmpcount, npass, nekid, gid
    iBase_TagHandle nekGidTag
    IBASE_HANDLE_T tag_ptr
    pointer(tag_ptr, nekid(1))

! check whether we have nek gid tag yet
    nekGidTag = 0
    call iMesh_getTagHandle(%VAL(imeshh), "NEK_GID", nekGidTag, ierr)

    if (ierr == iBase_SUCCESS) then
    ! initialize global numbers
        call izero(lglel, nelgt)
        call izero(gllel, nelgt)
        call izero(gllnid, nelgt)
                 
    ! - foreach fluid/solid set
        mcount = 0
        do i = 1, numflu+numoth
            atend = 0
            count = 0
            if (iecount(i) /= 0) then
                call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
                ierr)
                IMESH_ASSERT
            else
                atend = 1
            endif

            do while (atend == 0)
            !   . get ptr to nekgid tag storage
                call iMesh_tagIterate(%VAL(imeshh), %VAL(nekGidTag), &
                %VAL(ieiter(i)), tag_ptr, tmpcount, ierr)
            !   . foreach elem in set
                do k = 1, tmpcount
                !     - check that gid is <= nelt and, if fluid, is <= nelv
                    gid = nekid(k)
                    if (gid > nelgt .OR. &
                    (i < numflu .AND. gid > nelgv)) then
                        call exitti('Global id is too large for elem ', &
                        count+k)
                    endif
                !     - set lglel, gllel, gllnid
                    lglel(mcount+count+k) = gid
                    gllel(gid) = mcount+count + k
                    gllnid(gid) = nid_
                enddo

                call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
                %VAL(tmpcount), atend, ierr)
                IMESH_ASSERT
                count = count + tmpcount
                mcount = mcount + tmpcount
            enddo
        enddo
    else
        ietmp1(1) = 0
        do i = 1, numflu
            ietmp1(1) = ietmp1(1) + iecount(i)
        enddo
        ietmp1(2) = 0
        do i = numflu+1, numflu+numoth
            ietmp1(2) = ietmp1(2) + iecount(i)
        enddo
        #ifdef MPI
        call mpi_scan(ietmp1, ietmp2, 2, MPI_INTEGER, &
        MPI_SUM, nekcomm, ierr)
        if (ierr /= MPI_SUCCESS) ierr = iBase_FAILURE
        IMESH_ASSERT
        #else
        ietmp2(1) = ietmp1(1)
        ietmp2(2) = ietmp1(2)
        #endif
    !     set returned nums to exclusive start - 1
        ietmp2(1) = ietmp2(1) - ietmp1(1)
        ietmp2(2) = ietmp2(2) - ietmp1(2)
    !     set gids for local fluid, other elems
        call izero(lglel, nelgt)
        call izero(gllel, nelgt)
        call izero(gllnid, nelgt)
        do i = 1, nelv
            lglel(i) = ietmp2(1)+i
            gllel(lglel(i)) = i
            gllnid(lglel(i)) = nid_
        enddo
        do i = nelv+1, nelt
            lglel(i) = nelgv+ietmp2(2)-nelv+i
            gllel(lglel(i)) = i
            gllnid(lglel(i)) = nid_
        enddo

    endif

!     now communicate to other procs (taken from set_proc_map in map2.f)
    npass = 1 + nelgt/lelt
    k=1
    do ipass = 1,npass
        m = nelgt - k + 1
        m = min(m,lelt)
        if (m > 0) call igop(gllel(k),iwork,'+  ',m)
        k = k+m
    enddo

    call igop(GLLNID, iwork, 'M  ', NELGT)

!     set the global id tag on elements

!      print *, 'Local to global'
!      print *, lglel(:)

! create the nekgidtag if it wasn\'t there before
    if (0 == nekGidTag) then
        call nekMOAB_create_find_tag("NEK_GID", &
        " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1", &
        &         1, iBase_INTEGER, nekGidTag, .TRUE. , ierr)
        IMESH_ASSERT
    endif

    do i = 1, numflu+numoth
        atend = 0
        count = 0
        if (iecount(i) /= 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            ierr)
            IMESH_ASSERT
        else
            atend = 1
        endif

        do while (atend == 0)
        !            print *, 'Estart:', iestart(i),
        !     $           'local-to-global:',lglel(iestart(i)+count)

        ! use the same iterator for all variables, since the elems are the same
            call nekMOAB_set_int_tag(ieiter(i), globalIdTag, 1, &
            tmpcount, lglel(iestart(i)+count))
            call nekMOAB_set_int_tag(ieiter(i), nekGidTag, 1, &
            tmpcount, lglel(iestart(i)+count))
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT
            count = count + tmpcount
        enddo
    enddo

    return
    end subroutine nekMOAB_get_gids
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_set_infos(imeshh, fileset, nekcomm, &
    matsetTag, numflu, numoth, numsts, iestart, iecount, ieiter, &
    matsets, matids, nelgv, nelgt)

    implicit none

    #ifdef PTRSIZE8
    #define POINTER_SIZE 8
    #else
    #define POINTER_SIZE 4
    #endif
    #ifdef MPI
    #include "iMeshP_f.h"
    #else
    #include "iMesh_f.h"
    #endif

    integer :: iestart(*), iecount(*), nelgv, nelgt, numsts, &
    numflu, numoth, nekcomm, matids(*)
    iBase_EntityArrIterator ieiter(*)
    iBase_TagHandle matsetTag
    iBase_EntitySetHandle fileset, matsets(*)
    iMesh_Instance imeshh

    include 'SIZE'

    integer :: dumval, dumsize, tmp, j, dumalloc, dumnum, i, ierr, ilast
    integer :: iglsum
    iBase_EntitySetHandle dumsets(numsts)
    IBASE_HANDLE_T valptr, setsptr

! get fluid, other material sets, and count elements in them
    valptr = loc(dumval)
    setsptr = loc(dumsets(1))
    dumalloc = numsts
    ilast = 0
    tmp = 1
    do i = 1, numflu+numoth
        dumval = matids(i)
        dumsize = numsts
    ! get the set by matset number
        call iMesh_getEntSetsByTagsRec(%VAL(imeshh), %VAL(fileset), &
        matsetTag, valptr, %VAL(1), %VAL(1), &
        setsptr, dumsize, dumsize, ierr)
        if (dumsize > 1) then
            call exitti('More than one material set with id ', dumval)
        endif
        IMESH_ASSERT
    ! get the number of hexes
        if (dumsize > 0) then
            call iMesh_getNumOfTopoRec(%VAL(imeshh), %VAL(dumsets(1)), &
            %VAL(iMesh_HEXAHEDRON), %VAL(1), dumnum, ierr)
            IMESH_ASSERT
            matsets(i) = dumsets(1)
            iestart(i) = ilast + 1
            ilast = ilast + dumnum
            iecount(i) = dumnum
        !     get an iterator for this set, used later
            call iMesh_initEntArrIterRec(%VAL(imeshh), %VAL(dumsets(1)), &
            %VAL(iBase_REGION), %VAL(iMesh_HEXAHEDRON), &
            %VAL(dumnum), %VAL(0), %VAL(1), ieiter(i), ierr)
        else
            matsets(i) = 0
            iestart(i) = ilast + 1
            iecount(i) = 0
        endif

        tmp = ilast+1
        if (i == numflu) then
            nelv = ilast
        endif
    ! this is if, not elseif, to handle numoth=0
        if (i == numflu+numoth) then
            nelt = ilast
        endif
                    
    enddo

! set remaining values to default values
    do i = numflu+numoth+1, numsts
        iecount(i) = -1
        iestart(i) = -1
        ieiter(i) = 0
        matsets(i) = 0
    enddo

! check local size
    if(nelt <= 0 .OR. nelv > lelv) then
        print *, 'ABORT: nelv is invalid in nekmoab_proc_map'
        print *, 'nelv, lelv = ', nelv, lelv
        call exitt
    endif

! reduce to get global numbers of fluid, other elements, and check size
    nelgv = iglsum(nelv,1)
    nelgt = iglsum(nelt,1)
    if (NELGT > LELG) then
        print *, 'ABORT; increase lelg ',nelgv,lelg
        call exitt
    endif

    return
    end subroutine nekMOAB_get_set_infos
!-----------------------------------------------------------------------
    subroutine nekMOAB_gllnid(matset, iter, count)
    implicit none

! initialize global ids for hexes in this set
    iBase_EntitySetHandle matset
    iBase_EntityArrIterator iter
    integer :: itmp, i, j, ierr, gid, count, atend
    IBASE_HANDLE_T tag_ptr
    pointer(tag_ptr, gid(1))

    #include "NEKMOAB"
    include 'PARALLEL'

    call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(iter), ierr)
    IMESH_ASSERT
    i = 0
    atend = 0
    do while (atend == 0)
        call iMesh_tagIterate(%VAL(imeshh), %VAL(globalIdTag), &
        %VAL(iter), tag_ptr, itmp, ierr)
        IMESH_ASSERT

    ! set the global ids to this proc
        do j = 1, itmp
            if (gid(j) > nelgt) &
            call exitti('Global id greater than NELGT', gid(i))
            gllnid(gid(j)) = nid
        enddo
    ! step the iterator
        call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(iter), %VAL(itmp), &
        atend, ierr)
        IMESH_ASSERT
        i = i + itmp
    enddo

! assert we got the right number of elements
    if (i < count) call exitti( &
    'Wrong number of entities in region iterator', i)

    return
    end subroutine nekMOAB_gllnid
!-----------------------------------------------------------------------
    subroutine nekMOAB_loadConn(vertex, nelgt, ncrnr)

!     Fill the vertex array with connectivity data from imesh

!     vertex(ncrnr, nelt): int array, global id of vertices in element nelt
!     nelgt: int, global number of elements
!     ncrnr: int, number of corner vertices per element (should be 8)

    implicit none
    #include "NEKMOAB"

    integer :: vertex(ncrnr, *), i


! get corner vertex gids
    integer :: e_in_set, eid, j, k, nv, ierr, e_in_chunk, v_per_e
    integer :: gids(27)
    iBase_EntityArrIterator iter
    IBASE_HANDLE_T connect_ptr
    iBase_EntityHandle connect
    pointer (connect_ptr, connect(0:1))

    integer :: l2c(8)
    save    l2c
    data    l2c / 1, 2, 4, 3, 5, 6, 8, 7 /

    do i = 1, numflu+numoth
        nv = 8
        e_in_set = 0
        eid = iestart(i)
        do while (e_in_set < iecount(i))
            if (e_in_set == 0) then
                call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
                ierr)
                IMESH_ASSERT
            endif

        !     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), &
            connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

        !     for each element
            do j = 0, e_in_chunk-1
            !     get vertex gids for this e
                call iMesh_getIntArrData(%VAL(imeshh),  & !iMesh_Instance instance,
                connect(j*v_per_e), %VAL(8), %VAL(globalIdTag), &
                loc(gids), nv, nv, ierr)
                IMESH_ASSERT
            !     permute into vertex array
                do k=1, 8
                    vertex(k, eid) = gids(l2c(k))
                enddo
                eid = eid + 1
            enddo

            e_in_set = e_in_set + e_in_chunk
        enddo
    enddo

    return
    end subroutine nekMOAB_loadConn
!-----------------------------------------------------------------------------
    subroutine nekMOAB_loadCoord(x27, y27, z27)

!     stuff the xyz coords of the 27 verts of each local element --
!     shared vertex coords are stored redundantly

    implicit none
    #include "NEKMOAB"
    real :: x27(27,*), y27(27,*), z27(27,*)
    IBASE_HANDLE_T connect_ptr
    iBase_EntityHandle connect
    pointer(connect_ptr, connect(0:1))
    integer :: i, j, k, ierr, e_in_chunk, e_in_set, v_per_e, e_tot

    e_tot = 0
    do i = 1, numflu+numoth
        e_in_set = 0
        do while (e_in_set < iecount(i))
            if (e_in_set == 0) then
                call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
                ierr)
                IMESH_ASSERT
            endif

        !     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), &
            connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

        !     for each element
            do j = 0, e_in_chunk-1
            !     get vertex gids for this e
                do k = 1, TWENTYSEVEN
                    call iMesh_getVtxCoord(%VAL(imeshh), &
                    %VAL(connect(j*v_per_e+k-1)), &
                    x27(k,e_tot+j+1), y27(k,e_tot+j+1), &
                    z27(k,e_tot+j+1), ierr)
                    IMESH_ASSERT
                enddo
            enddo
            e_tot = e_tot + e_in_chunk
            e_in_set = e_in_set + e_in_chunk
        enddo
    enddo

    return
    end subroutine nekMOAB_loadCoord
!-----------------------------------------------------------------------------
    subroutine nekMOAB_BC()

!     Mark the "side" (hex/face) with the index of the corresponding set in ibcsts array


    implicit none
    #include "NEKMOAB"
    common /mbc/ moabbc(6,lelt,ldimt1)
    integer :: moabbc

    IBASE_HANDLE_T hentSet(*)
    pointer (rpentSet, hentSet)
    integer :: entSetSize

    integer :: ierr, i, j, set_idx, set_ids(numsts)

    call ifill(moabbc, -1, 6*lelt*ldimt1)

! idesets in cubit come in as entity sets with the NEUMANN_SET -- see sample file
    call iMesh_getTagHandle(%VAL(imeshh), &
    "NEUMANN_SET", neuSetTag, ierr)
    IMESH_ASSERT

    rpentSet = IMESH_NULL
    entSetSize      = 0
    call iMesh_getEntSetsByTagsRec(%VAL(imeshh), &
    %VAL(fileset), neuSetTag, %VAL(IMESH_NULL), %VAL(1), %VAL(0), &
    rpentSet, entSetSize, entSetSize, ierr)
    IMESH_ASSERT

!     get the set ids
    if (entSetSize > numsts) then
    !     too many sets, need to bail
        ierr = iBase_FAILURE
        IMESH_ASSERT
    endif
    do i = 1, entSetSize
        call iMesh_getEntSetIntData(%VAL(imeshh), &
        %VAL(hentSet(i)), %VAL(neuSetTag), set_ids(i), ierr)
        IMESH_ASSERT
    enddo

!     loop over local neusets, will be <= numsts
    do i=1, entSetSize
    !     find the right set
        set_idx = -1
        do j = 1, numsts
            if (set_ids(i) == ibcsts(j)) then
                set_idx = j
                if (set_idx /= -1) then
                    call nekMOAB_intBC(moabbc, hentSet(i), set_ids(i), &
                    bcf(set_idx))
                endif
            endif
        enddo
    enddo
          
    call iMesh_freeMemory(%VAL(imeshh), rpentSet)

    return
    end subroutine nekMOAB_BC
!-----------------------------------------------------------------------
    subroutine nekMOAB_intBC(bcdata, setHandle, setId, field)

!     Internal function, don\'t call directly


    implicit none
    #include "NEKMOAB"

    integer :: bcdata(6,lelt,ldimt1), field

    iBase_EntitySetHandle setHandle

    IBASE_HANDLE_T hset
    integer :: setId !coming from cubit

    IBASE_HANDLE_T faces(*)
    pointer (rpfaces, faces)
    integer :: facesSize

    IBASE_HANDLE_T ahex(*)
    pointer (rpahex, ahex)
    integer :: ahexSize

    IBASE_HANDLE_T fvtx(*)
    pointer (rpfvtx, fvtx)
    integer :: fvtxSize

    IBASE_HANDLE_T hvtx(*)
    pointer (rphvtx, hvtx)
    integer :: hvtxSize

    integer :: ierr, i, j, elno

    integer :: side_no, side_offset, side_sense
    integer :: hexMBCNType
    integer :: num_sides
    integer :: elnos(2)
    integer :: numids

    numids = 0

    call iMesh_MBCNType(%VAL(iMesh_HEXAHEDRON), hexMBCNType)

! hat faces are in this set?
    facesSize = 0
    rpfaces = IMESH_NULL
    call iMesh_getEntitiesRec(%VAL(imeshh), &
    %VAL(setHandle), %VAL(iBase_FACE), %VAL(iMesh_QUADRILATERAL), &
    %VAL(1), rpfaces, facesSize, facesSize, ierr)
    IMESH_ASSERT

    num_sides = 0

    do i=1, facesSize
    ! et vertices defining the face
        fvtxSize = 0
        rpfvtx = IMESH_NULL
        call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), &
        %VAL(iBase_VERTEX), rpfvtx, fvtxSize, fvtxSize, ierr)
        IMESH_ASSERT

    ! et hexes adjacent to the face (1 or 2)
    ! -- hopefully only returns hexes on local proc, but untested
        ahexSize = 0
        rpahex = IMESH_NULL
        call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), &
        %VAL(iBase_REGION), rpahex, ahexSize, ahexSize, ierr)
        IMESH_ASSERT

        do j=1, ahexSize
        ! et verts adjacent to the hex
            hvtxSize = 0
            rphvtx = IMESH_NULL
            call iMesh_getEntAdj(%VAL(imeshh), &
            %VAL(ahex(j)), %VAL(iBase_VERTEX), &
            rphvtx, hvtxSize, hvtxSize, ierr)
            IMESH_ASSERT

        ! et the side number
            call MBCN_SideNumberUlong(%VAL(rphvtx), &
            %VAL(hexMBCNType), %VAL(rpfvtx), &
            %VAL(4), %VAL(2),      & !4 vertices, want 2D side number
            side_no, side_sense, side_offset)
            if ((side_no < 0) .OR. (side_no > 5)) ierr = 1
            IMESH_ASSERT

            side_no = side_no + 1 !moab side number is zero based
            num_sides = num_sides + 1

        ! all nekMOAB_getGlobElNo(ahex(j), elno)
        ! rint *, 'INFO: ', elno, side_no, setId

            call nekMOAB_getElNo(ahex(j), elno)
            if (ahexSize == 2) elnos(j) = elno

            if (bcdata(side_no, elno, field) /= -1) &
            print *, 'Warning: resetting BC, bcno, elno, sideno = ', &
            setId, elno, side_no
            bcdata(side_no, elno, field) = setId
            numids = numids + 1

            call iMesh_freeMemory(%VAL(imeshh), rphvtx)
        enddo

    !		 Check if we share a side and warn iff not solving a
    !		 conjugate-heat transfer problem i.e., only if elements of
    !		 non-fluid type are present.
        if (ahexSize == 2 .AND. numoth == 0 ) &
        print *, 'Warning: face shared by 2 hexes: ', elnos(1), &
        elnos(2)

        call iMesh_freeMemory(%VAL(imeshh), rpahex)
        call iMesh_freeMemory(%VAL(imeshh), rpfvtx)

    enddo

    call iMesh_freeMemory(%VAL(imeshh), rpfaces)

!     Should probably print the following only if verbose output needed
    print *, 'Setid, numids = ', setId, numids

    return
    end subroutine nekMOAB_intBC
!-----------------------------------------------------------------------
    subroutine nekMOAB_getElNo(handle, elno)

!     Given the imesh handle of an element (hex),
!     return its local nek element number in elno
!     Cannot be used until the GLLEL array has been set (in map2.f)

    implicit none
    #include "NEKMOAB"
    include 'PARALLEL'

    IBASE_HANDLE_T handle
    integer :: elno

! irst get global id
    call nekMOAB_getGlobElNo(handle, elno)
! hen convert to local id
    elno = GLLEL(elno)

    end subroutine nekMOAB_getElNo
!-----------------------------------------------------------------------
    subroutine nekMOAB_getGlobElNo(handle, elno)

    implicit none
    #include "NEKMOAB"
    include 'PARALLEL'

    IBASE_HANDLE_T handle
    integer :: elno
          
    integer :: ierr

    call iMesh_getIntData(%VAL(imeshh), %VAL(handle), &
    %VAL(globalIdTag), elno, ierr)
    IMESH_ASSERT
          
    if (GLLNID(elno) /= NID) then
        print *, 'Wrong proc for element; gid, proc, gllnid = ', &
        elno, nid, gllnid(elno)
        call exitt
    endif
                              
    end subroutine nekMOAB_getGlobElNo
!-----------------------------------------------------------------------
    subroutine nekMOAB_bcs() ! fill the nek cbc arrays
    implicit none

    #include "NEKMOAB"
    common /mbc/ moabbc(6,lelt,ldimt1)
    integer :: moabbc

    integer :: e,f,l, j
    character(3) :: cbi

    integer :: ibcs(3), i, lcbc, nface
    data ibcs / 0, 0, 0 /

    lcbc=18*lelt*(ldimt1 + 1)
    call blank(cbc,lcbc)

    nface = 2*ndim
    do l=1,nfield
        do e=1,nelt
            do f=1,nface
                if (moabbc(f,e,l) /= -1) then
                !     moabbc stores local set index, retrieve the character bc type
                    do j = 1, numsts
                        if (moabbc(f,e,l) == ibcsts(j) &
                         .AND. bcf(j) == l) then
                            cbc(f, e, l) = bctyps(j)
                        endif
                    enddo
                else
                    cbc(f, e, l) = 'E  '
                endif
            enddo
        enddo
    enddo
    return
    end subroutine nekMOAB_bcs
!-----------------------------------------------------------------------
    subroutine moab_geometry (xmlo,ymlo,zmlo)

    implicit none
    #include "NEKMOAB"

    real ::      xmlo(nx1*ny1*nz1,1) &
    , ymlo(nx1*ny1*nz1,1) &
    , zmlo(nx1*ny1*nz1,1)
    integer ::   e, nmoab

    common /tcrmg/ x27(27,lelt), y27(27,lelt), z27(27,lelt)
    real :: x27, y27, z27

    call nekMOAB_loadCoord(x27, y27, z27) !     Get coords from moab

!     call outmat(x27,27,8,'x27dat',nelt)
!     call outmat(y27,27,8,'y27dat',nelt)
!     call outmat(z27,27,8,'z27dat',nelt)

    nmoab = 3 !not used
    do e=1,nelt   !  Interpolate for each element
        call permute_and_map(xmlo(1,e),x27(1,e),nx1,nmoab,e)
        call permute_and_map(ymlo(1,e),y27(1,e),ny1,nmoab,e)
        if (if3d) call permute_and_map(zmlo(1,e),z27(1,e),nz1,nmoab,e)
    enddo

!     param(66) = 0
!     ifxyo = .true.
!     ifvo  = .true.
!     call outpost(xm1,ym1,zm1,pr,t,'   ')
!     write(6,*) 'DONE PERMUTE; ABORT'
!     call exitt

    return
    end subroutine moab_geometry
!-----------------------------------------------------------------------
    subroutine xml2xc
          
    implicit none
    #include "NEKMOAB"
    include 'GEOM'
    integer :: i, j, k, l

    integer :: e

    l = 0
    if (if3d) then
        do e=1,nelt
            do k=1,nz1,nz1-1
                do j=1,ny1,ny1-1
                    do i=1,nx1,nx1-1
                        l = l+1
                        xc(l,1) = xm1(i,j,k,e)
                        yc(l,1) = ym1(i,j,k,e)
                        zc(l,1) = zm1(i,j,k,e)
                    enddo
                enddo
            enddo
        enddo


    else  ! 2D
        do e=1,nelt
            do j=1,ny1,ny1-1
                do i=1,nx1,nx1-1
                    l = l+1
                    xc(l,1) = xm1(i,j,1,e)
                    yc(l,1) = ym1(i,j,1,e)
                enddo
            enddo
        enddo
    endif

    do e=1,nelt ! flip corners back to pre-proc notation
        dtmp = xc(3,e)
        xc(3,e) = xc(4,e)
        xc(4,e) = dtmp
        dtmp = yc(3,e)
        yc(3,e) = yc(4,e)
        yc(4,e) = dtmp
        dtmp = zc(3,e)
        zc(3,e) = zc(4,e)
        zc(4,e) = dtmp

        if(if3d) then
            dtmp = xc(7,e)
            xc(7,e) = xc(8,e)
            xc(8,e) = dtmp
            dtmp = yc(7,e)
            yc(7,e) = yc(8,e)
            yc(8,e) = dtmp
            dtmp = zc(7,e)
            zc(7,e) = zc(8,e)
            zc(8,e) = dtmp
        endif
    enddo

    return
    end subroutine xml2xc
!-----------------------------------------------------------------------
    subroutine permute_and_map(x,x27,nx,nmoab,e)

    implicit none
    #include "NEKMOAB"
    integer :: nx, nmoab, e

    real :: x(1),x27(0:1)

    common /ctmp0/ z3(3),zpt(lx1) &
    , xt(3,3,3) &
    , wk(3*lx1*ly1*lz1) &
    , xw(3*lx1*ly1*lz1)
    real :: z3, zpt, xt, wk, xw

    real :: interp(lx1*3),interpt(lx1*3),zl(lx1*3)
    save interp,interpt

    integer :: moabmap(27)
    save    moabmap
    data    moabmap &
    /  0,  8,  1, 11, 24,  9,  3, 10,  2 &
    , 12, 20, 13, 23, 26, 21, 15, 22, 14 &
    ,  4, 16,  5, 19, 25, 17,  7, 18,  6  /

    integer :: nmlast,nxlast
    save    nmlast,nxlast
    data    nmlast,nxlast / 0,0 /

    real :: gh_edge(3,3,3),gh_vtx(3,3,3),zgh(3)
    save zgh
    data zgh / -1., 0., 1. /
    integer :: gh_type, i, j, ldw

    if (nx /= nxlast) then ! define interp. op
        nxlast = nx
        call zwgll (zl,interp,nx)
        call igllm (interp,interpt,zgh,zl,3,nx,3,nx)
    endif

    do i=1,3**ndim     ! currently support only 3x3x3 in moab
        j = moabmap(i)
        xt(i,1,1) = x27(j)
    enddo

    gh_type = 1 ! vertex only extension
    gh_type = 2 ! edge extension
    call gh_face_extend(xt,zgh,3,gh_type,gh_edge,gh_vtx)

!     Interpolate from 3x3x3 to (nx1 x ny1 x nz1) SEM mesh
    ldw = 3*lx1*ly1*lz1
    call map_to_crs(x,nx1,xt,3,if3d,wk,ldw)

    return
    end subroutine permute_and_map
!-----------------------------------------------------------------------
    subroutine nekMOAB_export_vars
    implicit none
    #include "NEKMOAB"
    include 'GEOM'
    include 'SOLN'
    include 'PARALLEL'
    integer :: nekcomm, nekgroup, nekreal, nid_, np_
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
          
    integer :: i, j, ierr, ntot, tmpcount, count, atend
    integer :: junk, numelems_iterated
    real :: tag_ptr(1)
    real :: density
    common /moabvp/  density (lx1,ly1,lz1,lelt)             ! fuel density for MOAB
    double precision :: start_t, end_t
    double precision :: wtime

    ntot = nx1*ny1*nz1
    numelems_iterated = 0

!     TODO: Modify/combine the two loops below over fluid/solid
!     materials for element/vertex based tags into one loop

!     do element based tags only
    do i = 1, numflu+numoth
        atend = 0
        count = 0
        if (iecount(i) /= 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            ierr)
            IMESH_ASSERT
        else
            atend = 1
        endif

    !     use the same iterator for all variables, since the elems are the same
        do while (atend == 0)

            if (xm1Tag /= 0) then
                call nekMOAB_set_tag(ieiter(i), xm1Tag, ntot, tmpcount, &
                xm1)
            endif
            if (ym1Tag /= 0) then
                call nekMOAB_set_tag(ieiter(i), ym1Tag, ntot, tmpcount, &
                ym1)
            endif
            if (zm1Tag /= 0) then
                call nekMOAB_set_tag(ieiter(i), zm1Tag, ntot, tmpcount, &
                zm1)
            endif

            if (vxTag /= 0) then
                call nekMOAB_set_tag(ieiter(i), vxTag, ntot, tmpcount, &
                vx)
            endif
            if (vyTag /= 0) then
                call nekMOAB_set_tag(ieiter(i), vyTag, ntot, tmpcount, &
                vy)
            endif
            if (vzTag /= 0) then
                call nekMOAB_set_tag(ieiter(i), vzTag, ntot, tmpcount, &
                vz)
            endif

            if (tTag /= 0) then
                tmpcount = numelems_iterated
                call nekMOAB_set_tag_avg(ieiter(i), tTag, ntot, tmpcount, &
                t)
            endif

            if (pTag /= 0 .AND. nx2 == nx1 .AND. ny2 == ny1 .AND. &
            nz2 == nz1) then
                tmpcount = numelems_iterated
                call nekMOAB_set_tag_avg(ieiter(i), pTag, ntot, tmpcount, &
                pr)
            endif

            if (dTag /= 0) then
                tmpcount = numelems_iterated
                call nekMOAB_set_tag_avg(ieiter(i), dTag, ntot, tmpcount, &
                density)
            !     $           vtrans)
            endif

        !     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT

            count = count + tmpcount
            numelems_iterated = numelems_iterated + tmpcount

        enddo

    !     double-check the total number of elements in this set
        if (count /= iecount(i)) then
            call exitti('Wrong no of elems iterating over matset ', &
            matids(i))
        endif
    enddo

!     now do vertex based tags only
    count = 0
    numelems_iterated = 0
    do i = 1, numflu+numoth
        atend = 0
        if (iecount(i) /= 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            ierr)
            IMESH_ASSERT
        else
            atend = 1
        endif

        do while (atend == 0)

        !     vertex-based temperature
            if (vtTag /= 0) then
                tmpcount = count
                call nekMOAB_set_vertex_tag(ieiter(i), vtTag, tmpcount, &
                t, ierr)
            endif
        !     vertex-based density
            if (vdTag /= 0) then
                tmpcount = count
                call nekMOAB_set_vertex_tag(ieiter(i), vdTag, tmpcount, &
                density, ierr)
            endif
        !     vertex-based pressure, but only if its there
            if (vpTag /= 0 .AND. &
            nx2 == nx1 .AND. ny2 == ny1 .AND. nz2 == nz1) then
                tmpcount = count
                call nekMOAB_set_vertex_tag(ieiter(i), vpTag, &
                tmpcount, pr, ierr)
            endif

        !     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT

            count = count + tmpcount

        enddo

        numelems_iterated = numelems_iterated + iecount(i)

    !     double-check the total number of elements in this set
        if (count > numelems_iterated) then
            call exitti('Wrong no of elems iterating over matset ', &
            matids(i))
        endif
    enddo

!     write (6,*) 'Printing temperature vertex values in tags'
!     call print_tag_values(vtTag, 99, 1)
!     write (6,*) 'Printing density vertex values in tags'
!     call print_tag_values(vdTag, 99, 1)
!     write (6,*) 'Printing temperature element values in tags'
!     call print_tag_values(tTag, 99, 0)

    call mpi_barrier(nekcomm, ierr)

    return
    end subroutine nekMOAB_export_vars
!-----------------------------------------------------------------------
    subroutine nekMOAB_import_tagvalues(tagname, is_v, field)
    implicit none
    #include "NEKMOAB"

    character(*) tagname
    integer :: is_v, ierr
    real :: field(lx1,ly1,lz1,lelt)

    iBase_TagHandle tagh

    call iMesh_getTagHandle(%VAL(imeshh), tagname, tagh, ierr)
          
    call nekMOAB_import_vars(tagh, is_v, field)

    return
    end subroutine nekMOAB_import_tagvalues
!-----------------------------------------------------------------------
    subroutine nekMOAB_import_vars(tagh, is_v, field)
    implicit none
    #include "NEKMOAB"

    real :: field(lx1,ly1,lz1,lelt)

    include 'GEOM'
    include 'SOLN'
    integer :: i, j, ierr, ntot, tmpct, count, atend, is_v
    integer :: numelems_iterated
    iBase_TagHandle tagh

    ntot = nx1*ny1*nz1
    numelems_iterated = 0
    count = 0
    if (tagh == 0) return

!      call print_tag_values(tagh, 5, is_v)

    do i = 1, numflu+numoth

        atend = 0
        if (iecount(i) /= 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            ierr)
            IMESH_ASSERT
        else
            atend = 1
        endif

        do while (atend == 0)
            tmpct = count
        ! use the same iterator for all variables, since the elems are the same
            if (is_v == 1) then
                call nekMOAB_get_vertex_tag(ieiter(i), tagh, tmpct, &
                field, ierr)
            else
                call nekMOAB_get_tag_avg(ieiter(i), tagh, ntot, tmpct, &
                field, ierr)
            endif
            IMESH_ASSERT
        !     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), &
            %VAL(tmpct), atend, ierr)
            IMESH_ASSERT

            count = count + tmpct
        enddo

        numelems_iterated = numelems_iterated + iecount(i)

    ! double-check the total number of elements in this set
        if (count > numelems_iterated) then
            call exitti('Wrong no of elems iterating over matset ', &
            matids(i))
        endif
    enddo

    return
    end subroutine nekMOAB_import_vars
!-----------------------------------------------------------------------
    subroutine nekMOAB_set_tag_avg(iter, tagh, size, count, vals)
    implicit none

    #include "NEKMOAB"
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, j, size, count, tmpcount, offset
    real :: vals(*), tag_vals, avg_vals
    pointer(tag_ptr, tag_vals(1))

    call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), &
    %VAL(iter), tag_ptr, tmpcount, ierr)
    IMESH_ASSERT

! set the tag vals by looping over each element and then updating the average value
    do j = 1, tmpcount
        avg_vals = 0.0
        offset = (count+j-1)*size
        do i = 1, size
            avg_vals = avg_vals + vals(offset+i)
        enddo
        avg_vals = avg_vals/size
        tag_vals(j) = avg_vals
    enddo

    count = tmpcount

    return
    end subroutine nekMOAB_set_tag_avg
!-----------------------------------------------------------------------
    subroutine nekMOAB_set_tag(iter, tagh, size, count, vals)
    implicit none

    #include "NEKMOAB"
    include 'PARALLEL'
          
    integer :: nekcomm, nekgroup, nekreal, nid_, np_
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, ivals, size, count!, tmpcount, offset
    real :: vals(*), tag_vals
    pointer(tag_ptr, tag_vals(1))

    call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), &
    %VAL(iter), tag_ptr, count, ierr)
    IMESH_ASSERT

! set the tag vals
    ivals = size * count
    do i = 1, ivals
        tag_vals(i) = vals(i)
    enddo

    return
    end subroutine nekMOAB_set_tag
!-----------------------------------------------------------------------
    subroutine nekMOAB_set_int_tag(iter, tagh, size, count, vals)
    implicit none

    #include "NEKMOAB"
    include 'PARALLEL'
    integer :: nekcomm, nekgroup, nekreal, nid_, np_
    common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
       
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, ivals, size, count
    integer :: vals(*), tag_vals
    pointer(tag_ptr, tag_vals(1))

    call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), &
    %VAL(iter), tag_ptr, count, ierr)
    IMESH_ASSERT

! set the tag vals
    ivals = size * count
    do i = 1, ivals
        tag_vals(i) = vals(i)
    enddo

    return
    end subroutine nekMOAB_set_int_tag
!-----------------------------------------------------------------------
    subroutine nekMOAB_set_vertex_tag(iter, tagh, count, vals, ierr)
    implicit none

    #include "NEKMOAB"
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, ic, j, ivals, size, count, v_per_e, ntot
    integer :: loccount, nv
    real :: vals(*), tag_vals(27)
    integer :: gids(27)
    real :: avg_vals
    iBase_EntityHandle connect
    pointer (connect_ptr, connect(0:1))

    integer :: l2c(8)
    save    l2c
    data    l2c / 1, 2, 4, 3, 5, 6, 8, 7 /
      
    integer :: jj

    call iMesh_connectIterate(%VAL(imeshh), %VAL(iter), &
    connect_ptr, v_per_e, loccount, ierr)
    IMESH_ASSERT

! only works if nx, ny, nz are equal, and if v_per_e is 27
    if (nx1 /= ny1 .OR. nx1 /= nz1 .OR. v_per_e /= 27) then
        ierr = iBase_FAILURE
        IMESH_ASSERT
    endif

    ntot = nx1 * ny1 * nz1

! set the tag vals
    nv = 8
    do ic = 0, loccount-1
        do jj = 1,v_per_e
            tag_vals(jj) = 0.0
        enddo

    !        write(*,*) '--', ic, '--'
    !     permute into vertex array
        avg_vals = 0
        do j=1, nv
            tag_vals(j) = vals((count+ic)*ntot+l2c(j))
            avg_vals = avg_vals + tag_vals(j)
        !          print *, ((count+ic)*ntot+l2c(j)), tag_vals(j)
        enddo
    !        write(*,*) '--'

        avg_vals = avg_vals/nv
        do j = nv+1, v_per_e
            tag_vals(j) = avg_vals
        enddo

        call iMesh_setDblArrData(%VAL(imeshh), &
        connect(ic*v_per_e), %VAL(v_per_e), %VAL(tagh), &
        tag_vals(1), %VAL(v_per_e), ierr)
        IMESH_ASSERT
    enddo

    count = loccount

    return
    end subroutine nekMOAB_set_vertex_tag
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_tag_avg(iter,tagh,size,count,vals,ierr)
    implicit none

    #include "NEKMOAB"
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, j, size, count, offset, tmpcount
    real :: vals(*), tag_vals
    pointer(tag_ptr, tag_vals(1))

    call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), &
    %VAL(iter), tag_ptr, tmpcount, ierr)
! assert and break if there is a problem
    IMESH_ASSERT

    offset = count*size

! set the tag vals
    do i = 1, tmpcount
        do j = 1, size
            vals(offset+(i-1)*size+j) = tag_vals(i)
        enddo
    enddo

    count = tmpcount

    return
    end subroutine nekMOAB_get_tag_avg
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_tag(iter, tagh, size, count, vals, ierr)
    implicit none

    #include "NEKMOAB"
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, ivals, size, count
    real :: vals(*), tag_vals
    pointer(tag_ptr, tag_vals(1))

    call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), &
    %VAL(iter), tag_ptr, count, ierr)
! assert and break if there is a problem
    IMESH_ASSERT

! set the tag vals
    ivals = size * count
    do i = 1, ivals
        vals(i) = tag_vals(i)
    enddo

    return
    end subroutine nekMOAB_get_tag
!-----------------------------------------------------------------------
    subroutine nekMOAB_get_vertex_tag(iter, tagh, count, vals, ierr)
    implicit none

    #include "NEKMOAB"
    iBase_EntityArrIterator iter
    iBase_TagHandle tagh
    integer :: ierr, i, j, ivals, size, count, v_per_e, ntot
    integer :: loccount, tmpind
    integer :: tagv_size
    real :: vals(lx1*ly1*lz1*lelt), tag_vals
    IBASE_HANDLE_T tagv_ptr
    iBase_EntityHandle connect
    pointer (connect_ptr, connect(1))
    pointer (tagv_ptr, tag_vals(1))
    real :: avg

! since v_per_e is required to be 27, use the MOAB transformation
! map explicitly to set the MOAB tag values appropriately
    integer :: moabmap(27)
    save    moabmap
    data    moabmap &
    /  0,  8,  1, 11, 24,  9,  3, 10,  2 &
    , 12, 20, 13, 23, 26, 21, 15, 22, 14 &
    ,  4, 16,  5, 19, 25, 17,  7, 18,  6  /

    call iMesh_connectIterate(%VAL(imeshh), %VAL(iter), &
    connect_ptr, v_per_e, loccount, ierr)
! don't assert here, just return
    if (iBase_SUCCESS /= ierr) return

! only works if nx1, ny1, nz1 are equal, and if v_per_e is 27
    if (nx1 /= ny1 .OR. nx1 /= nz1 .OR. v_per_e /= 27) then
        ierr = iBase_FAILURE
        return
    endif

    ntot = nx1 * ny1 * nz1

! set the tag vals
    ivals = 1
    do i = 1, loccount
        tmpind = (i+count-1)*ntot
        tagv_size = 0
    !       transfer spectral variable from vertex variable corners only
        call iMesh_getDblArrData(%VAL(imeshh), &
        connect(ivals), %VAL(v_per_e), %VAL(tagh), tagv_ptr, &
        tagv_size, tagv_size, ierr)
    ! don't assert here, just return
        if (iBase_SUCCESS /= ierr) return
        if (tagv_size /= v_per_e) then
            ierr = 5
            return
        endif

    !        Overwrite values based on indices directly
        if (lx1 == 3) then
        ! Use the MOABMAP directly only for lx1=3
            do j = 1, v_per_e
                vals(tmpind+j) = tag_vals(moabmap(j)+1)
            enddo
        else
            avg = 0.0d0
            do j = 1, v_per_e
                avg = avg + tag_vals(j)
            enddo
            avg = avg/v_per_e
        ! set all values to the average computed
            do j = 1, ntot
                vals(tmpind+j) = avg
            enddo

        ! TODO: the optimal and accurate way to perform this
        ! operation is to do a local L_2 projection of the
        ! quadratic solution defined on MOAB mesh to the
        ! NEK's GLL mesh using a proper mass matrix for
        ! the linear transformation.
        ! The hacky way done below is to use injection of
        ! the corner vertices alone since the rest of the
        ! unknowns have been provided the average value.
        ! This will probably not preserve the averages and
        ! should be removed in the future.
        ! For coupled runs, it is recommended to use
        ! lx1=3 in which case, the 'if' block will be
        ! executed.

            vals(tmpind+1) = tag_vals(moabmap(1)+1) ! bottom left
            vals(tmpind+lx1) = tag_vals(moabmap(3)+1) ! bottom right
            vals(tmpind+lx1*ly1)=tag_vals(moabmap(9)+1) ! top right
            vals(tmpind+1+lx1*(ly1-1)) = tag_vals(moabmap(7)+1) ! top left

            vals(tmpind+1+lx1*ly1*(lz1-1)) = tag_vals(moabmap(19)+1)
            vals(tmpind+lx1+lx1*ly1*(lz1-1)) = tag_vals(moabmap(21)+1)
            vals(tmpind+lx1+lx1*(ly1-1)+lx1*ly1*(lz1-1))= &
            tag_vals(moabmap(27)+1)
            vals(tmpind+1+lx1*(ly1-1)+lx1*ly1*(lz1-1)) = &
            tag_vals(moabmap(25)+1)
        endif

                 
    ! update ivals by the necessary offset to get the next element
        ivals = ivals + v_per_e
    enddo
          
    count = loccount

    return
    end subroutine nekMOAB_get_vertex_tag
!-----------------------------------------------------------------------
    subroutine nekMOAB_compute_diagnostics()
    implicit none

    #include "NEKMOAB"
    include 'GEOM'

    integer :: i, j, k, l

    integer :: e, intv
    real :: avg(3)

    intv = nelt / 10
    e = 1
    do while (e < nelt)
        avg(1) = 0.0
        avg(2) = 0.0
        avg(3) = 0.0
        do k=1,nz1
            do j=1,ny1
                do i=1,nx1
                    avg(1) = avg(1) + xm1(i,j,k,e)
                    avg(2) = avg(2) + ym1(i,j,k,e)
                    if (if3d) &
                    avg(3) = avg(3) + zm1(i,j,k,e)
                enddo
            enddo
        enddo
        avg(1) = avg(1) / (nx1*ny1*nz1)
        avg(2) = avg(2) / (nx1*ny1*nz1)
        avg(3) = avg(3) / (nx1*ny1*nz1)
        print *, "Average for e is ", e, avg(1), avg(2), avg(3)

        e = e + intv
    enddo

    return
    end subroutine nekMOAB_compute_diagnostics

!-----------------------------------------------------------------------
    block data nekMOABdata
    #include "MOABCORE"
    data imeshh/0/, hPartn/0/, fileset/0/, &
!	  The following variables have been removed from MOABCORE interface
!     $     rpParts/0/, rpHexes/0/,
!     $     rpxm1/0/, rpym1/0/, rpzm1/0/,
!     $     rpvx/0/, rpvy/0/, rpvz/0/, rpt/0/, rpp/0/,
    globalIdTag/0/, matsetTag/0/, neusetTag/0/, &
    matsets/numsts*0/, ieiter/numsts*0/, &
    xm1Tag/0/, ym1Tag/0/, zm1Tag/0/, vxTag/0/, vyTag/0/, &
    vzTag/0/, tTag/0/, &
    pTag/0/, dTag/0/, powTag/0/, vtTag/0/, vpTag/0/, vdTag/0/, &
    vpowTag/0/, senseTag/0/, &
    iCreatedImesh/0/, iCreatedPartn/0/, iCreatedFileset/0/, &
    iestart/numsts*0/, iecount/numsts*0/
!	  The following variables have been removed from MOABCORE interface
!     $     partsSize/0/, hexesSize/0/
    END PROGRAM

