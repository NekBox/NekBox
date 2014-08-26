!===============================================================================
!     pff@cfm.brown.edu   3/19/96


!     This  is a suite of routines for solving overlapping subdomain
!     problems with finite elements on distributed memory machines.

!     The overall hierarchy of data structures is as follows:

!         System        -  index set denoted by       _glob

!            Processor  -  index set denoted by       _pglob

!              .Domain  -  index set denoted by       _dloc  (1,2,3,...,n_dloc)

!              .Sp.Elem -  index set denoted by       _sloc  (1,2,3,...,n_sloc)


!     A critical component of the parallel DD solver is the gather-scatter
!     operation.   As there may be more than one domain on a processor, and
!     communication must be minimized, it is critical that communication be
!     processor oriented, not domain oriented.  Hence domains will access data
!     via the dloc_to_pglob interface, and the pglob indexed variables will
!     be accessed via a gather-scatter which interacts via the pglob_glob
!     interface.   Noticed that, in a uni-processor application, the pglob_glob
!     map will be simply the identity.

!===============================================================================

    subroutine set_overlap

!     Set up arrays for overlapping Schwartz algorithm *for pressure solver*

    use size_m
    use domain
    use esolv
    use input
    use tstep

    REAL*8 :: dnekclock,t0

    parameter (          n_tri = 7*ltotd )
    common /scrns/  tri (n_tri)
    integer ::         tri,elem

    common /screv/ x(2*ltotd)
    common /scrvh/ y(2*ltotd)
    common /scrch/ z(2*ltotd)

    common /ctmp0/ nv_to_t(2*ltotd)

    parameter (lia = ltotd - 2 - 2*lelt)
    common /scrcg/ ntri(lelt+1),nmask(lelt+1) &
    , ia(lia)

    common /scruz/ color   (4*ltotd)
    common /scrmg/ ddmask  (4*ltotd)
    common /ctmp1/ mask    (4*ltotd)

    parameter(lxx=lx1*lx1, levb=lelv+lbelv)
    common /fastd/  df(lx1*ly1*lz1,levb) &
    ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)


    integer :: e

    if (lx1 == 2) param(43)=1.
    if (lx1 == 2 .AND. nid == 0) write(6,*) 'No mgrid for lx1=2!'

    if (ifaxis) ifmgrid = .FALSE. 
    if (param(43) /= 0) ifmgrid = .FALSE. 

    npass = 1
    if (ifmhd) npass = 2
    do ipass=1,npass
        ifield = 1

        if (ifsplit .AND. ifmgrid) then

            if (ipass > 1) ifield = ifldmhd

            call swap_lengths
            call gen_fast_spacing(x,y,z)
             
            call hsmg_setup
            call h1mg_setup

        elseif ( .NOT. ifsplit) then ! Pn-Pn-2
#if 0
            if (ipass > 1) ifield = ifldmhd

            if (param(44) == 1) then !  Set up local overlapping solves
                call set_fem_data_l2(nel_proc,ndom,n_o,x,y,z,tri)
            else
                call swap_lengths
            endif
             
            e = 1
            if (ifield > 1) e = nelv+1

            call gen_fast_spacing(x,y,z)
            call gen_fast(df(1,e),sr(1,e),ss(1,e),st(1,e),x,y,z)

            call init_weight_op
            if (param(43) == 0) call hsmg_setup
#endif
        endif

        call set_up_h1_crs

    enddo
     
    return
    end subroutine set_overlap
!-----------------------------------------------------------------------
#if 0
    subroutine map_one_face12(x2,x1,iface,i12,i12t,w1,w2)
    use size_m
    use input

!     Interpolate iface of x1 (GLL pts) onto interior of face of x2 (GL pts).
!     Work arrays should be of size nx1*nx1 each.

    real :: x2(nx1,ny1,nz1)
    real :: x1(nx1,ny1,nz1)
    real :: w1(1),w2(1)
    real :: i12(1),i12t(1)


!     Extract surface data from x1

    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
    i=0
    do iz=kz1,kz2
        do iy=ky1,ky2
            do ix=kx1,kx2
                i     = i+1
                w1(i) = x1(ix,iy,iz)
            enddo
        enddo
    enddo

!     Interpolate from mesh 1 to 2

    if (if3d) then
        call mxm(i12 ,nx2,w1,nx1,w2,nx1)
        call mxm(w2,nx2,i12t,nx1,w1,nx2)
    else
        call mxm(i12 ,nx2,w1,nx1,w2,  1)
        call copy(w1,w2,nx2)
    endif

!     Write to interior points on face of x2

    kx1=min(kx1+1,nx1,kx2)
    kx2=max(kx2-1,  1,kx1)
    ky1=min(ky1+1,ny1,ky2)
    ky2=max(ky2-1,  1,ky1)
    kz1=min(kz1+1,nz1,kz2)
    kz2=max(kz2-1,  1,kz1)

    i=0
    do iz=kz1,kz2
        do iy=ky1,ky2
            do ix=kx1,kx2
                i     = i+1
                x2(ix,iy,iz) = w1(i)
            enddo
        enddo
    enddo

    return
    end subroutine map_one_face12
#endif
!-----------------------------------------------------------------------
