!-----------------------------------------------------------------------
    subroutine gen_gtp_vertex (vertex,ncrnr)

    use size_m
    use input
    use parallel
    use zper

!     integer vertex((2**ldim)*lelt)    ! local -- 2D for now !?  long?
    integer :: vertex(2**ldim,1)         ! local -- 2D for now !?  long?
    integer :: e,eg,ex,ey,ez

    common /ctmp1/ indx(lelx+1) &
    , indy(lely+1) &
    , indz(lelz+1)

    nptsx = nelx + 1
    nptsy = nely + 1
    nptsz = nelz + 1

    call jjnt(indx,nptsx)
    call jjnt(indy,nptsy)
    call jjnt(indz,nptsz)

    ifld = 2                                          ! ifield > 2?
    if (ifflow) ifld = 1

    nbc = 2*ndim
    write(6,6) nid,ifld,(gtp_cbc(k,ifld),k=1,nbc)
    6 format(2i4,'  GTP BC:',6(2x,a3))

    if (gtp_cbc(1,ifld) == 'P  ') indx(nptsx) = 1  ! Set last point to "1"
    if (gtp_cbc(3,ifld) == 'P  ') indy(nptsy) = 1
    if (gtp_cbc(5,ifld) == 'P  ') indz(nptsz) = 1

    if (gtp_cbc(1,ifld) == 'P  ') nptsx = nelx
    if (gtp_cbc(3,ifld) == 'P  ') nptsy = nely
    if (gtp_cbc(5,ifld) == 'P  ') nptsz = nelz


    do e=1,nelv                                       ! nelt/nel?
        eg = lglel(e)
        call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
        if (if3d) then
            l = 0
            do k=0,1
                do j=0,1
                    do i=0,1
                        l = l+1
                        i1 = i+ex
                        j1 = j+ey
                        k1 = k+ez
                        vertex(l,e) =  indx(i1) &
                        + (indy(j1)-1)*nptsx &
                        + (indz(k1)-1)*nptsx*nptsy
                    enddo
                enddo
            enddo
        else
            l = 0
            do j=0,1
                do i=0,1
                    l = l+1
                    i1 = i+ex
                    j1 = j+ey
                    vertex(l,e) =  indx(i1) &
                    + (indy(j1)-1)*nptsx
                enddo
            enddo
        endif
    enddo

    return
    end subroutine gen_gtp_vertex
!-----------------------------------------------------------------------
    subroutine jjnt(x,n)
    integer :: x(1)
    do i=1,n
        x(i) = i
    enddo
    return
    end subroutine jjnt
!-----------------------------------------------------------------------
