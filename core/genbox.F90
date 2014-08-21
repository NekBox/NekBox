!-----------------------------------------------------------------------
#if 0
    subroutine genbox

    use size_m
    use input
    include 'ZPER'

    character(132) :: string
    character(1) ::  string1(132)
    equivalence (string,string1)


    call blank (string,132)
!     call gets  (string,132,iend,9) ! read name of "box" file from .rea file
!     read (9,132) string       ! *** MESH?
    132 format(a132)

    if (nid == 0) then
    
    !        Read in name of previously generated NEKTON data set.
    !        (Must have same dimension and number of fields as current run)
    
        write(6,*) 'opening file:'
    !        write(6,*) string
        lfname = ltrunc(reafle,132) - 4
        call blank (string,132)
        call chcopy(string,reafle,lfname)
        call chcopy(string1(lfname+1),'.box',4)
        open (unit=7,file=string,status='old',err=999)
        call gets(string,132,iend,7)  ! This is a dummy read, for compatibility
    
    !        here is where the 2d/3d determination is made....
    
        call geti1(ndim,iend,7)
        ndim = iabs(ndim)
        call getrv(rfld,1,iend,7)
        nfld = int(rfld)                  ! Determine number of fields
        if(rfld /= nfld) nfld = nfld +1
        call getbox(xgtp,ygtp,zgtp,nfld)  ! Read in the .box file
        close (unit=7)
    
    endif

    call bcpbox(nfld)

    call mapbox(nelx,nely,nelz)  ! Map elements to processors
    call makebox(xgtp,ygtp,zgtp,nfld) ! Fill up local elements on _this_ processor

!     call outbox_mesh

    ifgfdm = .TRUE. 
    param(116) = nelx
    param(117) = nely
    param(118) = nelz

    if3d = .FALSE. 
    if (ndim == 3) if3d= .TRUE. 

    return

    999 continue
    if (nid == 0) write(6,*) 'ABORT: Could not find box file ',string
    call exitt

    end subroutine genbox
#endif
!-----------------------------------------------------------------------
    subroutine gen_gtp_vertex (vertex,ncrnr)

    use size_m
    use input
    use parallel
    include 'ZPER'

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
#if 0
    subroutine getbox(x,y,z,nfld)

    use size_m
    use input
    include 'ZPER'

    real :: x(0:1),y(0:1),z(0:1)

    logical :: ifevenx,ifeveny,ifevenz,if_multi_seg
    character(1) :: boxcirc

    integer :: nelxyz(3)


    ifevenx = .FALSE. 
    ifeveny = .FALSE. 
    ifevenz = .FALSE. 

    call gets(boxcirc,1,iend,7)
    if (iend == 1) goto 99
    if (boxcirc == 'c' .OR. boxcirc == 'C') &
    call getr2(xc,yc,iend,7)
    if_multi_seg = .FALSE. 
    if (boxcirc == 'm' .OR. boxcirc == 'M') &
    if_multi_seg = .TRUE. 

    if (if_multi_seg) then
        maxx = 1
        call get_multi_seg &
        (nelxyz,x(0),y(0),z(0),maxx,if3d)
        nelx = nelxyz(1)
        nely = nelxyz(2)
        nelz = 1
        if (if3d) nelz = nelxyz(3)
        ifld0 = 2
        if (ifflow) ifld0 = 1
        ifld1 = ifld0-1 + nfld
        nbc = 2*ndim
        do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,6,iend,7)
            write(6,*) 'CBC1',(gtp_cbc(k,ifld),k=1,nbc),ifld
        enddo
    elseif (if3d) then
        call geti3(nelx,nely,nelz,iend,7)
    
        if (nelx < 0) ifevenx = .TRUE. 
        if (nely < 0) ifeveny = .TRUE. 
        if (nelz < 0) ifevenz = .TRUE. 
        nelx = abs(nelx)
        nely = abs(nely)
        nelz = abs(nelz)
    
        if (nelx > lelx .OR. nely > lely .OR. nelz > lelz) then
            write(6,*) 'Abort, increase lelx,..lelz and recompile', &
            nelx,nely,nelz,lelx,lely,lelz
            call exitt
        endif
    
        nery = nely
        if (boxcirc == 'C') nery = 1
        if (iend == 1) goto 99
        ninbox = nelx*nely*nelz
        write(6,6) ninbox,nelx,nely,nelz,nfld
    
        if (ifevenx) then
            call getr3(x0,x1,ratio,iend,7)
            if (ratio <= 0) ratio=1.
            dx = (x1-x0)/nelx
            x(0) = 0.
            do k=1,nelx
                x(k) = x(k-1) + dx
                dx        = dx*ratio
            enddo
            scale = (x1-x0)/x(nelx)
            do k=0,nelx
                x(k) = x0 + scale*x(k)
            enddo
        else
            call getrv(x(0),nelx+1,iend,7)
        endif
    
        if (ifeveny) then
            call getr3(y0,y1,ratio,iend,7)
            if (ratio <= 0) ratio=1.
            dy = (y1-y0)/nely
            y(0) = 0.
            do k=1,nely
                y(k) = y(k-1) + dy
                dy        = dy*ratio
            enddo
            scale = (y1-y0)/y(nely)
            do k=0,nely
                y(k) = y0 + scale*y(k)
            enddo
        else
            call getrv(y(0),nely+1,iend,7)
        endif
    
        if (ifevenz) then
            call getr3(z0,z1,ratio,iend,7)
            if (ratio <= 0) ratio=1.
            dz = (z1-z0)/nelz
            z(0) = 0.
            do k=1,nelz
                z(k) = z(k-1) + dz
                dz        = dz*ratio
            enddo
            scale = (z1-z0)/z(nelz)
            do k=0,nelz
                z(k) = z0 + scale*z(k)
            enddo
        else
            call getrv(z(0),nelz+1,iend,7)
        endif
    
        ifld0 = 2
        if (ifflow) ifld0 = 1
        ifld1 = ifld0-1 + nfld
        nbc   = 2*ndim
        write(6,*) 'this is ifld01:',ifld0,ifld1
        do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,6,iend,7)
            write(6,*) 'CBC2',(gtp_cbc(k,ifld),k=1,nbc),ifld
        enddo
    else
        nelz = 1
        call geti2(nelx,nely,iend,7)
    
        if (nelx < 0) ifevenx = .TRUE. 
        if (nely < 0) ifeveny = .TRUE. 
        nelx = abs(nelx)
        nely = abs(nely)
        if (nelx > lelx .OR. nely > lely) then
            write(6,*) 'Abort, increase lelx,lely and recompile', &
            nelx,nely,lelx,lely
            call exitt
        endif
    
        nery = nely
        if (boxcirc == 'C') nery = 1
        if (iend == 1) goto 99
        ninbox = nelx*nely*nelz
        write(6,6) ninbox,nelx,nely,nery,nfld
    
        if (ifevenx) then
            call getr3(x0,x1,ratio,iend,7)
            if (ratio <= 0) ratio=1.
            dx = (x1-x0)/nelx
            x(0) = 0.
            do k=1,nelx
                x(k) = x(k-1) + dx
                dx        = dx*ratio
            enddo
            scale = (x1-x0)/x(nelx)
            do k=0,nelx
                x(k) = x0 + scale*x(k)
            enddo
        else
            call getrv(x(0),nelx+1,iend,7)
        endif
    
        if (ifeveny) then
            call getr3(y0,y1,ratio,iend,7)
            if (ratio <= 0) ratio=1.
            dy = (y1-y0)/nely
            y(0) = 0.
            do k=1,nely
                y(k) = y(k-1) + dy
                dy        = dy*ratio
            enddo
            scale = (y1-y0)/y(nely)
            do k=0,nely
                y(k) = y0 + scale*y(k)
            enddo
        else
            call getrv(y(0),nely+1,iend,7)
        endif
    
        ifld0 = 2
        if (ifflow) ifld0 = 1
        ifld1 = ifld0-1 + nfld
        nbc   = 2*ndim
        do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,4,iend,7)
            write(6,*) 'CBC3',(gtp_cbc(k,ifld),k=1,nbc),ifld
        enddo
    endif
    nel = nelx*nely*nelz
    6 format('Reading',i12,' =',3i6,' elements for box',i3,'.')
    99 continue

    if (nid == 0) then
        write(6,*) 'Beginning construction of new mesh.' &
        ,nel,nelx,nely,nelz
    endif

    return
    end subroutine getbox
!-----------------------------------------------------------------------
    subroutine mapbox(melx,mely,melz)

!     Given melx  x  mely  x  melz elements, distribute to np processors

    use size_m
    use parallel
    include 'ZPER'

    nelx = melx
    nely = mely
    nelz = melz

    nelgv = nelx*nely*nelz
    nelgt = nelx*nely*nelz

    if (nelgt > lelg .OR. nelgv > lelg) then
        write(6,*) nid,' NELGT too large:',nelgt,nelgv,lelg
        call exitt
    endif

    call mapelpr

    if (nelt > lelt .OR. nelv > lelv) then
        write(6,*) nid,' NEL too large:',nelt,nelv,lelv,lelt
        call exitt
    endif


    return
    end subroutine mapbox
!-----------------------------------------------------------------------
    subroutine makebox(x,y,z,nfld)

    use size_m
    use input
    use parallel
    include 'ZPER'

    real :: x(0:1),y(0:1),z(0:1)
    integer :: e,eg,ex,ey,ez

    character(3) :: cbc1,   cbc2,   cbc3,   cbc4,   cbc5,   cbc6
    real ::        rbc1(5),rbc2(5),rbc3(5),rbc4(5),rbc5(5),rbc6(5)

    integer :: eface(6)
    save    eface
    data    eface  / 4,2,1,3,5,6 /

!     Set up TEMPORARY value for NFIELD - NFLDT

    nfldt = 1
    if (ifheat) nfldt=2+npscal
    if (ifmhd ) nfldt=2+npscal+1
    nbcs      = nfldt
    ibcs      = 2
    if (ifflow) ibcs = 1

    nbc = 5*6*lelt*(ldimt1+1)
    call rzero(bc,nbc)

    eg = 0
    do ez=1,nelz
        do ey=1,nely
            do ex=1,nelx
                eg   = eg+1
                jnid = gllnid(eg)
                if (jnid == nid) then
                    e = gllel(eg)
                
                    xc(1,e) = x(ex-1)
                    xc(2,e) = x(ex  )
                    xc(3,e) = x(ex  )
                    xc(4,e) = x(ex-1)
                
                    yc(1,e) = y(ey-1)
                    yc(2,e) = y(ey-1)
                    yc(3,e) = y(ey  )
                    yc(4,e) = y(ey  )
                
                    if (if3d) then
                    
                        xc(5,e) = x(ex-1)
                        xc(6,e) = x(ex  )
                        xc(7,e) = x(ex  )
                        xc(8,e) = x(ex-1)
                    
                        yc(5,e) = y(ey-1)
                        yc(6,e) = y(ey-1)
                        yc(7,e) = y(ey  )
                        yc(8,e) = y(ey  )
                    
                        zc(1,e) = z(ez-1)
                        zc(2,e) = z(ez-1)
                        zc(3,e) = z(ez-1)
                        zc(4,e) = z(ez-1)
                    
                        zc(5,e) = z(ez  )
                        zc(6,e) = z(ez  )
                        zc(7,e) = z(ez  )
                        zc(8,e) = z(ez  )
                    
                    endif
                
                    do ifld=ibcs,nbcs
                    
                        ilftz = mod1(ez-1+nelz,nelz)
                        irgtz = mod1(ez+1+nelz,nelz)
                        ilfty = mod1(ey-1+nely,nely)
                        irgty = mod1(ey+1+nely,nely)
                        ilftx = mod1(ex-1+nelx,nelx)
                        irgtx = mod1(ex+1+nelx,nelx)
                    
                        rbc1(1) = nely*nelx*(ez-1)+nelx*(ey-1) + ilftx
                        rbc1(2) = eface(2)
                        rbc2(1) = nely*nelx*(ez-1)+nelx*(ey-1) + irgtx
                        rbc2(2) = eface(1)
                    
                        rbc3(1) = nely*nelx*(ez-1)+nelx*(ilfty-1) + ex
                        rbc3(2) = eface(4)
                        rbc4(1) = nely*nelx*(ez-1)+nelx*(irgty-1) + ex
                        rbc4(2) = eface(3)
                    
                        rbc5(1) = nely*nelx*(ilftz-1)+nelx*(ey-1) + ex
                        rbc5(2) = eface(6)
                        rbc6(1) = nely*nelx*(irgtz-1)+nelx*(ey-1) + ex
                        rbc6(2) = eface(5)
                    
                        if (ex == 1) then
                            cbc1=gtp_cbc(1,ifld)
                            if (cbc1 /= 'P  ') call rzero(rbc1,5)
                        else
                            cbc1='E  '
                        endif
                    
                        if (ex == nelx) then
                            cbc2=gtp_cbc(2,ifld)
                            if (cbc2 /= 'P  ') call rzero(rbc2,5)
                        else
                            cbc2='E  '
                        endif
                    
                        if (ey == 1) then
                            cbc3=gtp_cbc(3,ifld)
                            if (cbc3 /= 'P  ') call rzero(rbc3,5)
                        else
                            cbc3='E  '
                        endif
                    
                        if (ey == nely) then
                            cbc4=gtp_cbc(4,ifld)
                            if (cbc4 /= 'P  ') call rzero(rbc4,5)
                        else
                            cbc4='E  '
                        endif
                    
                        if (ez == 1) then
                            cbc5=gtp_cbc(5,ifld)
                            if (cbc5 /= 'P  ') call rzero(rbc5,5)
                        else
                            cbc5='E  '
                        endif
                    
                        if (ez == nelz) then
                            cbc6=gtp_cbc(6,ifld)
                            if (cbc6 /= 'P  ') call rzero(rbc6,5)
                        else
                            cbc6='E  '
                        endif
                    
                        cbc(1,e,ifld) = cbc3
                        cbc(2,e,ifld) = cbc2
                        cbc(3,e,ifld) = cbc4
                        cbc(4,e,ifld) = cbc1
                        cbc(5,e,ifld) = cbc5
                        cbc(6,e,ifld) = cbc6
                    
                        call copy(bc(1,1,e,ifld),rbc3,5)
                        call copy(bc(1,2,e,ifld),rbc2,5)
                        call copy(bc(1,3,e,ifld),rbc4,5)
                        call copy(bc(1,4,e,ifld),rbc1,5)
                        call copy(bc(1,5,e,ifld),rbc5,5)
                        call copy(bc(1,6,e,ifld),rbc6,5)
                    
                    enddo
                endif
            
            enddo
        enddo
    enddo

    return
    end subroutine makebox
!-----------------------------------------------------------------------
    subroutine geti1(i,iend,io)

!     Get an integer from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,*) i
    close(unit=99)
    return
    end subroutine geti1
!-----------------------------------------------------------------------
    subroutine geti2(i1,i2,iend,io)

!     Get an integer from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,*) i1,i2
    close(unit=99)
    return
    end subroutine geti2
!-----------------------------------------------------------------------
    subroutine geti3(i1,i2,i3,iend,io)

!     Get an integer from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,*) i1,i2,i3
    close(unit=99)
    return
    end subroutine geti3
!-----------------------------------------------------------------------
    subroutine getr2(r1,r2,iend,io)

!     Get two reals from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,*) r1,r2
    close(unit=99)
    return
    end subroutine getr2
!-----------------------------------------------------------------------
    subroutine getr3(r1,r2,r3,iend,io)

!     Get two reals from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,*) r1,r2,r3
    close(unit=99)
    return
    end subroutine getr3
!-----------------------------------------------------------------------
    subroutine getrv(r,n,iend,io)
    real :: r(n)
    parameter (big=1.e29)

!     Get reals from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return

!     serious hack for f90
    call cfill(r,big,n)

    open(unit=99,file='box.tmp')
    read(99,*,end=1) (r(k),k=1,n)
    close(unit=99)
    return
    1 continue
    close(unit=99)
    do k=1,n
        if (r(k) == big) goto 15
    enddo
    15 continue
    write(6,*) 'this is k:',k,n
    read( 7,*,end=2) (r(j),j=k,n)
    return
    2 continue
    iend=1
    return
    end subroutine getrv
!-----------------------------------------------------------------------
    subroutine getiv(r,n,iend,io)
    integer :: r(n)
    parameter (ibig=99999999)

!     Get integers from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return

!     serious hack for f90
    call ifill(r,ibig,n)

    open(unit=99,file='box.tmp')
    read(99,*,end=1) (r(k),k=1,n)
    close(unit=99)
    return
    1 continue
    close(unit=99)
    do k=1,n
        if (r(k) == ibig) goto 15
    enddo
    15 continue
    write(6,*) 'this is k:',k,n
    read( 7,*,end=2) (r(j),j=k,n)
    return
    2 continue
    iend=1
    return
    end subroutine getiv
!-----------------------------------------------------------------------
    subroutine getcv0(c,m,n,iend,io)
    character(1) :: c(m,n)
    character(1) :: adum(132)

!     Get character strings, with no seperator, from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,1,end=2) (adum(k),k=1,132)
    1 format(132a1)
    2 continue

    i = 0
    do l=1,n
        do k=1,m
            i = i+1
            c(k,l) = adum(i)
        enddo
    enddo
    do j=1,n
        write(6,3) m,n,(c(i,j),i=1,m)
        3 format(2i4,'getcv0:',132a1)
    enddo

    close(unit=99)
    return
    end subroutine getcv0
!-----------------------------------------------------------------------
    subroutine getcv(c,m,n,iend,io)
    character(1) :: c(m,n)
    character(1) :: adum(132)

!     Get character strings, with single seperator, from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,1,end=2) (adum(k),k=1,132)
    1 format(132a1)
    2 continue

    i = 0
    do l=1,n
        do k=1,m
            i = i+1
            c(k,l) = adum(i)
        enddo
    
    !        bump pointer for "," in input string
        i = i+1
    enddo
    do j=1,n
        write(6,3) m,n,(c(i,j),i=1,m)
        3 format(2i4,'getcv:',132a1)
    enddo

    close(unit=99)
    return
    end subroutine getcv
!-----------------------------------------------------------------------
    subroutine gets(c,n,iend,io)
    character(1) :: c(n)

!     Get character string from first uncommented line

    call scannocom(iend,io)
    if (iend /= 0) return
    open(unit=99,file='box.tmp')
    read(99,1) (c(k),k=1,n)
    1 format(132a1)
    close(unit=99)
    return
    end subroutine gets
!-----------------------------------------------------------------------
    subroutine get_multi_seg(nelxyz,x,y,z,m,if3d)

    integer :: nelxyz(3)
    real :: x(0:m),y(0:m),z(0:m)
    logical :: if3d

    integer :: nels (1000)
    real ::    gains(1000)
    real ::    xs   (0:1000)

!     This routine allows the user to input a complex sequence of
!     segments for each of the x, y and z directions, so that
!     genbox will generate the desired tensor product array.

!     Input format for each coordinat is

!     nsegs
!     nel_1    nel_2 ...  nel_nsegs
!     gain_1  gain_2 ... gain_nsegs
!     x_0, x_1, x_2, ...    x_nsegs

!     The output will be:

!     x(0)...x(nel_1)...x_(nel_2)...x(nel_segs),

!     where the segment between x(e-1) and x(e) is filled with
!     a distribution determined by gain_e.   Uniform spacing
!     corresponds to gain_e = 1, otherwise, a geometric sequence
!     is generated, with dx_i+1 = dx_i * gain_e

    ndim = 2
    if (if3d) ndim=3

    do id=1,ndim
    
        call geti1(nseg,iend,7)
        call getiv(nels ,nseg  ,iend,7)
        call getrv(xs   ,nseg+1,iend,7)
        call getrv(gains,nseg  ,iend,7)
        write(6,*) 'NSEG: ',id,nseg,(nels(k),k=1,nseg)
    
        nelxyz(id) = ilsum(nels,nseg)
        if (nelxyz(id) > m) then
            write(6,*) 'NEL too large:',nelxyz(id),m,id
            write(6,*) 'Abort in get_multi_seg.'
            call exitt
        endif
    
        k = 0
        do is=1,nseg
            nel = nels(is)
            x0  = xs(is-1)
            x1  = xs(is)
            if (id == 1) call geometric_x(x(k),nel,x0,x1,gains(is))
            if (id == 2) call geometric_x(y(k),nel,x0,x1,gains(is))
            if (id == 3) call geometric_x(z(k),nel,x0,x1,gains(is))
            k = k + nel
        enddo
    enddo

    return
    end subroutine get_multi_seg
!-----------------------------------------------------------------------
    subroutine geometric_x(x,n,x0,x1,gain)

    real :: x(0:n)

    x(0) = 0.
    dx = (x1-x0)/n
    do i=1,n
        x(i) = x(i-1) + dx
        dx   = dx*gain
    enddo

    scale = (x1-x0)/(x(n)-x(0))
    do i=0,n
        x(i) = x0 + scale*x(i)
    enddo

    return
    end subroutine geometric_x
!-----------------------------------------------------------------------
    subroutine get_xyz_distribution (x,nelx)
    real :: x(0:2)
    if (nelx >= 0) return ! else, uniform or geometric spacing

    x0 = x(0)
    x1 = x(1)
    r  = x(2)
    write(6,*) 'x0:',x0,x1,r,nelx

    nelx = abs(nelx)
    if (r == 1.0) then
        dx = (x1-x0)/nelx
        do i=1,nelx
            x(i) = x0 + i*dx
        enddo
        x(nelx) = x1
    else
        dx = (x1-x0)/nelx
        x(0) = 0.
        do i=1,nelx
            x(i) = x(i-1) + dx
            dx = r*dx
        enddo
        scale = (x1-x0)/x(nelx)
        call cmult(x,scale,nelx+1)
        call cadd (x,x0   ,nelx+1)
        x(nelx) = x1
    endif

    return
    end subroutine get_xyz_distribution
!-----------------------------------------------------------------------
    subroutine scannocom(iend,infile)

!     scan through infile until "no comment" is found

    character(132) :: string
    character(1) ::  string1(132)
    equivalence (string1,string)

    character(1) :: comment
    save        comment
    data        comment /'#'/

    iend = 0
    do line=1,100000
        call blank(string,132)
        read (infile ,80,end=100,err=100) string
    
        write(6,80) string
        if   (indx1(string,comment,1) /= 1) then
            open(unit=99,file='box.tmp')
            len = ltrunc(string,132)
            write(99,81) (string1(k),k=1,len)
            close(unit=99)
            return
        else
            i1 = indx1(string,comment,1)
        !             write(6,*) i1,' i1', comment
        endif
    
    enddo
    80 format(a132)
    81 format(132a1)

    100 continue
    iend = 1
    return
    end subroutine scannocom
!-----------------------------------------------------------------------
    function ilsum(x,n)
    integer :: x(1)
    ilsum = 0
    do i=1,n
        ilsum = ilsum + x(i)
    enddo
    return
    end function ilsum
#endif
!-----------------------------------------------------------------------
    subroutine jjnt(x,n)
    integer :: x(1)
    do i=1,n
        x(i) = i
    enddo
    return
    end subroutine jjnt
!-----------------------------------------------------------------------
#if 0
    subroutine bcpbox(nfld)

    use size_m
    use input
    include 'ZPER'
    use parallel

    call bcast(ndim,isize)
    call bcast(nfld,isize)
    call bcast(nelx,isize)
    call bcast(nely,isize)
    call bcast(nelz,isize)

    lnx = wdsize*(nelx+1)
    lny = wdsize*(nely+1)
    lnz = wdsize*(nelz+1)
    call bcast(xgtp,lnx)
    call bcast(ygtp,lny)
    call bcast(zgtp,lnz)

    ifld0 = 2
    if (ifflow) ifld0 = 1
    ifld1 = ifld0-1 + nfld
    do ifld=ifld0,ifld1
        call bcast(gtp_cbc(1,ifld),18)
    enddo

    return
    end subroutine bcpbox
!-----------------------------------------------------------------------
    subroutine outbox_mesh
    use size_m
    use input
    include 'ZPER'
    use parallel
    integer :: e

    if (np > 1) return
    if (ndim == 2) return

    do e=1,nelv
        write(88,1) e
        write(88,2) (xc(k,e),k=1,4)
        write(88,2) (yc(k,e),k=1,4)
        write(88,2) (zc(k,e),k=1,4)
        write(88,2) (xc(k,e),k=5,8)
        write(88,2) (yc(k,e),k=5,8)
        write(88,2) (zc(k,e),k=5,8)
    enddo
    1 format('            ELEMENT',i12,' [    1a]    GROUP     0')
    2 format(1p4e15.7)

    return
    end subroutine outbox_mesh
#endif
!-----------------------------------------------------------------------
