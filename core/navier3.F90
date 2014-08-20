    SUBROUTINE EPREC2(Z2,R2)
!----------------------------------------------------------------

!     Precondition the explicit pressure operator (E) with
!     a Neumann type (H1) Laplace operator: JT*A*J.
!     Invert A by conjugate gradient iteration or multigrid.

!     NOTE: SCRNS is used.

!----------------------------------------------------------------
    use size_m
    use geom
    INCLUDE 'INPUT'
    INCLUDE 'MASS'
    INCLUDE 'PARALLEL'
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    REAL ::           Z2   (LX2,LY2,LZ2,LELV)
    REAL ::           R2   (LX2,LY2,LZ2,LELV)
    COMMON /SCRNS/ MASK (LX1,LY1,LZ1,LELV) &
    ,R1   (LX1,LY1,LZ1,LELV) &
    ,X1   (LX1,LY1,LZ1,LELV) &
    ,W2   (LX2,LY2,LZ2,LELV) &
    ,H1   (LX1,LY1,LZ1,LELV) &
    ,H2   (LX1,LY1,LZ1,LELV)
    REAL ::    MASK

    integer :: icalld
    save    icalld
    data    icalld/0/
    icalld=icalld+1

    ntot2  = nx2*ny2*nz2*nelv
    call rzero(z2,ntot2)




!  Both local and global solver...
    call dd_solver (z2,r2)



!  Local solver only
!      call local_solves_fdm (z2,r2)



    return
    END SUBROUTINE EPREC2
!-----------------------------------------------------------------------
    subroutine dd_solver(u,v)

    use ctimer
    use size_m
    use domain
    include 'INPUT'
    include 'PARALLEL'
    include 'SOLN'

    real :: u(1),v(1)
    common /scrprc/ uc(lx1*ly1*lz1*lelt)

    if (icalld == 0) then
        tddsl=0.0
        tcrsl=0.0
        nddsl=0
        ncrsl=0
    endif
    icalld = icalld + 1
    nddsl  = nddsl  + 1
    ncrsl  = ncrsl  + 1

    ntot  = nx2*ny2*nz2*nelv
    call rzero(u,ntot)

    etime1=dnekclock()
    call local_solves_fdm    (u,v)
    tddsl=tddsl+dnekclock()-etime1

    etime1=dnekclock()
    call crs_solve_l2 (uc,v)
    tcrsl=tcrsl+dnekclock()-etime1

    alpha = 10.
!     if (param(89).ne.0.) alpha = abs(param(89))
    call add2s2(u,uc,alpha,ntot)

    return
    end subroutine dd_solver
!-----------------------------------------------------------------------
    subroutine rar2_out(x,name13)
    use size_m

    real :: x(lx2,ly2,lz2,lelt)
    character(13) :: name13

    if (nelv > 20) return
    write(6,*)
    write(6,1) name13
    1 format(a13)
    if (nelv > 2) then
        write(6,*)
        do j=ny2,1,-1
            write(6,6) (x(k,j,1,3),k=1,nx2),(x(k,j,1,4),k=1,nx2)
        enddo
        write(6,*)
        write(6,*)
    endif

    do j=ny2,1,-1
        write(6,6) (x(k,j,1,1),k=1,nx2),(x(k,j,1,2),k=1,nx2)
    enddo
    write(6,*)
    6 format(3f8.4,5x,3f8.4)
    return
    end subroutine rar2_out
!-----------------------------------------------------------------------
    subroutine rarr_out2(x,name13)
    use size_m
    include 'INPUT'

    real :: x(lx2,ly2,lz2,lelt)
    character(13) :: name13

    if (nelv > 20) return
    write(6,*)
    write(6,1) name13
    1 format('rarr2',3x,a13)

!     3 D

    if (if3d) then
        do iz=1,nz1
            write(6,*)
            do j=ny1,1,-1
                write(6,3) (x(k,j,iz,1),k=1,nx2),(x(k,j,iz,2),k=1,nx2)
            enddo
        enddo
        write(6,*)
        write(6,*)
        return
    endif

!     2 D

    if (nelv > 2) then
        write(6,*)
        do j=ny2,1,-1
            write(6,6) (x(k,j,1,3),k=1,nx2),(x(k,j,1,4),k=1,nx2)
        enddo
        write(6,*)
        write(6,*)
    endif

    do j=ny2,1,-1
        write(6,6) (x(k,j,1,1),k=1,nx2),(x(k,j,1,2),k=1,nx2)
    enddo
    write(6,*)
    3 format(4f6.2,5x,4f6.2)
    6 format(4f8.5,5x,4f8.5)
    return
    end subroutine rarr_out2
!-----------------------------------------------------------------------
