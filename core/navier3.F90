    SUBROUTINE EPREC2(Z2,R2)
!----------------------------------------------------------------

!     Precondition the explicit pressure operator (E) with
!     a Neumann type (H1) Laplace operator: JT*A*J.
!     Invert A by conjugate gradient iteration or multigrid.

!     NOTE: SCRNS is used.

!----------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use parallel
    use soln
    use tstep
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
    use input
    use parallel
    use soln

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
