!-----------------------------------------------------------------------
!     To do:
!        Differing BC's imposed for ophinv, incomprn, etc.
!        1-shot Fast solver for Helmholtz and pressure
!-----------------------------------------------------------------------

    subroutine compute_cfl(cfl,u,v,w,dt)

!     Given velocity field (u,v,w) and dt, compute current CFL number.

    use size_m
    use geom
    use input
    use soln
    use wz_m

    real :: u(nx1,ny1,nz1,nelv),v(nx1,ny1,nz1,nelv),w(nx1,ny1,nz1,nelv)

!     Store the inverse jacobian to speed up this operation

    common /cfldx/ dri(lx1),dsi(ly1),dti(lz1)

    integer :: e

    integer :: icalld
    save    icalld
    data    icalld /0/

    if (icalld == 0) then
        icalld=1
        call getdr(dri,zgm1(1,1),nx1)
        call getdr(dsi,zgm1(1,2),ny1)
        if (if3d) call getdr(dti,zgm1(1,3),nz1)
    endif

    cfl = 0.
    l   = 0

    if (if3d) then
        nxyz = nx1*ny1*nz1
        do e=1,nelv
            do k=1,nz1
                do j=1,ny1
                    do i=1,nx1
                        l = l+1
                        ur = ( u(i,j,k,e)*rxm1(i,j,k,e) &
                        +   v(i,j,k,e)*rym1(i,j,k,e) &
                        +   w(i,j,k,e)*rzm1(i,j,k,e) ) * jacmi(l,1)
                        us = ( u(i,j,k,e)*sxm1(i,j,k,e) &
                        +   v(i,j,k,e)*sym1(i,j,k,e) &
                        +   w(i,j,k,e)*szm1(i,j,k,e) ) * jacmi(l,1)
                        ut = ( u(i,j,k,e)*txm1(i,j,k,e) &
                        +   v(i,j,k,e)*tym1(i,j,k,e) &
                        +   w(i,j,k,e)*tzm1(i,j,k,e) ) * jacmi(l,1)
                         
                        cflr = abs(dt*ur*dri(i))
                        cfls = abs(dt*us*dsi(j))
                        cflt = abs(dt*ut*dti(k))
                         
                        cflm = cflr + cfls + cflt
                        cfl  = max(cfl,cflm)

                        cflf(i,j,k,e) = cflm
                         
                    enddo
                enddo
            enddo
        enddo
    else
        nxyz = nx1*ny1
        do e=1,nelv
            do j=1,ny1
                do i=1,nx1
                    l = l+1
                    ur = ( u(i,j,1,e)*rxm1(i,j,1,e) &
                    +   v(i,j,1,e)*rym1(i,j,1,e) ) * jacmi(l,1)
                    us = ( u(i,j,1,e)*sxm1(i,j,1,e) &
                    +   v(i,j,1,e)*sym1(i,j,1,e) ) * jacmi(l,1)

                    cflr = abs(dt*ur*dri(i))
                    cfls = abs(dt*us*dsi(j))

                    cflm = cflr + cfls
                    cfl  = max(cfl,cflm)

                    cflf(i,j,1,e) = cflm

                enddo
            enddo
        enddo
    endif

    cfl = glmax(cfl,1)

    return
    end subroutine compute_cfl
!-----------------------------------------------------------------------
    subroutine getdr(dri,zgm1,nx1)
    real :: dri(nx1),zgm1(nx1)

    dri(1) = zgm1(2) - zgm1(1)   !  Compute 1/Dx
    do i=2,nx1-1
        dri(i) = 0.5*( zgm1(i+1) - zgm1(i-1) )
    enddo
    dri(nx1) = zgm1(nx1) - zgm1(nx1-1)

    call invcol1(dri,nx1)

    return
    end subroutine getdr
!-----------------------------------------------------------------------
    subroutine set_dealias_rx

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

    use size_m
    use geom
    use input
    use tstep ! for istep

    common /dealias1/ zd(lxd),wd(lxd)
    integer :: e

    integer :: ilstep
    save    ilstep
    data    ilstep /-1/

    if ( .NOT. ifgeom .AND. ilstep > 1) return  ! already computed
    if (ifgeom .AND. ilstep == istep)  return  ! already computed
    ilstep = istep

    nxyz1 = nx1*ny1*nz1
    nxyzd = nxd*nyd*nzd

    call zwgl (zd,wd,nxd)  ! zwgl -- NOT zwgll!

    if (if3d) then
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),rzm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,5,e),sym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,6,e),szm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,7,e),txm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,8,e),tym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,9,e),tzm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd

            l = 0
            do k=1,nzd
                do j=1,nyd
                    do i=1,nxd
                        l = l+1
                        w = wd(i)*wd(j)*wd(k)
                        do ii=1,9
                            rx(l,ii,e) = w*rx(l,ii,e)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    else ! 2D
    
        do e=1,nelv

        !           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),sxm1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sym1(1,1,1,e),nx1,nxd,if3d,0) ! 0 --> fwd

            l = 0
            do j=1,nyd
                do i=1,nxd
                    l = l+1
                    w = wd(i)*wd(j)
                    do ii=1,4
                        rx(l,ii,e) = w*rx(l,ii,e)
                    enddo
                enddo
            enddo
        enddo

    endif

    return
    end subroutine set_dealias_rx
!-----------------------------------------------------------------------
