!-----------------------------------------------------------------------
    subroutine tensr3(v,nv,u,nu,A,Bt,Ct,w)

!     -  Tensor product application of v = (C x B x A) u
!        NOTE -- the transpose of B & C must be input, rather than B & C.

!     -  scratch arrays: w(nu*nu*nv)


    use size_m
    use input
    real :: v(nv,nv,nv),u(nu,nu,nu)
    real :: A(1),Bt(1),Ct(1)
    real :: w(1)

    if (nu > nv) then
        write(6,*) nid,nu,nv,' ERROR in tensr3. Contact P.Fischer.'
        write(6,*) nid,nu,nv,' Memory problem.'
        call exitt
    endif

    if (if3d) then
        nuv = nu*nv
        nvv = nv*nv
        call mxm(A,nv,u,nu,v,nu*nu)
        k=1
        l=1
        do iz=1,nu
            call mxm(v(k,1,1),nv,Bt,nu,w(l),nv)
            k=k+nuv
            l=l+nvv
        enddo
        call mxm(w,nvv,Ct,nu,v,nv)
    else
        call mxm(A,nv,u,nu,w,nu)
        call mxm(w,nv,Bt,nu,v,nv)
    endif
    return
    end subroutine tensr3
!-----------------------------------------------------------------------
    subroutine s_face_to_int(x,c)

!     Scale face and add to interior of element

    use size_m
    use input
    real :: x(nx1,ny1,nz1,1)

    do ie=1,nelv
    
        if (if3d) then
        
            do iz=2,nz1-1
                do ix=2,nx1-1
                    x(ix,2    ,iz,ie) = c*x(ix,1  ,iz,ie) + x(ix,2    ,iz,ie)
                    x(ix,ny1-1,iz,ie) = c*x(ix,ny1,iz,ie) + x(ix,ny1-1,iz,ie)
                enddo
            enddo
        
            do iz=2,nz1-1
                do iy=2,ny1-1
                    x(2    ,iy,iz,ie) = c*x(1  ,iy,iz,ie) + x(2    ,iy,iz,ie)
                    x(nx1-1,iy,iz,ie) = c*x(nx1,iy,iz,ie) + x(nx1-1,iy,iz,ie)
                enddo
            enddo
        
            do iy=2,ny1-1
                do ix=2,nx1-1
                    x(ix,iy,2    ,ie) = c*x(ix,iy,1  ,ie) + x(ix,iy,2    ,ie)
                    x(ix,iy,nz1-1,ie) = c*x(ix,iy,nz1,ie) + x(ix,iy,nz1-1,ie)
                enddo
            enddo
        
        else
        !         2D
            do ix=2,nx1-1
                x(ix,2    ,1,ie) = c*x(ix,1  ,1,ie) + x(ix,2    ,1,ie)
                x(ix,ny1-1,1,ie) = c*x(ix,ny1,1,ie) + x(ix,ny1-1,1,ie)
            enddo
            do iy=2,ny1-1
                x(2    ,iy,1,ie) = c*x(1  ,iy,1,ie) + x(2    ,iy,1,ie)
                x(nx1-1,iy,1,ie) = c*x(nx1,iy,1,ie) + x(nx1-1,iy,1,ie)
            enddo
        endif
    enddo
    return
    end subroutine s_face_to_int
!-----------------------------------------------------------------------
    subroutine dface_ext(x)
!     Extend interior to face of element

    use size_m
    use input
    real :: x(nx1,ny1,nz1,1)

    do ie=1,nelv
    
        if (if3d) then
        
            do iz=2,nz1-1
                do ix=2,nx1-1
                    x(ix,1  ,iz,ie) = x(ix,2    ,iz,ie)
                    x(ix,ny1,iz,ie) = x(ix,ny1-1,iz,ie)
                enddo
            enddo
        
            do iz=2,nz1-1
                do iy=2,ny1-1
                    x(1  ,iy,iz,ie) = x(2    ,iy,iz,ie)
                    x(nx1,iy,iz,ie) = x(nx1-1,iy,iz,ie)
                enddo
            enddo
        
            do iy=2,ny1-1
                do ix=2,nx1-1
                    x(ix,iy,1  ,ie) = x(ix,iy,2    ,ie)
                    x(ix,iy,nz1,ie) = x(ix,iy,nz1-1,ie)
                enddo
            enddo
        
        else
        
            do ix=2,nx1-1
                x(ix,1  ,1,ie) = x(ix,2    ,1,ie)
                x(ix,ny1,1,ie) = x(ix,ny1-1,1,ie)
            enddo
            do iy=2,ny1-1
                x(1  ,iy,1,ie) = x(2    ,iy,1,ie)
                x(nx1,iy,1,ie) = x(nx1-1,iy,1,ie)
            enddo
        
        endif
    enddo
    return
    end subroutine dface_ext
!-----------------------------------------------------------------------
    subroutine dface_add1si(x,c)
!     Scale interior and add to face of element

    use size_m
    use input
    real :: x(nx1,ny1,nz1,1)

    do ie=1,nelv
    
        if (if3d) then
        
            do iz=2,nz1-1
                do ix=2,nx1-1
                    x(ix,1  ,iz,ie) = x(ix,1  ,iz,ie) + c*x(ix,2    ,iz,ie)
                    x(ix,ny1,iz,ie) = x(ix,ny1,iz,ie) + c*x(ix,ny1-1,iz,ie)
                enddo
            enddo
        
            do iz=2,nz1-1
                do iy=2,ny1-1
                    x(1  ,iy,iz,ie) = x(1  ,iy,iz,ie) + c*x(2    ,iy,iz,ie)
                    x(nx1,iy,iz,ie) = x(nx1,iy,iz,ie) + c*x(nx1-1,iy,iz,ie)
                enddo
            enddo
        
            do iy=2,ny1-1
                do ix=2,nx1-1
                    x(ix,iy,1  ,ie) = x(ix,iy,1  ,ie) + c*x(ix,iy,2    ,ie)
                    x(ix,iy,nz1,ie) = x(ix,iy,nz1,ie) + c*x(ix,iy,nz1-1,ie)
                enddo
            enddo
        
        else
        
        !         2D
        
            do ix=2,nx1-1
                x(ix,1  ,1,ie) = x(ix,1  ,1,ie) + c*x(ix,2    ,1,ie)
                x(ix,ny1,1,ie) = x(ix,ny1,1,ie) + c*x(ix,ny1-1,1,ie)
            enddo
            do iy=2,ny1-1
                x(1  ,iy,1,ie) = x(1  ,iy,1,ie) + c*x(2    ,iy,1,ie)
                x(nx1,iy,1,ie) = x(nx1,iy,1,ie) + c*x(nx1-1,iy,1,ie)
            enddo
        
        endif
    enddo
    return
    end subroutine dface_add1si
!-----------------------------------------------------------------------
    subroutine do_weight_op(x)
    use size_m
    use input
    use tstep
    parameter(levb=lelv+lbelv)
    common /weightop/ w(lx2,lz2,2,3,levb)
    real :: w

    real :: x(0:nx1-1,0:ny1-1,0:nz1-1,1)
    integer :: e,e0,eb

    e0  = 0
    if (ifield > 1) e0 = nelv

    if (if3d) then
        do e=1,nelv
            eb = e0 + e
            do k=1,nz2
                do j=1,ny2
                    x(  1,j,k,e)=w(j,k,1,1,eb)*x(  1,j,k,e)
                    x(nx2,j,k,e)=w(j,k,2,1,eb)*x(nx2,j,k,e)
                enddo
            enddo
            do k=1,nz2
                do i=2,nx2-1
                    x(i,  1,k,e)=w(i,k,1,2,eb)*x(i,  1,k,e)
                    x(i,ny2,k,e)=w(i,k,2,2,eb)*x(i,ny2,k,e)
                enddo
            enddo
            do j=2,ny2-1
                do i=2,nx2-1
                    x(i,j,  1,e)=w(i,j,1,3,eb)*x(i,j,  1,e)
                    x(i,j,nz2,e)=w(i,j,2,3,eb)*x(i,j,nz2,e)
                enddo
            enddo
        enddo
    else
        do e=1,nelv
            eb = e0 + e
            do j=1,ny2
                x(  1,j,0,e)=w(j,1,1,1,eb)*x(  1,j,0,e)
                x(nx2,j,0,e)=w(j,1,2,1,eb)*x(nx2,j,0,e)
            enddo
            do i=2,nx2-1
                x(i,  1,0,e)=w(i,1,1,2,eb)*x(i,  1,0,e)
                x(i,ny2,0,e)=w(i,1,2,2,eb)*x(i,ny2,0,e)
            enddo
        enddo
    endif
    end subroutine do_weight_op
