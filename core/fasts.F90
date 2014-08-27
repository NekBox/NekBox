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
