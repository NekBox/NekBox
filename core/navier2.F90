!-----------------------------------------------------------------------
    subroutine aspect_ratios(ar)

!     7/6/96

!     This routine returns the aspect ratio of a
!     conglomerate of a set of simplices defined by (tri,nt)

    use size_m
    use geom
    include 'INPUT'
    include 'MASS'

    real :: ar(1)
    real :: xx(9)

    nxyz = nx1*ny1*nz1

    if (if3d) then
        do ie=1,nelv
            vol = vlsum(bm1(1,1,1,ie),nxyz)
            x0  = vlsc2(xm1(1,1,1,ie),bm1(1,1,1,ie),nxyz)/vol
            y0  = vlsc2(ym1(1,1,1,ie),bm1(1,1,1,ie),nxyz)/vol
            z0  = vlsc2(zm1(1,1,1,ie),bm1(1,1,1,ie),nxyz)/vol
        
            call rzero(xx,9)
            do i=1,nxyz
                x10 = xm1(i,1,1,ie) - x0
                y10 = ym1(i,1,1,ie) - y0
                z10 = zm1(i,1,1,ie) - z0
            
                xx(1) = xx(1) + bm1(i,1,1,ie)*x10*x10
                xx(2) = xx(2) + bm1(i,1,1,ie)*x10*y10
                xx(3) = xx(3) + bm1(i,1,1,ie)*x10*z10
                xx(4) = xx(4) + bm1(i,1,1,ie)*y10*x10
                xx(5) = xx(5) + bm1(i,1,1,ie)*y10*y10
                xx(6) = xx(6) + bm1(i,1,1,ie)*y10*z10
                xx(7) = xx(7) + bm1(i,1,1,ie)*z10*x10
                xx(8) = xx(8) + bm1(i,1,1,ie)*z10*y10
                xx(9) = xx(9) + bm1(i,1,1,ie)*z10*z10
            enddo
            vi = 1./vol
            call cmult(xx,vi,9)
        !         call eig3(xx,eign,eig1)
            call eig2(xx,eign,eig1)
            ar(ie) = sqrt(eign/eig1)
        enddo
    else
        do ie=1,nelv
            vol = vlsum(bm1(1,1,1,ie),nxyz)
            x0  = vlsc2(xm1(1,1,1,ie),bm1(1,1,1,ie),nxyz)/vol
            y0  = vlsc2(ym1(1,1,1,ie),bm1(1,1,1,ie),nxyz)/vol
        
            call rzero(xx,4)
            do i=1,nxyz
                x10 = xm1(i,1,1,ie) - x0
                y10 = ym1(i,1,1,ie) - y0
                xx(1) = xx(1) + bm1(i,1,1,ie)*x10*x10
                xx(2) = xx(2) + bm1(i,1,1,ie)*x10*y10
                xx(3) = xx(3) + bm1(i,1,1,ie)*y10*x10
                xx(4) = xx(4) + bm1(i,1,1,ie)*y10*y10
            enddo
            vi = 1./vol
            call cmult(xx,vi,4)
            call eig2(xx,eign,eig1)
        !         write(6,6) ie,vol,eign,eig1
        !  6      format(i5,' veig:',1p3e16.6)
            ar(ie) = sqrt(eign/eig1)
        enddo
    endif

!     do ie=1,nelv
!        write(6,*) ar(ie),ie,' aspect'
!     enddo

    return
    end subroutine aspect_ratios
!-----------------------------------------------------------------------
    subroutine eig2(AA,eign,eig1)

!     return max and min eigenvalues of a 2x2 matrix

    real :: aa(2,2)

    a = aa(1,1)
    b = aa(1,2)
    c = aa(2,1)
    d = aa(2,2)

    qa = 1.
    qb = -(a+d)
    qc = (a*d - c*b)
    call quadratic_h(eign,eig1,qa,qb,qc)

    return
    end subroutine eig2
!-----------------------------------------------------------------------
    subroutine quadratic_h(x1,x2,a,b,c)

!     Stable routine for computation of real roots of quadratic

!     The "_h" denotes the hack below so we don't need to worry
!     about complex arithmetic in the near-double root case. pff 1/22/97

!     Upon return,  | x1 | >= | x2 |

    x1 = 0.
    x2 = 0.

    if (a == 0.) then
        if (b == 0) then
            write(6,10) x1,x2,a,b,c
            return
        endif
        x1 = -c/b
        write(6,11) x1,a,b,c
        return
    endif

    d = b*b - 4.*a*c
    if (d < 0) then
    !        write(6,12) a,b,c,d
    !        hack, for this application we'll assume d<0 by epsilon, and just set
        d = 0
    endif
    if (d > 0) d = sqrt(d)

    q = -0.5 * (b + b/abs(b) * d)
    x1 = q/a
    x2 = c/q

    if (abs(x2) > abs(x1)) then
        tmp = x2
        x2  = x1
        x1  = tmp
    endif

    10 format('ERROR: Both a & b zero in routine quadratic NO ROOTS.' &
    ,1p5e12.4)
    11 format('ERROR: a = 0 in routine quadratic, only one root.' &
    ,1p5e12.4)
    12 format('ERROR: negative discriminate in routine quadratic.' &
    ,1p5e12.4)

    return
    end subroutine quadratic_h
!-----------------------------------------------------------------------
    subroutine out_sem(iel)

    use size_m
    use geom
    include 'INPUT'


    open(unit=33,file='v1')
    if (iel == 0) then
        do ie=1,nelv
            write(33,33) xm1(  1,  1,1,ie),ym1(  1,  1,1,ie)
            write(33,33) xm1(nx1,  1,1,ie),ym1(nx1,  1,1,ie)
            write(33,33) xm1(nx1,ny1,1,ie),ym1(nx1,ny1,1,ie)
            write(33,33) xm1(  1,ny1,1,ie),ym1(  1,ny1,1,ie)
        enddo
    else
        ie = iel
        write(33,33) xm1(  1,  1,1,ie),ym1(  1,  1,1,ie)
        write(33,33) xm1(nx1,  1,1,ie),ym1(nx1,  1,1,ie)
        write(33,33) xm1(nx1,ny1,1,ie),ym1(nx1,ny1,1,ie)
        write(33,33) xm1(  1,ny1,1,ie),ym1(  1,ny1,1,ie)
    endif
    33 format(f14.6)
    close(unit=33)
    return
    end subroutine out_sem
!-----------------------------------------------------------------------
    subroutine gradm11(ux,uy,uz,u,e)

!     Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)

!     Single element case

    use size_m
    use dxyz
    use geom
    include 'INPUT'
    include 'TSTEP'

    parameter (lxyz=lx1*ly1*lz1)
    real :: ux(lxyz),uy(lxyz),uz(lxyz),u(lxyz,1)

    common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)

    integer :: e


    N = nx1-1
    if (if3d) then
        call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
        do i=1,lxyz
            ux(i) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e) &
            + us(i)*sxm1(i,1,1,e) &
            + ut(i)*txm1(i,1,1,e) )
            uy(i) = jacmi(i,e)*(ur(i)*rym1(i,1,1,e) &
            + us(i)*sym1(i,1,1,e) &
            + ut(i)*tym1(i,1,1,e) )
            uz(i) = jacmi(i,e)*(ur(i)*rzm1(i,1,1,e) &
            + us(i)*szm1(i,1,1,e) &
            + ut(i)*tzm1(i,1,1,e) )
        enddo
    else
        call local_grad2(ur,us,u,N,e,dxm1,dxtm1)
        do i=1,lxyz
            ux(i) =jacmi(i,e)*(ur(i)*rxm1(i,1,1,e) &
            + us(i)*sxm1(i,1,1,e) )
            uy(i) =jacmi(i,e)*(ur(i)*rym1(i,1,1,e) &
            + us(i)*sym1(i,1,1,e) )
        enddo
    endif

    return
    end subroutine gradm11
!-----------------------------------------------------------------------
    subroutine gradm11ts(u,ux,uy,uz,e)
!                 T
!     Compute grad  of (ux,uy,uz) and sum to u.

!     Single element case

    use size_m
    include 'TOTAL'

    real :: u(1),ux(1),uy(1),uz(1)
    integer :: e

    parameter (lxyz=lx1*ly1*lz1)
    common /ctmp1/ v(lxyz),w(lxyz)

    if (if3d) then
        call col3    (w,ux,rxm1(1,1,1,e),lxyz)
        call addcol3 (w,uy,rym1(1,1,1,e),lxyz)
        call addcol3 (w,uz,rzm1(1,1,1,e),lxyz)
        call col2    (w,w3m1,lxyz)
        call mxm     (dxtm1,lx1,w,lx1,v,ly1*lz1)
        call add2    (u,v,lxyz)
    
        call col3    (w,ux,sxm1(1,1,1,e),lxyz)
        call addcol3 (w,uy,sym1(1,1,1,e),lxyz)
        call addcol3 (w,uz,szm1(1,1,1,e),lxyz)
        call col2    (w,w3m1,lxyz)
        l=1
        do k=1,lz1
            call mxm  (w(l),lx1,dym1,ly1,v(l),ly1)
            l = l+lx1*ly1
        enddo
        call add2    (u,v,lxyz)
    
        call col3    (w,ux,txm1(1,1,1,e),lxyz)
        call addcol3 (w,uy,tym1(1,1,1,e),lxyz)
        call addcol3 (w,uz,tzm1(1,1,1,e),lxyz)
        call col2    (w,w3m1,lxyz)
        call mxm     (w,lx1*ly1,dzm1,lz1,v,lz1)
        call add2    (u,v,lxyz)
    else
        call col3    (w,ux,rxm1(1,1,1,e),lxyz)
        call addcol3 (w,uy,rym1(1,1,1,e),lxyz)
        call col2    (w,jacmi(1,e),lxyz)
        call mxm     (dxtm1,lx1,w,lx1,v,ly1)
        call add2    (u,v,lxyz)
    
        call col3    (w,ux,sxm1(1,1,1,e),lxyz)
        call addcol3 (w,uy,sym1(1,1,1,e),lxyz)
        call col2    (w,jacmi(1,e),lxyz)
        call mxm     (w,lx1,dym1,ly1,v,ly1)
        call add2    (u,v,lxyz)
    endif

    return
    end subroutine gradm11ts
!-----------------------------------------------------------------------
    subroutine makemsf(afx,afy,afz)

!     This updates the RHS for the modified stress tensor.

!     NOTE:  We assume that vdiff() contains *alpha* times the desired
!            diffusion coefficient, because this is what is required on
!            the lhs of the equation.

!            For the incompressible NS (as in Julie Mullen Thesis, 99),
!            alpha = 2.   For the anelastic case, alpha=4/3.

    Include 'SIZE'
    include 'TOTAL'

    parameter(lxyz = lx1*ly1*lz1)
    real ::  afx(lxyz,lelv),afy(lxyz,lelv),afz(lxyz,lelv)

    common /scrns/ w(lxyz,3,3)
    common /scruz/ v(lxyz,3,3)

    integer :: e

    one = 1.
    two = 2.
    thr = 3.

    c1  = -one/3.
    c2  =  two/3.


    do e=1,nelv
    
        call gradm11(w(1,1,1),w(1,2,1),w(1,3,1),vx,e)
        call gradm11(w(1,1,2),w(1,2,2),w(1,3,2),vy,e)
    
        if (if3d) then
            call gradm11(w(1,1,3),w(1,2,3),w(1,3,3),vz,e)
            do i=1,lxyz
                tmp      = vdiff(i,1,1,e,1)*bm1(i,1,1,e)
                q        = w(i,1,1) + w(i,2,2) + w(i,3,3)
            
                v(i,1,1) = tmp * c2*q
                v(i,2,1) = tmp * (w(i,2,1)-w(i,1,2))
                v(i,3,1) = tmp * (w(i,3,1)-w(i,1,3))
            
                v(i,1,2) = tmp * (w(i,1,2)-w(i,2,1))
                v(i,2,2) = tmp * c2*q
                v(i,3,2) = tmp * (w(i,3,2)-w(i,2,3))
            
                v(i,1,3) = tmp * (w(i,1,3)-w(i,3,1))
                v(i,2,3) = tmp * (w(i,2,3)-w(i,3,2))
                v(i,3,3) = tmp * c2*q
            enddo
            call gradm11ts(afz(1,e),v(1,1,3),v(1,2,3),v(1,3,3),e)
        
        else
            do i=1,lxyz
                tmp      = vdiff(i,1,1,e,1)*bm1(i,1,1,e)
                q        = w(i,1,1) + w(i,2,2)
                v(i,1,1) = tmp * c2*q
                v(i,2,1) = tmp * (w(i,2,1)-w(i,1,2))
                v(i,1,2) = tmp * (w(i,1,2)-w(i,2,1))
                v(i,2,2) = tmp * c2*q
            !              write(6,1) e,i,q,w(i,1,1),w(i,2,2)
            !    1         format(2i5,1p3e15.7)
            enddo
        endif
    
        call gradm11ts(afy(1,e),v(1,1,2),v(1,2,2),v(1,3,2),e)
        call gradm11ts(afx(1,e),v(1,1,1),v(1,2,1),v(1,3,1),e)
    
    enddo

    return
    end subroutine makemsf
!-----------------------------------------------------------------------
