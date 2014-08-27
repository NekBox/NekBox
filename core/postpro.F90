!-----------------------------------------------------------------------
    subroutine gen_int_gz(j,jt,g,n,z,m)

!     Generate interpolater from m z points to n g points

!        j   = interpolation matrix, mapping from z to g
!        jt  = transpose of interpolation matrix
!        m   = number of points on z grid
!        n   = number of points on g grid

    real :: j(n,m),jt(m,n),g(n),z(m)

    mpoly  = m-1
    do i=1,n
        call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
    enddo

    call transpose(j,n,jt,m)

    return
    end subroutine gen_int_gz
!-----------------------------------------------------------------------
    subroutine zuni(z,np)

!     Generate equaly spaced np points on the interval [-1:1]

    real :: z(1)

    dz = 2./(np-1)
    z(1) = -1.
    do i = 2,np-1
        z(i) = z(i-1) + dz
    enddo
    z(np) = 1.

    return
    end subroutine zuni
!-----------------------------------------------------------------------
