!-----------------------------------------------------------------------
    subroutine calcz(d,e,n,dmax,dmin,z,ierr)

!     Num. Rec. 2nd Ed., p.473

!     Note:  d(1:n) is the diagonal of the sym. tridiagonal matrix
!            e(1:n) is the upper diagonal of the tridiagonal matrix,
!                   WITH e(n) ARBITRARY (a slight shift from Num.Rec.)

!            z(n:n) is the packed array of eigenvectors

    real ::  d(n),e(n),z(n,n)
    real ::  smalln,small

    call ident(z,n)
    one = 1.

!     Find smallest number  (pff mod to Num. Rec. 2nd Ed., p.473)

    small = 0.5
    do i = 1,100
        smalln = small * small
        if (smalln == 0) then
            do j=1,1000
                smalln = 0.5*small
                if (smalln == 0) goto 10
                small = smalln
            enddo
        endif
        small = smalln
    enddo
    10 continue
    small = 10.*small
    small = max(small,1e-99)
          
!     write(6,*) 'this is small:',small

    do 15 l=1,n
        iter = 0
    
        1 do 12 m=l,n-1
            dd = abs( d(m) ) + abs( d(m+1) )
            de = e(m) + dd
            df = abs(dd - de)
        !           write(6,112) iter,m,'de:',dd,de,df,small
            if ( df <= small ) goto 2
        12 END DO
        112 format(i3,i9,1x,a3,1p4e16.8)
        m = n
    
        2 if ( m /= l ) then
            if ( iter == 600 ) then
                write (6,*) 'too many iterations in calc'
            !              n10 = min(n,10)
            !              do i=1,n
            !                 write(6,9) d(i),(z(i,j),j=1,n10)
            !              enddo
            !   9          format(1pe12.4,' e:',1p10e12.4)
            !              call exitt
                ierr=1
                return
            endif
        
            iter = iter + 1
            g = ( d(l+1) - d(l) ) / ( 2.0 * e(l) )
            r = pythag(g,one)
            g = d(m) - d(l) + e(l)/(g+sign(r,g))
            s = 1.0
            c = 1.0
            p = 0.0
        
            do 14 i = m-1,l,-1
                f = s * e(i)
                b = c * e(i)
                r = pythag(f,g)
                e(i+1) = r
                if ( abs(r) <= small ) then
                    d(i+1) = d(i+1) - p
                    e(m)   = 0.
                    goto 1
                endif
                s = f/r
                c = g/r
                g = d(i+1) - p
                r = ( d(i)-g )*s + 2.*c*b
                p = s*r
                d(i+1) = g + p
                g = c*r - b
            !      ...     find eigenvectors ... (new, 11/19/94, pff, p.363. Num.Rec.I.)
                do 13 k=1,n
                    f = z(k,i+1)
                    z(k,i+1) = s*z(k,i)+c*f
                    z(k,i  ) = c*z(k,i)-s*f
                13 END DO
            !      ...     end of eigenvector section ...
            14 END DO
        
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0
            goto 1
        endif
    
    15 END DO

!     write (8,8) (d(j),j=1,n)
!   8 format('eig:',8f10.4)

    dmin = d(1)
    dmax = d(1)
    do 40 i = 1 , n
        dmin = min( d(i) , dmin )
        dmax = max( d(i) , dmax )
    40 END DO

!     Output eigenvectors

!     n10 = min(n,10)
!     do i=1,n
!        write(6,9) d(i),(z(i,j),j=1,n10)
!     enddo
!   9 format(1pe12.4,' e:',1p10e12.4)

    ierr=0
    return
    end subroutine calcz
!-----------------------------------------------------------------------
    function pythag(a,b)

!     Compute sqrt(a*a + b*b) w/o under- or over-flow.

    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
        pythag = absa*sqrt(1. + (absb/absa)**2 )
    else
        if (absb == 0.) then
            pythag = 0.
        else
            pythag = absb*sqrt(1. + (absa/absb)**2 )
        endif
    endif
    return
    end function pythag
!-----------------------------------------------------------------------
    subroutine ident(a,n)
    real ::  a(n,n)
    call rzero(a,n*n)
    do i=1,n
        a(i,i) = 1.0
    enddo
    return
    end subroutine ident
!-----------------------------------------------------------------------
