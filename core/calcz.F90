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
