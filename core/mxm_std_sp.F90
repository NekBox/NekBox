!>     unrolled loop version
subroutine mxmf2_sp(A,N1,B,N2,C,N3)
  use kinds, only : SP
  integer :: n1, n2, n3
  real(SP) :: a(n1,n2),b(n2,n3),c(n1,n3)

  select case (n2)
    case (1 : 8)
      select case (n2)
        case (8)
          call mxf8_sp(a,n1,b,n2,c,n3)
        case (1)
          call mxf1_sp(a,n1,b,n2,c,n3)
        case (2)
          call mxf2_sp(a,n1,b,n2,c,n3)
        case (3)
          call mxf3_sp(a,n1,b,n2,c,n3)
        case (4)
          call mxf4_sp(a,n1,b,n2,c,n3)
        case (5)
          call mxf5_sp(a,n1,b,n2,c,n3)
        case (6)
          call mxf6_sp(a,n1,b,n2,c,n3)
        case (7)
          call mxf7_sp(a,n1,b,n2,c,n3)
      end select
    case (9 : 16)
      select case (n2)
        case (12)
          call mxf12_sp(a,n1,b,n2,c,n3)
        case (9)
          call mxf9_sp(a,n1,b,n2,c,n3)
        case (10)
          call mxf10_sp(a,n1,b,n2,c,n3)
        case (11)
          call mxf11_sp(a,n1,b,n2,c,n3)
        case (13)
          call mxf13_sp(a,n1,b,n2,c,n3)
        case (14)
          call mxf14_sp(a,n1,b,n2,c,n3)
        case (15)
          call mxf15_sp(a,n1,b,n2,c,n3)
        case (16)
          call mxf16_sp(a,n1,b,n2,c,n3)
      end select
    case (17 : 24)
      select case (n2)
        case (17)
          call mxf17_sp(a,n1,b,n2,c,n3)
        case (18)
          call mxf18_sp(a,n1,b,n2,c,n3)
        case (19)
          call mxf19_sp(a,n1,b,n2,c,n3)
        case (20)
          call mxf20_sp(a,n1,b,n2,c,n3)
        case (21)
          call mxf21_sp(a,n1,b,n2,c,n3)
        case (22)
          call mxf22_sp(a,n1,b,n2,c,n3)
        case (23)
          call mxf23_sp(a,n1,b,n2,c,n3)
        case (24)
          call mxf24_sp(a,n1,b,n2,c,n3)
      end select
    case default 
      call mxm44_0_sp(a,n1,b,n2,c,n3)
  end select 

  return
end subroutine mxmf2_sp

!-----------------------------------------------------------------------
    subroutine mxf1_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,1),b(1,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
        enddo
    enddo
    return
    end subroutine mxf1_sp
!-----------------------------------------------------------------------
    subroutine mxf2_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,2),b(2,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j)
        enddo
    enddo
    return
    end subroutine mxf2_sp
!-----------------------------------------------------------------------
    subroutine mxf3_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,3),b(3,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j)
        enddo
    enddo
    return
    end subroutine mxf3_sp
!-----------------------------------------------------------------------
    subroutine mxf4_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,4),b(4,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j)
        enddo
    enddo
    return
    end subroutine mxf4_sp
!-----------------------------------------------------------------------
    subroutine mxf5_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,5),b(5,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j)
        enddo
    enddo
    return
    end subroutine mxf5_sp
!-----------------------------------------------------------------------
    subroutine mxf6_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,6),b(6,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j)
        enddo
    enddo
    return
    end subroutine mxf6_sp
!-----------------------------------------------------------------------
    subroutine mxf7_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,7),b(7,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j)
        enddo
    enddo
    return
    end subroutine mxf7_sp
!-----------------------------------------------------------------------
    subroutine mxf8_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,8),b(8,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j)
        enddo
    enddo
    return
    end subroutine mxf8_sp
!-----------------------------------------------------------------------
    subroutine mxf9_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,9),b(9,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j)
        enddo
    enddo
    return
    end subroutine mxf9_sp
!-----------------------------------------------------------------------
    subroutine mxf10_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,10),b(10,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j)
        enddo
    enddo
    return
    end subroutine mxf10_sp
!-----------------------------------------------------------------------
    subroutine mxf11_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,11),b(11,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j)
        enddo
    enddo
    return
    end subroutine mxf11_sp
!-----------------------------------------------------------------------
    subroutine mxf12_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,12),b(12,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j)
        enddo
    enddo
    return
    end subroutine mxf12_sp
!-----------------------------------------------------------------------
    subroutine mxf13_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,13),b(13,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j)
        enddo
    enddo
    return
    end subroutine mxf13_sp
!-----------------------------------------------------------------------
    subroutine mxf14_sp(a,n1,b,n2,c,n3)
use kinds, only : SP

    real(SP) :: a(n1,14),b(14,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j)
        enddo
    enddo
    return
    end subroutine mxf14_sp
!-----------------------------------------------------------------------
    subroutine mxf15_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,15),b(15,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j)
        enddo
    enddo
    return
    end subroutine mxf15_sp
!-----------------------------------------------------------------------
    subroutine mxf16_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,16),b(16,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
#if 0
          c(i,j) = 0
          do k=1,16
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
#else
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j)
#endif
        enddo
    enddo
    return
    end subroutine mxf16_sp
!-----------------------------------------------------------------------
    subroutine mxf17_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,17),b(17,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j)
        enddo
    enddo
    return
    end subroutine mxf17_sp
!-----------------------------------------------------------------------
    subroutine mxf18_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,18),b(18,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j)
        enddo
    enddo
    return
    end subroutine mxf18_sp
!-----------------------------------------------------------------------
    subroutine mxf19_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,19),b(19,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j)
        enddo
    enddo
    return
    end subroutine mxf19_sp
!-----------------------------------------------------------------------
    subroutine mxf20_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,20),b(20,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j) &
            + a(i,20)*b(20,j)
        enddo
    enddo
    return
    end subroutine mxf20_sp
!-----------------------------------------------------------------------
    subroutine mxf21_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,21),b(21,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j) &
            + a(i,20)*b(20,j) &
            + a(i,21)*b(21,j)
        enddo
    enddo
    return
    end subroutine mxf21_sp
!-----------------------------------------------------------------------
    subroutine mxf22_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,22),b(22,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j) &
            + a(i,20)*b(20,j) &
            + a(i,21)*b(21,j) &
            + a(i,22)*b(22,j)
        enddo
    enddo
    return
    end subroutine mxf22_sp
!-----------------------------------------------------------------------
    subroutine mxf23_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,23),b(23,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j) &
            + a(i,20)*b(20,j) &
            + a(i,21)*b(21,j) &
            + a(i,22)*b(22,j) &
            + a(i,23)*b(23,j)
        enddo
    enddo
    return
    end subroutine mxf23_sp
!-----------------------------------------------------------------------
    subroutine mxf24_sp(a,n1,b,n2,c,n3)

    use kinds, only : SP
    real(SP) :: a(n1,24),b(24,n3),c(n1,n3)

    do j=1,n3
        do i=1,n1
            c(i,j) = a(i,1)*b(1,j) &
            + a(i,2)*b(2,j) &
            + a(i,3)*b(3,j) &
            + a(i,4)*b(4,j) &
            + a(i,5)*b(5,j) &
            + a(i,6)*b(6,j) &
            + a(i,7)*b(7,j) &
            + a(i,8)*b(8,j) &
            + a(i,9)*b(9,j) &
            + a(i,10)*b(10,j) &
            + a(i,11)*b(11,j) &
            + a(i,12)*b(12,j) &
            + a(i,13)*b(13,j) &
            + a(i,14)*b(14,j) &
            + a(i,15)*b(15,j) &
            + a(i,16)*b(16,j) &
            + a(i,17)*b(17,j) &
            + a(i,18)*b(18,j) &
            + a(i,19)*b(19,j) &
            + a(i,20)*b(20,j) &
            + a(i,21)*b(21,j) &
            + a(i,22)*b(22,j) &
            + a(i,23)*b(23,j) &
            + a(i,24)*b(24,j)
        enddo
    enddo
    return
    end subroutine mxf24_sp
!-----------------------------------------------------------------------
    subroutine mxm44_0_sp(a, m, b, k, c, n)

! matrix multiply with a 4x4 pencil

    use kinds, only : SP
    real(SP) :: a(m,k), b(k,n), c(m,n)
    real(SP) :: s11, s12, s13, s14, s21, s22, s23, s24
    real(SP) :: s31, s32, s33, s34, s41, s42, s43, s44

    mresid = iand(m,3)
    nresid = iand(n,3)
    m1 = m - mresid + 1
    n1 = n - nresid + 1

    do i=1,m-mresid,4
        do j=1,n-nresid,4
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            s14 = 0.0d0
            s24 = 0.0d0
            s34 = 0.0d0
            s44 = 0.0d0
            do l=1,k
                s11 = s11 + a(i,l)*b(l,j)
                s12 = s12 + a(i,l)*b(l,j+1)
                s13 = s13 + a(i,l)*b(l,j+2)
                s14 = s14 + a(i,l)*b(l,j+3)

                s21 = s21 + a(i+1,l)*b(l,j)
                s22 = s22 + a(i+1,l)*b(l,j+1)
                s23 = s23 + a(i+1,l)*b(l,j+2)
                s24 = s24 + a(i+1,l)*b(l,j+3)

                s31 = s31 + a(i+2,l)*b(l,j)
                s32 = s32 + a(i+2,l)*b(l,j+1)
                s33 = s33 + a(i+2,l)*b(l,j+2)
                s34 = s34 + a(i+2,l)*b(l,j+3)

                s41 = s41 + a(i+3,l)*b(l,j)
                s42 = s42 + a(i+3,l)*b(l,j+1)
                s43 = s43 + a(i+3,l)*b(l,j+2)
                s44 = s44 + a(i+3,l)*b(l,j+3)
            enddo
            c(i,j)     = s11
            c(i,j+1)   = s12
            c(i,j+2)   = s13
            c(i,j+3)   = s14

            c(i+1,j)   = s21
            c(i+2,j)   = s31
            c(i+3,j)   = s41

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42

            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43

            c(i+1,j+3) = s24
            c(i+2,j+3) = s34
            c(i+3,j+3) = s44
        enddo
    ! Residual when n is not multiple of 4
        if (nresid /= 0) then
            if (nresid == 1) then
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,n)
                    s21 = s21 + a(i+1,l)*b(l,n)
                    s31 = s31 + a(i+2,l)*b(l,n)
                    s41 = s41 + a(i+3,l)*b(l,n)
                enddo
                c(i,n)     = s11
                c(i+1,n)   = s21
                c(i+2,n)   = s31
                c(i+3,n)   = s41
            elseif (nresid == 2) then
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                s12 = 0.0d0
                s22 = 0.0d0
                s32 = 0.0d0
                s42 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,j)
                    s12 = s12 + a(i,l)*b(l,j+1)

                    s21 = s21 + a(i+1,l)*b(l,j)
                    s22 = s22 + a(i+1,l)*b(l,j+1)

                    s31 = s31 + a(i+2,l)*b(l,j)
                    s32 = s32 + a(i+2,l)*b(l,j+1)

                    s41 = s41 + a(i+3,l)*b(l,j)
                    s42 = s42 + a(i+3,l)*b(l,j+1)
                enddo
                c(i,j)     = s11
                c(i,j+1)   = s12

                c(i+1,j)   = s21
                c(i+2,j)   = s31
                c(i+3,j)   = s41

                c(i+1,j+1) = s22
                c(i+2,j+1) = s32
                c(i+3,j+1) = s42
            else
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                s12 = 0.0d0
                s22 = 0.0d0
                s32 = 0.0d0
                s42 = 0.0d0
                s13 = 0.0d0
                s23 = 0.0d0
                s33 = 0.0d0
                s43 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,j)
                    s12 = s12 + a(i,l)*b(l,j+1)
                    s13 = s13 + a(i,l)*b(l,j+2)

                    s21 = s21 + a(i+1,l)*b(l,j)
                    s22 = s22 + a(i+1,l)*b(l,j+1)
                    s23 = s23 + a(i+1,l)*b(l,j+2)

                    s31 = s31 + a(i+2,l)*b(l,j)
                    s32 = s32 + a(i+2,l)*b(l,j+1)
                    s33 = s33 + a(i+2,l)*b(l,j+2)

                    s41 = s41 + a(i+3,l)*b(l,j)
                    s42 = s42 + a(i+3,l)*b(l,j+1)
                    s43 = s43 + a(i+3,l)*b(l,j+2)
                enddo
                c(i,j)     = s11
                c(i+1,j)   = s21
                c(i+2,j)   = s31
                c(i+3,j)   = s41
                c(i,j+1)   = s12
                c(i+1,j+1) = s22
                c(i+2,j+1) = s32
                c(i+3,j+1) = s42
                c(i,j+2)   = s13
                c(i+1,j+2) = s23
                c(i+2,j+2) = s33
                c(i+3,j+2) = s43
            endif
        endif
    enddo

! Residual when m is not multiple of 4
    if (mresid == 0) then
        return
    elseif (mresid == 1) then
        do j=1,n-nresid,4
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,j)
                s12 = s12 + a(m,l)*b(l,j+1)
                s13 = s13 + a(m,l)*b(l,j+2)
                s14 = s14 + a(m,l)*b(l,j+3)
            enddo
            c(m,j)     = s11
            c(m,j+1)   = s12
            c(m,j+2)   = s13
            c(m,j+3)   = s14
        enddo
    ! mresid is 1, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n)
            enddo
            c(m,n) = s11
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s12 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n-1)
                s12 = s12 + a(m,l)*b(l,n)
            enddo
            c(m,n-1) = s11
            c(m,n) = s12
            return
        else
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n-2)
                s12 = s12 + a(m,l)*b(l,n-1)
                s13 = s13 + a(m,l)*b(l,n)
            enddo
            c(m,n-2) = s11
            c(m,n-1) = s12
            c(m,n) = s13
            return
        endif
    elseif (mresid == 2) then
        do j=1,n-nresid,4
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s21 = 0.0d0
            s22 = 0.0d0
            s23 = 0.0d0
            s24 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,j)
                s12 = s12 + a(m-1,l)*b(l,j+1)
                s13 = s13 + a(m-1,l)*b(l,j+2)
                s14 = s14 + a(m-1,l)*b(l,j+3)

                s21 = s21 + a(m,l)*b(l,j)
                s22 = s22 + a(m,l)*b(l,j+1)
                s23 = s23 + a(m,l)*b(l,j+2)
                s24 = s24 + a(m,l)*b(l,j+3)
            enddo
            c(m-1,j)   = s11
            c(m-1,j+1) = s12
            c(m-1,j+2) = s13
            c(m-1,j+3) = s14
            c(m,j)     = s21
            c(m,j+1)   = s22
            c(m,j+2)   = s23
            c(m,j+3)   = s24
        enddo
    ! mresid is 2, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n)
            enddo
            c(m-1,n) = s11
            c(m,n) = s21
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n-1)
                s12 = s12 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n-1)
                s22 = s22 + a(m,l)*b(l,n)
            enddo
            c(m-1,n-1) = s11
            c(m-1,n)   = s12
            c(m,n-1)   = s21
            c(m,n)     = s22
            return
        else
            s11 = 0.0d0
            s21 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n-2)
                s12 = s12 + a(m-1,l)*b(l,n-1)
                s13 = s13 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n-2)
                s22 = s22 + a(m,l)*b(l,n-1)
                s23 = s23 + a(m,l)*b(l,n)
            enddo
            c(m-1,n-2) = s11
            c(m-1,n-1) = s12
            c(m-1,n)   = s13
            c(m,n-2)   = s21
            c(m,n-1)   = s22
            c(m,n)     = s23
            return
        endif
    else
    ! mresid is 3
        do j=1,n-nresid,4
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0

            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0

            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0

            s14 = 0.0d0
            s24 = 0.0d0
            s34 = 0.0d0

            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,j)
                s12 = s12 + a(m-2,l)*b(l,j+1)
                s13 = s13 + a(m-2,l)*b(l,j+2)
                s14 = s14 + a(m-2,l)*b(l,j+3)

                s21 = s21 + a(m-1,l)*b(l,j)
                s22 = s22 + a(m-1,l)*b(l,j+1)
                s23 = s23 + a(m-1,l)*b(l,j+2)
                s24 = s24 + a(m-1,l)*b(l,j+3)

                s31 = s31 + a(m,l)*b(l,j)
                s32 = s32 + a(m,l)*b(l,j+1)
                s33 = s33 + a(m,l)*b(l,j+2)
                s34 = s34 + a(m,l)*b(l,j+3)
            enddo
            c(m-2,j)   = s11
            c(m-2,j+1) = s12
            c(m-2,j+2) = s13
            c(m-2,j+3) = s14

            c(m-1,j)   = s21
            c(m-1,j+1) = s22
            c(m-1,j+2) = s23
            c(m-1,j+3) = s24

            c(m,j)     = s31
            c(m,j+1)   = s32
            c(m,j+2)   = s33
            c(m,j+3)   = s34
        enddo
    ! mresid is 3, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n)
            enddo
            c(m-2,n) = s11
            c(m-1,n) = s21
            c(m,n) = s31
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n-1)
                s12 = s12 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n-1)
                s22 = s22 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n-1)
                s32 = s32 + a(m,l)*b(l,n)
            enddo
            c(m-2,n-1) = s11
            c(m-2,n)   = s12
            c(m-1,n-1) = s21
            c(m-1,n)   = s22
            c(m,n-1)   = s31
            c(m,n)     = s32
            return
        else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n-2)
                s12 = s12 + a(m-2,l)*b(l,n-1)
                s13 = s13 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n-2)
                s22 = s22 + a(m-1,l)*b(l,n-1)
                s23 = s23 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n-2)
                s32 = s32 + a(m,l)*b(l,n-1)
                s33 = s33 + a(m,l)*b(l,n)
            enddo
            c(m-2,n-2) = s11
            c(m-2,n-1) = s12
            c(m-2,n)   = s13
            c(m-1,n-2) = s21
            c(m-1,n-1) = s22
            c(m-1,n)   = s23
            c(m,n-2)   = s31
            c(m,n-1)   = s32
            c(m,n)     = s33
            return
        endif
    endif

    return
    end subroutine mxm44_0_sp
!-----------------------------------------------------------------------
