!> \brief Arrays for direct stiffness summation.
!! cleaned
module topol
  use size_m
  implicit none

  integer :: NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6), EFACE(6),EFACE1(6)

  integer :: ESKIP(-12:12,3), NEDG(3), NCMP &
    ,IXCN(8),NOFFST(3,0:LDIMT1) &
    ,MAXMLT,NSPMAX(0:LDIMT1) &
    ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1), NGCOMM(2,0:LDIMT1)

  integer, allocatable :: NUMSCN(:,:),NUMSED(:,:) 
  integer, allocatable :: GCNNUM(:,:,:), LCNNUM(:,:,:)
  integer, allocatable :: GEDNUM(:,:,:), LEDNUM(:,:,:)
  integer, allocatable :: GEDTYP(:,:,:) 

  integer :: IEDGE(20),IEDGEF(2,4,6,0:1) &
    ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10) &
    ,INDX(8),INVEDG(27)

    DATA    IEDGFC /  5,7,9,11,  6,8,10,12, &
    &                   1,3,9,10,  2,4,11,12, &
    &                   1,2,5,6,   3,4,7,8    /
    DATA    ICEDG / 1,2,1,   3,4,1,   5,6,1,   7,8,1, &
    &                 1,3,2,   2,4,2,   5,7,2,   6,8,2, &
    &                 1,5,3,   2,6,3,   3,7,3,   4,8,3, &
!      -2D-
    &                 1,2,1,   3,4,1,   1,3,2,   2,4,2 /
    DATA    ICFACE/ 1,3,5,7, 2,4,6,8, &
    &                 1,2,5,6, 3,4,7,8, &
    &                 1,2,3,4, 5,6,7,8, &
!      -2D-
    &                 1,3,0,0, 2,4,0,0, &
    &                 1,2,0,0, 3,4,0,0  /



  contains

  subroutine init_topol()
    use size_m
    implicit none
#if 0
     allocate(NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1) &
    ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1) &
    ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1) &
    ,GEDTYP(12,LELT,0:LDIMT1))
#endif

  end subroutine init_topol

end module topol
