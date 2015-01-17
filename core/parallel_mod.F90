!> cleaned
module parallel
!     Communication information
!     NOTE: NID is stored in 'SIZE' for greater accessibility
  use size_m, only : nid
  implicit none

  INTEGER ::        NODE,PID,NP,NULLPID,NODE0

!     Maximum number of elements (limited to 2**31/12, at least for now)
  integer, parameter :: NELGT_MAX = 178956970

  integer, allocatable :: nelg(:), lglel(:)
  integer :: nvtot, nelgv, nelgt

  LOGICAL :: IFGPRNT
  INTEGER :: WDSIZE,ISIZE,LSIZE,CSIZE,WDSIZI
  LOGICAL :: IFDBLAS

!     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solve)

  integer :: cr_h, gsh
  integer, allocatable :: gsh_fld(:), xxth(:)

  integer :: nekcomm, nekgroup, nekreal

  ! map information
  integer :: queue_dim(30)
  integer :: queue_fac(30)
  integer :: queue_div(30)
  integer :: num_queue
  integer :: proc_shape(3)


  contains

  subroutine init_parallel()
    use size_m
    implicit none

    allocate(NELG(0:LDIMT1), LGLEL(LELT))
    allocate(gsh_fld(0:ldimt3), xxth(ldimt3))

  end subroutine init_parallel

  subroutine init_gllnid()
    use mesh, only : shape_x
    implicit none

    integer :: i
    integer :: np_targ
    integer :: my_shape(3)
    integer :: factors(30), num_fac
    integer :: largest_idx

    factors = -1
    np_targ = np
    num_fac = 0
    do while (np_targ > 1)
      do i = 2, np_targ
        if ((np_targ / i) * i == np_targ) then
          num_fac = num_fac + 1
          factors(num_fac) = i
          np_targ = np_targ / i
          exit
        endif
      enddo
    enddo

    my_shape = shape_x
    num_queue = 0
    do while (num_fac > 0) 
      largest_idx = 1
      if (my_shape(2) > my_shape(largest_idx)) largest_idx = 2
      if (my_shape(3) > my_shape(largest_idx)) largest_idx = 3

      if ((my_shape(largest_idx) / factors(num_fac)) * factors(num_fac) /= my_shape(largest_idx)) then
          if (nid == 0) write(*,*) "Largest dimension isn't divisible by largest factor"
      endif

      my_shape(largest_idx) = my_shape(largest_idx) / factors(num_fac)
      num_queue = num_queue + 1
      queue_dim(num_queue) = largest_idx
      queue_fac(num_queue) = factors(num_fac)
      queue_div(num_queue) = my_shape(largest_idx)
      num_fac = num_fac - 1
    enddo
    proc_shape = my_shape

  end subroutine init_gllnid

  integer function gllnid(ieg)
    use mesh, only : ieg_to_xyz
    implicit none
    integer, intent(in) :: ieg

    integer :: queue_pos
    integer :: ix(3) 

    ix = ieg_to_xyz(ieg)
    gllnid = 0
    do queue_pos = 1, num_queue
      gllnid = queue_fac(queue_pos) * gllnid + ix(queue_dim(queue_pos)) / queue_div(queue_pos)
      ix(queue_dim(queue_pos)) = mod(ix(queue_dim(queue_pos)), queue_div(queue_pos))
    enddo

    return
  end function gllnid

  integer function gllel(ieg)
    use mesh, only : ieg_to_xyz
    implicit none
    integer, intent(in) :: ieg

    integer :: ix(3) 
    ix = ieg_to_xyz(ieg)
    ix(1) = mod(ix(1), proc_shape(1))
    ix(2) = mod(ix(2), proc_shape(2))
    ix(3) = mod(ix(3), proc_shape(3))
    gllel = 1 + ix(1) + ix(2)*proc_shape(1) + ix(3)*proc_shape(1)*proc_shape(2)
    return
  end function gllel

end module parallel

