!> cleaned
module parallel
!     Communication information
!     NOTE: NID is stored in 'SIZE' for greater accessibility
  use size_m, only : nid
  implicit none

  INTEGER ::        NODE,PID,NP,NULLPID,NODE0

!     Maximum number of elements (limited to 2**31/12, at least for now)
  integer, parameter :: NELGT_MAX = 178956970

  integer, allocatable :: nelg(:) 
  integer :: nvtot, nelgv, nelgt

  LOGICAL :: IFGPRNT
  INTEGER :: WDSIZE,ISIZE,LSIZE,CSIZE,WDSIZI
  LOGICAL :: IFDBLAS

!     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solve)

  integer :: cr_h, gsh
  integer, allocatable :: gsh_fld(:), xxth(:)

  integer :: nekcomm, nekgroup, nekreal

  ! map information
  integer, private :: queue_dim(30)
  integer, private :: queue_fac(30)
  integer, private :: queue_div(30)
  integer, private :: num_queue
  integer :: proc_pos(3) 
  integer :: proc_shape(3)


  contains

  subroutine init_parallel()
    use size_m
    implicit none

    allocate(NELG(0:LDIMT1))
    allocate(gsh_fld(0:ldimt3), xxth(ldimt3))

  end subroutine init_parallel

  subroutine init_gllnid()
    use mesh, only : shape_x
    use size_m, only : lelt
    implicit none

    integer :: i, l, iel, ieg
    integer :: np_targ
    integer :: my_shape(3), my_nid
    integer :: factors(30), num_fac
    integer :: largest_idx
    integer :: queue_pos

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
      largest_idx = 3
      if (my_shape(2) > my_shape(largest_idx)) largest_idx = 2
      if (my_shape(1) > my_shape(largest_idx)) largest_idx = 1

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

    proc_pos = 0
    my_nid = nid
    do queue_pos = num_queue, 1, -1
      l = queue_dim(queue_pos) 
      proc_pos(l) = proc_pos(l) + mod(my_nid, queue_fac(queue_pos)) * my_shape(l)
      my_shape(l) = my_shape(l) * queue_fac(queue_pos)
      my_nid = my_nid / queue_fac(queue_pos)
    enddo

    do iel = 1, lelt
      ieg = lglel(iel)
      if (gllnid(ieg) /= nid .or. gllel(ieg) /= iel) then
        write(*,*) "LGL/GLL mismatch", nid, gllnid(ieg), iel, gllel(ieg) 
      endif
    enddo
    if (nid == 0) write(*,*) "LGL/GLL checks out"

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

  integer function lglel(iel)
    use mesh, only : xyz_to_ieg
    implicit none
    integer, intent(in) :: iel

    integer :: my_pos(3)

    my_pos(1) = proc_pos(1) + mod((iel - 1),                               proc_shape(1))
    my_pos(2) = proc_pos(2) + mod((iel - 1)/(proc_shape(1)),               proc_shape(2))
    my_pos(3) = proc_pos(3) + mod((iel - 1)/(proc_shape(1)*proc_shape(2)), proc_shape(3))

    lglel = xyz_to_ieg(my_pos)
    return

  end function lglel

end module parallel

