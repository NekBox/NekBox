!==============================================================================
!> \file poisson_mod.F90
!! \brief Spectral coarse solver for poisson equation 
!! \date November 2014
!! \author Max Hutchinson
!!
!! This module implements a coarse solve (preconditioner) for the pressure
!! Poisson equation.
module poisson
  use, intrinsic :: iso_c_binding
  use kinds, only : DP
  implicit none

  public spectral_solve
  private

  integer :: comm_xy, comm_yz
  logical :: interface_initialized = .false.
  logical :: mesh_to_grid_initialized = .false.

  integer :: alloc_local_xy, nin_local_xy, nout_local_xy, idx_in_local_xy, idx_out_local_xy
  integer :: alloc_local_yz, nin_local_yz, nout_local_yz, idx_in_local_yz, idx_out_local_yz

  type real_p
    real(DP), allocatable :: p(:)
  end type real_p
  type int_p
    integer, allocatable :: p(:)
  end type int_p

  type(real_p), allocatable :: send_buffers(:)
  type(real_p), allocatable :: rec_buffers(:)
  integer, allocatable :: dest_pids(:)
  integer, allocatable :: dest_slots(:)
  integer, allocatable :: dest_indexes(:)
  integer, allocatable :: dest_lengths(:)

  integer, allocatable :: src_pids(:)
  integer, allocatable :: src_lengths(:)
  integer, allocatable :: src_slots(:,:,:)
  integer, allocatable :: src_indexes(:,:,:)

  integer :: comm_size

contains

!> \brief 
subroutine spectral_solve(u,rhs)!,h1,mask,mult,imsh,isd)
  use kinds, only : DP
  use geom, only : bm1, binvm1, volvm1
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nekcomm, nid, lglel
  use soln, only : vmult
  use size_m, only : nx1, ny1, nz1, nelv
  use ctimer, only : dnekclock

  use fft, only : P_FORWARD, P_BACKWARD, W_FORWARD, W_BACKWARD
  use fft, only : fft_r2r, transpose_grid

  REAL(DP), intent(out)   :: U    (:)
  REAL(DP), intent(inout) :: RHS  (:)

  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  real(DP), allocatable :: tmp_fine(:,:,:,:)
  integer :: nelm
  integer :: i, j
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: rescale
  real(DP) :: h2(1,1,1,1)
  integer :: ix(3)
  real(DP), save :: thistime, tottime = 0._dp

  thistime = - dnekclock()

  nelm = size(rhs) / 8

  if (.not. interface_initialized) then
    call init_comm_infrastructure(nekcomm, shape_x)
  endif

  ! map onto fine mesh to use dssum
  !> \todo replace this with coarse grid dssum  
  allocate(tmp_fine(nx1, ny1, nz1, nelm))
  tmp_fine = 0._dp
  forall (i = 1: nelm)
    tmp_fine(1,  1,  1,   i) = rhs(1 + (i-1)*8)
    tmp_fine(nx1,1,  1,   i) = rhs(2 + (i-1)*8)
    tmp_fine(1,  ny1,1,   i) = rhs(3 + (i-1)*8)
    tmp_fine(nx1,ny1,1,   i) = rhs(4 + (i-1)*8)
    tmp_fine(1,  1,  nz1, i) = rhs(5 + (i-1)*8)
    tmp_fine(nx1,1,  nz1, i) = rhs(6 + (i-1)*8)
    tmp_fine(1,  ny1,nz1, i) = rhs(7 + (i-1)*8)
    tmp_fine(nx1,ny1,nz1, i) = rhs(8 + (i-1)*8)
  end forall
  call dssum(tmp_fine)

  ! convert RHS to coarse mesh
  allocate(rhs_coarse(nelm))
  forall(i = 1 : nelm) rhs_coarse(i) = (tmp_fine(1,1,1,i) + tmp_fine(1,1,nz1,i))/2._dp
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_FORWARD, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_FORWARD, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_FORWARD, rescale)

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_BACKWARD, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_BACKWARD, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_BACKWARD, rescale)

  ! normalize the FFTs
  rescale = 1._dp / (rescale * sum(bm1(:,:,:,1)))
  plane_xy = plane_xy * rescale

  ! reorder to local elements
  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)

  ! populate U
  tmp_fine = 0._dp
  forall (i = 1: nelm)
    tmp_fine(1,  1,  1,   i)   = 4.*soln_coarse(i) 
    tmp_fine(1,  1,  nz1,   i) = 4.*soln_coarse(i) 
  end forall
  call dssum(tmp_fine)
  tmp_fine = tmp_fine * vmult

  ! extract coarse values
  u = 0._dp
  forall (i = 1: nelm)
    u(1+(i-1)*8) = tmp_fine(1,  1,  1,   i) 
    u(2+(i-1)*8) = tmp_fine(nx1,1,  1,   i) 
    u(3+(i-1)*8) = tmp_fine(1,  ny1,1,   i) 
    u(4+(i-1)*8) = tmp_fine(nx1,ny1,1,   i) 
    u(5+(i-1)*8) = tmp_fine(1,  1,  nz1, i) 
    u(6+(i-1)*8) = tmp_fine(nx1,1,  nz1, i) 
    u(7+(i-1)*8) = tmp_fine(1,  ny1,nz1, i) 
    u(8+(i-1)*8) = tmp_fine(nx1,ny1,nz1, i) 
  end forall

  thistime = thistime + dnekclock()
  tottime = tottime + thistime
  if (nid == 0) write(*,*) "SCPS timers:", thistime, tottime

  return
 
end subroutine spectral_solve

!> \brief one-time setup of communication infrastructure for poisson_mod
subroutine init_comm_infrastructure(comm_world, shape_x)
  use fftw3, only : FFTW_MPI_DEFAULT_BLOCK
  use fftw3, only : fftw_mpi_local_size_many_transposed
  use fftw3, only : fftw_mpi_init
  use mpif, only : MPI_UNDEFINED
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: shape_x(3) !>!< Shape of mesh

  integer(C_INTPTR_T) :: shape_c(3)
  integer(C_INTPTR_T), parameter :: one = 1
  integer :: nxy, nyz, ixy, iyz
  integer :: nid, ierr, i

  call MPI_Comm_rank(comm_world, nid, ierr) 
  call MPI_Comm_size(comm_world, comm_size, ierr) 

  !call fftw_mpi_init()

  comm_size = min(comm_size, 4096)
  !comm_size = min(comm_size, 1024)
  !comm_size = min(comm_size, 256)
  !comm_size = min(comm_size, 64)

  nxy =  int(2**int(log(real(comm_size))/log(2.) / 2))
  comm_size = nxy*nxy
  if (nid < comm_size) then
    ixy = ((nxy*nid/comm_size) * shape_x(3)) / nxy
  else
    ixy = comm_size + 1
  endif
  call MPI_Comm_split(comm_world, ixy, 0, comm_xy, ierr)
  if (ierr /= 0) write(*,*) "Comm split xy failed", nid
  if (nid >= comm_size) then
    call MPI_Comm_free(comm_xy, ierr)
    if (ierr /= 0) write(*,*) "Comm free xy failed", nid
  endif

  nyz =  comm_size/nxy
  if (nid < comm_size) then
    iyz = (mod(nid,nyz) * shape_x(2)) / nyz
  else
    iyz = comm_size + 1
  endif
  call MPI_Comm_split(comm_world, iyz, 0, comm_yz, ierr)
  if (ierr /= 0) write(*,*) "Comm split yz failed", nid
  if (nid >= comm_size) then
    call MPI_Comm_free(comm_yz, ierr)
    if (ierr /= 0) write(*,*) "Comm free yz failed", nid
  endif

  if (nid < comm_size) then
    shape_c = shape_x
#if 1
    idx_in_local_yz = ixy; idx_out_local_yz = ixy
    idx_in_local_xy = iyz; idx_out_local_xy = iyz
    nin_local_yz = shape_x(3) / nxy; nout_local_yz = nin_local_yz
    nin_local_xy = shape_x(2) / nyz; nout_local_xy = nin_local_xy
#else
    alloc_local_xy = fftw_mpi_local_size_many_transposed( 2, &
                  (/shape_c(2), shape_c(1)/), &
                  one, shape_c(2)/nxy, shape_c(1)/nxy, &
                  comm_xy, &
                  nin_local_xy, idx_in_local_xy, nout_local_xy, idx_out_local_xy)
    alloc_local_yz = fftw_mpi_local_size_many_transposed( 2, &
                  (/shape_c(3), shape_c(2)/), &
                  one, shape_c(3)/nyz, shape_c(2)/nyz, &
                  comm_yz, &
                  nin_local_yz, idx_in_local_yz, nout_local_yz, idx_out_local_yz)
#endif

    if (ixy /= idx_in_local_yz .or. ixy /= idx_out_local_yz) write(*,*) "fail 1", nid, ixy, idx_in_local_yz, idx_out_local_yz
    if (iyz /= idx_in_local_xy .or. iyz /= idx_out_local_xy) write(*,*) "fail 2", nid, iyz, idx_in_local_xy, idx_out_local_xy
    !if (nxy /= 8 .or. nyz /= 8) write(*,*) "fail 3", nid, nxy, nyz
    !if (nxy /= 16 .or. nyz /= 16) write(*,*) "fail 3", nid, nxy, nyz
    !if (nxy /= 32 .or. nyz /= 32) write(*,*) "fail 3", nid, nxy, nyz
    !if (nxy /= 64 .or. nyz /= 64) write(*,*) "fail 3", nid, nxy, nyz
    !if (shape_c(1) /= alloc_local_xy .or. shape_c(2) /= alloc_local_yz) write(*,*) "fail 3", nid
  else
    nin_local_xy = 0; nout_local_xy = 0
    nin_local_yz = 0; nout_local_yz = 0
    idx_in_local_xy = -1; idx_out_local_xy = -1
    idx_in_local_yz = -1; idx_out_local_yz = -1
    alloc_local_xy = 0; alloc_local_yz = 0
  endif

#if 1
  do i = 0, 65
    call nekgsync()
    if (nid == i) write(*,'(A,15(I5))') "MAX:", nid, alloc_local_xy, nin_local_xy, nout_local_xy, &
      alloc_local_yz, nin_local_yz, nout_local_yz, &
      idx_in_local_xy, idx_out_local_xy, idx_in_local_yz, idx_out_local_yz, &
      nxy, nyz, ixy, iyz
    !if (nid == i) write(*,'(A,6(I5))') "MAX:", nid, alloc_local_xy, nin_local_xy, nout_local_xy, idx_in_local_xy, idx_out_local_xy
  enddo
  call nekgsync()
!  write(*,'(A,6(I5))') "MAX:", nid, alloc_local_yz, nin_local_yz, idx_in_local_yz, nout_local_yz, idx_out_local_yz
#endif

  call nekgsync()

  call transpose_test()
  call shuffle_test()
  call cos_test()

  interface_initialized = .true.

end subroutine init_comm_infrastructure

function ieg_to_xyz(ieg, shape_x) result(xyz)
  integer, intent(in) :: ieg
  integer, intent(in) :: shape_x(3)
  integer :: xyz(3)

  xyz(1) = mod(ieg - 1, shape_x(1))
  xyz(2) = mod((ieg-1)/shape_x(1), shape_x(2))
  xyz(3) = mod((ieg-1)/(shape_x(1)*shape_x(2)), shape_x(3))

end function ieg_to_xyz

integer function xyz_to_pid(ix, iy, iz, shape_x, shape_p)
  integer, intent(in) :: ix, iy, iz
  integer, intent(in) :: shape_x(3)
  integer, intent(in) :: shape_p(2)

  xyz_to_pid = (iz * shape_p(2) / shape_x(3)) * shape_p(1) + (iy * shape_p(1) / shape_x(2))

end function

subroutine init_mesh_to_grid(nelm, shape_x)
  use parallel, only : lglel, gllnid, nid
  integer, intent(in) :: nelm
  integer, intent(in) :: shape_x(3)

  integer :: i, j, ieg, src_pid, dest_pid, npids, slot
  integer :: idx, idy, idz
  integer :: ix(3)
  integer :: shape_p(2)

  shape_p(1) = int(sqrt(real(comm_size)))
  shape_p(2) = int(sqrt(real(comm_size)))

  allocate(dest_pids(nelm))
  npids = 0
  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg, shape_x)
    dest_pid = xyz_to_pid(ix(1), ix(2), ix(3), shape_x, shape_p)
    slot =  -1
    do j = 1, npids
      if (dest_pid == dest_pids(j)) then
        slot = j
        exit
      endif
    enddo
    if (slot == -1) then
      npids = npids + 1
      slot = npids
      dest_pids(slot) = dest_pid
    endif 
  enddo
  deallocate(dest_pids)
  allocate(dest_pids(npids))
  allocate(dest_lengths(npids))
  allocate(dest_slots(nelm))
  allocate(dest_indexes(nelm))
  npids = 0
  dest_lengths = 0
  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg, shape_x)
    dest_pid = xyz_to_pid(ix(1), ix(2), ix(3), shape_x, shape_p)
    slot =  -1
    do j = 1, npids
      if (dest_pid == dest_pids(j)) then
        slot = j
        exit
      endif
    enddo
    if (slot == -1) then
      npids = npids + 1
      slot = npids
      dest_pids(slot) = dest_pid
    endif 
    dest_slots(i) = slot
    dest_lengths(slot) = dest_lengths(slot) + 1
    dest_indexes(i) = dest_lengths(slot) 
  enddo
 
  allocate(src_pids(shape_x(1) * nin_local_xy * nin_local_yz))
  npids = 0
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
  do idx = 0, shape_x(1)-1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        src_pid = gllnid(IEG)
        slot =  -1
        do j = 1, npids
          if (src_pid == src_pids(j)) then
            slot = j
            exit
          endif
        enddo
        if (slot == -1) then
          npids = npids + 1
          slot = npids
          src_pids(slot) = src_pid
        endif 
      enddo
    enddo
  enddo
  deallocate(src_pids)

  allocate(src_pids(npids))
  allocate(src_lengths(npids))
  allocate(src_slots(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  allocate(src_indexes(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  npids = 0
  src_lengths = 0
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
  do idx = 0, shape_x(1)-1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        src_pid = gllnid(IEG)
        slot =  -1
        do j = 1, npids
          if (src_pid == src_pids(j)) then
            slot = j
            exit
          endif
        enddo
        if (slot == -1) then
          npids = npids + 1
          slot = npids
          src_pids(slot) = src_pid
        endif
        src_slots(idx,idy-idx_in_local_xy, idz-idx_in_local_yz) = slot
        src_lengths(slot) = src_lengths(slot) + 1
        src_indexes(idx,idy-idx_in_local_xy, idz-idx_in_local_yz) = src_lengths(slot) 
      enddo
    enddo
  enddo

  allocate(send_buffers(size(dest_pids)))
  do i = 1, size(dest_lengths)
    allocate(send_buffers(i)%p(dest_lengths(i)))
  enddo
  allocate(rec_buffers(size(src_pids)))
  do i = 1, size(src_lengths)
    allocate(rec_buffers(i)%p(src_lengths(i)))
  enddo

  call nekgsync()
  if (nid == 0) write(*,*) "Finished init", nelm

end subroutine init_mesh_to_grid

subroutine mesh_to_grid(mesh, grid, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nid, nekreal
  use parallel, only : lglel, gllnid
  use mpif, only : MPI_STATUS_IGNORE

  real(DP), intent(in) :: mesh(:)
  real(DP), intent(out) :: grid(0:,0:,0:)
  integer, intent(in) :: shape_x(3)

  integer, allocatable :: mpi_reqs(:)
  integer :: n_mpi_reqs
  integer :: dest_pid, src_pid, i, idx, idy, idz, ieg, ierr
  integer :: ix(3)
  integer :: nelm
  integer :: slot, index_in_slot
  nelm = size(mesh)

  if (.not. mesh_to_grid_initialized) then
    call init_mesh_to_grid(nelm, shape_x)
    mesh_to_grid_initialized = .true.
  endif

  ! go through our stuff 
  do i = 1, nelm
    send_buffers(dest_slots(i))%p(dest_indexes(i)) = mesh(i)
  enddo

  allocate(mpi_reqs(size(src_pids)+size(dest_pids))); n_mpi_reqs = 0
  do i = 1, size(dest_pids)
    n_mpi_reqs = n_mpi_reqs + 1
    call MPI_Isend(send_buffers(i)%p, dest_lengths(i), &
                   nekreal, dest_pids(i), nid+dest_pids(i)*comm_size, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
  enddo

  do i = 1, size(src_pids)
    n_mpi_reqs = n_mpi_reqs + 1
    call MPI_Irecv(rec_buffers(i)%p, src_lengths(i), &
                   nekreal, src_pids(i), src_pids(i)+nid*comm_size, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
  enddo

  do i = 1, n_mpi_reqs
    call MPI_Wait(mpi_reqs(i), MPI_STATUS_IGNORE, ierr)        
  enddo


  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        slot = src_slots(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        index_in_slot = src_indexes(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz) = &
          rec_buffers(slot)%p(index_in_slot)
      enddo
    enddo
  enddo

end subroutine mesh_to_grid

subroutine grid_to_mesh(grid, mesh, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nid, nekreal
  use parallel, only : lglel, gllel, gllnid
  use mpif, only : MPI_STATUS_IGNORE

  real(DP), intent(in) :: grid(0:,0:,0:)
  real(DP), intent(out) :: mesh(:)
  integer, intent(in) :: shape_x(3)

  integer, allocatable :: mpi_reqs(:)
  integer :: n_mpi_reqs
  integer :: dest_pid, src_pid, i, idx, idy, idz, ieg, ierr
  integer :: ix(3)
  integer :: slot, index_in_slot
  integer :: nelm
  nelm = size(mesh)

  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        slot = src_slots(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        index_in_slot = src_indexes(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        rec_buffers(slot)%p(index_in_slot) = grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
      enddo
    enddo
  enddo

  allocate(mpi_reqs(size(src_pids)+size(dest_pids))); n_mpi_reqs = 0
  do i = 1, size(src_pids)
    n_mpi_reqs = n_mpi_reqs + 1
    call MPI_Isend(rec_buffers(i)%p, src_lengths(i), &
                   nekreal, src_pids(i), src_pids(i)+nid*1024, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
  enddo


  do i = 1, size(dest_pids)
    n_mpi_reqs = n_mpi_reqs + 1
    call MPI_Irecv(send_buffers(i)%p, dest_lengths(i), &
                   nekreal, dest_pids(i), nid+dest_pids(i)*1024, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
  enddo

  do i = 1, n_mpi_reqs
    call MPI_Wait(mpi_reqs(i), MPI_STATUS_IGNORE, ierr)        
  enddo


  ! go through our stuff 
  do i = 1, nelm
    mesh(i) = send_buffers(dest_slots(i))%p(dest_indexes(i))
  enddo

end subroutine grid_to_mesh

subroutine poisson_kernel(grid, shape_x, start_x, end_x)
  use kinds, only : DP
  use tstep, only : pi 
  use fft, only : wavenumber, P_FORWARD, W_FORWARD

  real(DP), intent(inout) :: grid(0:,0:,0:)
  integer,  intent(in) :: shape_x(3)
  real(DP), intent(in) :: start_x(3)
  real(DP), intent(in) :: end_x(3)
  real(DP) :: kx, ky, kz
  real(DP), allocatable, save :: ks(:,:,:)
  integer :: idz, idy, idx

  ! if we don't have the kernel weights, generate them 
  if (.not. allocated(ks)) then
    allocate(ks(0:ubound(grid,1), 0:ubound(grid,2), 0:ubound(grid,3))) 
    do idz = 0, shape_x(3) - 1
      do idy = 0, nout_local_yz - 1
        do idx = 0, nout_local_xy - 1
          kx = wavenumber(idx + idx_out_local_xy, shape_x(1), &
                          end_x(1)-start_x(1), P_FORWARD)

          ky = wavenumber(idy + idx_out_local_yz, shape_x(2), &
                          end_x(2)-start_x(2), P_FORWARD)
 
          kz = wavenumber(idz, shape_x(3), &
                          end_x(3)-start_x(3), W_FORWARD)
 
          if (kx**2. + ky**2. + kz**2. < 1.e-9_dp) then
            ks(idz,idx,idy) = 0._dp
          else
            ks(idz, idx, idy) = 1._dp / ( &
              (kz)**2._dp + &
              (ky)**2._dp + &
              (kx)**2._dp)
          endif
        enddo
      enddo
    enddo
  endif

  ! vector-vector rescale
  grid = grid * ks

end subroutine poisson_kernel

subroutine shuffle_test()
  use kinds, only : DP
  use size_m, only : nelv
  use mesh, only : shape_x
  use parallel, only : nid
  use parallel, only : lglel

  use fft, only :P_FORWARD, P_BACKWARD, W_FORWARD, W_BACKWARD
  use fft, only : fft_r2r, transpose_grid

  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  integer :: nelm
  integer :: i
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: err
  real(DP) :: rescale


  ! convert RHS to coarse mesh
  nelm = nelv
  allocate(rhs_coarse(nelm))
  do i = 1, nelm
    rhs_coarse(i) = lglel(i)
  enddo
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )
  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_FORWARD, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_FORWARD, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_FORWARD, rescale)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_BACKWARD, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_BACKWARD, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_BACKWARD, rescale)

  plane_xy = plane_xy * (1._dp/ rescale)

  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  ! reorder to local elements
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)

  do i = 1, nelm
    err = abs(soln_coarse(i) - lglel(i))
    if (err > 0.001) then
      write(*,*) "WARNING: shuffle not working", nid, i, nid*nelm + i, err
      exit
    endif
  enddo
  if (nid == 0) write(*,*) "Passed shuffle test" 
end subroutine shuffle_test

subroutine transpose_test()
  use kinds, only : DP
  use size_m, only : nelv
  use mesh, only : shape_x
  use parallel, only : nid, lglel

  use fft, only : transpose_grid

  real(DP), allocatable :: rhs_coarse(:)
  integer :: nelm
  integer :: i, ieg
  integer :: idx, idy, idz
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: err
  real(DP) :: rescale


  ! convert RHS to coarse mesh
  nelm = nelv
  allocate(rhs_coarse(nelm))

  do i = 1, nelm
    rhs_coarse(i) = lglel(i)
  enddo
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        err = abs(plane_xy(idx, idy-idx_in_local_xy, idz-idx_in_local_yz) - ieg)
        if (err > 0.001) then
          write(*,'(A,6(I6))') "WARNING: confused about k after init", nid, idx, idy, idz, ieg, int(plane_xy(idy,idx,idz))
          return
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed init"

  ! forward FFT
  rescale = 1._dp
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  do idx = 0, nout_local_xy - 1
    do idy = 0, shape_x(2) - 1
      do idz = 0, nin_local_yz - 1
        ieg = 1 + idx + idx_out_local_xy + shape_x(1)*idy &
                + shape_x(1) * shape_x(2) * (idz + idx_in_local_yz)
        err = abs(plane_yx(idy,idx,idz) - ieg)
        if (err > 0.001) then
          write(*,'(A,6(I6))') "WARNING: confused about k after xy", nid, idx, idy, idz, ieg, int(plane_yx(idy,idx,idz))
          return
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed xy transpose"

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  do idx = 0, nout_local_xy - 1
    do idy = 0, nout_local_yz - 1
      do idz = 0, shape_x(3) - 1
        ieg = 1 + idx + idx_out_local_xy + shape_x(1)*(idy + idx_out_local_yz) &
                + shape_x(1) * shape_x(2) * idz
        err = abs(plane_zy(idz,idx,idy) - ieg)
        if (err > 0.001) then
          write(*,'(A,6(I6))') "WARNING: confused about k after yz", nid, idx, idy, idz, ieg, int(plane_zy(idz,idy,idx))
          return
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed yz transpose"
 
end subroutine transpose_test

subroutine cos_test()
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nelv
  use geom, only : bm1
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nid, lglel
  use tstep, only : PI, imesh
  use soln, only : pmask

  use fft, only : P_FORWARD, P_BACKWARD
  use fft, only : W_FORWARD, W_BACKWARD
  use fft, only : fft_r2r, transpose_grid
  use semhat, only : zh

  real(DP), allocatable :: rhs_fine(:,:,:,:), soln_fine(:,:,:,:)
  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  real(DP), allocatable :: tmp_fine(:,:,:,:)
  integer :: nelm
  integer :: i, j
  integer :: ix(3)
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: rescale
  real(DP) :: c
  real(DP) :: ax, bx, x
  real(DP) :: ay, by, y
  real(DP) :: az, bz, z
  real(DP), allocatable :: h1(:,:,:,:), h2(:,:,:,:)

  ! convert RHS to coarse mesh
  nelm = nelv
  allocate(rhs_fine(lx1,ly1,lz1,nelm))
  allocate(rhs_coarse(nelm))
  do i = 1, nelm
    ix = ieg_to_xyz(lglel(i), shape_x)
    rhs_coarse(i) =  cos(8.*pi * (ix(3) / shape_x(3))) & 
                   + cos(2.*pi * (ix(3) / shape_x(3)))
  enddo
  if (nid == 0) write(*,*) "COS TEST: rhs ", sqrt(sum(rhs_coarse*rhs_coarse))
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_FORWARD, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_FORWARD, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_FORWARD, rescale)

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_BACKWARD, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_BACKWARD, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_BACKWARD, rescale)
  plane_xy = plane_xy * (1._dp/ rescale)


  ! reorder to local elements
  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)
  if (nid == 0) write(*,*) "COS TEST: u_c ", &
    sqrt(sum(soln_coarse*soln_coarse)) * (2._dp * pi /(end_x(1)-start_x(1)))**2._dp

  return
 
end subroutine cos_test


end module poisson
