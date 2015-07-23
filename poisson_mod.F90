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
  logical, save :: interface_initialized = .false.
  logical, save :: mesh_to_grid_initialized = .false.

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
  real(DP), allocatable :: buffer(:)

  integer :: comm_size
  integer :: mesh_to_grid_handle

contains

!> \brief 
subroutine spectral_solve(u,rhs)!,h1,mask,mult,imsh,isd)
  use kinds, only : DP
  use geom, only : bm1
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nekcomm, lglel
  use soln, only : vmult
  use size_m, only : nx1, ny1, nz1
  use ctimer, only : nscps, tscps, dnekclock

  use fft, only : P_FORWARD, P_BACKWARD, W_FORWARD, W_BACKWARD
  use fft, only : fft_r2r, transpose_grid
  use mesh, only : boundaries

  REAL(DP), intent(out)   :: U    (:)
  REAL(DP), intent(inout) :: RHS  (:)

  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  real(DP), allocatable :: tmp_fine(:,:,:,:)
  integer :: nelm
  integer :: i
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: rescale
  real(DP) :: etime


  nelm = size(rhs) / 8

  if (.not. interface_initialized) then
    call init_comm_infrastructure(nekcomm, shape_x)
  endif

  nscps = nscps + 1

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
  !forall(i = 1 : nelm) rhs_coarse(i) = (tmp_fine(1,1,1,i) + tmp_fine(1,1,nz1,i))/2._dp
  forall(i = 1 : nelm) rhs_coarse(i) = ( &
                                      + tmp_fine(1,1,1,i) &
                                      + tmp_fine(nx1,1,1,i) &
                                      + tmp_fine(1,ny1,1,i) &
                                      + tmp_fine(nx1,ny1,1,i) &
                                      + tmp_fine(1,1,nz1,i) &
                                      + tmp_fine(nx1,1,nz1,i) &
                                      + tmp_fine(1,ny1,nz1,i) &
                                      + tmp_fine(nx1,ny1,nz1,i) &
                                       )/8._dp
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)
  etime = dnekclock()

  ! forward FFT
  rescale = 1._dp
  if (boundaries(2) == 'P  ') then
    call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_FORWARD, rescale)
  else
    call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), W_FORWARD, rescale)
  endif
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  if (boundaries(1) == 'P  ') then
    call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_FORWARD, rescale)
  else
    call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), W_FORWARD, rescale)
  endif

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  if (boundaries(5) == 'P  ') then
    call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), P_FORWARD, rescale)
  else
    call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_FORWARD, rescale)
  endif

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x, boundaries)

  ! reverse FFT
  if (boundaries(5) == 'P  ') then
    call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_BACKWARD, rescale)
  endif

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  if (boundaries(1) == 'P  ') then
    call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), W_BACKWARD, rescale)
  endif

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  if (boundaries(2) == 'P  ') then
    call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), W_BACKWARD, rescale)
  endif

  ! normalize the FFTs
  rescale = 1._dp / (rescale * sum(bm1(:,:,:,1)))
  plane_xy = plane_xy * rescale

  tscps = tscps + (dnekclock() - etime) 
  ! reorder to local elements
  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)

  ! populate U
  tmp_fine = 0._dp
  forall (i = 1: nelm)
    tmp_fine(1,   1,   1,   i)   = soln_coarse(i) 
    tmp_fine(nx1, 1,   1,   i)   = soln_coarse(i) 
    tmp_fine(1  , ny1, 1,   i)   = soln_coarse(i) 
    tmp_fine(nx1, ny1, 1,   i)   = soln_coarse(i) 
    tmp_fine(1,   1,   nz1, i)   = soln_coarse(i) 
    tmp_fine(nx1, 1,   nz1, i)   = soln_coarse(i) 
    tmp_fine(1  , ny1, nz1, i)   = soln_coarse(i) 
    tmp_fine(nx1, ny1, nz1, i)   = soln_coarse(i) 
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


  return
 
end subroutine spectral_solve

!> \brief one-time setup of communication infrastructure for poisson_mod
subroutine init_comm_infrastructure(comm_world, shape_x)
  use input, only : param
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: shape_x(3) !>!< Shape of mesh

  integer(C_INTPTR_T) :: shape_c(3)
  integer :: nxy, nyz, ixy, iyz
  integer :: nid, ierr, i

  call MPI_Comm_rank(comm_world, nid, ierr) 
  call MPI_Comm_size(comm_world, comm_size, ierr) 
  !nid = nid/2 + mod(nid, 2) * comm_size / 2

  if (param(49) >= 1) then
    comm_size = min(comm_size, int(param(49)))
  endif
  comm_size = min(comm_size, shape_x(1) * shape_x(2))
  comm_size = min(comm_size, shape_x(2) * shape_x(3))
  comm_size = min(comm_size, shape_x(1) * shape_x(3))

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
    idx_in_local_xy = (mod(nid,nxy) * shape_x(2)) / nxy
    idx_out_local_xy = (mod(nid,nxy) * shape_x(1)) / nxy
    idx_in_local_yz = ((nxy*nid/comm_size) * shape_x(3)) / nyz
    idx_out_local_yz = ((nxy*nid/comm_size) * shape_x(2)) / nyz
    nin_local_yz = shape_x(3) / nxy; nout_local_yz = shape_x(2)/nxy
    nin_local_xy = shape_x(2) / nyz; nout_local_xy = shape_x(1)/nyz

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

  if (param(75) < 1) then
    call transpose_test()
    call shuffle_test()
    call cos_test()
  endif

  interface_initialized = .true.

end subroutine init_comm_infrastructure

integer function xyz_to_pid(ix, iy, iz, shape_x, shape_p)
  integer, intent(in) :: ix, iy, iz
  integer, intent(in) :: shape_x(3)
  integer, intent(in) :: shape_p(2)

  xyz_to_pid = (iz * shape_p(2) / shape_x(3)) * shape_p(1) + (iy * shape_p(1) / shape_x(2))
  !xyz_to_pid = xyz_to_pid * 2

end function


integer function xyz_to_glo(ix, iy, iz, shape_x)
  integer, intent(in) :: ix, iy, iz
  integer, intent(in) :: shape_x(3)
  xyz_to_glo = ix + (iy + shape_x(2) * iz) * shape_x(1) + 1
end function xyz_to_glo

subroutine init_mesh_to_grid(nelm, shape_x, comm_world)
  use kinds, only : i8
  use parallel, only : lglel, gllnid, nid
  use parallel, only : np
  use mesh, only : ieg_to_xyz
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: nelm
  integer, intent(in) :: shape_x(3)

  integer, external :: iglmax

  integer :: ix(3), i, j
  integer(i8), allocatable :: glo_num(:)
  integer :: nxy_max, nyz_max
  integer :: idx, idy, idz, ieg

  nxy_max = iglmax(nin_local_xy,1)
  nyz_max = iglmax(nin_local_yz,1)

  allocate(glo_num(nelm + shape_x(1) * nxy_max * nyz_max))
  glo_num = 0

  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg)
    glo_num(i) = xyz_to_glo(ix(1), ix(2), ix(3), shape_x)
  enddo

  i = nelm + 1
  do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idx = 0, shape_x(1)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  allocate(buffer(i))
  call gs_setup(mesh_to_grid_handle,glo_num,nelm+ shape_x(1) * nxy_max * nyz_max,comm_world,np)

  call nekgsync()
  if (nid == 0) write(*,*) "Finished init", nelm

end subroutine init_mesh_to_grid

subroutine mesh_to_grid(mesh, grid, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nekreal
  use parallel, only : lglel, gllnid
  use mpif, only : MPI_STATUS_IGNORE

  real(DP), intent(in) :: mesh(:)
  real(DP), intent(out) :: grid(0:,0:,0:)
  integer, intent(in) :: shape_x(3)

  integer, allocatable :: mpi_reqs(:)
  integer :: n_mpi_reqs
  integer :: i, idx, idy, idz, ierr
  integer :: nelm
  integer :: slot, index_in_slot
  nelm = size(mesh)

  if (.not. mesh_to_grid_initialized) then
    call init_mesh_to_grid(nelm, shape_x, nekcomm)
    mesh_to_grid_initialized = .true.
  endif

  ! go through our stuff 
  buffer = 0._dp
  buffer(1:nelm) = mesh(1:nelm)

  call gs_op(mesh_to_grid_handle,buffer,1,1,0)

  i = nelm + 1
  do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idx = 0, shape_x(1)-1
        grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz) = buffer(i)
        i = i + 1
      enddo
    enddo
  enddo

end subroutine mesh_to_grid

subroutine grid_to_mesh(grid, mesh, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nekreal
  use parallel, only : lglel, gllel, gllnid
  use mpif, only : MPI_STATUS_IGNORE

  real(DP), intent(in) :: grid(0:,0:,0:)
  real(DP), intent(out) :: mesh(:)
  integer, intent(in) :: shape_x(3)

  integer, allocatable :: mpi_reqs(:)
  integer :: n_mpi_reqs
  integer :: i, idx, idy, idz, ierr
  integer :: slot, index_in_slot
  integer :: nelm
  nelm = size(mesh)

  buffer = 0._dp
  i = nelm + 1
  do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idx = 0, shape_x(1)-1
        buffer(i) = grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        i = i + 1
      enddo
    enddo
  enddo

  call gs_op(mesh_to_grid_handle,buffer,1,1,1)

  mesh(1:nelm) = buffer(1:nelm)
 
end subroutine grid_to_mesh

subroutine poisson_kernel(grid, shape_x, start_x, end_x, boundaries)
  use kinds, only : DP
  use fft, only : wavenumber, P_FORWARD, W_FORWARD

  real(DP), intent(inout) :: grid(0:,0:,0:)
  integer,  intent(in) :: shape_x(3)
  real(DP), intent(in) :: start_x(3)
  real(DP), intent(in) :: end_x(3)
  character(3), intent(in) :: boundaries(6)
  real(DP) :: kx, ky, kz
  real(DP), allocatable, save :: ks(:,:,:)
  integer :: idz, idy, idx

  ! if we don't have the kernel weights, generate them 
  if (.not. allocated(ks)) then
    allocate(ks(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1)) 
    do idz = 0, shape_x(3) - 1
      do idy = 0, nout_local_yz - 1
        do idx = 0, nout_local_xy - 1
          if (boundaries(2) == 'P  ') then
            kx = wavenumber(idx + idx_out_local_xy, shape_x(1), &
                            end_x(1)-start_x(1), P_FORWARD)
          else
            kx = wavenumber(idx + idx_out_local_xy, shape_x(1), &
                            end_x(1)-start_x(1), W_FORWARD)
          endif

          if (boundaries(1) == 'P  ') then
            ky = wavenumber(idy + idx_out_local_yz, shape_x(2), &
                            end_x(2)-start_x(2), P_FORWARD)
          else
            ky = wavenumber(idy + idx_out_local_yz, shape_x(2), &
                            end_x(2)-start_x(2), W_FORWARD)
          endif

          if (boundaries(5) == 'P  ') then
            kz = wavenumber(idz, shape_x(3), &
                            end_x(3)-start_x(3), P_FORWARD)
          else
            kz = wavenumber(idz, shape_x(3), &
                            end_x(3)-start_x(3), W_FORWARD)
          endif
 
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
          write(*,'(A,6(I6))') "WARNING: confused about k after init", nid, idx, idy, idz, ieg, &
          int(plane_xy(idx,idy-idx_in_local_xy,idz-idx_in_local_yz))
          !return
        endif
      enddo
    enddo
  enddo
  return
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
          write(*,'(A,6(I6))') "WARNING: confused about k after yz", nid, idx, idy, idz, ieg, int(plane_zy(idz,idx,idy))
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
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nid, lglel
  use tstep, only : PI
  use mesh, only : ieg_to_xyz
  use mesh, only : boundaries

  use fft, only : P_FORWARD, P_BACKWARD
  use fft, only : W_FORWARD, W_BACKWARD
  use fft, only : fft_r2r, transpose_grid

  real(DP), allocatable :: rhs_fine(:,:,:,:) 
  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  integer :: nelm
  integer :: i
  integer :: ix(3)
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: rescale

  ! convert RHS to coarse mesh
  nelm = nelv
  allocate(rhs_fine(lx1,ly1,lz1,nelm))
  allocate(rhs_coarse(nelm))
  do i = 1, nelm
    ix = ieg_to_xyz(lglel(i))
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
  call poisson_kernel(plane_zy, shape_x, start_x, end_x, boundaries)

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
