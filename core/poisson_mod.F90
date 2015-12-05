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

  logical, save :: interface_initialized = .false.
  logical, save :: mesh_to_grid_initialized = .false.

!  integer :: nin_local_xy, nout_local_xy, idx_in_local_xy, idx_out_local_xy
!  integer :: nin_local_yz, nout_local_yz, idx_in_local_yz, idx_out_local_yz

  integer :: n_x(2), idx_x(2)
  integer :: n_y(2), idx_y(2)
  integer :: n_z(2), idx_z(2)

  real(DP), allocatable :: buffer(:)

  integer :: mesh_to_grid_handle
  integer :: transpose_xy_handle, transpose_yx_handle
  integer :: transpose_yz_handle, transpose_zy_handle

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
  use fft, only : fft_r2r
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
  allocate(plane_xy(0:shape_x(1)-1, 0:n_x(1)-1, 0:n_x(2)-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)
  etime = dnekclock()

  ! forward FFT
  rescale = 1._dp
  if (boundaries(2) == 'P  ') then
    call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), P_FORWARD, rescale)
  else
    call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), W_FORWARD, rescale)
  endif
  
  allocate(plane_yx(0:shape_x(2)-1, 0:n_y(1)-1, 0:n_y(2)-1) )
  call transpose_grid(plane_xy, plane_yx, 1)
  deallocate(plane_xy)

  if (boundaries(1) == 'P  ') then
    call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), P_FORWARD, rescale)
  else
    call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), W_FORWARD, rescale)
  endif

  allocate(plane_zy(0:shape_x(3)-1, 0:n_z(1)-1, 0:n_z(2)-1) )
  call transpose_grid(plane_yx, plane_zy, 2)
  deallocate(plane_yx)

  if (boundaries(5) == 'P  ') then
    call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), P_FORWARD, rescale)
  else
    call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), W_FORWARD, rescale)
  endif

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x, boundaries)

  ! reverse FFT
  if (boundaries(5) == 'P  ') then
    call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), W_BACKWARD, rescale)
  endif

  allocate(plane_yx(0:shape_x(2)-1, 0:n_y(2)-1, 0:n_y(2)-1) )
  call transpose_grid(plane_zy, plane_yx, -2)
  deallocate(plane_zy)

  if (boundaries(1) == 'P  ') then
    call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), W_BACKWARD, rescale)
  endif

  allocate(plane_xy(0:shape_x(1)-1, 0:n_x(1)-1, 0:n_x(2)-1))
  call transpose_grid(plane_yx, plane_xy, -1)
  deallocate(plane_yx)

  if (boundaries(2) == 'P  ') then
    call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), P_BACKWARD, rescale)
  else
    call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), W_BACKWARD, rescale)
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
  use parallel, only : proc_pos, proc_shape
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: shape_x(3) !>!< Shape of mesh

  integer(C_INTPTR_T) :: shape_c(3)
  integer :: nxy, nyz
  integer :: nid, ierr, i
  integer :: comm_size
  integer :: proc_n(3)

  integer :: n_in, p_in

  call MPI_Comm_rank(comm_world, nid, ierr) 
  call MPI_Comm_size(comm_world, comm_size, ierr) 
  !nid = nid/2 + mod(nid, 2) * comm_size / 2

  comm_size = min(comm_size, shape_x(1)*shape_x(2))
  comm_size = min(comm_size, shape_x(1)*shape_x(3))
  comm_size = min(comm_size, shape_x(2)*shape_x(3))

  p_in = comm_size
  do while ((shape_x(2) * shape_x(3) / p_in) * p_in /= shape_x(2) * shape_x(3))
    p_in = p_in - 1
  enddo
  n_in = (shape_x(2)*shape_x(3) / p_in)
  if (nid < p_in) then
    n_x(1) = min(int(sqrt(float(n_in))), shape_x(2))
    do while ((n_in / n_x(1)) * n_x(1) /= n_in)
      n_x(1) = n_x(1) - 1
    enddo
    n_x(2) = n_in / n_x(1)
    idx_x(1) = mod(nid*n_x(1), shape_x(2)) 
    idx_x(2) = n_x(2)*(nid*n_x(1)/shape_x(2))
  else
    n_x = 0
    idx_x = -1
  endif

  p_in = comm_size
  do while ((shape_x(1) * shape_x(3) / p_in) * p_in /= shape_x(1) * shape_x(3))
    p_in = p_in - 1
  enddo
  n_in = (shape_x(1)*shape_x(3) / p_in)
  if (nid < p_in) then
    n_y(1) = min(int(sqrt(float(n_in))), shape_x(1))
    do while ((n_in / n_y(1)) * n_y(1) /= n_in)
      n_y(1) = n_y(1) - 1
    enddo
    n_y(2) = n_in / n_y(1)
    idx_y(1) = mod(nid*n_y(1), shape_x(1)) 
    idx_y(2) = n_y(2)*(nid*n_y(1)/shape_x(1))
  else
    n_y = 0
    idx_y = -1
  endif

  p_in = comm_size
  do while ((shape_x(2) * shape_x(1) / p_in) * p_in /= shape_x(2) * shape_x(1))
    p_in = p_in - 1
  enddo
  n_in = (shape_x(2)*shape_x(1) / p_in)
  if (nid < p_in) then
    n_z(1) = min(int(sqrt(float(n_in))), shape_x(1))
    do while ((n_in / n_z(1)) * n_z(1) /= n_in)
      n_z(1) = n_z(1) - 1
    enddo
    n_z(2) = n_in / n_z(1)
    idx_z(1) = mod(nid*n_z(1), shape_x(1)) 
    idx_z(2) = n_z(2)*(nid*n_z(1)/shape_x(1))
  else
    n_z = 0
    idx_z = -1
  endif
 
  if (param(75) < 1) then
    call transpose_test()
    call shuffle_test()
    !call cos_test()
  endif

  interface_initialized = .true.

end subroutine init_comm_infrastructure

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
  integer :: nx_max, ny_max, nz_max
  integer :: idx, idy, idz, ieg
  integer :: buff_size

  nx_max = iglmax(n_x(1)*n_x(2),1)
  ny_max = iglmax(n_y(1)*n_y(2),1)
  nz_max = iglmax(n_z(1)*n_z(2),1)  

  buff_size = max(nelm + shape_x(1) * nx_max, shape_x(1)*nx_max + shape_x(2)*ny_max)
  buff_size = max(buff_size, shape_x(2)*ny_max + shape_x(3) * nz_max)
  allocate(glo_num(buff_size))
  glo_num = 0

  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg)
    glo_num(i) = xyz_to_glo(ix(1), ix(2), ix(3), shape_x)
  enddo

  i = nelm + 1
  do idz = idx_x(2), idx_x(2) + n_x(2) - 1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idx = 0, shape_x(1)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  call gs_setup(mesh_to_grid_handle,glo_num,nelm+ shape_x(1) * nx_max,comm_world,np)
  
  glo_num = 0
  i = 1
  do idz = idx_x(2), idx_x(2) + n_x(2) - 1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idx = 0, shape_x(1)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  do idz = idx_y(2), idx_y(2) + n_y(2) - 1
    do idx = idx_y(1), idx_y(1) + n_y(1) - 1
      do idy = 0, shape_x(2)-1
        glo_num(i) = xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  call gs_setup(transpose_xy_handle,glo_num, shape_x(1)*nx_max + shape_x(2)*ny_max ,comm_world,np)

  glo_num = 0
  i = 1
  do idz = idx_y(2), idx_y(2) + n_y(2) - 1
    do idx = idx_y(1), idx_y(1) + n_y(1) - 1
      do idy = 0, shape_x(2)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  do idz = idx_x(2), idx_x(2) + n_x(2) - 1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idx = 0, shape_x(1)-1
        glo_num(i) = xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  call gs_setup(transpose_yx_handle,glo_num,shape_x(1)*nx_max + shape_x(2)*ny_max,comm_world,np)

  glo_num = 0
  i = 1
  do idz = idx_y(2), idx_y(2) + n_y(2) - 1
    do idx = idx_y(1), idx_y(1) + n_y(1) - 1
      do idy = 0, shape_x(2)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  do idy = idx_z(2), idx_z(2) + n_z(2) - 1
    do idx = idx_z(1), idx_z(1) + n_z(1) - 1
      do idz = 0, shape_x(3)-1
        glo_num(i) = xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  call gs_setup(transpose_yz_handle,glo_num, shape_x(2)*ny_max + shape_x(3)*nz_max,comm_world,np)

  glo_num = 0
  i = 1
  do idy = idx_z(2), idx_z(2) + n_z(2) - 1
    do idx = idx_z(1), idx_z(1) + n_z(1) - 1
      do idz = 0, shape_x(3)-1
        glo_num(i) = -xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  do idz = idx_y(2), idx_y(2) + n_y(2) - 1
    do idx = idx_y(1), idx_y(1) + n_y(1) - 1
      do idy = 0, shape_x(2)-1
        glo_num(i) = xyz_to_glo(idx, idy, idz, shape_x)
        i = i + 1
      enddo
    enddo
  enddo
  call gs_setup(transpose_zy_handle,glo_num,shape_x(2)*ny_max + shape_x(3)*nz_max,comm_world,np)

  deallocate(glo_num)
  allocate(buffer(buff_size))

  call nekgsync()
  if (nid == 0) write(*,*) "Finished init", nelm

end subroutine init_mesh_to_grid

subroutine mesh_to_grid(mesh, grid, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nekreal
  use parallel, only : lglel, gllnid

  real(DP), intent(in) :: mesh(:)
  real(DP), intent(out) :: grid(0:,0:,0:)
  integer, intent(in) :: shape_x(3)

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
  do idz = idx_x(2), idx_x(2) + n_x(2) - 1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idx = 0, shape_x(1)-1
        grid(idx, idy-idx_x(1), idz-idx_x(2)) = buffer(i)
        i = i + 1
      enddo
    enddo
  enddo

end subroutine mesh_to_grid

subroutine grid_to_mesh(grid, mesh, shape_x)
  use kinds, only : DP
  use parallel, only : nekcomm, nekreal
  use parallel, only : lglel, gllel, gllnid

  real(DP), intent(in) :: grid(0:,0:,0:)
  real(DP), intent(out) :: mesh(:)
  integer, intent(in) :: shape_x(3)

  integer :: i, idx, idy, idz, ierr
  integer :: slot, index_in_slot
  integer :: nelm
  nelm = size(mesh)

  buffer = 0._dp
  i = nelm + 1
  do idz = idx_x(2), idx_x(2) + n_x(2) - 1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idx = 0, shape_x(1)-1
        buffer(i) = grid(idx, idy-idx_x(1), idz-idx_x(2))
        i = i + 1
      enddo
    enddo
  enddo

  call gs_op(mesh_to_grid_handle,buffer,1,1,1)

  mesh(1:nelm) = buffer(1:nelm)
 
end subroutine grid_to_mesh

subroutine transpose_grid(plane_xy, plane_yx, dir)
  use kinds, only : DP
  use mesh, only : shape_x
  implicit none
  real(DP), intent(in) :: plane_xy(*)
  real(DP), intent(out) :: plane_yx(*)
  integer :: dir

  integer :: n, m
  if (dir == 1) then
    n = shape_x(1) * n_x(1) * n_x(2) 
  else if (dir == 2 .or. dir == -1) then
    n = shape_x(2) * n_y(1) * n_y(2) 
  else if (dir == -2) then
    n = shape_x(3) * n_z(1) * n_z(2) 
  endif
  
  buffer = 0._dp
  buffer(1:n) = plane_xy(1:n)
  if (dir == 1) then
    call gs_op(transpose_xy_handle, buffer, 1, 1, 1)
  else  if (dir == 2) then
    call gs_op(transpose_yz_handle, buffer, 1, 1, 1)
  else if (dir == -1) then    
    call gs_op(transpose_yx_handle, buffer, 1, 1, 1)
  else if (dir == -2) then
    call gs_op(transpose_zy_handle, buffer, 1, 1, 1)
  endif

  if (dir == 1 .or. dir == -2) then
    m = shape_x(2) * n_y(1) * n_y(2) 
  else if (dir == 2) then
    m = shape_x(3) * n_z(1) * n_z(2) 
  else if (dir == -1) then
    m = shape_x(1) * n_x(1) * n_x(2)
  endif

  plane_yx(1:m) = buffer(n+1:n+m)
  return

end subroutine transpose_grid


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
    allocate(ks(0:shape_x(3)-1, 0:n_z(1)-1, 0:n_z(2)-1)) 
    do idz = 0, shape_x(3) - 1
      do idy = 0, n_z(2) - 1
        do idx = 0, n_z(1) - 1
          if (boundaries(2) == 'P  ') then
            kx = wavenumber(idx + idx_z(1), shape_x(1), &
                            end_x(1)-start_x(1), P_FORWARD)
          else
            kx = wavenumber(idx + idx_z(1), shape_x(1), &
                            end_x(1)-start_x(1), W_FORWARD)
          endif

          if (boundaries(1) == 'P  ') then
            ky = wavenumber(idy + idx_z(2), shape_x(2), &
                            end_x(2)-start_x(2), P_FORWARD)
          else
            ky = wavenumber(idy + idx_z(2), shape_x(2), &
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
  use fft, only : fft_r2r 

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
  allocate(plane_xy(0:shape_x(1)-1, 0:n_x(1)-1, 0:n_x(2)-1) )
  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), P_FORWARD, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:n_y(1)-1, 0:n_y(2)-1) )
  call transpose_grid(plane_xy, plane_yx, 1) 
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), P_FORWARD, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:n_z(1)-1, 0:n_z(2)-1) )
  call transpose_grid(plane_yx, plane_zy, 2) 
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), W_FORWARD, rescale)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(n_z(1)*n_z(2)), W_BACKWARD, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:n_y(1)-1, 0:n_y(2)-1) )
  call transpose_grid(plane_zy, plane_yx, -2) 
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(n_y(1)*n_y(2)), P_BACKWARD, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:n_x(1)-1, 0:n_x(2)-1))
  call transpose_grid(plane_yx, plane_xy, -1) 
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(n_x(1)*n_x(2)), P_BACKWARD, rescale)

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
  allocate(plane_xy(0:shape_x(1)-1, 0:n_x(1)-1, 0:n_x(2)-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  do idx = 0, shape_x(1)-1
    do idy = idx_x(1), idx_x(1) + n_x(1) - 1
      do idz = idx_x(2), idx_x(2) + n_x(2) - 1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        err = abs(plane_xy(idx, idy-idx_x(1), idz-idx_x(2)) - ieg)
        if (err > 0.001) then
          write(*,'(A,6(I6))') "WARNING: confused about k after init", nid, idx, idy, idz, ieg, &
          int(plane_xy(idx,idy-idx_x(1),idz-idx_x(2)))
          return
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed init"

  ! forward FFT
  rescale = 1._dp
  
  allocate(plane_yx(0:shape_x(2)-1, 0:n_y(1)-1, 0:n_y(2)-1) )
  call transpose_grid(plane_xy, plane_yx, 1)
  deallocate(plane_xy)

  do idx = 0, n_y(1) - 1
    do idy = 0, shape_x(2) - 1
      do idz = 0, n_y(2) - 1
        ieg = 1 + idx + idx_y(1) + shape_x(1)*idy &
                + shape_x(1) * shape_x(2) * (idz + idx_y(2))
        err = abs(plane_yx(idy,idx,idz) - ieg)
        if (err > 0.001) then
          write(*,'(A,6(I6))') "WARNING: confused about k after xy", nid, idx, idy, idz, ieg, int(plane_yx(idy,idx,idz))
          return
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed xy transpose"

  allocate(plane_zy(0:shape_x(3)-1, 0:n_z(1)-1, 0:n_z(2)-1) )
  call transpose_grid(plane_yx, plane_zy, 2)
  deallocate(plane_yx)

  do idx = 0, n_z(1) - 1
    do idy = 0, n_z(2) - 1
      do idz = 0, shape_x(3) - 1
        ieg = 1 + idx + idx_z(1) + shape_x(1)*(idy + idx_z(2)) &
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


#if 0
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
  use fft, only : fft_r2r

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
  call transpose_grid(plane_xy, plane_yx, 1)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_FORWARD, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy,2)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_FORWARD, rescale)

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x, boundaries)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), W_BACKWARD, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx,-2)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), P_BACKWARD, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, -1)
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
#endif

end module poisson
