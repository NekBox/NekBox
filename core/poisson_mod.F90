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
  implicit none

  public spectral_solve
  private

  integer :: comm_xy, comm_yz
  logical :: interface_initialized = .false.

  integer(C_INTPTR_T) :: alloc_local_xy, nin_local_xy, nout_local_xy, idx_in_local_xy, idx_out_local_xy
  integer(C_INTPTR_T) :: alloc_local_yz, nin_local_yz, nout_local_yz, idx_in_local_yz, idx_out_local_yz

contains

!> \brief 
subroutine spectral_solve(u,rhs)!,h1,mask,mult,imsh,isd)
  use kinds, only : DP
  use geom, only : bm1, binvm1
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nekcomm, nid, lglel
  use soln, only : vmult
  use size_m, only : nx1, ny1, nz1, nelv

  use fftw3, only : FFTW_R2HC, FFTW_HC2R, FFTW_REDFT10, FFTW_REDFT01
  use fft, only : fft_r2r, transpose_grid

  REAL(DP), intent(out)   :: U    (:)
  !REAL(DP), intent(out)   :: U    (:,:,:,:)
  REAL(DP), intent(inout) :: RHS  (:)
  !REAL(DP), intent(inout) :: RHS  (:,:,:,:)
!  REAL(DP), intent(in)  :: H1   (:,:,:,:)
!  REAL(DP), intent(in)  :: MASK (:,:,:,:)
!  REAL(DP), intent(in)  :: MULT (:,:,:,:)
!  integer,  intent(in)  :: imsh
!  integer,  intent(in)  :: isd

  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  real(DP), allocatable :: tmp_fine(:,:,:,:)
  integer :: nelm
  integer :: i, j
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  real(DP) :: rescale
  real(DP) :: h2(1,1,1,1)
  integer :: ix(3)

  nelm = size(rhs) / 8

  if (.not. interface_initialized) then
    call init_comm_infrastructure(nekcomm, shape_x)
  endif
  
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
  forall(i = 1 : nelm) rhs_coarse(i) = tmp_fine(1,1,1,i)
  !if (nid == 0) write(*,*) "RHS Coarse", sqrt(sum(rhs_coarse * rhs_coarse)/512)
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_R2HC, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_R2HC, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT10, rescale)

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT01, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_HC2R, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_HC2R, rescale)

  ! normalize the FFTs
  plane_xy = plane_xy * (1._dp/ rescale)

  ! reorder to local elements
  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)

  ! populate U
  forall(i = 1:nelm) soln_coarse(i) = soln_coarse(i) / sum(bm1(:,:,:,i))
  tmp_fine = 0._dp
  forall (i = 1: nelm)
    tmp_fine(1,  1,  1,   i) = soln_coarse(i)
  end forall
  call dssum(tmp_fine)

  !call interpolate_into_element(soln_coarse, u)
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

  ! update residual
#if 0
  allocate(tmp_fine(size(u,1), size(u,2), size(u,3), size(u,4)))
  h2 = 0._dp
  call axhelm (tmp_fine, u, h1, h2, imsh, isd)
  tmp_fine = tmp_fine * mask
  call dssum   (tmp_fine)

  if (nid == 0) write(*,*) "RHS before: ", sqrt(sum(rhs * rhs))
  RHS = RHS - tmp_fine 
  if (nid == 0) write(*,*) "RHS after : ", sqrt(sum(rhs * rhs))
#endif

  return
 
end subroutine spectral_solve

!> \brief one-time setup of communication infrastructure for poisson_mod
subroutine init_comm_infrastructure(comm_world, shape_x)
  use fftw3, only : FFTW_MPI_DEFAULT_BLOCK
  use fftw3, only : fftw_mpi_local_size_many_transposed
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: shape_x(3) !>!< Shape of mesh

  integer(C_INTPTR_T) :: shape_c(3)
  integer(C_INTPTR_T), parameter :: one = 1
  integer :: nxy, nyz, ixy, iyz, offset_xy
  integer :: nid, comm_size, ierr

  call MPI_Comm_rank(comm_world, nid, ierr) 
  call MPI_Comm_size(comm_world, comm_size, ierr) 

  nxy =  2; ixy = 2*nid/comm_size; offset_xy = comm_size / nxy
  call MPI_Comm_split(comm_world, ixy, 0, comm_xy, ierr)
  nyz =  comm_size/nxy; iyz = mod(nid,nyz)
  call MPI_Comm_split(comm_world, iyz, 0, comm_yz, ierr)

  shape_c = shape_x
  alloc_local_xy = fftw_mpi_local_size_many_transposed( 2, &
                (/shape_c(2), shape_c(1)/), &
                one, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                comm_xy, &
                nin_local_xy, idx_in_local_xy, nout_local_xy, idx_out_local_xy)
  alloc_local_yz = fftw_mpi_local_size_many_transposed( 2, &
                (/shape_c(3), shape_c(2)/), &
                one, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                comm_yz, &
                nin_local_yz, idx_in_local_yz, nout_local_yz, idx_out_local_yz)

  write(*,'(A,6(I5))') "MAX:", nid, alloc_local_xy, nin_local_xy, idx_in_local_xy, nout_local_xy, idx_out_local_xy
  call nekgsync()
  write(*,'(A,6(I5))') "MAX:", nid, alloc_local_yz, nin_local_yz, idx_in_local_yz, nout_local_yz, idx_out_local_yz

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

integer function xyz_to_pid(ix, iy, iz, shape_x)
  integer, intent(in) :: ix, iy, iz
  integer, intent(in) :: shape_x(3)

  xyz_to_pid = (iz/nin_local_yz) * (shape_x(2)/nin_local_xy) + (iy/nin_local_xy)

end function

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
  nelm = size(mesh)

  ! go through our stuff
  allocate(mpi_reqs(2*nelm)); n_mpi_reqs = 0
  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg, shape_x)
    dest_pid = xyz_to_pid(ix(1), ix(2), ix(3), shape_x)
    if (dest_pid == nid) then
      ! store locally
      grid(ix(1), ix(2)-idx_in_local_xy, ix(3)-idx_in_local_yz) = mesh(i)
    else
      ! send somewhere
      n_mpi_reqs = n_mpi_reqs + 1
      call MPI_Isend(mesh(i), 1, nekreal, dest_pid, ieg, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
    endif
  enddo

  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        src_pid = gllnid(IEG)
        if (src_pid /= nid) then
          n_mpi_reqs = n_mpi_reqs + 1
          call MPI_Irecv(grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz),&
              1, nekreal, src_pid, ieg, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
        endif
      enddo
    enddo
  enddo
         
  do i = 1, n_mpi_reqs
    call MPI_Wait(mpi_reqs(i), MPI_STATUS_IGNORE, ierr)        
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
  integer :: nelm
  nelm = size(mesh)

  allocate(mpi_reqs(2*nelm)); n_mpi_reqs = 0
  n_mpi_reqs = 0
  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1 
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        dest_pid = gllnid(IEG)
        if (dest_pid /= nid) then
          ! send somewhere
          n_mpi_reqs = n_mpi_reqs + 1
          call MPI_Isend(grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz),&
              1, nekreal, dest_pid, ieg, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
        else
          mesh(gllel(ieg)) = grid(idx, idy-idx_in_local_xy, idz-idx_in_local_yz)
        endif
      enddo
    enddo
  enddo

  do i = 1, nelm
    ieg = lglel(i)
    ix = ieg_to_xyz(ieg, shape_x)
    src_pid = xyz_to_pid(ix(1), ix(2), ix(3), shape_x)
    if (src_pid /= nid) then
      n_mpi_reqs = n_mpi_reqs + 1
      call MPI_Irecv(mesh(i), 1, nekreal, src_pid, ieg, nekcomm, mpi_reqs(n_mpi_reqs), ierr)
    endif
  enddo

  do i = 1, n_mpi_reqs
    call MPI_Wait(mpi_reqs(i), MPI_STATUS_IGNORE, ierr)        
  enddo
end subroutine grid_to_mesh

subroutine poisson_kernel(grid, shape_x, start_x, end_x)
  use kinds, only : DP
  use tstep, only : pi 

  real(DP), intent(inout) :: grid(0:,0:,0:)
  integer,  intent(in) :: shape_x(3)
  real(DP), intent(in) :: start_x(3)
  real(DP), intent(in) :: end_x(3)
  real(DP) :: kx, ky, kz

  integer :: idx, idy, idz
  do idz = 0, shape_x(3) - 1
    do idy = 0, nout_local_yz - 1
      do idx = 0, nout_local_xy - 1
        if (idx + idx_out_local_xy <= shape_x(1) / 2) then
          kx = 2*pi*(idx +idx_out_local_xy)/(end_x(1)-start_x(1)) 
        else
          kx = 2*pi*(shape_x(1) - idx - idx_out_local_xy)/(end_x(1)-start_x(1)) 
        endif

        if (idy + idx_out_local_yz <= shape_x(2) / 2) then
          ky = 2*pi*(idy +idx_out_local_yz)/(end_x(2)-start_x(2)) 
        else
          ky = 2*pi*(shape_x(2) - idy - idx_out_local_yz)/(end_x(2)-start_x(2)) 
        endif

        kz = pi*(idz)/(end_x(3)-start_x(3)) 

        if (kx**2. + ky**2. + kz**2. < 1.e-9_dp) then
          grid(idz,idx,idy) = 0._dp
        else
          grid(idz, idx, idy) = grid(idz, idx, idy) / ( &
            (kz)**2._dp + &
            (ky)**2._dp + &
            (kx)**2._dp)
        endif
      enddo
    enddo
  enddo
end subroutine poisson_kernel

subroutine interpolate_into_element(soln_coarse, soln_fine)
  use kinds, only : DP 
  use semhat, only : zh
  use parallel, only : lglel
  use mesh, only : shape_x

  real(DP), intent(in) :: soln_coarse(:)
  real(DP), intent(out) :: soln_fine(:,:,:,:)

  real(DP) :: c
  real(DP) :: ax, bx, cx, dx, x
  real(DP) :: ay, by, cy, dy, y
  real(DP) :: az, bz, cz, dz, z
  real(DP) :: forward, backward
  integer :: i, j, ix(3)

  forall (i = 1: size(soln_coarse)) soln_fine(:,:,:,i) = soln_coarse(i)
  call dssum(soln_fine)

#if 0
  do i = 1, size(soln_coarse)
    soln_fine(:,:,:,i) = soln_coarse
  enddo
#else
  do i = 1, size(soln_coarse)
    c = soln_coarse(i) 
    forward = soln_fine(8,2,2,i) -c 
    backward = soln_fine(1,2,2,i) -c 
    ax = 0!.25_dp * (backward - forward - backward + forward) / 2.
    bx = .25_dp * (forward + backward - 2._dp * c) / 2._dp
    cx = .25_dp * (backward + forward )
    dx = .25_dp*(backward/2. + forward/2. + 3*c)

    forward = soln_fine(2,8,2,i) - c
    backward = soln_fine(2,1,2,i) - c 
    ay = 0!.25_dp * (backward - forward - backward + forward) / 2.
    by = .25_dp * (forward + backward - 2._dp * c) / 2._dp
    cy = .25_dp * (backward + forward )
    dy = .25_dp*(backward/2. + forward/2. + 3*c)

    ix = ieg_to_xyz(lglel(i), shape_x)
    if (ix(3) == 0) then
      forward = soln_fine(2,2,8,i) -c
      backward = c 
    else if (ix(3) == shape_x(3)) then
      forward = c 
      backward = soln_fine(2,2,1,i) - c 
    else
      forward = soln_fine(2,2,8,i) - c 
      backward = soln_fine(2,2,1,i) -c
    endif
    az = 0!.25_dp * (backward - forward - backward + forward) / 2.
    bz = .25_dp * (forward + backward - 2._dp * c) / 2._dp
    cz = .25_dp * (backward + forward )
    dz = .25_dp*(backward/2. + forward/2. + 3*c)                      !ay = 0._dp; by = 0._dp

    !az = 0._dp; bz = 0._dp

    do j = 1, 8
      x = zh(j)/2._dp
      soln_fine(j,:,:,i) = bx*x*x + cx*x + dx/3.
    enddo
    do j = 1, 8
      y = zh(j)/2._dp
      soln_fine(:,j,:,i) = by*y*y + cy*y + dy/3. + soln_fine(:,j,:,i)
    enddo
    do j = 1, 8
      z = zh(j)/2._dp
      soln_fine(:,:,j,i) = bz*z*z + cz*z + dz/3. + soln_fine(:,:,j,i)
    enddo
  enddo
#endif

end subroutine interpolate_into_element

subroutine shuffle_test()
  use kinds, only : DP
  use size_m, only : nelv
  use mesh, only : shape_x
  use parallel, only : nid
  use parallel, only : lglel

  use fftw3, only :FFTW_R2HC, FFTW_HC2R, FFTW_REDFT10, FFTW_REDFT01
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
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_R2HC, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_R2HC, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT10, rescale)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT01, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_HC2R, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_HC2R, rescale)

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
          write(*,'(A,6(I6))') "WARNING: confused about k after init", nid, idx, idy, idz, ieg, int(plane_yx(idy,idx,idz))
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

  use fftw3, only : FFTW_R2HC, FFTW_HC2R
  use fftw3, only : FFTW_REDFT10, FFTW_REDFT01
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
    do j = 1, 8 
      rhs_fine(:,:,j,i) = bm1(:,:,j,i) * (&
                          cos(8.*pi * (ix(3)+(zh(j)+1)/2._dp) / shape_x(3)) & 
                        + cos(2.*pi * (ix(3)+(zh(j)+1)/2._dp) / shape_x(3)) )
    enddo
  enddo
  forall(i = 1 : nelm) rhs_coarse(i) = sum(rhs_fine(:,:,:,i))
  if (nid == 0) write(*,*) "COS TEST: rhs ", sqrt(sum(rhs_coarse*rhs_coarse))
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

  ! forward FFT
  rescale = 1._dp
  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_R2HC, rescale)
  
  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_xy, plane_yx, shape_x, 1, 2, comm_xy)
  deallocate(plane_xy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_R2HC, rescale)

  allocate(plane_zy(0:shape_x(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  call transpose_grid(plane_yx, plane_zy, shape_x, 2, 3, comm_yz)
  deallocate(plane_yx)

  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT10, rescale)

  ! Poisson kernel
  call poisson_kernel(plane_zy, shape_x, start_x, end_x)

  ! reverse FFT
  call fft_r2r(plane_zy, shape_x(3), int(nout_local_xy * nout_local_yz), FFTW_REDFT01, rescale)

  allocate(plane_yx(0:shape_x(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  call transpose_grid(plane_zy, plane_yx, shape_x, 3, 2, comm_yz)
  deallocate(plane_zy)

  call fft_r2r(plane_yx, shape_x(2), int(nout_local_xy * nin_local_yz), FFTW_HC2R, rescale)

  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  call transpose_grid(plane_yx, plane_xy, shape_x, 2, 1, comm_xy)
  deallocate(plane_yx)

  call fft_r2r(plane_xy, shape_x(1), int(nin_local_xy * nin_local_yz), FFTW_HC2R, rescale)
  plane_xy = plane_xy * (1._dp/ rescale)


  ! reorder to local elements
  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)
  if (nid == 0) write(*,*) "COS TEST: u_c ", &
    sqrt(sum(soln_coarse*soln_coarse)) * (2._dp * pi /(end_x(1)-start_x(1)))**2._dp

  forall(i = 1 : nelm) soln_coarse(i)  = soln_coarse(i)  / sum(bm1(:,:,:,i))
  allocate(soln_fine(lx1,ly1,lz1,nelm))
  !forall(i = 1 : nelm) soln_fine(:,:,:,i) = binvm1(:,:,:,i) * soln_coarse(i) / 512._dp
#if 1
  call interpolate_into_element(soln_coarse, soln_fine)
#else
#if 0
  do i = 1, nelm
    c = soln_coarse(i)
    bx = (front_coarse(i) - back_coarse(i)) / 2._dp
    ax = front_coarse(i) - bx - c
    do j = 1, 8
      x = zh(j)/2._dp
      soln_fine(j,:,:,i) = c + bx * x + ax*x*x
    enddo
  enddo
#else
  forall(i = 1:nelm) soln_fine(:,:,:,i) = soln_coarse(i)
  call dssum(soln_fine)
  do i = 1, nelm
    c = soln_coarse(i)
    bx = (soln_fine(8,2,2,i) - soln_fine(1,2,2,i))/2._dp
    ax = soln_fine(8,2,2,i) - 2*c - bx
    by = (soln_fine(2,8,2,i) - soln_fine(2,1,2,i))/2._dp
    ay = soln_fine(2,8,2,i) - 2*c - by
    bz = (soln_fine(2,2,8,i) - soln_fine(2,2,1,i))/2._dp
    az = soln_fine(2,2,8,i) - 2*c - bz

    ix = ieg_to_xyz(lglel(i), shape_x)
    if (ix(3) == 0) then
      bz = (soln_fine(2,2,8,i) - 2.*soln_fine(2,2,1,i))/2._dp
      az = soln_fine(2,2,8,i) - 2*c - bz
    else if (ix(3) == shape_x(3)) then
      bz = (2.*soln_fine(2,2,8,i) - soln_fine(2,2,1,i))/2._dp
      az = 2.*soln_fine(2,2,8,i) - 2*c - bz
    endif

    do j = 1, 8
      x = zh(j)/2._dp
      soln_fine(j,:,:,i) = c + bx * x + ax*x*x
    enddo
    do j = 1, 8
      y = zh(j)/2._dp
      soln_fine(:,j,:,i) = soln_fine(:,j,:,i)  + by * y + ay*y*y
    enddo
    do j = 1, 8
      z = zh(j)/2._dp
      soln_fine(:,:,j,i) = soln_fine(:,:,j,i)  + bz * z + az*z*z
    enddo
  enddo
#endif

#if 0
  forall(i = 1:nelm) soln_fine(:,:,:,i) = soln_coarse(i)
  call dssum(soln_fine)
  soln_fine = soln_fine * vmult
#endif
#endif

  ! update residual
  allocate(tmp_fine(lx1,ly1,lz1,nelm))
  allocate(h1(lx1,ly1,lz1,nelm))
  allocate(h2(lx1,ly1,lz1,nelm))
  h1 = 1._dp; h2 = 0._dp
  !call axhelm (tmp_fine, soln_fine, h1, h2, imesh, 1)
  call axhelm (tmp_fine, soln_fine, h1, h2, imesh, 1)
  tmp_fine = tmp_fine * pmask
  call dssum   (tmp_fine)

  !rhs_fine = rhs_fine - tmp_fine 
  forall(i = 1 : nelm) rhs_coarse(i) = sum(tmp_fine(:,:,:,i))
  if (nid == 0) write(*,*) "COS TEST: Au ", sqrt(sum(rhs_coarse*rhs_coarse)) 
  rhs_fine = rhs_fine - tmp_fine 
  forall(i = 1 : nelm) rhs_coarse(i) = sum(rhs_fine(:,:,:,i))
  if (nid == 0) write(*,*) "COS TEST: r ", sqrt(sum(rhs_coarse*rhs_coarse)) 

  return
 
end subroutine cos_test


end module poisson
