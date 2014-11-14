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
  integer :: nxy, nyz, ixy, iyz, offset_xy
  logical :: interface_initialized = .false.

  integer(C_INTPTR_T) :: alloc_local_xy, nin_local_xy, nout_local_xy, idx_in_local_xy, idx_out_local_xy
  integer(C_INTPTR_T) :: alloc_local_yz, nin_local_yz, nout_local_yz, idx_in_local_yz, idx_out_local_yz

contains

!#define UNITARY_TEST
#define SHUFFLE_TEST
#define NO_FFT
!> \brief 
subroutine spectral_solve(u,rhs,h1,mask,mult,imsh,isd)
  use kinds, only : DP
  use geom, only : bm1, binvm1
  use mesh, only : shape_x, start_x, end_x
  use parallel, only : nekcomm, nid, nekreal
  use parallel, only : lglel, gllel, gllnid
  use tstep, only : PI

  use fftw3, only : fftw_mpi_local_size_many_transposed
  use fftw3, only : fftw_plan_many_r2r
  use fftw3, only : fftw_mpi_plan_many_transpose
  use fftw3, only : FFTW_MPI_DEFAULT_BLOCK
  use fftw3, only : FFTW_ESTIMATE, FFTW_R2HC, FFTW_HC2R, FFTW_REDFT00
  use fftw3, only : fftw_execute_r2r, fftw_mpi_execute_r2r
  use mpif, only : MPI_STATUS_IGNORE

  REAL(DP), intent(out)   :: U    (:,:,:,:)
  REAL(DP), intent(inout) :: RHS  (:,:,:,:)
  REAL(DP), intent(in)  :: H1   (:,:,:,:)
  REAL(DP), intent(in)  :: MASK (:,:,:,:)
  REAL(DP), intent(in)  :: MULT (:,:,:,:)
  integer,  intent(in)  :: imsh
  integer,  intent(in)  :: isd

  real(DP), allocatable :: rhs_coarse(:), soln_coarse(:)
  real(DP), allocatable :: tmp_fine(:,:,:,:)
  integer :: nelm
  integer :: i, ieg, ierr
  type(C_PTR) :: transpose_plan
  type(C_PTR) :: dft_plan
  integer(C_INTPTR_T) :: shape_c(3), dest_pid, src_pid
  integer(C_INTPTR_T), parameter :: one = 1
  integer :: ix(3), idx, idy, idz
  real(DP), allocatable :: plane_xy(:,:,:), plane_yx(:,:,:), plane_zy(:,:,:)
  integer, allocatable :: mpi_reqs(:)
  integer :: n_mpi_reqs
  real(DP) :: err
  real(DP) :: rescale
  real(DP) :: h2(1,1,1,1)
  real(DP) :: kx, ky, kz

  if (.not. interface_initialized) then
    call init_comm_infrastructure(nekcomm, shape_x)
  endif
  shape_c = shape_x

  ! convert RHS to coarse mesh
  nelm = size(rhs, 4)
  allocate(rhs_coarse(nelm))
  forall(i = 1 : nelm) rhs_coarse(i) = sum(bm1(:,:,:,i) * rhs(:,:,:,i))

#if defined(SHUFFLE_TEST) || defined(UNITARY_TEST) || defined(NO_FFT)
  do i = 1, nelm
    rhs_coarse(i) = lglel(i)
  enddo
#endif
 
  ! reorder onto sticks
  allocate(plane_xy(0:shape_x(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1) )

  call mesh_to_grid(rhs_coarse, plane_xy, shape_x)

#if defined(UNITARY_TESTT) || defined(NO_FFTT)
  do idx = 0, shape_c(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        ieg = 1 + idx + idy * shape_c(1) + idz * shape_c(1) * shape_c(2)
        plane_xy(idx, idy-idx_in_local_xy, idz-idx_in_local_yz) = ieg
      enddo
    enddo
  enddo
#endif

  ! forward FFT
  rescale = 1._dp
#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(1), int(nin_local_xy * nin_local_yz), &
                                plane_xy, shape_x(1), 1, shape_x(1), &
                                plane_xy, shape_x(1), 1, shape_x(1), &
                                (/FFTW_R2HC/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_xy, plane_xy)
  rescale = rescale * sqrt(real(shape_x(1), kind=DP)) 
#endif
  
  allocate(plane_yx(0:shape_c(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  do i = 0, nin_local_yz - 1
    transpose_plan = fftw_mpi_plan_many_transpose( &
                      shape_c(2), shape_c(1), one, &
                      nin_local_xy, nout_local_xy, &
                      plane_xy(:,:,i), plane_yx(:,:,i), comm_xy, FFTW_ESTIMATE)
    call fftw_mpi_execute_r2r(transpose_plan, plane_xy(:,:,i), plane_yx(:,:,i))
  enddo
  deallocate(plane_xy)

#ifdef NO_FFT
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
#endif

#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(2), int(nout_local_xy * nin_local_yz), &
                                plane_yx, shape_x(2), 1, shape_x(2), &
                                plane_yx, shape_x(2), 1, shape_x(2), &
                                (/FFTW_R2HC/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_yx, plane_yx)
  rescale = rescale * sqrt(real(shape_x(2), kind=DP)) 
#endif

  allocate(plane_zy(0:shape_c(3)-1, 0:nout_local_xy-1, 0:nout_local_yz-1) )
  do i = 0, nout_local_xy-1
    transpose_plan = fftw_mpi_plan_many_transpose( &
                      shape_c(3), shape_c(2), one, &
                      FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                      plane_yx(:,i,:), plane_zy(:,i,:), comm_yz, FFTW_ESTIMATE)
    call fftw_mpi_execute_r2r(transpose_plan, plane_yx(:,i,:), plane_zy(:,i,:))
  enddo
  deallocate(plane_yx)

#ifdef NO_FFT
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
#endif

#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(3), int(nout_local_xy * nout_local_yz), &
                                plane_zy, shape_x(3), 1, shape_x(3), &
                                plane_zy, shape_x(3), 1, shape_x(3), &
                                (/FFTW_REDFT00/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_zy, plane_zy)
  rescale = rescale * sqrt(2.*real(shape_x(3)-1, kind=DP)) 
#endif

  ! Poisson kernel
#if !(defined(UNITARY_TEST) || defined(SHUFFLE_TEST) || defined(NO_FFT))
  do idz = 0, shape_c(3)/nxy - 1
    do idy = 0, nout_local_yz - 1
      do idx = 0, nout_local_xy - 1
        if (idx + idx_out_local_xy <= shape_x(1) / 2) then
          kx = 2*pi*(idx +idx_out_local_xy)/(end_x(1)-start_x(1)) / shape_x(1)
        else
          kx = 2*pi*(shape_x(1) - idx - idx_out_local_xy)/(end_x(1)-start_x(1)) / shape_x(1)
        endif

        if (idy + idx_in_local_xy <= shape_x(2) / 2) then
          ky = 2*pi*(idy +idx_in_local_xy)/(end_x(2)-start_x(2)) / shape_x(2)
        else
          ky = 2*pi*(shape_x(2) - idy - idx_in_local_xy)/(end_x(2)-start_x(2)) / shape_x(2)
        endif

        kz = pi*(idz)/(end_x(3)-start_x(3)) / (shape_x(3) - 1)

        if (kx**2. + ky**2. + kz**2. < 1.e-6_dp) then
          plane_zy(idz,idx,idy) = 0._dp
        else
          plane_zy(idz, idx, idy) = plane_zy(idz, idx, idy) / ( &
            (kz)**2._dp + &
            (ky)**2._dp + &
            (kx)**2._dp)
        endif
      enddo
    enddo
  enddo
#endif


  ! reverse FFT
#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(3), int(nout_local_xy * nout_local_yz), &
                                plane_zy, shape_x(3), 1, shape_x(3), &
                                plane_zy, shape_x(3), 1, shape_x(3), &
                                (/FFTW_REDFT00/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_zy, plane_zy)
  rescale = rescale * sqrt(2.*real(shape_x(3)-1, kind=DP)) 
#endif

  allocate(plane_yx(0:shape_c(2)-1, 0:nout_local_xy-1, 0:nin_local_yz-1) )
  do i = 0, nout_local_xy - 1
  transpose_plan = fftw_mpi_plan_many_transpose( &
                    shape_c(2), shape_c(3), one, &
                    FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                    plane_zy(:,i,:), plane_yx(:,i,:), comm_yz, FFTW_ESTIMATE)
  call fftw_mpi_execute_r2r(transpose_plan, plane_zy(:,i,:), plane_yx(:,i,:))
  enddo
  deallocate(plane_zy)

#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(2), int(nout_local_xy * nin_local_yz), &
                                plane_yx, shape_x(2), 1, shape_x(2), &
                                plane_yx, shape_x(2), 1, shape_x(2), &
                                (/FFTW_HC2R/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_yx, plane_yx)
  rescale = rescale * sqrt(real(shape_x(2), kind=DP)) 
#endif

  allocate(plane_xy(0:shape_c(1)-1, 0:nin_local_xy-1, 0:nin_local_yz-1))
  do i = 0, nin_local_yz - 1
    transpose_plan = fftw_mpi_plan_many_transpose( &
                      shape_c(1), shape_c(2), one, &
                      FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                      plane_yx(:,:,i), plane_xy(:,:,i), comm_xy, FFTW_ESTIMATE)
    call fftw_mpi_execute_r2r(transpose_plan, plane_yx(:,:,i), plane_xy(:,:,i))
  enddo

#ifndef NO_FFT
  dft_plan = fftw_plan_many_r2r(1, shape_x(1), int(nin_local_xy * nin_local_yz), &
                                plane_xy, shape_x(1), 1, shape_x(1), &
                                plane_xy, shape_x(1), 1, shape_x(1), &
                                (/FFTW_HC2R/), FFTW_ESTIMATE)
  call fftw_execute_r2r(dft_plan, plane_xy, plane_xy)
  rescale = rescale * sqrt(real(shape_x(1), kind=DP)) 
#endif


  plane_xy = plane_xy * (1._dp/ rescale)
#ifdef UNITARY_TEST
  do idx = 0, shape_x(1)-1
    do idy = idx_in_local_xy, idx_in_local_xy + nin_local_xy - 1
      do idz = idx_in_local_yz, idx_in_local_yz + nin_local_yz - 1
        ieg = 1 + idx + idy * shape_x(1) + idz * shape_x(1) * shape_x(2)
        err = abs(plane_xy(idx, idy-idx_in_local_xy, idz - idx_in_local_yz) - ieg)
        if (err > 0.001) then
          write(*,*) "WARNING: spectral transform not unitary", nid, idx, idy, idz, ieg, err
          exit
        endif
      enddo
    enddo
  enddo
  if (nid == 0) write(*,*) "Passed unitary test"
#endif

  allocate(soln_coarse(nelm)); soln_coarse = 0._dp
  ! reorder to local elements
  call grid_to_mesh(plane_xy, soln_coarse, shape_x)

#ifdef SHUFFLE_TEST
  do i = 1, nelm
    err = abs(soln_coarse(i) - lglel(i))
    if (err > 0.001) then
      write(*,*) "WARNING: shuffle not working", nid, i, nid*nelm + i, err
      exit
    endif
  enddo
  if (nid == 0) write(*,*) "Passed shuffle test"
#endif

  ! populate U
  forall(i = 1 : nelm) u(:,:,:,i) = binvm1(:,:,:,i) * soln_coarse(i)
#if defined(UNITARY_TEST) || defined(SHUFFLE_TEST) || defined(NO_FFT)
  u = 0._dp
#endif

  ! update residual
  allocate(tmp_fine(size(u,1), size(u,2), size(u,3), size(u,4)))
  h2 = 0._dp
  call axhelm (tmp_fine, u, h1, h2, imsh, isd)
  tmp_fine = tmp_fine * mask
  call dssum   (tmp_fine)

  if (nid == 0) write(*,*) "RHS before: ", sqrt(sum(rhs * rhs))
  RHS = RHS - tmp_fine 
  if (nid == 0) write(*,*) "RHS after : ", sqrt(sum(rhs * rhs))
 
 
end subroutine spectral_solve

!> \brief one-time setup of communication infrastructure for poisson_mod
subroutine init_comm_infrastructure(comm_world, shape_x)
  use fftw3, only : FFTW_MPI_DEFAULT_BLOCK
  use fftw3, only : fftw_mpi_local_size_many_transposed
  integer, intent(in) :: comm_world !>!< Communicator in which to setup solver
  integer, intent(in) :: shape_x(3) !>!< Shape of mesh

  integer(C_INTPTR_T) :: shape_c(3)
  integer(C_INTPTR_T), parameter :: one = 1
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
  use parallel, only : lglel, gllel, gllnid
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

end module poisson
