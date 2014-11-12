!==============================================================================
!> \file poisson_mod.F90
!! \brief Spectral coarse solver for poisson equation 
!! \date November 2014
!! \author Max Hutchinson
!!
!! This module implements a coarse solve (preconditioner) for the pressure
!! Poisson equation.
module poisson
  implicit none

  public spectral_solve
  private

contains

!> \brief 
subroutine spectral_solve(u,rhs,h1,mask,mult,imsh,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, nelv, nelt, ndim, nid
  use ctimer, only : icalld, thmhz, nhmhz, etime1, dnekclock
  use fdmh1, only : kfldfdm
  use input, only : ifsplit, param
  use geom, only : bm1, binvm1
  use tstep, only : istep, nelfld, ifield

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
  integer :: i

  ! convert RHS to coarse mesh
  nelm = size(rhs, 4)
  allocate(rhs_coarse(nelm))
  forall(i = 1 : nelm) rhs_coarse(i) = sum(bm1(:,:,:,i) * rhs(:,:,:,i))
 
  ! reorder onto sticks

  ! forward FFT

  ! Poisson kernel
  
  ! reverse FFT

  ! reorder to local elements
  allocate(soln_coarse(nelm))

  ! populate U
  forall(i = 1 : nelm) u(:,:,:,i) = binvm1(:,:,:,i) * soln_coarse(i)
  allocate(tmp_fine(size(u,1), size(u,2), size(u,3), size(u,4)))
  call axhelm (tmp_fine, u, h1, 0, imsh, isd)
  RHS = RHS - tmp_fine 
 
 
end subroutine spectral_solve

end module poisson
