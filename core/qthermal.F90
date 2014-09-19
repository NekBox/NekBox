!> \file qthermal.F90 \copybrief qthermal

!-------------------------------------------------------------------------
!> \brief Compute the thermal divergence QTL
!! QTL := div(v) = -1/rho * Drho/Dt
!! If we use the ideal gas law and assume
!! that p,R is const we end up with
!! QTL = 1/(rho*cp) rho*cp*DT/Dt
!! where rho*cp*DT/Dt represents the RHS of the
!! energy equation expressed in terms of temperature.
subroutine qthermal()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, lx1, ly1, lz1, lelt
  use input, only : iflomach
  use geom, only : bm1, binvm1
  use soln, only : t, vdiff, vtrans, qtl
  use tstep, only : ifield
  implicit none

  real(DP), allocatable :: w2(:,:,:,:), tx(:,:,:,:), ty(:,:,:,:), tz(:,:,:,:)

  integer :: ntot, ifld_save

  ntot = nx1*ny1*nz1*nelv

  if ( .NOT. iflomach) then
      qtl = 0._dp
      return
  endif

  ifld_save = ifield

! - - Assemble RHS of T-eqn
  ifield=2
  call setqvol (QTL) ! volumetric heating source
  qtl = qtl * bm1

  ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
  allocate(tx(LX1,LY1,LZ1,LELT), ty(LX1,LY1,LZ1,LELT), tz(LX1,LY1,LZ1,LELT))
  call opgrad  (tx,ty,tz,T)
  call opdssum (tx,ty,tz)
  tx = tx * binvm1 * vdiff(:,:,:,:,2)
  ty = ty * binvm1 * vdiff(:,:,:,:,2)
  tz = tz * binvm1 * vdiff(:,:,:,:,2)

  allocate(w2(LX1,LY1,LZ1,LELT))
  call opdiv   (w2,tx,ty,tz)

  qtl = qtl + w2
  deallocate(w2,tx,ty,tz)

  qtl = qtl / vtrans(:,:,:,:,2) * T(:,:,:,:,1)

  call dssum   (QTL)
  qtl = qtl * binvm1

  ifield = ifld_save

  return
end subroutine qthermal
