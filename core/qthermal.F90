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
  use mass, only : bm1, binvm1
  use soln, only : t, vdiff, vtrans, qtl
  use tstep, only : ifield
  implicit none

  real(DP) :: w2(LX1,LY1,LZ1,LELT)
  real(DP) :: tx(LX1,LY1,LZ1,LELT), ty(LX1,LY1,LZ1,LELT), tz(LX1,LY1,LZ1,LELT)

  integer :: ntot, ifld_save

  ntot = nx1*ny1*nz1*nelv

  if ( .NOT. iflomach) then
      call rzero(qtl,ntot)
      return
  endif

  ifld_save = ifield

! - - Assemble RHS of T-eqn
  ifield=2
  call setqvol (QTL) ! volumetric heating source
  call col2    (QTL,BM1,ntot)

  ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
  call opgrad  (tx,ty,tz,T)
  call opdssum (tx,ty,tz)
  call opcolv  (tx,ty,tz,binvm1)
  call opcolv  (tx,ty,tz,vdiff(1,1,1,1,2))
  call opdiv   (w2,tx,ty,tz)

  call add2    (QTL,w2,ntot)

  call col3    (w2,vtrans(1,1,1,1,2),T,ntot)
  call invcol2 (QTL,w2,ntot)

  call dssum   (QTL,nx1,ny1,nz1)
  call col2    (QTL,binvm1,ntot)

  ifield = ifld_save

  return
end subroutine qthermal
