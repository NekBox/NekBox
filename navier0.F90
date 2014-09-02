!---------------------------------------------------------------------------
!> \brief Decide strategy for E-solver
!---------------------------------------------------------------------------
SUBROUTINE ESTRAT
  use size_m, only : nid
  use esolv, only : iesolv
  use input, only : ifsplit, solver_type, param
  use tstep, only : istep
  implicit none

  IESOLV = 1
  if (ifsplit) iesolv=0

  solver_type='itr'
  if (param(116) /= 0) solver_type='fdm'
!   if (param(90).ne.0)  solver_type='itn'

!   The following change recognizes that geometry is logically
!   tensor-product, but deformed:  pdm = Preconditioner is fdm

  if (param(59) /= 0 .AND. solver_type == 'fdm') solver_type='pdm'

  if (istep < 2 .AND. nid == 0) write(6,10) iesolv,solver_type
  10 format(2X,'E-solver strategy: ',I2,1X,A)

  RETURN
END SUBROUTINE ESTRAT
