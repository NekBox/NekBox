    SUBROUTINE ESTRAT
!---------------------------------------------------------------------------

!     Decide strategy for E-solver

!---------------------------------------------------------------------------
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    IESOLV = 1
    if (ifsplit) iesolv=0

    solver_type='itr'
    if (param(116) /= 0) solver_type='fdm'
!     if (param(90).ne.0)  solver_type='itn'

!     The following change recognizes that geometry is logically
!     tensor-product, but deformed:  pdm = Preconditioner is fdm

    if (param(59) /= 0 .AND. solver_type == 'fdm') solver_type='pdm'

    if (istep < 2 .AND. nid == 0) write(6,10) iesolv,solver_type
    10 format(2X,'E-solver strategy: ',I2,1X,A)

    RETURN
    END SUBROUTINE ESTRAT

!-----------------------------------------------------------------------------
!> \brief Initialize E-solver
!-----------------------------------------------------------------------------
SUBROUTINE EINIT
  implicit none
  CALL ESTRAT
  RETURN
END SUBROUTINE EINIT
!-----------------------------------------------------------------------
