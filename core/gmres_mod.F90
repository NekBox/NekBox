!> cleaned
module gmres
  implicit none

    real, allocatable :: x(:),r(:),w(:)
    real, allocatable :: h(:,:),gamma(:)
    real, allocatable :: c(:),s(:) ! store the Givens rotations
    real, allocatable :: v(:,:) ! stores the orthogonal Krylov subspace basis
    real, allocatable :: z(:,:) ! Z = M**(-1) V

    real, allocatable :: ml(:), mu(:)
              
  contains

  subroutine init_gmres()
    use size_m
    implicit none

    allocate(x(lx2*ly2*lz2*lelv) &
    , r(lx2*ly2*lz2*lelv), w(lx2*ly2*lz2*lelv) &
    , h(lgmres,lgmres), gamma(lgmres+1) &
    , c(lgmres), s(lgmres))

    allocate(v(lx2*ly2*lz2*lelv,lgmres)) ! verified
    allocate(z(lx2*ly2*lz2*lelv,lgmres)) ! verified

    allocate(ml(lx2*ly2*lz2*lelv), mu(lx2*ly2*lz2*lelv)) ! verified

  end subroutine init_gmres

end module gmres
