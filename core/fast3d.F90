!-----------------------------------------------------------------------
!> \brief Generate fast diagonalization matrices for each element
subroutine gen_fast_spacing()
  use input, only : param
  use speclib, only : zwgll, zwgl
  implicit none

#if 0
  real :: x(nx1,ny1,nz1,nelv)
  real :: y(nx1,ny1,nz1,nelv)
  real :: z(nx1,ny1,nz1,nelv)

  use size_m
  use parallel
  use soln
  use wz_m

  integer, parameter(lxx=lx1*lx1)

  common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4) &
  , llr(lelt),lls(lelt),llt(lelt) &
  , lmr(lelt),lms(lelt),lmt(lelt) &
  , lrr(lelt),lrs(lelt),lrt(lelt)
  real :: lr ,ls ,lt
  real :: llr,lls,llt
  real :: lmr,lms,lmt
  real :: lrr,lrs,lrt

  integer :: lbr,rbr,lbs,rbs,lbt,rbt,e

  real :: axwt(lx2)

  ierr = 0
#endif
  if (param(44) == 1) then
    write(*,*) "Oops: param(44)"
#if 0
  !                                    __ __ __
  !        Now, for each element, compute lr,ls,lt between specified planes
  
      n1 = nx2
      n2 = nx2+1
      nz0 = 1
      nzn = 1
      if (if3d) then
          nz0= 0
          nzn=n2
      endif
      eps = 1.e-7
      if (wdsize == 8)  eps = 1.e-14
  
  !        Find mean spacing between "left-most" planes
      call plane_space2(llr,lls,llt, 0,wxm2,x,y,z,n1,n2,nz0,nzn)
  
  !        Find mean spacing between "middle" planes
      call plane_space (lmr,lms,lmt, 1,n1,wxm2,x,y,z,n1,n2,nz0,nzn)
  
  !        Find mean spacing between "right-most" planes
      call plane_space2(lrr,lrs,lrt,n2,wxm2,x,y,z,n1,n2,nz0,nzn)
#endif 
  else
      call load_semhat_weighted    !   Fills the SEMHAT arrays
  endif

  return
end subroutine gen_fast_spacing

!-----------------------------------------------------------------------
!> \brief Here, spacing is based on harmonic mean.  pff 2/10/07
subroutine plane_space(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
  use kinds, only : DP
  use size_m, only : nelv
  use input, only : if3d
  implicit none

  integer :: nx, nz0, nxn, nzn, i1, i2
  real(DP) :: w(*),lr(*),ls(*),lt(*)
  real(DP) :: x(0:nxn,0:nxn,nz0:nzn,*)
  real(DP) :: y(0:nxn,0:nxn,nz0:nzn,*)
  real(DP) :: z(0:nxn,0:nxn,nz0:nzn,*)
  real(DP) :: lr2,ls2,lt2

  integer :: ny, nz, j1, k1, j2, k2, ie, k, j, i
  real(DP) :: wsum, weight

!   Now, for each element, compute lr,ls,lt between specified planes
  ny = nx
  nz = nx
  j1 = i1
  k1 = i1
  j2 = i2
  k2 = i2

  do ie=1,nelv
  
      if (if3d) then
          lr2  = 0.
          wsum = 0.
          do k=1,nz
              do j=1,ny
                  weight = w(j)*w(k)
              !              lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
              !    $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
              !    $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
              !    $                     *   weight
                  lr2  = lr2  +   weight / &
                  ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2 &
                  +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2 &
                  +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
                  wsum = wsum + weight
              enddo
          enddo
          lr2     = lr2/wsum
          lr(ie)  = 1./sqrt(lr2)
      
          ls2 = 0.
          wsum = 0.
          do k=1,nz
              do i=1,nx
                  weight = w(i)*w(k)
              !              ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
              !    $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
              !    $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
              !    $                     *   weight
                  ls2  = ls2  +   weight / &
                  ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2 &
                  +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2 &
                  +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
                  wsum = wsum + weight
              enddo
          enddo
          ls2     = ls2/wsum
          ls(ie)  = 1./sqrt(ls2)
      
          lt2 = 0.
          wsum = 0.
          do j=1,ny
              do i=1,nx
                  weight = w(i)*w(j)
              !              lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
              !    $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
              !    $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
              !    $                     *   weight
                  lt2  = lt2  +   weight / &
                  ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2 &
                  +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2 &
                  +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
                  wsum = wsum + weight
              enddo
          enddo
          lt2     = lt2/wsum
          lt(ie)  = 1./sqrt(lt2)
      
      else              ! 2D
          lr2 = 0.
          wsum = 0.
          do j=1,ny
              weight = w(j)
          !              lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
          !    $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
          !    $                     *   weight
              lr2  = lr2  + weight / &
              ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2 &
              + (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
              wsum = wsum + weight
          enddo
          lr2     = lr2/wsum
          lr(ie)  = 1./sqrt(lr2)
      
          ls2 = 0.
          wsum = 0.
          do i=1,nx
              weight = w(i)
          !              ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
          !    $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
          !    $                     *   weight
              ls2  = ls2  + weight / &
              ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2 &
              +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
              wsum = wsum + weight
          enddo
          ls2     = ls2/wsum
          ls(ie)  = 1./sqrt(ls2)
      !           write(6,*) 'lrls',ie,lr(ie),ls(ie)
      endif
  enddo
  return
end subroutine plane_space

!-----------------------------------------------------------------------
subroutine get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ierr)
  use size_m, only : ndim
  use input, only : cbc
  use parallel, only : lglel
  use topol, only : eface
  use tstep, only : ifield
  implicit none

  integer ::                lbr,rbr,lbs,rbs,lbt,rbt,e,bsym
  integer :: fbc(6)
  integer :: iface, ied, ibc, ierr

!   ibc = 0  <==>  Dirichlet
!   ibc = 1  <==>  Dirichlet, outflow (no extension)
!   ibc = 2  <==>  Neumann,

  do iface=1,2*ndim
      ied = eface(iface)
      ibc = -1

!max        if (ifmhd) call mhd_bc_dn(ibc,iface,e) ! can be overwritten by 'mvn'

      if (cbc(ied,e,ifield) == '   ') ibc = 0
      if (cbc(ied,e,ifield) == 'E  ') ibc = 0
      if (cbc(ied,e,ifield) == 'msi') ibc = 0
      if (cbc(ied,e,ifield) == 'MSI') ibc = 0
      if (cbc(ied,e,ifield) == 'P  ') ibc = 0
      if (cbc(ied,e,ifield) == 'p  ') ibc = 0
      if (cbc(ied,e,ifield) == 'O  ') ibc = 1
      if (cbc(ied,e,ifield) == 'ON ') ibc = 1
      if (cbc(ied,e,ifield) == 'o  ') ibc = 1
      if (cbc(ied,e,ifield) == 'on ') ibc = 1
      if (cbc(ied,e,ifield) == 'MS ') ibc = 1
      if (cbc(ied,e,ifield) == 'ms ') ibc = 1
      if (cbc(ied,e,ifield) == 'MM ') ibc = 1
      if (cbc(ied,e,ifield) == 'mm ') ibc = 1
      if (cbc(ied,e,ifield) == 'mv ') ibc = 2
      if (cbc(ied,e,ifield) == 'mvn') ibc = 2
      if (cbc(ied,e,ifield) == 'v  ') ibc = 2
      if (cbc(ied,e,ifield) == 'V  ') ibc = 2
      if (cbc(ied,e,ifield) == 'W  ') ibc = 2
      if (cbc(ied,e,ifield) == 'SYM') ibc = bsym
      if (cbc(ied,e,ifield) == 'SL ') ibc = 2
      if (cbc(ied,e,ifield) == 'sl ') ibc = 2
      if (cbc(ied,e,ifield) == 'SHL') ibc = 2
      if (cbc(ied,e,ifield) == 'shl') ibc = 2
      if (cbc(ied,e,ifield) == 'A  ') ibc = 2
      if (cbc(ied,e,ifield) == 'S  ') ibc = 2
      if (cbc(ied,e,ifield) == 's  ') ibc = 2
      if (cbc(ied,e,ifield) == 'J  ') ibc = 0
      if (cbc(ied,e,ifield) == 'SP ') ibc = 0

      fbc(iface) = ibc

      if (ierr == -1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
      1 format(2i3,i8,i3,2x,a3,'  get_fast_bc_error')

  enddo

  if (ierr == -1) call exitti('Error A get_fast_bc$',e)

  lbr = fbc(1)
  rbr = fbc(2)
  lbs = fbc(3)
  rbs = fbc(4)
  lbt = fbc(5)
  rbt = fbc(6)

  ierr = 0
  if (ibc < 0) ierr = lglel(e)

!   write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
! 6 format(i5,2x,4i3,3x,4(1x,a3),'  get_fast_bc')

  return
  end subroutine get_fast_bc

!-----------------------------------------------------------------------
!> \brief Note that this routine performs the following matrix multiplies
!!  after getting the matrices back from semhat:
!!  dgl = bgl dgl
!!  jgl = bgl jgl
subroutine load_semhat_weighted()    !   Fills the SEMHAT arrays
  use size_m, only : nx1
  use semhat, only : ah, bh, ch, dh, zh, dph, jph, bgl, zgl, dgl, jgl, wh
  implicit none
  integer :: nr

  nr = nx1-1
  call generate_semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,nr,wh)
  call do_semhat_weight(jgl,dgl,bgl,nr)

  return
end subroutine load_semhat_weighted

!----------------------------------------------------------------------
!> \brief re-weight jgl and dgl with bgl
subroutine do_semhat_weight(jgl,dgl,bgl,n)
  use kinds, only : DP
  implicit none
  integer :: n
  real(DP) :: bgl(1:n-1),jgl(1:n-1,0:n),dgl(1:n-1,0:n)
  integer :: i, j
  do j=0,n
      do i=1,n-1
          jgl(i,j)=bgl(i)*jgl(i,j)
      enddo
  enddo
  do j=0,n
      do i=1,n-1
          dgl(i,j)=bgl(i)*dgl(i,j)
      enddo
  enddo
  return
end subroutine do_semhat_weight

!-----------------------------------------------------------------------
!> \brief Generate matrices for single element, 1D operators:
!!
!!    a    = Laplacian
!!    b    = diagonal mass matrix
!!    c    = convection operator b*d
!!    d    = derivative matrix
!!    dgll = derivative matrix,    mapping from pressure nodes to velocity
!!    jgll = interpolation matrix, mapping from pressure nodes to velocity
!!    z    = GLL points
!!
!!    zgl  = GL points
!!    bgl  = diagonal mass matrix on GL
!!    dgl  = derivative matrix,    mapping from velocity nodes to pressure
!!    jgl  = interpolation matrix, mapping from velocity nodes to pressure
!!
!!    n    = polynomial degree (velocity space)
!!    w    = work array of size 2*n+2
!!
!! Currently, this is set up for pressure nodes on the interior GLL pts.
subroutine generate_semhat(a,b,c,d,z,dgll,jgll,bgl,zgl,dgl,jgl,n,w)
  use kinds, only : DP
  use speclib, only : zwgll, zwgl
  implicit none

  integer :: n
  real(DP) :: a(0:n,0:n),b(0:n),c(0:n,0:n),d(0:n,0:n),z(0:n)
  real(DP) :: dgll(0:n,1:n-1),jgll(0:n,1:n-1)

  real(DP) :: bgl(1:n-1),zgl(1:n-1)
  real(DP) :: dgl(1:n-1,0:n),jgl(1:n-1,0:n)

  real(DP) :: w(0:(n+1)*2)

  integer :: np, nm, n2
  integer :: i, j, k

  np = n+1
  nm = n-1
  n2 = n-2

  call zwgll (z,b,np)

  do i=0,n
      call fd_weights_full(z(i),z,n,1,w)
      do j=0,n
          d(i,j) = w(j+np)                   !  Derivative matrix
      enddo
  enddo

  if (n == 1) return                       !  No interpolation for n=1

  do i=0,n
      call fd_weights_full(z(i),z(1),n2,1,w(1))
      do j=1,nm
          jgll(i,j) = w(j   )                  !  Interpolation matrix
          dgll(i,j) = w(j+nm)                  !  Derivative    matrix
      enddo
  enddo

  a = 0._dp
  do j=0,n
      do i=0,n
          do k=0,n
              a(i,j) = a(i,j) + d(k,i)*b(k)*d(k,j)
          enddo
          c(i,j) = b(i)*d(i,j)
      enddo
  enddo

  call zwgl (zgl,bgl,nm)

  do i=1,n-1
      call fd_weights_full(zgl(i),z,n,1,w)
      do j=0,n
          jgl(i,j) = w(j   )                  !  Interpolation matrix
          dgl(i,j) = w(j+np)                  !  Derivative    matrix
      enddo
  enddo

  return
end subroutine generate_semhat

!-----------------------------------------------------------------------
!> \brief This routine evaluates the derivative based on all points
!!   in the stencils.  It is more memory efficient than "fd_weights"
!!
!!   This set of routines comes from the appendix of
!!   A Practical Guide to Pseudospectral Methods, B. Fornberg
!!   Cambridge Univ. Press, 1996.   (pff)
!!
!!   Input parameters:
!!     xx -- point at wich the approximations are to be accurate
!!     x  -- array of x-ordinates:   x(0:n)
!!     n  -- polynomial degree of interpolant (# of points := n+1)
!!     m  -- highest order of derivative to be approxxmated at xi
!!   Output:
!!     c  -- set of coefficients c(0:n,0:m).
!!           c(j,k) is to be applied at x(j) when
!!           the kth derivative is approxxmated by a
!!           stencil extending over x(0),x(1),...x(n).
subroutine fd_weights_full(xx,x,n,m,c)
  use kinds, only : DP
  implicit none

  integer :: n, m
  real(DP) :: xx, x(0:n),c(0:n,0:m)

  integer :: i, j, k, mn
  real(DP) :: c1, c2, c3, c4, c5

  c1       = 1.
  c4       = x(0) - xx

  do k=0,m
      do j=0,n
          c(j,k) = 0.
      enddo
  enddo
  c(0,0) = 1.

  do i=1,n
      mn = min(i,m)
      c2 = 1.
      c5 = c4
      c4 = x(i)-xx
      do j=0,i-1
          c3 = x(i)-x(j)
          c2 = c2*c3
          do k=mn,1,-1
              c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
          enddo
          c(i,0) = -c1*c5*c(i-1,0)/c2
          do k=mn,1,-1
              c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
          enddo
          c(j,0) = c4*c(j,0)/c3
      enddo
      c1 = c2
  enddo
!   call outmat(c,n+1,m+1,'fdw',n)
  return
end subroutine fd_weights_full

!-----------------------------------------------------------------------
subroutine swap_lengths()
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv
  use size_m, only : nx1, nelv
  use geom, only : xm1, ym1, zm1
  use hsmg, only : llr, lrr, lmr, lls, lms, lrs, llt, lmt, lrt
  use input, only : if3d
  use wz_m, only : wxm1
  implicit none

  real(DP), allocatable :: l(:,:,:,:)
  integer :: e, n2, nz0, nzn, nx, n, j, k

  allocate(l(lx1, ly1, lz1, lelv))

  n2 = nx1-1
  nz0 = 1
  nzn = 1
  nx  = nx1-2
  if (if3d) then
      nz0 = 0
      nzn = n2
  endif
  call plane_space(lmr,lms,lmt,0,n2,wxm1,xm1,ym1,zm1,nx,n2,nz0,nzn)

  n=n2+1
  if (if3d) then
      do e=1,nelv
          do j=2,n2
              do k=2,n2
                  l(1,k,j,e) = lmr(e)
                  l(n,k,j,e) = lmr(e)
                  l(k,1,j,e) = lms(e)
                  l(k,n,j,e) = lms(e)
                  l(k,j,1,e) = lmt(e)
                  l(k,j,n,e) = lmt(e)
              enddo
          enddo
      enddo
      call dssum(l)
      do e=1,nelv
          llr(e) = l(1,2,2,e)-lmr(e)
          lrr(e) = l(n,2,2,e)-lmr(e)
          lls(e) = l(2,1,2,e)-lms(e)
          lrs(e) = l(2,n,2,e)-lms(e)
          llt(e) = l(2,2,1,e)-lmt(e)
          lrt(e) = l(2,2,n,e)-lmt(e)
      enddo
  else
      do e=1,nelv
          do j=2,n2
              l(1,j,1,e) = lmr(e)
              l(n,j,1,e) = lmr(e)
              l(j,1,1,e) = lms(e)
              l(j,n,1,e) = lms(e)
          !           call outmat(l(1,1,1,e),n,n,' L    ',e)
          enddo
      enddo
  !        call outmat(l(1,1,1,25),n,n,' L    ',25)
      call dssum(l)
  !        call outmat(l(1,1,1,25),n,n,' L    ',25)
      do e=1,nelv
      !           call outmat(l(1,1,1,e),n,n,' L    ',e)
          llr(e) = l(1,2,1,e)-lmr(e)
          lrr(e) = l(n,2,1,e)-lmr(e)
          lls(e) = l(2,1,1,e)-lms(e)
          lrs(e) = l(2,n,1,e)-lms(e)
      enddo
  endif
  return
end subroutine swap_lengths

!----------------------------------------------------------------------
!> \brief zero the eth row of a
subroutine row_zero(a,m,n,e)
  use kinds, only : DP
  implicit none

  integer :: m,n,e
  real(DP) :: a(m,n)
  integer :: j
  do j=1,n
      a(e,j)=0.
  enddo
  return
end subroutine row_zero
!-----------------------------------------------------------------------
!> \brief Reorder vector using temporary buffer
subroutine swap(b,ind,n,temp)
  use kinds, only : DP
  implicit none
  real(DP) :: B(*),TEMP(*)
  integer :: n, IND(*)
  integer :: i, jj

!***
!***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
!***  INTO ITEM(I), WHERE JJ=IND(I).
!***
  DO I=1,N
      JJ=IND(I)
      TEMP(I)=B(JJ)
  END DO
  DO I=1,N
      B(I)=TEMP(I)
  END DO
  RETURN
end subroutine swap

