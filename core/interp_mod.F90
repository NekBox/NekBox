module interp
  use kinds, only : DP
  use size_m, only : lxd, ldim
  implicit none

  integer, parameter :: ldg=lxd**3, lwkd=4*lxd*lxd
  integer, parameter :: ld=2*lxd
  integer, parameter :: ldw = 2*(ld**ldim)

  real(DP), save :: jgl(ldg), jgt(ldg)
  real(DP), save :: dgl(ldg), dgt(ldg)
  integer, save :: pjgl(0:ld*ld) = 0
  integer, save ::  pdg(0:ld*ld) = 0

  real(DP) :: w(ld**ldim,2)
  real(DP) :: wkd(lwkd)

contains

!-----------------------------------------------------------------------
!> \brief Get pointer to jgl() for interpolation pair (mx,md)
!!
!! The interpolation matrices jgl, jgt are being memoized.
!!  pjgl is a map from (mx,md) pair to (ip) index of jgl, jgt
subroutine get_int_ptr (ip, mx,md) ! GLL-->GL pointer
  use kinds, only : DP
  use size_m
  implicit none

  integer, intent(out) :: ip
  integer, intent(in) :: mx, md

  integer :: ij, nstore, nwrkd

  ij = md + ld*(mx-1)
  ip = pjgl(ij)

  if (ip == 0) then

      nstore   = pjgl(0)
      pjgl(ij) = nstore+1
      nstore   = nstore + md*mx
      pjgl(0)  = nstore
      ip       = pjgl(ij)
  
      nwrkd = mx + md
      call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
      call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
  
      call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
  endif

  return
end subroutine get_int_ptr

!-----------------------------------------------------------------------
!> \brief Get pointer to GL-GL interpolation dgl() for pair (mx,md)
subroutine get_dgl_ptr (ip, mx,md)
  use kinds, only : DP
  use size_m, only : lxd 
  implicit none


  integer, intent(out) :: ip
  integer, intent(in) :: mx, md


  integer :: ij, nstore, nwrkd

  ij = md + ld*(mx-1)
  ip = pdg (ij)

  if (ip == 0) then

      nstore   = pdg (0)
      pdg (ij) = nstore+1
      nstore   = nstore + md*mx
      pdg (0)  = nstore
      ip       = pdg (ij)
  
      nwrkd = mx + md
      call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
      call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
  
      call gen_dgl(dgl(ip),dgt(ip),md,mx,wkd)
  endif

  return
end subroutine get_dgl_ptr


end module 

!-----------------------------------------------------------------------
!> \brief Generate interpolation from np GLL points to mp GL points
!!   jgl  = interpolation matrix, mapping from velocity nodes to pressure
!!   jgt  = transpose of interpolation matrix
!!   w    = work array of size (np+mp)
!!   np   = number of points on GLL grid
!!   mp   = number of points on GL  grid
subroutine gen_int(jgl,jgt,mp,np,w)
  use kinds, only : DP
  use speclib, only : zwgl, zwgll
  implicit none

  integer, intent(in) :: mp, np
  real(DP) :: jgl(mp,np),jgt(np*mp),w(*)

  integer :: iz, id, n, i, j

  iz = 1
  id = iz + np

  call zwgll (w(iz),jgt,np)
  call zwgl  (w(id),jgt,mp)

  n  = np-1
  do i=1,mp
      call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
      do j=1,np
          jgl(i,j) = jgt(j)                  !  Interpolation matrix
      enddo
  enddo

  call transpose(jgt,np,jgl,mp)

  return
end subroutine gen_int


!-----------------------------------------------------------------------
!> \brief Generate derivative from np GL points onto mp GL points
!!  dgl  = interpolation matrix, mapping from velocity nodes to pressure
!!  dgt  = transpose of interpolation matrix
!!  w    = work array of size (3*np+mp)
!!  np   = number of points on GLL grid
!!  mp   = number of points on GL  grid
subroutine gen_dgl(dgl,dgt,mp,np,w)
  use kinds, only : DP
  use speclib, only : zwgl
  implicit none

  integer, intent(in) :: mp, np
  real(DP) :: dgl(mp,np),dgt(np*mp),w(*)

  integer :: iz, id, ndgt, ldgt, n, i, j

  iz = 1
  id = iz + np

  call zwgl  (w(iz),dgt,np)  ! GL points
  call zwgl  (w(id),dgt,mp)  ! GL points

  ndgt = 2*np
  ldgt = mp*np
  call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')

  n  = np-1
  do i=1,mp
      call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
      do j=1,np
          dgl(i,j) = dgt(np+j)                       ! Derivative matrix
      enddo
  enddo

  call transpose(dgt,np,dgl,mp)

  return
end subroutine gen_dgl
!-----------------------------------------------------------------------
!> Check array limits
subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
  use size_m, only : nid            ! need nid
  implicit none
  character(5), intent(in) ::  avar5,lvar5
  character(10), intent(in) :: sub_name10
  integer, intent(in) :: n, m

  if (n > m) then
      write(6,1) nid,n,m,avar5,lvar5,sub_name10
      1 format(i8,' ERROR: :',2i12,2(1x,a5),1x,a10)
      call exitti('lim_chk problem. $',n)
  endif

  return
end subroutine lim_chk
!-----------------------------------------------------------------------
