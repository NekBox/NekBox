!-----------------------------------------------------------------------
!> \brief  Tensor product application of v = (C x B x A) u .
!!  NOTE -- the transpose of B & C must be input, rather than B & C.
!!  -  scratch arrays: w(nu*nu*nv)
subroutine tensr3(v,nv,u,nu,A,Bt,Ct,w)
  use kinds, only : DP
  use size_m, only : nid
  use input, only : if3d
  implicit none

  integer :: nv, nu
  real(DP) :: v(1),u(1)
  real(DP) :: A(1),Bt(1),Ct(1)
  real(DP) :: w(1)

  integer :: nuv, nvv, k, l, iz

  if (nu > nv) then
      write(6,*) nid,nu,nv,' ERROR in tensr3. Contact P.Fischer.'
      write(6,*) nid,nu,nv,' Memory problem.'
      call exitt
  endif

  if (if3d) then
      nuv = nu*nv
      nvv = nv*nv
      call mxm(A,nv,u,nu,v,nu*nu)
      k=1
      l=1
      do iz=1,nu
          call mxm(v(k),nv,Bt,nu,w(l),nv)
          k=k+nuv
          l=l+nvv
      enddo
      call mxm(w,nvv,Ct,nu,v,nv)
  else
      call mxm(A,nv,u,nu,w,nu)
      call mxm(w,nv,Bt,nu,v,nv)
  endif
  return
end subroutine tensr3

!-----------------------------------------------------------------------
