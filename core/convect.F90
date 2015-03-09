!-----------------------------------------------------------------------
!    Stability limits:
!    AB3:    .7236                     w/safety (1.2):   .603
!    RK3:    1.73   (sqrt 3)           w/safety (1.2):   1.44
!    RK4:    2.828                     w/safety (1.2):   2.36
!    SEM Safety factor:  1.52 for N=3
!                     <  1.20 for N=16
!                     ~  1.16 for N=256
!-----------------------------------------------------------------------
subroutine setup_convect(igeom)
  use dealias, only : vxd, vyd, vzd
  use input, only : param, ifchar, ifcons, ifpert
  use soln, only : vx, vy, vz
  implicit none

  integer, intent(in) :: igeom

  if (igeom == 1) return
  if (param(99) < 0) return ! no dealiasing

  if (ifchar) then
    write(*,*) "Oops: ifchar"
#if 0
    nelc = nelv
    if (ifmhd) nelc = max(nelv,nelfld(ifldmhd))
    if (ifmhd) call exitti('no characteristics for mhd yet$',istep)

    ifnew = .TRUE. 
    if (igeom > 2) ifnew = .FALSE. 
    call set_conv_char(ct_vx,c_vx,vx,vy,vz,nelc,time,ifnew)
#endif

  else

    if ( .NOT. ifpert) then
        if (ifcons) then
!max                call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
        else
            call set_convect_new  (vxd,vyd,vzd,vx,vy,vz)
        endif
    endif

  endif

  return
end subroutine setup_convect
!-----------------------------------------------------------------------
!> \brief GLL interpolation from mx to md.
!! If idir ^= 0, then apply transpose operator  (md to mx)
subroutine intp_rstd(ju,u,mx,md,if3d,idir) ! GLL->GL interpolation
  use kinds, only : DP
  use size_m, only : lxd, ldim
  use ctimer, only : nintp, tintp, intp_flop, intp_mop, dnekclock
  implicit none

  real(DP), intent(out) :: ju(*)
  real(DP), intent(in)  :: u(*)
  integer,  intent(in)  :: mx, md, idir
  logical,  intent(in)  :: if3d

  integer, parameter :: ldg=lxd**3
  real(DP), save :: jgl(ldg), jgt(ldg)

  integer, parameter :: ld=2*lxd
  real(DP) :: w(ld**ldim,2), etime
  integer :: ldw, i

  call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
  call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')

  ldw = 2*(ld**ldim)

  call get_int_ptr (i, jgl, jgt, mx,md)

  nintp = nintp + 1
  etime = dnekclock() 
  intp_flop = intp_flop + 2*(mx*mx*mx*md + mx*mx*md*md + mx*md*md*md)
  intp_mop = intp_mop + mx*mx*mx + md*md*md
  if (idir == 0) then
      call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
  endif
  tintp = tintp + (dnekclock() - etime)

  return
end subroutine intp_rstd
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
!> \brief Get pointer to jgl() for interpolation pair (mx,md)
!!
!! The interpolation matrices jgl, jgt are being memoized.
!!  pjgl is a map from (mx,md) pair to (ip) index of jgl, jgt
subroutine get_int_ptr (ip, jgl, jgt, mx,md) ! GLL-->GL pointer
  use kinds, only : DP
  use size_m
  implicit none

  integer, parameter :: ldg=lxd**3, lwkd=4*lxd*lxd

  integer, intent(out) :: ip
  real(DP), intent(inout) :: jgl(ldg), jgt(ldg)
  integer, intent(in) :: mx, md

  real(DP) :: wkd(lwkd)

  integer, parameter :: ld=2*lxd
  integer, save :: pjgl(0:ld*ld) = 0

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
subroutine get_dgl_ptr (ip, dg, dgt, wkd, mx,md)
  use kinds, only : DP
  use size_m, only : lxd 
  implicit none

  integer, parameter :: ld=2*lxd
  integer, parameter :: ldg=lxd**3, lwkd=4*lxd*lxd

  integer, intent(out) :: ip
  real(DP), intent(out) :: dg(ldg), dgt(ldg), wkd(lwkd)
  integer, intent(in) :: mx, md

  integer, save ::  pdg   (0:ld*ld)

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
  
      call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)
  endif

  return
end subroutine get_dgl_ptr
!-----------------------------------------------------------------------
subroutine grad_rst(ur,us,ut,u,md,if3d) ! Gauss-->Gauss grad
  use kinds, only : DP
  use ctimer, only : ngrst, tgrst, grst_flop, grst_mop, dnekclock
  use size_m
  implicit none

  real(DP) ::    ur(1),us(1),ut(1),u(1)
  integer, intent(in) :: md
  logical :: if3d

  integer, parameter :: ldg=lxd**3, lwkd=4*lxd*lxd
  real(DP), save :: dg(ldg), dgt(ldg)
  real(DP) ::  wkd(lwkd)
  integer :: m0, ip
  real(DP) :: etime

  ngrst = ngrst + 1
  etime = dnekclock()

  m0 = md-1
  call get_dgl_ptr (ip, dg, dgt, wkd, md,md)
  if (if3d) then
      grst_flop = grst_flop + 3*((2*md-1)*md**3)
      grst_mop  = grst_mop  + 4*md**3
      call local_grad3(ur,us,ut,u,m0,1,dg(ip),dgt(ip))
  else
!max        call local_grad2(ur,us   ,u,m0,1,dg(ip),dgt(ip))
  endif
  tgrst = tgrst + (dnekclock() - etime)
  return
end subroutine grad_rst
!-----------------------------------------------------------------------
!> \brief Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
subroutine convect_new(bdu,u,ifuf,cx,cy,cz,ifcf)
  use kinds, only : DP
  use size_m, only : nelv
  use size_m, only : lxd, lyd, lzd
  use size_m, only : nx1, ny1, nz1, nxd, nyd, nzd
  use input, only : if3d
  implicit none

  real(DP) :: bdu(*),u(*),cx(*),cy(*),cz(*)
  logical :: ifuf,ifcf            ! u and/or c already on fine mesh?

  integer, parameter :: ltd=lxd*lyd*lzd
  real(DP) :: ur(ltd), us(ltd), ut(ltd), tr(ltd,3), uf(ltd)

  integer :: e, iu, ic, ib, i
  integer :: nxyz1, nxyzd, nxyzu, nxyzc

!max  call set_dealias_rx()

  nxyz1 = nx1*ny1*nz1
  nxyzd = nxd*nyd*nzd

  nxyzu = nxyz1
  if (ifuf) nxyzu = nxyzd

  nxyzc = nxyz1
  if (ifcf) nxyzc = nxyzd

  iu = 1    ! pointer to scalar field u
  ic = 1    ! pointer to vector field C
  ib = 1    ! pointer to scalar field Bdu


  do e=1,nelv

      if (ifcf) then
          call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
          call copy(tr(1,2),cy(ic),nxyzd)
          if (if3d) call copy(tr(1,3),cz(ic),nxyzd)

      else  ! map coarse velocity to fine mesh (C-->F)
        write(*,*) "Oops: ifcf"
#if 0
          call intp_rstd(fx,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
          call intp_rstd(fy,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
          if (if3d) call intp_rstd(fz,cz(ic),nx1,nxd,if3d,0) ! 0 --> forward

          if (if3d) then  ! Convert convector F to r-s-t coordinates

              do i=1,nxyzd
                  tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
                  tr(i,2)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
                  tr(i,3)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
              enddo

          else

              do i=1,nxyzd
                  tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
                  tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
              enddo

          endif
#endif
      endif

      if (ifuf) then
          call grad_rst(ur,us,ut,u(iu),nxd,if3d)
      else
          call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward
          call grad_rst(ur,us,ut,uf,nxd,if3d)
      endif

      if (if3d) then
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
              uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
          enddo
      else
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
              uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
          enddo
      endif
      call intp_rstd(bdu(ib),uf,nx1,nxd,if3d,1) ! Project back to coarse

      ic = ic + nxyzc
      iu = iu + nxyzu
      ib = ib + nxyz1

  enddo

  return
  end subroutine convect_new
!-----------------------------------------------------------------------
!> \brief Put vxd,vyd,vzd into rst form on fine mesh
!! For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
subroutine set_convect_new(cr,cs,ct,ux,uy,uz)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lxd, lyd, lzd, ldim
  use size_m, only : nx1, ny1, nz1, nxd, nyd, nzd, nelv
  use geom, only : rx
  use input, only : if3d
  use mesh, only : if_ortho
  use ctimer, only : tscn, nscn, dnekclock
  implicit none

  integer, parameter :: lxy=lx1*ly1*lz1, ltd=lxd*lyd*lzd
  integer, parameter :: ld=2*lxd
  integer, parameter :: ldw=2*(ld**ldim)

  real(DP) :: cr(ltd,*),cs(ltd,*),ct(ltd,*)
  real(DP) :: ux(lxy,*),uy(lxy,*),uz(lxy,*)

  real(DP) :: fx(ltd), fy(ltd), fz(ltd)!, ur, us, ut, tr, uf
  real(DP) :: w(ld**ldim,2)
  real(DP) :: etime
  integer, parameter :: ldg=lxd**3
  !real(DP), save :: jgl(ldg), jgt(ldg)


  integer :: e, nxyz1, nxyzd, ic, i, j
  etime = dnekclock()
  nscn = nscn + 1
  call set_dealias_rx()

  nxyz1 = nx1*ny1*nz1
  nxyzd = nxd*nyd*nzd

  ic = 1    ! pointer to vector field C
  !call get_int_ptr (j, jgl, jgt, nx1,nxd) 

  do e=1,nelv

  !        Map coarse velocity to fine mesh (C-->F)

      !call specmpn(fx,nxd,ux(1,e),nx1,jgl(i),jgt(i),if3d,w,ldw)
      !call specmpn(fy,nxd,uy(1,e),nx1,jgl(i),jgt(i),if3d,w,ldw)
      !call specmpn(fz,nxd,uz(1,e),nx1,jgl(i),jgt(i),if3d,w,ldw)
      etime = etime - dnekclock()
      call intp_rstd(fx,ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
      call intp_rstd(fy,uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
      if (if3d) call intp_rstd(fz,uz(1,e),nx1,nxd,if3d,0) ! 0 --> forward
      etime = etime + dnekclock()

  !        Convert convector F to r-s-t coordinates

      if (if3d) then

          if (if_ortho) then
            do i=1,nxyzd
                cr(i,e)=rx(i,1,e)*fx(i)
                cs(i,e)=rx(i,2,e)*fy(i)
                ct(i,e)=rx(i,3,e)*fz(i)
            enddo
          else
            do i=1,nxyzd
                cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
                cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
                ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
            enddo
          endif
      else
#if 0
          do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
              cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
          enddo
#endif
      endif
  enddo
  tscn = tscn + (dnekclock() - etime)

  return
end subroutine set_convect_new
!-----------------------------------------------------------------------
