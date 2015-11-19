!-----------------------------------------------------------------------
!    Stability limits:
!    AB3:    .7236                     w/safety (1.2):   .603
!    RK3:    1.73   (sqrt 3)           w/safety (1.2):   1.44
!    RK4:    2.828                     w/safety (1.2):   2.36
!    SEM Safety factor:  1.52 for N=3
!                     <  1.20 for N=16
!                     ~  1.16 for N=256
!-----------------------------------------------------------------------
#define INLINE_INTP
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
!> \brief Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
subroutine convect_new(bdu,u,ifuf,cx,cy,cz,ifcf)
  use kinds, only : DP
  use size_m, only : nelv
  use size_m, only : lx1, lxd, lyd, lzd, ldim
  use size_m, only : nx1, ny1, nz1, nxd, nyd, nzd
  use ctimer, only : tscn, dnekclock
  use interp, only : jgl, jgt, dgl, dgt
  use interp, only : get_int_ptr, get_dgl_ptr
  implicit none

  real(DP) :: bdu(*),u(*),cx(*),cy(*),cz(*)
  logical :: ifuf,ifcf            ! u and/or c already on fine mesh?

  integer, parameter :: ltd=lxd*lyd*lzd
  real(DP) :: ur(ltd), us(ltd), ut(ltd), tr(ltd,3), uf(ltd)

  integer :: e, iu, ic, ib, i, iptr, iptr2
  integer :: nxyz1, nxyzd, nxyzu, nxyzc
  real(DP) :: etime
  real(DP) :: w((2*lxd)**ldim,2)
  real(DP) :: w1(lxd*lxd*lx1), w2(lxd*lx1*lx1)
  integer, parameter :: ldw = 2*(2*lxd)**ldim

  etime = dnekclock()
!max  call set_dealias_rx()

  if (ifuf .or. .not. ifcf) then
    write(*,*) "Oops: convect_new args unsupported"
    return
  endif

  nxyz1 = nx1*ny1*nz1
  nxyzd = nxd*nyd*nzd

  nxyzu = nxyz1

  nxyzc = nxyz1
  if (ifcf) nxyzc = nxyzd

  iu = 1    ! pointer to scalar field u
  ic = 1    ! pointer to vector field C
  ib = 1    ! pointer to scalar field Bdu

#ifdef INLINE_INTP
  call get_int_ptr (iptr,  nx1,nxd)
  call get_dgl_ptr (iptr2, nxd,nxd)
#endif

  do e=1,nelv
    call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
    call copy(tr(1,2),cy(ic),nxyzd)
    call copy(tr(1,3),cz(ic),nxyzd)

    call tensor_product_multiply(u(iu), nx1, uf, nxd, jgl(iptr), jgt(iptr), jgt(iptr), w2, w1)
    call local_grad3(ur,us,ut,uf,nxd-1,1,dgl(iptr2),dgt(iptr2))
    do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
        uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
    enddo
    call tensor_product_multiply(uf, nxd, bdu(ib), nx1, jgt(iptr), jgl(iptr), jgl(iptr), w1, w2)

    ic = ic + nxyzc
    iu = iu + nxyzu
    ib = ib + nxyz1

  enddo
  tscn = tscn + (dnekclock() - etime)

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
  use mesh, only : if_ortho
  use ctimer, only : tscn, nscn, dnekclock
  use interp, only : jgl, jgt
  use interp, only : get_int_ptr

  implicit none

  integer, parameter :: lxy=lx1*ly1*lz1, ltd=lxd*lyd*lzd

  real(DP) :: cr(ltd,*),cs(ltd,*),ct(ltd,*)
  real(DP) :: ux(lxy,*),uy(lxy,*),uz(lxy,*)

  real(DP) :: fx(ltd), fy(ltd), fz(ltd)!, ur, us, ut, tr, uf
  real(DP) :: w((2*lxd)**ldim,2)
  real(DP) :: w1(lxd*lxd*lx1), w2(lxd*lx1*lx1)
  integer, parameter :: ldw = 2*(2*lxd)**ldim

  real(DP) :: etime

  integer :: e, nxyz1, nxyzd, ic, i, j, iptr
  nscn = nscn + 1
  etime = dnekclock()
  call set_dealias_rx()

  nxyz1 = nx1*ny1*nz1
  nxyzd = nxd*nyd*nzd

  ic = 1    ! pointer to vector field C
#ifdef INLINE_INTP
  call get_int_ptr (iptr, nx1, nxd) 
#endif

  do e=1,nelv

  !      Map coarse velocity to fine mesh (C-->F)
#ifdef INLINE_INTP
    call tensor_product_multiply(ux(1,e), nx1, fx, nxd, jgl(iptr), jgt(iptr), jgt(iptr), w2, w1)
    call tensor_product_multiply(uy(1,e), nx1, fy, nxd, jgl(iptr), jgt(iptr), jgt(iptr), w2, w1)
    call tensor_product_multiply(uz(1,e), nx1, fz, nxd, jgl(iptr), jgt(iptr), jgt(iptr), w2, w1)
#else
    call intp_rstd(fx,ux(1,e),nx1,nxd,.true.,0) ! 0 --> forward
    call intp_rstd(fy,uy(1,e),nx1,nxd,.true.,0) ! 0 --> forward
    call intp_rstd(fz,uz(1,e),nx1,nxd,.true.,1) ! 0 --> forward
#endif
!      etime = etime + dnekclock()

  !        Convert convector F to r-s-t coordinates


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
  enddo
  tscn = tscn + (dnekclock() - etime)

  return
end subroutine set_convect_new
!-----------------------------------------------------------------------
