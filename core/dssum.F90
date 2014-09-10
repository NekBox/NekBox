subroutine setupds(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
  use kinds, only : DP, i8
  use size_m, only : nid
  use parallel, only : mp=>np, nekcomm
  implicit none

  integer :: gs_handle
  integer :: vertex(1)
  integer(i8) :: glo_num(1),ngv

  real(DP) :: t0, t1
  integer :: nx, nel, ntot, ny, nz, melg
  real(DP), external :: dnekclock

  t0 = dnekclock()

!   Global-to-local mapping for gs
  call set_vert(glo_num,ngv,nx,nel,vertex, .FALSE. )

!   Initialize gather-scatter code
  ntot      = nx*ny*nz*nel
  call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

!   call gs_chkr(glo_num)

  t1 = dnekclock() - t0
  if (nid == 0) then
      write(6,1) t1,gs_handle,nx,ngv,melg
      1 format('   setupds time',1pe11.4,' seconds ',2i3,2i12)
  endif

  return
end subroutine setupds
!-----------------------------------------------------------------------
subroutine dssum(u)
  use kinds, only : DP
  use ctimer, only : ifsync, icalld, tdsmx, tdsmn, etime1, dnekclock
  use ctimer, only : tdsum, ndsum
  use input, only : ifldmhd
  use parallel, only : gsh_fld
  use tstep, only : ifield
  implicit none

  real(DP) :: u(1)
  
  integer :: ifldt
  real(DP) :: timee

  ifldt = ifield
!   if (ifldt.eq.0)       ifldt = 1
  if (ifldt == ifldmhd) ifldt = 1
!   write(6,*) ifldt,ifield,gsh_fld(ifldt),imesh,' ifldt'

  if (ifsync) call nekgsync()

#ifndef NOTIMER
  if (icalld == 0) then
      tdsmx=0.
      tdsmn=0.
  endif
  icalld=icalld+1
  etime1=dnekclock()
#endif


!               T         ~  ~T  T
!   Implement QQ   :=   J Q  Q  J


!                T
!   This is the J  part,  translating child data
!    call apply_Jt(u,nx,ny,nz,nel)



!               ~ ~T
!   This is the Q Q  part

  call gs_op(gsh_fld(ifldt),u,1,1,0)  ! 1 ==> +



!   This is the J  part,  interpolating parent solution onto child
!    call apply_J(u,nx,ny,nz,nel)


#ifndef NOTIMER
  timee=(dnekclock()-etime1)
  tdsum=tdsum+timee
  ndsum=icalld
  tdsmx=max(timee,tdsmx)
  tdsmn=min(timee,tdsmn)
#endif

  return
end subroutine dssum

!-----------------------------------------------------------------------
subroutine dsop(u,op)
  use kinds, only : DP
  use ctimer, only : ifsync
  use input, only : ifldmhd
  use parallel, only : gsh_fld
  use tstep, only : ifield
  implicit none

  real(DP) :: u(1)
  character(3) :: op
  integer :: ifldt

!   o gs recognized operations:

!           o "+" ==> addition.
!           o "*" ==> multiplication.
!           o "M" ==> maximum.
!           o "m" ==> minimum.
!           o "A" ==> (fabs(x)>fabs(y)) ? (x) : (y), ident=0.0.
!           o "a" ==> (fabs(x)<fabs(y)) ? (x) : (y), ident=MAX_DBL
!           o "e" ==> ((x)==0.0) ? (y) : (x),        ident=0.0.

!           o note: a binary function pointer flavor exists.


!   o gs level:

!           o level=0 ==> pure tree
!           o level>=num_nodes-1 ==> pure pairwise
!           o level = 1,...num_nodes-2 ==> mix tree/pairwise.

  ifldt = ifield
!   if (ifldt.eq.0)       ifldt = 1
  if (ifldt == ifldmhd) ifldt = 1

!   if (nid.eq.0)
!  $   write(6,*) istep,' dsop: ',op,ifield,ifldt,gsh_fld(ifldt)

  if(ifsync) call nekgsync()

  if (op == '+  ') call gs_op(gsh_fld(ifldt),u,1,1,0)
  if (op == 'sum') call gs_op(gsh_fld(ifldt),u,1,1,0)
  if (op == 'SUM') call gs_op(gsh_fld(ifldt),u,1,1,0)

  if (op == '*  ') call gs_op(gsh_fld(ifldt),u,1,2,0)
  if (op == 'mul') call gs_op(gsh_fld(ifldt),u,1,2,0)
  if (op == 'MUL') call gs_op(gsh_fld(ifldt),u,1,2,0)

  if (op == 'm  ') call gs_op(gsh_fld(ifldt),u,1,3,0)
  if (op == 'min') call gs_op(gsh_fld(ifldt),u,1,3,0)
  if (op == 'mna') call gs_op(gsh_fld(ifldt),u,1,3,0)
  if (op == 'MIN') call gs_op(gsh_fld(ifldt),u,1,3,0)
  if (op == 'MNA') call gs_op(gsh_fld(ifldt),u,1,3,0)

  if (op == 'M  ') call gs_op(gsh_fld(ifldt),u,1,4,0)
  if (op == 'max') call gs_op(gsh_fld(ifldt),u,1,4,0)
  if (op == 'mxa') call gs_op(gsh_fld(ifldt),u,1,4,0)
  if (op == 'MAX') call gs_op(gsh_fld(ifldt),u,1,4,0)
  if (op == 'MXA') call gs_op(gsh_fld(ifldt),u,1,4,0)

  return
end subroutine dsop

!-----------------------------------------------------------------------
!> \brief Direct stiffness summation of the face data, for field U.
!!
!! Boundary condition data corresponds to component IFIELD of
!! the CBC array.
subroutine vec_dssum(u,v,w)
  use kinds, only : DP
  use size_m, only : ndim
  use ctimer, only : ifsync, icalld, tvdss, tgsum, nvdss, etime1, dnekclock
  use ctimer, only : tdsmx, tdsmn
  use input, only : ifldmhd
  use parallel, only : gsh_fld
  use tstep, only : ifield
  implicit none

  REAL(DP) :: U(1),V(1),W(1)
  integer :: ifldt
  real(DP) :: timee

  if(ifsync) call nekgsync()

#ifndef NOTIMER
  if (icalld == 0) tvdss=0.0d0
  if (icalld == 0) tgsum=0.0d0
  icalld=icalld+1
  nvdss=icalld
  etime1=dnekclock()
#endif

!============================================================================
!     execution phase
!============================================================================

  ifldt = ifield
!   if (ifldt.eq.0)       ifldt = 1
  if (ifldt == ifldmhd) ifldt = 1

  call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)

#ifndef NOTIMER
  timee=(dnekclock()-etime1)
  tvdss=tvdss+timee
  tdsmx=max(timee,tdsmx)
  tdsmn=min(timee,tdsmn)
#endif

  return
end subroutine vec_dssum

!-----------------------------------------------------------------------
!> \brief Direct stiffness summation of the face data, for field U.
!! Boundary condition data corresponds to component IFIELD of
!! the CBC array.
subroutine vec_dsop(u,v,w,op)
  use kinds, only : DP
  use size_m, only : ndim
  use ctimer, only : ifsync
  use input, only : ifldmhd
  use parallel, only : gsh_fld
  use tstep, only : ifield
  implicit none

  real(DP) :: u(1),v(1),w(1)
  character(3) :: op
  integer :: ifldt

!============================================================================
!     execution phase
!============================================================================

  ifldt = ifield
!   if (ifldt.eq.0)       ifldt = 1
  if (ifldt == ifldmhd) ifldt = 1

!   write(6,*) 'opdsop: ',op,ifldt,ifield
  if(ifsync) call nekgsync()

  if (op == '+  ' .OR. op == 'sum' .OR. op == 'SUM') &
  call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)


  if (op == '*  ' .OR. op == 'mul' .OR. op == 'MUL') &
  call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,2,0)


  if (op == 'm  ' .OR. op == 'min' .OR. op == 'mna' &
   .OR. op == 'MIN' .OR. op == 'MNA') &
  call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,3,0)


  if (op == 'M  ' .OR. op == 'max' .OR. op == 'mxa' &
   .OR. op == 'MAX' .OR. op == 'MXA') &
  call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,4,0)


  return
end subroutine vec_dsop

!-----------------------------------------------------------------------
