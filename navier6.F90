!===============================================================================
!     pff@cfm.brown.edu   3/19/96


!     This  is a suite of routines for solving overlapping subdomain
!     problems with finite elements on distributed memory machines.

!     The overall hierarchy of data structures is as follows:

!         System        -  index set denoted by       _glob

!            Processor  -  index set denoted by       _pglob

!              .Domain  -  index set denoted by       _dloc  (1,2,3,...,n_dloc)

!              .Sp.Elem -  index set denoted by       _sloc  (1,2,3,...,n_sloc)


!     A critical component of the parallel DD solver is the gather-scatter
!     operation.   As there may be more than one domain on a processor, and
!     communication must be minimized, it is critical that communication be
!     processor oriented, not domain oriented.  Hence domains will access data
!     via the dloc_to_pglob interface, and the pglob indexed variables will
!     be accessed via a gather-scatter which interacts via the pglob_glob
!     interface.   Noticed that, in a uni-processor application, the pglob_glob
!     map will be simply the identity.

!===============================================================================

!> \brief Set up arrays for overlapping Schwartz algorithm *for pressure solver*
subroutine set_overlap
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, lbelv, lelt, nid
  use domain, only : ltotd
  use input, only : param, ifaxis, ifmgrid, ifmhd, ifsplit, ifldmhd
  use tstep, only : ifield
  implicit none

  REAL*8 :: dnekclock,t0

  integer, parameter :: n_tri = 7*ltotd
  integer :: elem

  real(DP) :: x(2*ltotd), y(2*ltotd), z(2*ltotd)

  integer, parameter :: lia = ltotd - 2 - 2*lelt

  integer, parameter :: lxx=lx1*lx1, levb=lelv+lbelv
  integer :: e, npass, ipass

  if (lx1 == 2) param(43)=1.
  if (lx1 == 2 .AND. nid == 0) write(6,*) 'No mgrid for lx1=2!'

  if (ifaxis) ifmgrid = .FALSE. 
  if (param(43) /= 0) ifmgrid = .FALSE. 

  npass = 1
  if (ifmhd) npass = 2
  do ipass=1,npass
      ifield = 1

      if (ifsplit .AND. ifmgrid) then
          if (ipass > 1) ifield = ifldmhd

          call swap_lengths
          call gen_fast_spacing(x,y,z)
           
          call hsmg_setup
          call h1mg_setup

      elseif ( .NOT. ifsplit) then ! Pn-Pn-2
#if 0
          if (ipass > 1) ifield = ifldmhd

          if (param(44) == 1) then !  Set up local overlapping solves
              call set_fem_data_l2(nel_proc,ndom,n_o,x,y,z,tri)
          else
              call swap_lengths
          endif
           
          e = 1
          if (ifield > 1) e = nelv+1

          call gen_fast_spacing(x,y,z)
          call gen_fast(df(1,e),sr(1,e),ss(1,e),st(1,e),x,y,z)

          call init_weight_op
          if (param(43) == 0) call hsmg_setup
#endif
      endif

      call set_up_h1_crs

  enddo
   
  return
end subroutine set_overlap
!-----------------------------------------------------------------------
