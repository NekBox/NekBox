module ctimer
  use kinds, only : DP, i8
  implicit none

  REAL(DP) ::          tmxmf,tmxms,tdsum,tcopy,tinvc,tinv3
  REAL(DP) ::          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd &
  ,tpres,tgop ,tgop1,tdott,tbsol,tbso2 &
  ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn &
  ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee &
  ,tvdss,tspro,tgop_sync,tsyc &
  ,twal

  real(DP), save :: max_mops  = 0._dp
  real(DP), save :: max_flops = 0._dp
  real(DP), save :: allreduce_latency = 0._dp

  real(DP), save :: tproj  = 0._dp
  real(DP), save :: thconj = 0._dp
  real(DP), save :: taxhm  = 0._dp
  real(DP), save :: tcggo  = 0._dp
  real(DP), save :: tsetfast  = 0._dp
  real(DP), save :: tintp  = 0._dp
  real(DP), save :: tgrst  = 0._dp
  real(DP), save :: tdpc  = 0._dp
  real(DP), save :: tgmres  = 0._dp
  real(DP), save :: th1mg  = 0._dp
  real(DP), save :: tcps  = 0._dp
  real(DP), save :: tcrespsp  = 0._dp
  real(DP), save :: tcresvsp  = 0._dp
  real(DP), save :: theat2  = 0._dp
  real(DP), save :: tp4misc  = 0._dp
  real(DP), save :: tmakef  = 0._dp
  real(DP), save :: tmakeq  = 0._dp
  real(DP), save :: tprep  = 0._dp
  real(DP), save :: tscn  = 0._dp
  real(DP), save :: tmg_mask  = 0._dp
  real(DP), save :: tschw  = 0._dp
  real(DP), save :: tnmsc  = 0._dp
  real(DP), save :: tnmvc  = 0._dp
  real(DP), save :: tadvc  = 0._dp
  real(DP), save :: thmhz  = 0._dp
  real(DP), save :: tfoo   = 0._dp
  real(DP), save :: thslv  = 0._dp
  real(DP), save :: tsetn  = 0._dp
  real(DP), save :: tusrc  = 0._dp
  real(DP), save :: tqflt  = 0._dp

  integer, save :: nproj  = 0
  integer, save :: nhconj = 0
  integer, save :: naxhm  = 0
  integer, save :: ncggo  = 0
  integer, save :: nsetfast = 0
  integer, save :: nintp = 0
  integer, save :: ngrst = 0
  integer, save :: ndpc = 0
  integer, save :: ngmres = 0
  integer, save :: nh1mg = 0
  integer, save :: ncps = 0
  integer, save :: ncrespsp = 0
  integer, save :: ncresvsp = 0
  integer, save :: nheat2 = 0
  integer, save :: np4misc = 0
  integer, save :: nmakef = 0
  integer, save :: nmakeq = 0
  integer, save :: nprep = 0
  integer, save :: nscn = 0
  integer, save :: nmg_mask = 0
  integer, save :: nschw = 0
  integer, save :: nnmsc = 0
  integer, save :: nnmvc = 0
  integer, save :: nadvc = 0
  integer, save :: nhmhz = 0
  integer, save :: nfoo  = 0
  integer, save :: nhslv  = 0
  integer, save :: nsetn  = 0
  integer, save :: nusrc  = 0
  integer, save :: nqflt  = 0
  integer, save :: ndsum = 0

  integer :: nmxmf,nmxms,ncopy,ninvc,ninv3
  integer :: nsolv,ngsum,ndsnd,ndadd,ncdtp,nmltd &
  ,npres,ngop ,ngop1,ndott,nbsol,nbso2 &
  ,nsett,nslvb,nusbc,nddsl,ncrsl,ndsmx,ndsmn &
  ,ngsmn,ngsmx,neslv,nbbbb,ncccc,ndddd,neeee &
  ,nvdss,nspro,ngop_sync,nsyc,nwal

  REAL(DP) ::          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
  REAL(DP) ::          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd &
  ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2 &
  ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn &
  ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee &
  ,pvdss,pspro,pgop_sync,psyc,pwal

  REAL(DP) :: etime1,etime2,etime0,gtime1,tscrtch

  real(DP) ::          etimes,ttotal,tttstp,etims0,ttime

  integer, save :: icalld = 0

  logical ::         ifsync

  real(DP), save :: time_flop = 0._dp
  integer(i8), save :: total_flop, total_mop
  integer(i8), save :: axhelm_flop = 0, axhelm_mop = 0
  integer(i8), save :: h1mg_flop = 0, h1mg_mop = 0
  integer(i8), save :: schw_flop = 0, schw_mop = 0
  integer(i8), save :: hconj_flop = 0, hconj_mop = 0
  integer(i8), save :: proj_flop = 0, proj_mop = 0
  integer(i8), save :: cggo_flop = 0, cggo_mop = 0
  integer(i8), save :: gmres_flop = 0, gmres_mop = 0
  integer(i8), save :: conv_flop = 0, conv_mop = 0
  integer(i8), save :: grst_flop = 0, grst_mop = 0
  integer(i8), save :: othr_flop = 0, othr_mop = 0

contains

!-----------------------------------------------------------------------
real(DP) function dnekclock()
  use mpif, only : mpi_wtime
  implicit none

  dnekclock=mpi_wtime()

  return
END function

!-----------------------------------------------------------------------
real(DP) function dnekclock_sync()
  use mpif, only : mpi_wtime
  implicit none

  call nekgsync()
  dnekclock_sync=mpi_wtime()

  return
END function

!-----------------------------------------------------------------------
subroutine sum_flops()
  implicit none
  total_flop = axhelm_flop + proj_flop + hconj_flop + cggo_flop + gmres_flop &
             + conv_flop + grst_flop + h1mg_flop + schw_flop
  total_mop = axhelm_mop + proj_mop + hconj_mop + cggo_mop + gmres_mop &
            + conv_mop + grst_mop + h1mg_mop + schw_mop
  time_flop = taxhm + tproj + thconj + tcggo + tgmres + tintp + tgrst + th1mg + tschw
end subroutine sum_flops

!-----------------------------------------------------------------------
!> \brief compute the max bandwidth as in STREAM
!!
!! This is used when computing efficiencies in drive2
subroutine benchmark()
  use parallel, only : nid
  call benchmark_mxm()
  call benchmark_stream()
  call benchmark_mxm()
  call benchmark_stream()
  call benchmark_allreduce()
  if (nid == 0) write(*,'(2(A,F6.2),A,E11.3)') "GFLOP: ", max_flops/10.**9, &
                                        " GiB/s: ", max_mops * 8 / (2_8 **30), &
                                        " Allreduce: ", allreduce_latency
end subroutine

subroutine benchmark_allreduce()
  use kinds, only : DP, i8
  use parallel, only : nid

  integer(i8) :: i, n
  real(DP) :: a(1) 
  real(DP) :: work(10), etime_builtin, etime_gs, etime

  n = 100

  a = 1.
  do i = 1, n
    call nekgsync()
    call gop(a,work,'+  ',1)
    !call gop_gs(a,work,'+  ',1)
  enddo


  call nekgsync()
  etime = dnekclock()
  do i = 1, n
    call nekgsync()
    call gop(a,work,'+  ',1)
  enddo
  call nekgsync()
  etime_builtin = dnekclock() - etime 
  if (nid == 0) write(*,'(A,2E11.3)') "Builtin vs gs: ", etime_builtin / n

  call nekgsync()
  etime = dnekclock()
  do i = 1, n
    call nekgsync()
    call gop(a,work,'+  ',1)
  enddo
  call nekgsync()
  etime_gs = dnekclock() - etime 

  if (nid == 0) write(*,'(A,2E11.3)') "Builtin vs gs: ", etime_builtin / n, etime_gs / n

  allreduce_latency = etime_builtin / n

end subroutine benchmark_allreduce

subroutine benchmark_mxm()
  use kinds, only : DP
  use parallel, only : nid
  use size_m, only : lx1, ly1, lz1, lelt
#ifdef XSMM
  use iso_c_binding
  use LIBXSMM, only : LIBXSMM_DMMFUNCTION
  use LIBXSMM, only : libxsmm_dispatch, libxsmm_call
#endif

  real(DP), allocatable :: a(:,:,:), b(:,:,:), c(:,:)
  !DIR$ ATTRIBUTES ALIGN:64 :: a, b, c


  real(DP) :: etime
  integer(8) :: i, n, k, flops
  integer :: iz

#ifdef XSMM
  TYPE(LIBXSMM_DMMFUNCTION), save :: xmm1, xmm2, xmm3

  call libxsmm_dispatch(xmm1, lx1,     lx1*lx1, lx1, alpha=1.0_dp, beta=0.0_dp)
  call libxsmm_dispatch(xmm2, lx1,     lx1,     lx1, alpha=1.0_dp, beta=0.0_dp)
  call libxsmm_dispatch(xmm3, lx1*lx1, lx1,     lx1, alpha=1.0_dp, beta=0.0_dp)
#endif


  allocate(a(lx1, lx1, lx1), b(lx1, lx1, lx1), c(lx1, lx1))
  n = (2_8**40)/(lx1**4) / 4

  etime = dnekclock()
  do i = 1, n
#ifdef XSMM
  CALL libxsmm_call(xmm1, C_LOC(c), C_LOC(a(1,1,1)), C_LOC(b))
  do iz=1,lx1
      CALL libxsmm_call(xmm2, C_LOC(a(1,1,iz)), C_LOC(c), C_LOC(b(1,1,iz)))
  enddo
  CALL libxsmm_call(xmm3, C_LOC(a(1,1,1)), C_LOC(c), C_LOC(b))    
#else
  call mxm   (c,lx1,a(1,1,1),lx1,b,lx1*lx1)
  do iz=1,lx1
      call mxm   (a(1,1,iz),lx1,c,lx1,b(1,1,iz),lx1)
  END DO
  call mxm   (a(1,1,1),lx1*lx1,c,lx1,b,lx1)
#endif
  enddo
  etime = dnekclock() - etime
  deallocate(a,b,c)
  flops = lx1*lx1*lx1*n * 3 * ( 2*lx1 - 1)
  max_flops = flops / etime

  if (nid == 0) write(*,*) "Got max flop rate of ", max_flops, etime

end subroutine benchmark_mxm

!-----------------------------------------------------------------------
!> \brief compute the max bandwidth as in STREAM
!!
!! This is used when computing efficiencies in drive2
subroutine benchmark_stream()
  use kinds, only : DP
  use parallel, only : nid
  use size_m, only : lx1, ly1, lz1, lelt
  implicit none

  real(DP), allocatable :: a(:), b(:), c(:)
  !DIR$ ATTRIBUTES ALIGN:64 :: a, b, c

  integer(8) :: i, n, k
  real(DP) :: foo, etime, etime_total(4)


  ! just replicate STREAM
  n = 2**24 / 8 
  k = 8
  allocate(a(n), b(n), c(n)) 
  a = 2._dp; b = 0.5_dp; c = 0._dp

  etime = dnekclock_sync() 
  a = 0.5*a
  etime = dnekclock_sync() - etime
  foo = .5 * a(1)

  etime_total = 0 
  do i = 1, k
    etime = dnekclock_sync()
    a(1) = a(1) + etime
    !DEC$ vector aligned nontemporal
    c = a
    etime = dnekclock_sync() - etime
    c(n) = c(n) + etime
    etime_total(1) = etime_total(1) + etime

    etime = dnekclock_sync()
    c(1) = c(1) + etime
    !DEC$ vector aligned nontemporal
    b = foo * c
    etime = dnekclock_sync() - etime 
    b(n) = b(n) + etime
    etime_total(2) = etime_total(2) + etime

    etime = dnekclock_sync()
    a(1) = a(1) + etime
    !DEC$ vector aligned nontemporal
    c = a + b
    etime = dnekclock_sync() - etime
    c(n) = c(n) + etime
    etime_total(3) = etime_total(3) + etime

    etime = dnekclock_sync()
    b(1) = b(1) + etime
    !DEC$ vector aligned nontemporal
    a = b + foo*c
    etime = dnekclock_sync() - etime
    a(n) = a(n) + etime
    etime_total(4) = etime_total(4) + etime
 
  enddo
  deallocate(a,b,c)
  max_mops = real((n*k*10))/sum(etime_total)
  if (nid == 0) write(*,'(A,5E12.5)') "STREAM: ", real(n*k*8) * (/2, 2, 3, 3/) / etime_total, max_mops*8

end subroutine benchmark_stream
end module ctimer
