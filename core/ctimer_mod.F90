module ctimer
  use kinds, only : DP, i8
  implicit none

  REAL(DP) ::          tmxmf,tmxms,tdsum,tcopy,tinvc,tinv3
  REAL(DP) ::          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep &
  ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2 &
  ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn &
  ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee &
  ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc &
  ,twal
  real(DP), save :: tproj  = 0._dp
  real(DP), save :: thconj = 0._dp
  real(DP), save :: taxhm  = 0._dp
  real(DP), save :: tcggo  = 0._dp
  real(DP), save :: tsetfast  = 0._dp
  real(DP), save :: tintp  = 0._dp
  real(DP), save :: tdpc  = 0._dp
  real(DP), save :: tgmres  = 0._dp


  integer, save :: nproj  = 0
  integer, save :: nhconj = 0
  integer, save :: naxhm  = 0
  integer, save :: ncggo  = 0
  integer, save :: nsetfast = 0
  integer, save :: nintp = 0
  integer, save :: ndpc = 0
  integer, save :: ngmres = 0
  integer :: nmxmf,nmxms,ndsum,ncopy,ninvc,ninv3
  integer :: nsolv,ngsum,ndsnd,ndadd,ncdtp,nmltd,nprep &
  ,npres,nhmhz,ngop ,ngop1,ndott,nbsol,nbso2 &
  ,nsett,nslvb,nusbc,nddsl,ncrsl,ndsmx,ndsmn &
  ,ngsmn,ngsmx,neslv,nbbbb,ncccc,ndddd,neeee &
  ,nvdss,nadvc,nspro,ngop_sync,nsyc,nwal

  REAL(DP) ::          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
  REAL(DP) ::          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep &
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
  integer(i8), save :: hconj_flop = 0, hconj_mop = 0
  integer(i8), save :: proj_flop = 0, proj_mop = 0
  integer(i8), save :: cggo_flop = 0, cggo_mop = 0
  integer(i8), save :: intp_flop = 0, intp_mop = 0

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
  total_flop = axhelm_flop + proj_flop + hconj_flop + cggo_flop + intp_flop
  total_mop = axhelm_mop + proj_mop + hconj_mop + cggo_mop + intp_mop
  time_flop = taxhm + tproj + thconj + tcggo + tintp
end subroutine sum_flops


end module ctimer
