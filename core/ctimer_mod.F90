module ctimer
  use kinds, only : DP

  REAL(DP) ::          tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
  REAL(DP) ::          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep &
  ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2 &
  ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn &
  ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee &
  ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc &
  ,twal


  integer :: nmxmf,nmxms,ndsum,naxhm,ncopy,ninvc,ninv3
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

end module ctimer
