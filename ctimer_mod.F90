module ctimer

    REAL*8 ::          tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
    REAL*8 ::          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep &
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

    REAL*8 ::          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
    REAL*8 ::          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep &
    ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2 &
    ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn &
    ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee &
    ,pvdss,pspro,pgop_sync,psyc,pwal

    REAL*8 :: etime1,etime2,etime0,gtime1,tscrtch
    REAL*8, external :: dnekclock,dnekclock_sync

    real*8 ::          etimes,ttotal,tttstp,etims0,ttime

    integer :: icalld 
    save    icalld
    data    icalld /0/

    logical ::         ifsync

end module ctimer
