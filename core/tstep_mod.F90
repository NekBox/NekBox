!> cleaned
module tstep
  use kinds, only : DP
  use size_m
  implicit none

  real(DP) :: TIME,TIMEF,FINTIM,TIMEIO &
    ,DT,DTINIT,DTINVM,COURNO,CTARG &
    ,AB(10),BD(10),ABMSH(10) &
    ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1) &
    ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX &
    ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE &
    ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN &
    ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT) &
    ,TNRML8(LDIMT),TMEAN(LDIMT)


  real(DP) :: re_cell = 0._dp !>!< Cell Reynolds number

  integer :: IFIELD,IMESH,NSTEPS,IOSTEP,LASTEP,IOCOMM &
    ,INSTEP &
    ,NAB,NBDINP,NTAUBD &
    ,NMXH,NMXP,NMXE,NMXNL,NINTER &
    ,NELFLD(0:LDIMT1) &
    ,nconv,nconv_max, ntdump

  integer :: istep = 1, NBD=0
  real(DP) :: DTLAG(10) = 0._dp

  real(DP) :: PI, BETAG, GTHETA
  LOGICAL :: IFPRNT,if_full_pres

  real(DP) :: lyap(3,lpert)
  real(DP) :: mixing_alpha = 1.0_dp, mixing_beta = 1.0_dp

end module tstep
