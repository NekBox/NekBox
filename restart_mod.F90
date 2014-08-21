module restart
  use kinds, only : DP
  use size_m
  implicit none
!     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservative VERSION
  integer, parameter :: lelr=lelg              ! THIS IS THE MEMORY INTENSIVE VERSION

  real(DP) :: max_rst           ! for full restart

  integer ::  nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler &
    , nxo,nyo,nzo,nrg &
    , wdsizr,wdsizo &
    , nfileo,nproc_o,nfldr &
    , er(lelr),nelB,nelBr

  integer, parameter :: iHeaderSize=132

  real(DP) :: timer

    character(3) ::  ihdr
    character(10) :: rdcode
    character(80) :: mfi_fname
    character(1) ::  rdcode1(10)
!max    equivalence (rdcode,rdcode1)

    logical :: &
    ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps(ldimt1),ifgtim &
    ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr &
    ,if_byte_sw &
    ,ifgetz,ifgetw &
    ,ifdiro

    integer ::         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00

    integer ::          nekcomm_io,ifh_mbyte

end module restart
