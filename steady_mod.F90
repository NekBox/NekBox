module steady
  use kinds, only : DP
  use size_m
  implicit none

  real(DP) :: TAUSS(LDIMT1), TXNEXT(LDIMT1) 
  integer :: nsskip
  LOGICAL :: IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1) &
    ,                 IFEXVT, IFEXTR(LDIMT1)
  real(DP) :: DVNNH1, DVNNSM, DVNNL2, DVNNL8 &
    , DVDFH1, DVDFSM, DVDFL2, DVDFL8 &
    , DVPRH1, DVPRSM, DVPRL2, DVPRL8

end module steady
