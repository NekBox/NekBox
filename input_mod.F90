!> \brief Input parameters from preprocessors.
!!
!!     Note that in parallel implementations, we distinguish between
!!     distributed data (LELT) and uniformly distributed data.
!!     Input common block structure:
!!     INPUT1:  REAL            INPUT5: REAL      with LELT entries
!!     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
!!     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
!!     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
!> cleaned
module input


  use kinds, only : DP
  implicit none

  real(DP) :: RSTIM, VNEKTON
  real(DP) :: NKTONV
  integer :: IRSTV, IRSTT, IRSTIM, NOBJ, NGEOM
  integer :: nhis, ipscal, npscal, ipsco, ifldmhd

  real(DP), allocatable ::  PARAM(:), CPFLD(:,:), CPGRP(:,:,:), QINTEG(:,:)
  integer, allocatable :: matype(:,:), lochis(:,:), nmember(:)
  logical, allocatable :: IFTMSH(:),IFNONL(:),IFVARP(:),IFPSCO(:),IFPSO(:) 
  logical, allocatable, target :: IFADVC(:)
  
  logical :: IF3D &
    ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT &
    ,IFMGRID &
    ,IFMVBD,IFNATC,IFCHAR &
    ,IFVPS &
    ,IFMODEL,IFKEPS &
    ,IFINTQ,IFCONS &
    ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFFMTIN &
    ,IFBO &
    ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE &
    ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP, &
    IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO, &
    IFXYO_,ifaziv,IFNEKNEK
  logical :: ifprint

  LOGICAL, pointer :: IFNAV

    CHARACTER(1), allocatable :: HCODE(:,:)
    CHARACTER(2), allocatable ::     OCODE(:)
    CHARACTER(10), allocatable ::    DRIVC(:)
    CHARACTER(14) ::    RSTV,RSTT
    CHARACTER(40), allocatable, target ::    TEXTSW(:,:)
    CHARACTER(40), pointer ::    TURBMOD
    CHARACTER(132), allocatable ::    INITC(:)

    CHARACTER(132) ::   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
    CHARACTER(132) ::   SESSION,PATH,RE2FLE,H5MFLE

! proportional to LELT

  real(DP), allocatable :: XC(:,:),YC(:,:),ZC(:,:) &
    ,BC(:,:,:,:) &
    ,CURVE(:,:,:) &
    ,CERROR(:)

  integer, allocatable :: IGROUP(:),OBJECT(:,:,:)

  CHARACTER(1), allocatable :: CCURVE(:,:),CDOF(:,:)
  CHARACTER(3), allocatable :: CBC(:,:,:)
  character(3) :: solver_type

  integer, allocatable :: ieact(:)
  integer :: neact

! material set ids, BC set ids, materials (f=fluid, s=solid), bc types
  integer, parameter :: numsts = 50
  integer, allocatable :: matindx(:), matids(:), imatie(:), ibcsts(:), bcf(:)
  character(3), allocatable :: BCTYPS(:)
  integer :: numflu, numoth, numbcs

  contains

  subroutine init_input()
    use size_m
    implicit none
    allocate(param(200), CPFLD(LDIMT1,3), CPGRP(-5:10,LDIMT1,3), QINTEG(LDIMT3,MAXOBJ))

    allocate(MATYPE(-5:10,LDIMT1), LOCHIS(4,lhis), NMEMBER(MAXOBJ))

    allocate(IFADVC(LDIMT1),IFTMSH(0:LDIMT1), IFNONL(LDIMT1), IFVARP(LDIMT1), IFPSCO(LDIMT1))
    allocate(IFPSO(LDIMT1))
    ifnav => ifadvc(1)

    allocate(HCODE(11,lhis),OCODE(8),DRIVC(5), INITC(15),TEXTSW(100,2))
    turbmod => textsw(1,1)

    allocate(XC(8,LELT),YC(8,LELT),ZC(8,LELT), BC(5,6,LELT,0:LDIMT1))
    allocate(CURVE(6,12,LELT),CERROR(LELT))

    allocate(IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2))

    allocate(CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT), CDOF(6,LELT), IEACT(LELT))

    allocate(MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT), IBCSTS (NUMSTS))
  end subroutine init_input

end module input
