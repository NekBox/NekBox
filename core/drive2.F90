!> \file drive2.F90
!! \brief second level drivers called from drive1.F90

!-------------------------------------------------------------------
!> Transfer array dimensions to common
!-------------------------------------------------------------------
subroutine initdim
  use size_m
  use input
  implicit none

  NX1=LX1
  NY1=LY1
  NZ1=LZ1

  NX2=LX2
  NY2=LY2
  NZ2=LZ2

  NX3=LX3
  NY3=LY3
  NZ3=LZ3

  NXD=LXD
  NYD=LYD
  NZD=LZD


  NELT=LELT
  NELV=LELV
  NDIM=LDIM

  RETURN
end subroutine initdim

!--------------------------------------------------------------------
!>     Initialize and set default values.
!--------------------------------------------------------------------
subroutine initdat
  use kinds, only : DP
  use size_m,   only : lx1, lx2, lelt, nx1, ny1, nz1, nx2, ny2, nz2
  use input,    only : IFCVODE, IFEXPLVIS, ifsplit, param, ccurve, xc, yc, zc
  use parallel, only : ifgprnt 
  use soln,     only : abx1, abx2, aby1, aby2, abz1, abz2, vgradt1, vgradt2
  use tstep,    only : if_full_pres
  implicit none

  integer :: nel8, ntot

!     Set default logicals
  IFCVODE   = .FALSE. 
  IFEXPLVIS = .FALSE. 

  ifsplit = .FALSE. 
  if (lx1 == lx2) ifsplit= .TRUE. 

  if_full_pres = .FALSE. 

!     Turn off (on) diagnostics for communication
  IFGPRNT= .FALSE. 

  param = 0._dp

!     The initialization of CBC is done in READAT

!      LCBC = 3*6*LELT*(LDIMT1+1)
!      CALL BLANK(CBC,LCBC)

  CALL BLANK(CCURVE ,12*LELT)
  NEL8 = 8*LELT
  xc = 0._dp
  yc = 0._dp
  zc = 0._dp

  NTOT=NX1*NY1*NZ1*LELT
  abx1 = 0._dp
  abx2 = 0._dp
  aby1 = 0._dp
  aby2 = 0._dp
  abz1 = 0._dp
  abz2 = 0._dp
  vgradt1 = 0._dp
  vgradt2 = 0._dp

  NTOT=NX2*NY2*NZ2*LELT
!max  CALL RZERO(USRDIV,NTOT)

  RETURN
end subroutine initdat

!---------------------------------------------------------------------
!>     No need to comment !!
!---------------------------------------------------------------------
subroutine comment
  use kinds,  only : DP
  use ctimer, only : dnekclock, ttime
  use geom,   only : ifwcno
  use input,  only : ifadvc, iftran
  use tstep,  only : nid, istep, ifield, nfield, lastep, time, dt, courno
  implicit none

  LOGICAL ::  IFCOUR
  SAVE     IFCOUR
  REAL(DP) :: EETIME0,EETIME1,EETIME2, TTIME_STP
  SAVE   EETIME0,EETIME1,EETIME2
  DATA   EETIME0,EETIME1,EETIME2 /0.0, 0.0, 0.0/


!     Only node zero makes comments.
  IF (NID /= 0) RETURN


  IF (EETIME0 == 0.0 .AND. ISTEP == 1) EETIME0=DNEKCLOCK()
  EETIME1=EETIME2
  EETIME2=DNEKCLOCK()

  IF (ISTEP == 0) THEN
      IFCOUR  = .FALSE. 
      DO 10 IFIELD=1,NFIELD
          IF (IFADVC(IFIELD)) IFCOUR = .TRUE. 
      10 END DO
      IF (IFWCNO) IFCOUR = .TRUE. 
  ELSEIF (ISTEP > 0 .AND. LASTEP == 0 .AND. IFTRAN) THEN
      TTIME_STP = EETIME2-EETIME1   ! time per timestep
      TTIME     = EETIME2-EETIME0   ! sum of all timesteps
      IF(ISTEP == 1) THEN
          TTIME_STP = 0
          TTIME     = 0
      ENDIF
      IF (     IFCOUR) &
      WRITE (6,100) ISTEP,TIME,DT,COURNO,TTIME,TTIME_STP
      IF ( .NOT. IFCOUR) WRITE (6,101) ISTEP,TIME,DT
  ELSEIF (LASTEP == 1) THEN
      TTIME_STP = EETIME2-EETIME1   ! time per timestep
      TTIME     = EETIME2-EETIME0   ! sum of all timesteps
  ENDIF
  100 FORMAT('Step',I7,', t=',1pE14.7,', DT=',1pE14.7 &
  ,', C=',0pF7.3,2(1pE11.4))
  101 FORMAT('Step',I7,', time=',1pE12.5,', DT=',1pE11.3)

  RETURN
end subroutine comment

!------------------------------------------------------------------------
!>     Initialize variables
!------------------------------------------------------------------------
subroutine setvar
  use kinds, only : DP
  use size_m, only : nfield, nid, lorder, nelt, ldimt, lzd, lyd, lxd
  use size_m, only : nxd, nyd, nzd, nelv
  use geom, only : ifgmsh3
  use input, only : param, ifcvode, ifnav, ifnatc, iflomach
  use input, only : ifadvc, iftmsh, ifmvbd, ifmodel, ifmhd, iftran
  use input, only : npscal, ifstrs, ifflow, ifsplit, ngeom, ifheat
  use tstep, only : tolpdf, lastep, iostep, timeio, iocomm, nsteps, fintim
  use tstep, only : dtinit, dt, gtheta, betag, nmxnl, dtlag, nbd, ifield, tolnl
  use tstep, only : prelax, nbdinp, tolrel, tolhdf, pi, ctarg, tolabs, nmxe
  use tstep, only : nelfld, nmxh, nmxp

  implicit none

  integer :: NFLDTM, nfldt, mfield, NABMSH, IADV, IFLD1
  real(DP) :: TLFAC, one

!     Enforce splitting/Uzawa according to the way the code was compiled
  nxd = lxd
  nyd = lyd
  nzd = lzd

!     Geometry on Mesh 3 or 1?
  IFGMSH3 = .TRUE. 
  IF ( IFSTRS )           IFGMSH3 = .FALSE. 
  IF ( .NOT. IFFLOW)        IFGMSH3 = .FALSE. 
  IF ( IFSPLIT )          IFGMSH3 = .FALSE. 

  NGEOM  = 2

  NFIELD = 1
  IF (IFHEAT) THEN
      NFIELD = 2 + NPSCAL
      NFLDTM = 1 + NPSCAL
  ENDIF

  nfldt = nfield
  if (ifmhd) then
      nfldt  = nfield + 1
      nfldtm = nfldtm + 1
  endif


  IF (IFMODEL) write(*,*) "Oops: turb"
#if 0
  IF (IFMODEL) CALL SETTMC
  IF (IFMODEL .AND. IFKEPS) THEN
      write(*,*) "Oops: turb"
      NPSCAL = 1
      NFLDTM = NPSCAL + 1
      IF (LDIMT < NFLDTM) THEN
          WRITE (6,*) 'k-e turbulence model activated'
          WRITE (6,*) 'Insufficient number of field arrays'
          WRITE (6,*) 'Rerun through PRE or change SIZE file'
          call exitt
      ENDIF
      NFIELD = NFIELD + 2
      CALL SETTURB
  ENDIF
#endif
  MFIELD = 1
  IF (IFMVBD) MFIELD = 0

  DO 100 IFIELD=MFIELD,nfldt+(LDIMT-1 - NPSCAL)
      IF (IFTMSH(IFIELD)) THEN
          NELFLD(IFIELD) = NELT
      ELSE
          NELFLD(IFIELD) = NELV
      ENDIF
  100 END DO

  NMXH   = 1000
  if (iftran) NMXH   = 100
  NMXP   = 1000 !  (for testing) 100 !  2000
  NMXE   = 100 !  1000
  NMXNL  = 10  !  100

  PARAM(86) = 0 ! No skew-symm. convection for now

  BETAG  = 0 ! PARAM(3)
  GTHETA = 0 ! PARAM(4)
  DT     = abs(PARAM(12))
  DTINIT = DT
  FINTIM = PARAM(10)
  NSTEPS = int(PARAM(11))
  IOCOMM = int(PARAM(13))
  TIMEIO = PARAM(14)
  IOSTEP = int(PARAM(15))
  LASTEP = 0
  TOLPDF = abs(PARAM(21))
  TOLHDF = abs(PARAM(22))
  TOLREL = abs(PARAM(24))
  TOLABS = abs(PARAM(25))
  CTARG  = PARAM(26)
  NBDINP = int(PARAM(27))
  NABMSH = int(PARAM(28))

  if (nbdinp > lorder) then
      if (nid == 0) then
          write(6,*) 'ERROR: torder > lorder.',nbdinp,lorder
          write(6,*) 'Change SIZEu and recompile entire code.'
      endif
      call exitt
  endif

  if(abs(PARAM(16)) >= 2) IFCVODE = .TRUE. 


!     Check accuracy requested.

  IF (TOLREL <= 0.) TOLREL = 0.01

!     Relaxed pressure iteration; maximum decrease in the residual.

  PRELAX = 0.1*TOLREL
  IF ( .NOT. IFTRAN .AND. .NOT. IFNAV) PRELAX = 1.E-5

!     Tolerance for nonlinear iteration

  TOLNL  = 1.E-4

!     Fintim overrides nsteps
  IF (FINTIM /= 0.) NSTEPS = 1000000000
  IF ( .NOT. IFTRAN ) NSTEPS = 1

!     Print interval defaults to 1
  IF (IOCOMM == 0)  IOCOMM = nsteps+1


!     Set logical for Boussinesq approx (natural convection)
  IFNATC = .FALSE. 
  IF (BETAG > 0.) IFNATC= .TRUE. 
  IF(IFLOMACH) IFNATC = .FALSE. 

!     Set default for mesh integration scheme
  IF (NABMSH <= 0 .OR. NABMSH > 3) THEN
      NABMSH    = NBDINP
      PARAM(28) = (NABMSH)
  ENDIF

!     Set default for mixing length factor
  TLFAC = 0.14
!     IF (PARAM(49) .LE. 0.0) PARAM(49) = TLFAC

!     Courant number only applicable if convection in ANY field.
  IADV  = 0
  IFLD1 = 1
  IF ( .NOT. IFFLOW) IFLD1 = 2
  DO 200 IFIELD=IFLD1,nfldt
      IF (IFADVC(IFIELD)) IADV = 1
  200 END DO

!     If characteristics, need number of sub-timesteps (DT/DS).
!     Current sub-timeintegration scheme: RK4.
!     If not characteristics, i.e. standard semi-implicit scheme,
!     check user-defined Courant number.
  IF (IADV == 1) CALL SETCHAR

!     Initialize order of time-stepping scheme (BD)
!     Initialize time step array.
  NBD    = 0
  dtlag = 0._dp

!     Useful constants
  one = 1.
  PI  = 4.*ATAN(one)

  RETURN
end subroutine setvar

!> \brief Echo the nonzero parameters from the readfile to the logfile
subroutine echopar
  use kinds, only : DP
  use size_m, only : nid, ndim, ldim
  use input, only : reafle, vnekton, nktonv, param
  use string, only : ltrunc
  implicit none

  CHARACTER(132) :: tmp_string
  CHARACTER(1) ::  tmp_string1(132)
  EQUIVALENCE (tmp_string,tmp_string1)

  real(DP) :: vnekmin
  integer :: ls, nparam, j, i

  IF (nid /= 0) RETURN

  OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
  REWIND(UNIT=9)


  READ(9,*,ERR=400)
  READ(9,*,ERR=400) VNEKTON
  NKTONV=VNEKTON
  VNEKMIN=2.5
  IF(VNEKTON < VNEKMIN)THEN
      PRINT*,' Error: This NEKTON Solver Requires a .rea file'
      PRINT*,' from prenek version ',VNEKMIN,' or higher'
      PRINT*,' Please run the session through the preprocessor'
      PRINT*,' to bring the .rea file up to date.'
      call exitt
  ENDIF
  READ(9,*,ERR=400) NDIM
!   error check
  IF(NDIM /= LDIM)THEN
      WRITE(6,10) LDIM,NDIM
      10 FORMAT(//,2X,'Error: This NEKTON Solver has been compiled' &
      /,2X,'       for spatial dimension equal to',I2,'.' &
      /,2X,'       The data file has dimension',I2,'.')
      CALL exitt
  ENDIF

  CALL BLANK(tmp_string,132)
  CALL CHCOPY(tmp_string,REAFLE,132)
  Ls=LTRUNC(tmp_string,132)
  READ(9,*,ERR=400) NPARAM
  WRITE(6,82) NPARAM,(tmp_string1(j),j=1,Ls)

  DO 20 I=1,NPARAM
      CALL BLANK(tmp_string,132)
      READ(9,80,ERR=400) tmp_string
      Ls=LTRUNC(tmp_string,132)
      IF (PARAM(i) /= 0.0) WRITE(6,81) I,(tmp_string1(j),j=1,Ls)
  20 END DO
  80 FORMAT(A132)
  81 FORMAT(I4,3X,132A1)
  82 FORMAT(I4,3X,'Parameters from file:',132A1)
  CLOSE (UNIT=9)
  write(6,*) ' '

!    if(param(2).ne.param(8).and.nid.eq.0) then
!       write(6,*) 'Note VISCOS not equal to CONDUCT!'
!       write(6,*) 'Note VISCOS  =',PARAM(2)
!       write(6,*) 'Note CONDUCT =',PARAM(8)
!    endif

  if (param(62) > 0) then
      if(nid == 0) write(6,*) &
      'enable byte swap for output'
      call set_bytesw_write(1)
  endif

  return

!   Error handling:

  400 CONTINUE
  WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA' &
  ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
  CALL exitt

  WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA' &
  ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
  CALL exitt

  RETURN
end subroutine echopar

!----------------------------------------------------------------------
!> \brief Generate geometry data
!----------------------------------------------------------------------
subroutine gengeom (igeom)
  use size_m, only : nid
  use geom, only : ifgeom
  use tstep, only : istep
  implicit none

  integer, intent(in) :: igeom

  if (nid == 0 .AND. istep <= 1) write(6,*) 'generate geometry data'

  IF (IGEOM == 1) THEN
      RETURN
  ELSEIF (IGEOM == 2) THEN
      if (ifgeom) CALL LAGMASS
      IF (ISTEP == 0) CALL GENCOOR ()
!max      IF (ISTEP >= 1) CALL UPDCOOR
      CALL GEOM1 ()!XM3,YM3,ZM3)
      CALL GEOM2
!max      CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
      CALL SETDEF
      CALL SFASTAX
!max      IF (ISTEP >= 1) CALL EINIT
  ELSEIF (IGEOM == 3) THEN
  !        Take direct stiffness avg of mesh
    write(*,*) "Oops: igeom"
#if 0  
      ifieldo = ifield
      CALL GENCOOR (XM3,YM3,ZM3)
      if (ifheat) then
          ifield = 2
          CALL dssum(xm3,nx3,ny3,nz3)
          call col2 (xm3,tmult,ntot3)
          CALL dssum(ym3,nx3,ny3,nz3)
          call col2 (ym3,tmult,ntot3)
          if (if3d) then
              CALL dssum(xm3,nx3,ny3,nz3)
              call col2 (xm3,tmult,ntot3)
          endif
      else
          ifield = 1
          CALL dssum(xm3,nx3,ny3,nz3)
          call col2 (xm3,vmult,ntot3)
          CALL dssum(ym3,nx3,ny3,nz3)
          call col2 (ym3,vmult,ntot3)
          if (if3d) then
              CALL dssum(xm3,nx3,ny3,nz3)
              call col2 (xm3,vmult,ntot3)
          endif
      endif
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
!max      CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
      CALL SETDEF
      CALL SFASTAX
      ifield = ifieldo
#endif
  ENDIF

  if (nid == 0 .AND. istep <= 1) then
      write(6,*) 'done :: generate geometry data'
      write(6,*) ' '
  endif

  return
end subroutine gengeom

!-----------------------------------------------------------------------
!>    Defines machine specific input and output file names.
subroutine files
  use size_m, only : nid
  use input, only : session, path
  use input, only : reafle, re2fle, fldfle, hisfle, schfle
  use input, only : dmpfle, orefle, nrefle
  use string, only : ltrunc, indx1
  implicit none

  integer :: ls, lpp, lsp, l1, ln, len

  CHARACTER(132) :: NAME, slash
  CHARACTER(1) ::  DMP(4),FLD(4),REA(4),HIS(4),SCH(4) ,ORE(4), NRE(4)
  CHARACTER(1) ::  RE2(4)
  CHARACTER(4) ::  DMP4  ,FLD4  ,REA4  ,HIS4  ,SCH4   ,ORE4  , NRE4
  CHARACTER(4) ::  RE24
  EQUIVALENCE (DMP,DMP4), (FLD,FLD4), (REA,REA4), (HIS,HIS4) &
  , (SCH,SCH4), (ORE,ORE4), (NRE,NRE4) &
  , (RE2,RE24)
  DATA DMP4,FLD4,REA4 /'.dmp','.fld','.rea'/
  DATA HIS4,SCH4      /'.his','.sch'/
  DATA ORE4,NRE4      /'.ore','.nre'/
  DATA RE24           /'.re2'       /
  CHARACTER(78) ::  tmp_string

!    Find out the session name:

!     CALL BLANK(SESSION,132)
!     CALL BLANK(PATH   ,132)

!     ierr = 0
!     IF(NID.EQ.0) THEN
!       OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD',ERR=24)
!       READ(8,10) SESSION
!       READ(8,10) PATH
! 10      FORMAT(A132)
!       CLOSE(UNIT=8)
!       GOTO 23
! 24    ierr = 1
! 23  ENDIF
!     call err_chk(ierr,' Cannot open SESSION.NAME!$')

  slash = '/'

  len = ltrunc(path,132)
  if(indx1(path(len:len),slash,1) < 1) then
      call chcopy(path(len+1:len+1),slash,1)
  endif

!     call bcast(SESSION,132*CSIZE)
!     call bcast(PATH,132*CSIZE)

  CALL BLANK(REAFLE,132)
  CALL BLANK(RE2FLE,132)
  CALL BLANK(FLDFLE,132)
  CALL BLANK(HISFLE,132)
  CALL BLANK(SCHFLE,132)
  CALL BLANK(DMPFLE,132)
  CALL BLANK(OREFLE,132)
  CALL BLANK(NREFLE,132)
  CALL BLANK(NAME  ,132)

!    Construct file names containing full path to host:

  LS=LTRUNC(SESSION,132)
  LPP=LTRUNC(PATH,132)
  LSP=LS+LPP

  call chcopy(name(1:1),path,lpp)
  call chcopy(name(lpp+1:lpp+1),session,ls )
  l1 = lpp+ls+1
  ln = lpp+ls+4


!.rea file
  call chcopy(name(l1:l1),rea , 4)
  call chcopy(reafle    ,name,ln)
!     write(6,*) 'reafile:',reafle

!.re2 file
  call chcopy(name(l1:l1),re2 , 4)
  call chcopy(re2fle    ,name,ln)

!.fld file
  call chcopy(name(l1:l1),fld , 4)
  call chcopy(fldfle    ,name,ln)

!.his file
  call chcopy(name(l1:l1),his , 4)
  call chcopy(hisfle    ,name,ln)

!.sch file
  call chcopy(name(l1:l1),sch , 4)
  call chcopy(schfle    ,name,ln)


!.dmp file
  call chcopy(name(l1:l1),dmp , 4)
  call chcopy(dmpfle    ,name,ln)

!.ore file
  call chcopy(name(l1:l1),ore , 4)
  call chcopy(orefle    ,name,ln)

!.nre file
  call chcopy(name(l1:l1),nre , 4)
  call chcopy(nrefle    ,name,ln)

!    Write the name of the .rea file to the logfile.

  IF (NID == 0) THEN
      CALL CHCOPY(tmp_string,REAFLE,78)
      WRITE(6,1000) tmp_string
      WRITE(6,1001)
      1000 FORMAT(//,2X,'Beginning session:',/,2X,A78)
      1001 FORMAT(/,' ')
  ENDIF


  RETURN

end subroutine files

!----------------------------------------------------------------------
!> \brief setup time-stepping
!!
!! Store old time steps and compute new time step, time and timef.
!! Set time-dependent coefficients in time-stepping schemes.
!----------------------------------------------------------------------
subroutine settime
  use kinds, only : DP
  use geom,  only : ifsurt
  use input, only : ifmvbd, param, ifprint
  use tstep, only : dtlag, ab, abmsh, bd, dt, iocomm
  use tstep, only : istep, nab, nbd, nbdinp, time, timef
  implicit none

  integer :: ilag, irst, nabmsh, nbdmsh

  irst = int(param(46))

!   Set time step.

  DO 10 ILAG=10,2,-1
      DTLAG(ILAG) = DTLAG(ILAG-1)
  10 END DO
  CALL SETDT
  DTLAG(1) = DT
  IF (ISTEP == 1 .AND. irst <= 0) DTLAG(2) = DT

!   Set time.

  TIMEF    = TIME
  TIME     = TIME+DT

!   Set coefficients in AB/BD-schemes.

  CALL SETORDBD
  if (irst > 0) nbd = nbdinp
  bd = 0._dp
  CALL SETBD (BD,DTLAG,NBD)
  NAB = 3
  IF (ISTEP <= 2 .AND. irst <= 0) NAB = ISTEP
  ab = 0._dp
  CALL SETABBD (AB,DTLAG,NAB,NBD)
  IF (IFMVBD) THEN
      NBDMSH = 1
      NABMSH = int(PARAM(28))
      IF (NABMSH > ISTEP .AND. irst <= 0) NABMSH = ISTEP
      IF (IFSURT)          NABMSH = NBD
      abmsh = 0._dp
      CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)
  ENDIF


!   Set logical for printout to screen/log-file

  IFPRINT = .FALSE. 
  IF (IOCOMM > 0 .AND. MOD(ISTEP,IOCOMM) == 0) IFPRINT= .TRUE. 
  IF (ISTEP == 1  .OR. ISTEP == 0           ) IFPRINT= .TRUE. 

  RETURN
  end subroutine settime

!-----------------------------------------------------------------------
!> \brief Compute eigenvalues.
!!
!! Used for automatic setting of tolerances and to find critical
!! time step for explicit mode.
!! Currently eigenvalues are computed only for the velocity mesh.
!-----------------------------------------------------------------------
subroutine geneig (igeom)
  use eigen, only : ifaa, ifae, ifas, ifast, ifga, ifge, ifgs, ifgst
  use input, only : ifflow, ifheat
  use tstep, only : ifield, imesh, tolev, tolhdf, tolhe, tolhr, tolhs
  use tstep, only : tolpdf, tolps

  implicit none

  integer, intent(in) :: igeom

  IF (IGEOM == 1) RETURN

!   Decide which eigenvalues to be computed.

  IF (IFFLOW) THEN
  
      IFAA  = .FALSE. 
      IFAE  = .FALSE. 
      IFAS  = .FALSE. 
      IFAST = .FALSE. 
      IFGA  = .TRUE. 
      IFGE  = .FALSE. 
      IFGS  = .FALSE. 
      IFGST = .FALSE. 
  
  !        For now, only compute eigenvalues during initialization.
  !        For deforming geometries the eigenvalues should be
  !        computed every time step (based on old eigenvectors => more memory)
  
      IMESH  = 1
      IFIELD = 1
      TOLEV  = 1.E-3
      TOLHE  = TOLHDF
      TOLHR  = TOLHDF
      TOLHS  = TOLHDF
      TOLPS  = TOLPDF
      CALL EIGENV
      CALL ESTEIG
  
  ELSEIF (IFHEAT .AND. .NOT. IFFLOW) THEN
  
      CALL ESTEIG
  
  ENDIF

  RETURN
end subroutine geneig

!-----------------------------------------------------------------------
!> \brief Driver for solving the incompressible Navier-Stokes equations.
!!
!! Current version:
!! (1) Velocity/stress formulation.
!! (2) Constant/variable properties.
!! (3) Implicit/explicit time stepping.
!! (4) Automatic setting of tolerances .
!! (5) Lagrangian/"Eulerian"(operator splitting) modes
!-----------------------------------------------------------------------
subroutine fluid (igeom)
  use kinds, only : DP
  use size_m, only : nid
  use ctimer, only : dnekclock
  use input, only : ifnav, ifsplit, iftran
  use tstep, only : ifield, imesh, istep, time

  implicit none

  integer, intent(inout) :: igeom
 
  real(DP) :: ts
   
  ifield = 1
  imesh  = 1
  call unorm
  call settolv

  ts = dnekclock()

  if(nid == 0 .AND. igeom == 1) &
  write(6,*) 'Solving for fluid',ifsplit,iftran,ifnav

  if (ifsplit) then

  !        PLAN 4: TOMBO SPLITTING
  !                - Time-dependent Navier-Stokes calculation (Re>>1).
  !                - Same approximation spaces for pressure and velocity.
  !                - Incompressibe or Weakly compressible (div u .ne. 0).
       
      call plan4
      igeom = 2
#if 0
      call twalluz (igeom) ! Turbulence model
#endif
!max      call chkptol         ! check pressure tolerance
!max      if (param(55) /= 0) call vol_flow        ! check for fixed flow rate

  elseif (iftran) then
#if 0

  !        call plan1 (igeom)       !  Orig. NEKTON time stepper

      call plan3 (igeom)       !  Same as PLAN 1 w/o nested iteration
  !  Std. NEKTON time stepper  !
      if (ifmodel)    call twalluz (igeom) ! Turbulence model
      if (igeom >= 2) call chkptol         ! check pressure tolerance
!max      if (igeom >= 2 and param(55) /= 0) call vol_flow ! check for fixed flow rate

#endif
  else   !  steady Stokes, non-split

  !             - Steady/Unsteady Stokes/Navier-Stokes calculation.
  !             - Consistent approximation spaces for velocity and pressure.
  !             - Explicit treatment of the convection term.
  !             - Velocity/stress formulation.
#if 0
      call plan1 (igeom) ! The NEKTON "Classic".
#endif
  endif

  if(nid == 0 .AND. igeom >= 2) &
  write(*,'(4x,i7,1x,1p2e12.4,a)') &
  istep,time,dnekclock()-ts,' Fluid done'

  return
end subroutine fluid

!-----------------------------------------------------------------------
!> \brief Driver for temperature or passive scalar.
!!
!!  Current version:
!!  (1) Varaiable properties.
!!  (2) Implicit time stepping.
!!  (3) User specified tolerance for the Helmholtz solver
!!      (not based on eigenvalues).
!!  (4) A passive scalar can be defined on either the
!!      temperatur or the velocity mesh.
!!  (5) A passive scalar has its own multiplicity (B.C.).
!-----------------------------------------------------------------------
subroutine heat (igeom)
  use kinds, only : DP
  use size_m, only : nid, nfield
  use ctimer, only : dnekclock
  use input, only : ifcvode, ifsplit, iftmsh
  use tstep, only : ifield, imesh, istep, time

  implicit none

  integer, intent(inout) :: igeom
  real(DP) :: ts
  integer :: intype, igeo

  ts = dnekclock()

  if (nid == 0 .AND. igeom == 1) &
  write(*,'(13x,a)') 'Solving for heat'

  if (ifcvode) then
#if 0

      call cdscal_cvode(igeom)
      igeom = 2
#endif
  elseif (ifsplit) then

      do igeo=1,2
          do ifield=2,nfield
              intype        = -1
              if ( .NOT. iftmsh(ifield)) imesh = 1
              if (     iftmsh(ifield)) imesh = 2
              call unorm
              call settolt
              call cdscal (igeo)
          enddo
      enddo
      igeom = 2

  else  ! PN-PN-2

      do ifield=2,nfield
          intype        = -1
          if ( .NOT. iftmsh(ifield)) imesh = 1
          if (     iftmsh(ifield)) imesh = 2
          call unorm
          call settolt
          call cdscal (igeom)
      enddo

  endif

  if (nid == 0 .AND. igeom >= 2) &
  write(*,'(4x,i7,1x,1p2e12.4,a)') &
  istep,time,dnekclock()-ts,' Heat done'

  return
end subroutine heat

!-----------------------------------------------------------------------
!> \brief zero the ctimer
subroutine time00
  use ctimer
  implicit none

  nmxmf=0
  nmxms=0
  ndsum=0
  nvdss=0
  nsett=0
  ncdtp=0
  npres=0
  nmltd=0
  ngsum=0
  nprep=0
  ndsnd=0
  ndadd=0
  nhmhz=0
  naxhm=0
  ngop =0
  nusbc=0
  ncopy=0
  ninvc=0
  ninv3=0
  nsolv=0
  nslvb=0
  nddsl=0
  ncrsl=0
  ndott=0
  nbsol=0
  nadvc=0
  nspro=0

  tmxmf=0.0
  tmxms=0.0
  tdsum=0.0
  tvdss=0.0
  tvdss=0.0
  tdsmn=9.9e9
  tdsmx=0.0
  tsett=0.0
  tcdtp=0.0
  tpres=0.0
  teslv=0.0
  tmltd=0.0
  tgsum=0.0
  tgsmn=9.9e9
  tgsmx=0.0
  tprep=0.0
  tdsnd=0.0
  tdadd=0.0
  thmhz=0.0
  taxhm=0.0
  tgop =0.0
  tusbc=0.0
  tcopy=0.0
  tinvc=0.0
  tinv3=0.0
  tsolv=0.0
  tslvb=0.0
  tddsl=0.0
  tcrsl=0.0
  tdott=0.0
  tbsol=0.0
  tbso2=0.0
  tspro=0.0
  tadvc=0.0
  ttime=0.0

  return
end subroutine time00

!-----------------------------------------------------------------------
!> \brief print run statistics from ctimer
subroutine runstat
#ifndef NOTIMER
  use kinds, only : DP
  use ctimer, only : ifsync, nadvc, naxhm, ncdtp, ncrsl, ndadd, nddsl
  use ctimer, only : pinvc, pinv3, phmhz, peslv, pdsum, pddsl, pdadd
  use ctimer, only : pvdss, pusbc, pspro, psolv, ppres, pprep, pmltd, pcrsl
  use ctimer, only : pcdtp, paxhm
  use ctimer, only : tgop_sync, tgop, teslv, tdsum, tddsl, tdadd, tcrsl, tcdtp
  use ctimer, only : tspro, tsolv, tpres, tprep, tmltd, tinvc, tinv3, thmhz
  use ctimer, only : taxhm, tadvc, twal, tvdss, tusbc, tttstp, ttime, tsyc
  use ctimer, only : nwal, nvdss, nusbc, nsyc, nspro, nsolv, npres, nprep
  use ctimer, only : ninvc, ninv3, nhmhz, ngop, neslv, nmltd, ndsum
  use ctimer, only : dnekclock
  use size_m, only : nid
  use parallel, only : np
  implicit none

  real(DP) :: min_dsum, max_dsum, avg_dsum
  real(DP) :: min_vdss, max_vdss, avg_vdss
  real(DP) :: min_gop,  max_gop,  avg_gop
  real(DP) :: min_gop_sync,  max_gop_sync,  avg_gop_sync
  real(DP) :: min_crsl, max_crsl, avg_crsl
  real(DP) :: min_usbc, max_usbc, avg_usbc
  real(DP) :: min_syc, max_syc, avg_syc
  real(DP) :: min_wal, max_wal, avg_wal
  real(DP) :: min_irc, max_irc, avg_irc
  real(DP) :: min_isd, max_isd, avg_isd
  real(DP) :: min_comm, max_comm, avg_comm
  real(DP) :: tstop, tirc, tisd, trc, tsd, tcomm, wwork, padvc
  integer :: nirc, nisd

  real(DP) :: comm_timers(8)
  integer :: comm_counters(8)
  character(132) :: s132

  tstop=dnekclock()
  tttstp=ttime         ! sum over all timesteps

!      call opcount(3)      ! print op-counters

  call nek_comm_getstat(comm_timers,comm_counters)
  tgop      = comm_timers(1)
  tgop_sync = comm_timers(2)
  twal      = comm_timers(3)
  tsyc      = comm_timers(4)
  tirc      = comm_timers(5)
  tisd      = comm_timers(6)
  trc       = comm_timers(7)
  tsd       = comm_timers(8)
  ngop      = comm_counters(1)
  nwal      = comm_counters(3)
  nsyc      = comm_counters(4)
  nirc      = comm_counters(5)
  nisd      = comm_counters(6)

  tcomm  = tisd + tirc + tsyc + tgop + twal + trc + tsd
  min_comm = tcomm
  call gop(min_comm,wwork,'m  ',1)
  max_comm = tcomm
  call gop(max_comm,wwork,'M  ',1)
  avg_comm = tcomm
  call gop(avg_comm,wwork,'+  ',1)
  avg_comm = avg_comm/np

  min_isd = tisd
  call gop(min_isd,wwork,'m  ',1)
  max_isd = tisd
  call gop(max_isd,wwork,'M  ',1)
  avg_isd = tisd
  call gop(avg_isd,wwork,'+  ',1)
  avg_isd = avg_isd/np

  min_irc = tirc
  call gop(min_irc,wwork,'m  ',1)
  max_irc = tirc
  call gop(max_irc,wwork,'M  ',1)
  avg_irc = tirc
  call gop(avg_irc,wwork,'+  ',1)
  avg_irc = avg_irc/np

  min_syc = tsyc
  call gop(min_syc,wwork,'m  ',1)
  max_syc = tsyc
  call gop(max_syc,wwork,'M  ',1)
  avg_syc = tsyc
  call gop(avg_syc,wwork,'+  ',1)
  avg_syc = avg_syc/np

  min_wal = twal
  call gop(min_wal,wwork,'m  ',1)
  max_wal = twal
  call gop(max_wal,wwork,'M  ',1)
  avg_wal = twal
  call gop(avg_wal,wwork,'+  ',1)
  avg_wal = avg_wal/np

  min_gop = tgop
  call gop(min_gop,wwork,'m  ',1)
  max_gop = tgop
  call gop(max_gop,wwork,'M  ',1)
  avg_gop = tgop
  call gop(avg_gop,wwork,'+  ',1)
  avg_gop = avg_gop/np

  min_gop_sync = tgop_sync
  call gop(min_gop_sync,wwork,'m  ',1)
  max_gop_sync = tgop_sync
  call gop(max_gop_sync,wwork,'M  ',1)
  avg_gop_sync = tgop_sync
  call gop(avg_gop_sync,wwork,'+  ',1)
  avg_gop_sync = avg_gop_sync/np

  min_vdss = tvdss
  call gop(min_vdss,wwork,'m  ',1)
  max_vdss = tvdss
  call gop(max_vdss,wwork,'M  ',1)
  avg_vdss = tvdss
  call gop(avg_vdss,wwork,'+  ',1)
  avg_vdss = avg_vdss/np

  min_dsum = tdsum
  call gop(min_dsum,wwork,'m  ',1)
  max_dsum = tdsum
  call gop(max_dsum,wwork,'M  ',1)
  avg_dsum = tdsum
  call gop(avg_dsum,wwork,'+  ',1)
  avg_dsum = avg_dsum/np


  min_crsl = tcrsl
  call gop(min_crsl,wwork,'m  ',1)
  max_crsl = tcrsl
  call gop(max_crsl,wwork,'M  ',1)
  avg_crsl = tcrsl
  call gop(avg_crsl,wwork,'+  ',1)
  avg_crsl = avg_crsl/np

  min_usbc = tusbc
  call gop(min_usbc,wwork,'m  ',1)
  max_usbc = tusbc
  call gop(max_usbc,wwork,'M  ',1)
  avg_usbc = tusbc
  call gop(avg_usbc,wwork,'+  ',1)
  avg_usbc = avg_usbc/np

  tttstp = tttstp + 1e-7
  if (nid == 0) then
      write(6,'(A)') 'runtime statistics:'
      write(6,*) 'total time',tttstp

  !         pcopy=tcopy/tttstp
  !         write(6,*) 'copy time',ncopy,tcopy,pcopy
  !         pmxmf=tmxmf/tttstp
  !         write(6,*) 'mxmf time',nmxmf,tmxmf,pmxmf

      pinv3=tinv3/tttstp
      write(6,*) 'inv3 time',ninv3,tinv3,pinv3
      pinvc=tinvc/tttstp
      write(6,*) 'invc time',ninvc,tinvc,pinvc
      pmltd=tmltd/tttstp
      write(6,*) 'mltd time',nmltd,tmltd,pmltd
      pcdtp=tcdtp/tttstp
      write(6,*) 'cdtp time',ncdtp,tcdtp,pcdtp
      peslv=teslv/tttstp
      write(6,*) 'eslv time',neslv,teslv,peslv

  !        Pressure solver timings
      ppres=tpres/tttstp
      write(6,*) 'pres time',npres,tpres,ppres

  !        Coarse grid solver timings
      pcrsl=tcrsl/tttstp
      write(6,*) 'crsl time',ncrsl,tcrsl,pcrsl
      write(6,*) 'crsl min ',min_crsl
      write(6,*) 'crsl max ',max_crsl
      write(6,*) 'crsl avg ',avg_crsl

  !        Helmholz solver timings
      phmhz=thmhz/tttstp
      write(6,*) 'hmhz time',nhmhz,thmhz,phmhz

      pspro=tspro/tttstp
      write(6,*) 'spro time',nspro,tspro,pspro

  !        USERBC timings
      pusbc=tusbc/tttstp
      write(6,*) 'usbc time',nusbc,tusbc,pusbc
      write(6,*) 'usbc min ',min_usbc
      write(6,*) 'usbc max ',max_usbc
      write(6,*) 'usb  avg ',avg_usbc

  !        Axhelm timings
      paxhm=taxhm/tttstp
      write(6,*) 'axhm time',naxhm,taxhm,paxhm

  !        Convection timings
      padvc=tadvc/tttstp
      write(6,*) 'advc time',nadvc,tadvc,padvc

  !        Vector direct stiffness summuation timings
      pvdss=tvdss/tttstp
      write(6,*) 'vdss time',nvdss,tvdss,pvdss
      write(6,*) 'vdss min ',min_vdss
      write(6,*) 'vdss max ',max_vdss
      write(6,*) 'vdss avg ',avg_vdss

  !        Direct stiffness summuation timings
      pdsum=tdsum/tttstp
      write(6,*) 'dsum time',ndsum,tdsum,pdsum
      write(6,*) 'dsum min ',min_dsum
      write(6,*) 'dsum max ',max_dsum
      write(6,*) 'dsum avg ',avg_dsum

  !         pgsum=tgsum/tttstp
  !         write(6,*) 'gsum time',ngsum,tgsum,pgsum

  !         pdsnd=tdsnd/tttstp
  !         write(6,*) 'dsnd time',ndsnd,tdsnd,pdsnd

      pdadd=tdadd/tttstp
      write(6,*) 'dadd time',ndadd,tdadd,pdadd

  !         pdsmx=tdsmx/tttstp
  !         write(6,*) 'dsmx time',ndsmx,tdsmx,pdsmx
  !         pdsmn=tdsmn/tttstp
  !         write(6,*) 'dsmn time',ndsmn,tdsmn,pdsmn
  !         pslvb=tslvb/tttstp
  !         write(6,*) 'slvb time',nslvb,tslvb,pslvb
      pddsl=tddsl/tttstp
      write(6,*) 'ddsl time',nddsl,tddsl,pddsl
  
      psolv=tsolv/tttstp
      write(6,*) 'solv time',nsolv,tsolv,psolv

  !         psett=tsett/tttstp
  !         write(6,*) 'sett time',nsett,tsett,psett

      pprep=tprep/tttstp
      write(6,*) 'prep time',nprep,tprep,pprep
  !         pbsol=tbsol/tttstp
  !         write(6,*) 'bsol time',nbsol,tbsol,pbsol
  !         pbso2=tbso2/tttstp
  !         write(6,*) 'bso2 time',nbso2,tbso2,pbso2

#ifdef MPITIMER
      write(6,'(/,A)') 'MPI timings'
  !        MPI timings
      write(6,*) 'total comm time',tcomm, max_comm/ttime
      write(6,*) 'comm min ',min_comm
      write(6,*) 'comm max ',max_comm
      write(6,*) 'comm avg ',avg_comm

  !        MPI_Barrier timings
      psyc=tsyc/tcomm
      write(6,*) 'barrier time',nsyc,tsyc,psyc
      write(6,*) 'barrier min ',min_syc
      write(6,*) 'barrier max ',max_syc
      write(6,*) 'barrier avg ',avg_syc

  !        MPI_Waitall timings
      pwal=twal/tcomm
      write(6,*) 'waitall time',nwal,twal,pwal
      write(6,*) 'waitall min ',min_wal
      write(6,*) 'waitall max ',max_wal
      write(6,*) 'waitall avg ',avg_wal

  !        MPI_Allreduce timings
      pgop=tgop/tcomm
      write(6,*) 'allreduce  time',ngop,tgop,pgop
      write(6,*) 'allreduce  min ',min_gop
      write(6,*) 'allreduce  max ',max_gop
      write(6,*) 'allreduce  avg ',avg_gop

  !        MPI_Allreduce(sync) timings
      pgop_sync=tgop_sync/tcomm
      write(6,*) 'allreduce_sync  time',tgop_sync,pgop_sync
      write(6,*) 'allreduce_sync  min ',min_gop_sync
      write(6,*) 'allreduce_sync  max ',max_gop_sync
      write(6,*) 'allreduce_sync  avg ',avg_gop_sync
#endif
  endif

  if (nid == 0)   & ! header for timing
  write(6,1) 'tusbc','tdadd','tcrsl','tvdss','tdsum',' tgop',ifsync
  1 format(/,'#',2x,'nid',6(7x,a5),4x,'qqq',1x,l4)

  call blank(s132,132)
  write(s132,132) nid,tusbc,tdadd,tcrsl,tvdss,tdsum,tgop
  132 format(i12,1p6e12.4,' qqq')
  call pprint_all(s132,132,6)

#endif

  return
end subroutine runstat

!-----------------------------------------------------------------------
!> \brief pretty-print from runstat
subroutine pprint_all(s,n_in,io)
  use size_m, only : nid
  use parallel, only : np
  use string, only : ltrunc
  implicit none

  integer, intent(in) :: n_in, io
  character(n_in), intent(in) :: s
  character(132) :: w
  integer :: n, mtag, m, i, k, l

  n = min(132,n_in)

  mtag = 999
  m    = 1
  call nekgsync()

  if (nid == 0) then
      l = ltrunc(s,n)
      write(io,1) (s(k:k),k=1,l)
      1 format(132a1)

      do i=1,np-1
          call csend(mtag,s,1,i,0)    ! send handshake
          m = 132
          call blank(w,m)
          call crecv(i,w,m)
          if (m <= 132) then
              l = ltrunc(w,m)
              write(io,1) (w(k:k),k=1,l)
          else
              write(io,*) 'pprint long message: ',i,m
              l = ltrunc(w,132)
              write(io,1) (w(k:k),k=1,l)
          endif
      enddo
  else
      call crecv(mtag,w,m)          ! wait for handshake
      l = ltrunc(s,n)
      call csend(nid,s,l,0,0)       ! send data to node 0
  endif
  return
end subroutine pprint_all

!-----------------------------------------------------------------------
!> \brief init opcounter
subroutine opcount(ICALL)
  use opctr, only : nrout, dcount, ncall, dct, maxrts
  implicit none

  integer, intent(in) :: icall
  integer :: i

  if (icall == 1) then
      nrout=0
  endif
  if (icall == 1 .OR. icall == 2) then
      dcount = 0.0
      do 100 i=1,maxrts
          ncall(i) = 0
          dct(i)   = 0.0
      100 END DO
  endif
  if (icall == 3) then
    write(*,*) "Oops: icall = 3"
#if 0
  !        Sort and print out diagnostics
  
      if (nid == 0) then
          write(6,*) nid,' opcount',dcount
          do i = 1,np-1
              call csend(i,idum,4,i,0)
              call crecv(i,ddcount,wdsize)
              write(6,*) i,' opcount',ddcount
          enddo
      else
          call crecv (nid,idum,4)
          call csend (nid,dcount,wdsize,0,0)
      endif

      dhc = dcount
      call gop(dhc,dwork,'+  ',1)
      if (nid == 0) then
          write(6,*) ' TOTAL OPCOUNT',dhc
      endif
  
      CALL DRCOPY(rct,dct,nrout)
      CALL SORT(rct,ind,nrout)
      CALL CHSWAPR(rname,6,ind,nrout,sname)
      call iswap(ncall,ind,nrout,idum)
  
      if (nid == 0) then
          do 200 i=1,nrout
              write(6,201) nid,rname(i),rct(i),ncall(i)
          200 END DO
          201 format(2x,' opnode',i4,2x,a6,g18.7,i12)
      endif
#endif
  endif
  return
end subroutine opcount

!-----------------------------------------------------------------------
!> \brief count degrees of freedom
subroutine dofcnt
  use kinds, only : DP, i8
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use size_m, only : nx2, ny2, nz2
  use input, only : ifflow, ifsplit
  use parallel, only : nvtot
  use soln, only : vmult, tmult, tmask, v1mask, pmask
  implicit none

  real(DP) :: WORK(1)!LCTMP1)

  integer :: nxyz, nel, ntot1, ntot2
  real(DP) :: vpts, ppts
  real(DP), external :: glsum
  integer, external :: glsc2, i8glsum
  integer(i8) :: ntot,ntotp,ntotv

  nxyz  = nx1*ny1*nz1
  nel   = nelv

! unique points on v-mesh
  vpts = glsum(vmult,nel*nxyz) + .1
  nvtot=int(vpts)

! unique points on pressure mesh
  work(1)=nel*nxyz
  ppts = glsum(work,1) + .1
  ntot=int(ppts)

  if (nid == 0) write(6,'(A,2i13)') &
  'gridpoints unique/tot: ',nvtot,ntot

  ntot1=nx1*ny1*nz1*nelv
  ntot2=nx2*ny2*nz2*nelv

  ntotv = glsc2(tmult,tmask,ntot1)
  ntotp = i8glsum(ntot2,1)

  if (ifflow)  ntotv = glsc2(vmult,v1mask,ntot1)
  if (ifsplit) ntotp = glsc2(vmult,pmask ,ntot1)
  if (nid == 0) write(6,*) ' dofs:',ntotv,ntotp

  return
end subroutine dofcnt
!-----------------------------------------------------------------------
