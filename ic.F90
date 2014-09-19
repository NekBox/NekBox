!-----------------------------------------------------------------------
!> \brief Set initial conditions.
!-----------------------------------------------------------------------
subroutine setics
  use kinds, only : DP, i8
  use size_m, only : lx1, ly1, lz1, lelv, ldimt, ldimt1
  use size_m, only : nx1, ny1, nz1, nelt, nx2, ny2, nz2, nelv
  use size_m, only : nid, lpert, npert, nfield
  use geom, only : ifvcor, xm1, ym1, zm1
  use input, only : ifheat, ifsplit, ifflow, ifmhd, ifpert, ifmodel, ifkeps
  use input, only : param, ifmvbd, npscal
  use input, only : iftmsh, ifldmhd, ifadvc
  use mvgeom, only : wx, wy, wz
  use parallel, only : nelgv
  use soln, only : vx, vy, vz, pr, t, jp, vmult, bx, by, bz
  use soln, only : tmult
  use tstep, only : ifield, nbdinp, time, timeio, nelfld, ntdump
  implicit none
   
  logical  iffort(  ldimt1,0:lpert) &
  , ifrest(0:ldimt1,0:lpert) &
  , ifprsl(  ldimt1,0:lpert)
   
  LOGICAL ::  IFANYP
  real(DP), allocatable :: work(:,:,:,:)
  integer(i8) :: ntotg,nn

  real :: psmax(LDIMT)

  integer :: nxyz2, ntot2, nxyz1, ntott, ntotv, irst, maxfld, mfldt
  integer :: itest, nbdmax, nbdsav, i, ntot, ifldsave, nfiles
  real(DP) :: rdif, rtotg, vxmax, vymax, vzmax, prmax, ttmax, small
  real(DP) :: xxmax, yymax, zzmax
  real(DP), external :: glsum, glamax, glmin, glmax

  if(nid == 0) write(6,*) 'set initial conditions'

!   Initialize all fields:
  nxyz2=nx2*ny2*nz2
  ntot2=nxyz2*nelv
  nxyz1=nx1*ny1*nz1
  ntott=nelt*nxyz1
  ntotv=nelv*nxyz1

  vx = 0._dp; vy = 0._dp; vz = 0._dp; pr = 0._dp; t = 0._dp

  jp = 0                  ! set counter for perturbation analysis

  irst = int(param(46))        ! for lee's restart (rarely used)
  if (irst > 0)  call setup_convect(2)


!   If moving geometry then add a perturbation to the
!   mesh coordinates (see Subroutine INIGEOM)

!max    if (ifmvbd) call ptbgeom

!   Find out what type of i.c. is requested
!   Current options:

!   (1) - User specified fortran function (default is zero i.c.)
!   (2) - Restart from file(s)
!   (3) - Activate pre-solver => steady diffusion / steady Stokes

!   If option (2) is requested, also return with the name of the
!   restart file(s) together with the associated dump number

  call slogic (iffort,ifrest,ifprsl,nfiles)

!   Set up proper initial values for turbulence model arrays
#if 0
  IF (IFMODEL) CALL PRETMIC
#endif

!    ***** TEMPERATURE AND PASSIVE SCALARS ******

!   Check if any pre-solv necessary for temperature/passive scalars

  IFANYP = .FALSE. 
  DO 100 IFIELD=2,NFIELD
      IF (IFPRSL(IFIELD,jp)) THEN
          IF (NID == 0) WRITE(6,101) IFIELD
          IFANYP = .TRUE. 
      ENDIF
  100 END DO
  101 FORMAT(2X,'Using PRESOLVE option for field',I2,'.')

!   Fortran function initial conditions for temp/pass. scalars.
  maxfld = nfield
  if (ifmodel .AND. ifkeps) maxfld = nfield-2
  if (ifmhd) maxfld = npscal+3

!   Always call nekuic (pff, 12/7/11)
  do ifield=1,maxfld
      if (nid == 0) write(6,*) 'nekuic (1) for ifld ', ifield
      call nekuic
  enddo

!   If any pre-solv, do pre-solv for all temperatur/passive scalar fields
!max    if (ifanyp) call prsolvt

  jp = 0 ! jp=0 --> base field, not perturbation field
  do 200 ifield=2,maxfld
      if (iffort(ifield,jp)) then
          if (nid == 0) write(6,*) 'call nekuic for ifld ', ifield
          call nekuic
      endif
  200 END DO

  if (ifpert) then
      ifield=2
      do jp=1,npert
          if (nid == 0) write(6,*) 'nekuicP',ifield,jp,iffort(ifield,jp)
          if (iffort(ifield,jp)) call nekuic
      enddo
  endif
  jp = 0
       

  call nekgsync()
  call restart_driver(nfiles) !  Check restart files
  call nekgsync()


!    ***** VELOCITY ******
!   (If restarting for V, we're done,
!   ...else, do pre-solv for fluid if requested.)

  ifield = 1
!max    if (ifprsl(ifield,jp)) call prsolvv


!   Fortran function initial conditions for velocity.
  ifield = 1
  if (iffort(ifield,jp)) then
      if (nid == 0) write(6,*) 'call nekuic for vel  '
      call nekuic
  endif

  if (ifpert) then
      ifield=1
      do jp=1,npert
          if (iffort(ifield,jp)) call nekuic
          if (nid == 0) write(6,*) 'ic vel pert:',iffort(1,jp),jp
      enddo
  endif
  jp = 0

  ntotv = nx1*ny1*nz1*nelv

!   Fortran function initial conditions for turbulence k-e model
  if (ifmodel .AND. ifkeps) then
      mfldt = nfield - 1
      do 300 ifield=mfldt,nfield
          if (iffort(ifield,jp)) call nekuic
      300 END DO
  endif

!   Initial mesh velocities
  if (ifmvbd) call opcopy (wx,wy,wz,vx,vy,vz)
!max    if (ifmvbd .AND. .NOT. ifrest(0,jp)) call meshv (2)

!   Compute additional initial values for turbulence model arrays
!   based on I.C.
#if 0
  if (ifmodel) call postmic
#endif

!   If convection-diffusion of a passive scalar with a fixed velocity field,
!   make sure to fill up lagged arrays since this will not be done in
!   the time-stepping procedure (no flow calculation) (01/18/91 -EMR).

  if ( .NOT. ifflow .AND. ifheat) then
      ITEST=0
      DO 400 IFIELD=2,NFIELD
          IF (IFADVC(IFIELD)) ITEST=1
      400 END DO
      IF (ITEST == 1) THEN
          NBDMAX = 3
          NBDSAV = NBDINP
          NBDINP = NBDMAX
          DO 500 I=1,NBDMAX
              CALL LAGVEL
          500 END DO
          NBDINP = NBDSAV
      ENDIF
  ENDIF
       
!   Ensure that all processors have the same time as node 0.
  if (nid /= 0) time=0.0
  time=glsum(time,1)
  ntdump=0
  if (timeio /= 0.0) ntdump = int( time/timeio )

!   Ensure that initial field is continuous!

  nxyz1=nx1*ny1*nz1
  ntott=nelt*nxyz1
  ntotv=nelv*nxyz1
  nn = nxyz1
  ntotg=nelgv*nn

  ifield = 2
  if (ifflow) ifield = 1
  allocate(work(lx1,ly1,lz1,lelv))
  call rone(work,ntotv)
  ifield = 1
  call dssum(work)
  work = work * vmult
  rdif = glsum(work,ntotv)
  deallocate(work)
  rtotg = ntotg
  rdif = (rdif-rtotg)/rtotg
  if (abs(rdif) > 1e-14) then
      if (nid == 0) write(*,*) 'ERROR: dssum test has failed!',rdif
      call exitt
  endif

  vxmax = glamax(vx,ntotv)
  vymax = glamax(vy,ntotv)
  vzmax = glamax(vz,ntotv)
  prmax = glamax(pr,ntot2)

  ntot = nxyz1*nelfld(2)
  ttmax = glamax(t ,ntot)

  do i=1,NPSCAL
      ntot = nx1*ny1*nz1*nelfld(i+2)
      psmax(i) = glamax(T(1,1,1,1,i+1),ntot)
  enddo

  small=1.0E-20
  ifldsave = ifield
!  if (vxmax == 0.0) call perturb(vx,1,small)
!  if (vymax == 0.0) call perturb(vy,1,small)
!  if (vzmax == 0.0) call perturb(vz,1,small)
!  if (prmax == 0.0 .AND. ifsplit) call perturb(pr,1,small)
!  if (ttmax == 0.0) call perturb(t ,2,small)

  do i=1,npscal
      ntot = nxyz1*nelfld(i+2)
!      if(psmax(i) == 0) call perturb(t(1,1,1,1,1+i),i+2,small)
  enddo
  ifield = ifldsave
      
  if(ifflow) then
      ifield = 1
      call opdssum(vx,vy,vz)
      vx = vx * vmult; vy = vy * vmult; vz = vz * vmult
      if (ifsplit) call dsavg(pr)  ! continuous pressure
      if (ifvcor)  call ortho(pr)  ! remove any mean
  endif

  if (ifmhd) then
      ifield = ifldmhd
      call opdssum(bx,by,bz)
      bx = bx * vmult; by = by * vmult; bz = bz * vmult
  endif

  if (ifheat) then
      ifield = 2
      call dssum(t)
      t = t * tmult
      do ifield=3,nfield
          call dssum(t(1,1,1,1,i-1))
          if(iftmsh(ifield)) then
              t(:,:,:,:,i-1) = t(:,:,:,:,i-1) * tmult(:,:,:,:,1)
!              call col2 (t(1,1,1,1,i-1),tmult,ntott)
          else
              t(:,:,:,:,i-1) = t(:,:,:,:,i-1) * vmult
!              call col2 (t(1,1,1,1,i-1),vmult,ntotv)
          endif
      enddo
  endif

  if (ifpert) then
#if 0
      do jp=1,npert
          ifield = 1
          call opdssum(vxp(1,jp),vyp(1,jp),vzp(1,jp))
          call opcolv (vxp(1,jp),vyp(1,jp),vzp(1,jp),vmult)
          ifield = 2
          call dssum(tp(1,1,jp))
          call col2 (tp(1,1,jp),tmult,ntotv)
      !           note... must be updated for addl pass. scal's. pff 4/26/04
          vxmax = glamax(vxp(1,jp),ntotv)
          vymax = glamax(vyp(1,jp),ntotv)
          if (nid == 0) write(6,111) jp,vxmax,vymax
          111 format(i5,1p2e12.4,' max pert vel')
      enddo
#endif
  endif
  jp = 0

! print min values
  xxmax = glmin(xm1,ntott)
  yymax = glmin(ym1,ntott)
  zzmax = glmin(zm1,ntott)

  vxmax = glmin(vx,ntotv)
  vymax = glmin(vy,ntotv)
  vzmax = glmin(vz,ntotv)
  prmax = glmin(pr,ntot2)

  ntot = nxyz1*nelfld(2)
  ttmax = glmin(t ,ntott)

  do i=1,ldimt-1
      ntot = nxyz1*nelfld(i+2)
      psmax(i) = glmin(T(1,1,1,1,i+1),ntot)
  enddo

  if (nid == 0) then
      write(6,19) xxmax,yymax,zzmax
      19 format(' xyz min  ',5g13.5)
  endif
  if (nid == 0) then
      write(6,20) vxmax,vymax,vzmax,prmax,ttmax
      20 format(' uvwpt min',5g13.5)
  endif
  if (ldimt-1 > 0) then
      if (nid == 0) write(6,21) (psmax(i),i=1,LDIMT-1)
      21 format(' PS min   ',50g13.5)
  endif

! print max values
  xxmax = glmax(xm1,ntott)
  yymax = glmax(ym1,ntott)
  zzmax = glmax(zm1,ntott)

  vxmax = glmax(vx,ntotv)
  vymax = glmax(vy,ntotv)
  vzmax = glmax(vz,ntotv)
  prmax = glmax(pr,ntot2)

  ntot = nxyz1*nelfld(2)
  ttmax = glmax(t ,ntott)

  do i=1,ldimt-1
      ntot = nxyz1*nelfld(i+2)
      psmax(i) = glmax(T(1,1,1,1,i+1),ntot)
  enddo

  if (nid == 0) then
      write(6,16) xxmax,yymax,zzmax
      16 format(' xyz max  ',5g13.5)
  endif

  if (nid == 0) then
      write(6,17) vxmax,vymax,vzmax,prmax,ttmax
      17 format(' uvwpt max',5g13.5)
  endif

  if (ldimt-1 > 0) then
      if (nid == 0)  then
          write(6,18) (psmax(i),i=1,ldimt-1)
          18 format(' PS max   ',50g13.5)
      endif
  endif


  if (ifrest(0,jp)) then !  mesh has been read in.
      if (nid == 0) write(6,*) 'Restart: recompute geom. factors.'
      call geom_reset(1)  !  recompute geometric factors
  endif

!   ! save velocity on fine mesh for dealiasing
  call setup_convect(2)

!   call outpost(vx,vy,vz,pr,t,'   ')
!   call exitti('setic exit$',nelv)

  if(nid == 0) then
      write(6,*) 'done :: set initial conditions'
      write(6,*) ' '
  endif

  return
end subroutine setics

!---------------------------------------------------------------------
!> \brief Set up logicals for initial conditions.
!---------------------------------------------------------------------
subroutine slogic (iffort,ifrest,ifprsl,nfiles)
  use size_m, only : nfield, npert, nid, ldimt1, lpert
  use input, only : ifmhd, initc, npscal, ifpert
  use restart, only : ifgetx, ifgetu, ifgett, ifgtps
  use string, only : indx1, ltrunc, indx_cut, csplit, ljust, capit
  implicit none

  logical  iffort(  ldimt1,0:lpert) &
  , ifrest(0:ldimt1,0:lpert) &
  , ifprsl(  ldimt1,0:lpert)

  character(132) :: line,fname,cdum
  character(2) ::  s2
  character(1) ::  line1(132)
  equivalence (line1,line)

  integer :: ifield, ndumps, iline, ip, ll, l, nfldt, jp, ifld, nfiles

!   Default is user specified fortran function (=0 if not specified)

  nfldt = nfield
  if (ifmhd) nfldt = nfield+1

  do jp=0,npert
      ifrest(0,jp) = .FALSE. 
      do ifld=1,nfldt
          iffort(ifld,jp) = .TRUE. 
          ifrest(ifld,jp) = .FALSE. 
          ifprsl(ifld,jp) = .FALSE. 
      enddo
  enddo

  jp = 0
  nfiles=0

!   Check for Presolve options

  DO 1000 ILINE=1,15
      LINE=INITC(ILINE)
      CALL CAPIT(LINE,132)
      IF (INDX1(LINE,'PRESOLV',7) /= 0) THEN
      !           found a presolve request
          CALL BLANK(INITC(ILINE),132)
          CALL LJUST(LINE)
          CALL CSPLIT(CDUM,LINE,' ',1)
      
          IF (LTRUNC(LINE,132) == 0) THEN
              IF (NID == 0) WRITE(6,700)
              700 FORMAT(/,2X,'Presolve options: ALL')
          !              default - all fields are presolved.
              DO 800 IFIELD=1,nfldt
                  ifprsl(ifield,jp) = .TRUE. 
                  iffort(ifield,jp) = .FALSE. 
              800 END DO
          ELSE
          !           check line for arguments
          
              LL=LTRUNC(LINE,132)
              IF (NID == 0) WRITE(6,810) (LINE1(L),L=1,LL)
              810 FORMAT(/,2X,'Presolve options: ',132A1)
          
              IF (INDX_CUT(LINE,'U',1) /= 0) THEN
                  ifprsl(1,jp) = .TRUE. 
                  iffort(1,jp) = .FALSE. 
              ENDIF
          
              IF (INDX_CUT(LINE,'T',1) /= 0) THEN
                  ifprsl(2,jp) = .TRUE. 
                  iffort(2,jp) = .FALSE. 
              ENDIF
          
              DO 900 IFIELD=3,NPSCAL+2
                  IP=IFIELD-2
                  WRITE(S2,901) IP
                  IF (INDX_CUT(LINE,S2,2) /= 0) THEN
                      ifprsl(ifield,jp) = .TRUE. 
                      iffort(ifield,jp) = .FALSE. 
                  ENDIF
              900 END DO
              901 FORMAT('P',I1)
          ENDIF
      ENDIF
  1000 END DO

!   Check for restart options

  jp = 0
  DO 2000 ILINE=1,15
      if (ifpert) jp=iline-1
      LINE=INITC(ILINE)
      IF (LTRUNC(LINE,132) /= 0) THEN
      !           found a filename
          NFILES=NFILES+1
          INITC(NFILES)=LINE
      
          IF (NID == 0 .AND. NFILES == 1) WRITE(6,1010) LINE
          1010 FORMAT(1X,'Checking restart options: ',A132)
      !            IF (NID.EQ.0) WRITE(6,'(A132)') LINE
      
      !           Parse restart options
           
          call sioflag(ndumps,fname,line)

          if (ifgetx) then
              ifrest(0,jp) = .TRUE. 
          endif
          if (ifgetu) then
              iffort(1,jp) = .FALSE. 
              ifprsl(1,jp) = .FALSE. 
              ifrest(1,jp) = .TRUE. 
          endif
          if (ifgett) then
              iffort(2,jp) = .FALSE. 
              ifprsl(2,jp) = .FALSE. 
              ifrest(2,jp) = .TRUE. 
          endif
          do 1900 ifield=3,nfldt
          !              write(6,*) 'ifgetps:',(ifgtps(k),k=1,ldimt-1)
              if (ifgtps(ifield-2)) then
                  iffort(ifield,jp) = .FALSE. 
                  ifprsl(ifield,jp) = .FALSE. 
                  ifrest(ifield,jp) = .TRUE. 
              endif
          1900 END DO
      endif
  2000 END DO

  return
end subroutine slogic

!-----------------------------------------------------------------------
!> \brief driver for restarts
!!  (1) Open restart file(s)
!!  (2) Check previous spatial discretization
!!  (3) Map (K1,N1) => (K2,N2) if necessary
!!  nfiles > 1 has several implications:
!!  i.   For std. run, data is taken from last file in list, unless
!!       explicitly specified in argument list of filename
!!  ii.  For MHD and perturbation cases, 1st file is for U,P,T;
!!       subsequent files are for B-field or perturbation fields
!----------------------------------------------------------------------
subroutine restart_driver(nfiles)
  use kinds, only : DP, r4
  use size_m, only : lx1, ly1, lz1, lelt, ldimt, ldimt1
  use size_m, only : nid, nx1, ny1, nz1, nx2, ny2, nz2
  use restart, only : nxr, nyr, nzr
  use restart, only : ifgetx, ifgetz, ifgetu, ifgetw, ifgetp, ifgett, ifgtim
  use restart, only : ifgtps
  use geom, only : xm1, ym1, zm1
  use input, only : ifpert, ifmhd, if3d, npscal, param, initc
  use parallel, only : nelgt, isize, gllnid
  use soln, only : vx, vy, vz, t
  use string, only : ltrunc, i1_from_char
  use tstep, only : time
  implicit none

  integer, intent(in) :: nfiles

  integer :: nelrr

  integer, parameter :: LXR=LX1+6
  integer, parameter :: LYR=LY1+6
  integer, parameter :: LZR=LZ1+6
  integer, parameter :: LXYZR=LXR*LYR*LZR
  integer, parameter :: LXYZT=LX1*LY1*LZ1*LELT
  integer, parameter :: LPSC9=LDIMT+9

  real(DP), allocatable :: sdump(:,:)
  integer :: mesg(40)

  real(r4), allocatable :: tdump(:,:)

  REAL(DP), allocatable :: SDMP2(:,:)

!   cdump comes in via PARALLEL (->TOTAL)

  character(30) :: excoder
  character(1) ::  excoder1(30)
  equivalence (excoder,excoder1)

  character(132) :: fname
  character(1) ::  fname1(132)
  equivalence (fname1,fname)

  integer ::       hnami (30)
  character(132) :: hname

  CHARACTER(132) :: header

!   Local logical flags to determine whether to copy data or not.
  logical :: ifok,iffmat
  integer :: iposx,iposz,iposu,iposw,iposp,ipost,ipsps(ldimt1)

  logical :: ifbytsw, if_byte_swap_test
  real(r4) ::   bytetest

  real(DP) :: p67, rstime, cdump
  integer :: ifile, ndumps, ierr, len, idump, neltr, istepr, i, icase, ipass
  integer :: nxyz2, mid, ieg, nerr, ips, ii, nxyzr, ixyzz, nouts
  integer :: jxyz, ntotv, ntott
  integer :: nxyz1, lname, is, nps0, nps1, i1, iposv, iposy, nps

  allocate(TDUMP(LXYZR,LPSC9))

  ifok= .FALSE. 
  ifbytsw = .FALSE. 

  if(nfiles < 1) return
  if(nid == 0) write(6,*) 'Reading checkpoint data'

! use new reader (only binary support)
  p67 = abs(param(67))
  if (p67 == 6.0) then
    write(*,*) "Oops: p67"
#if 0
      do ifile=1,nfiles
          call sioflag(ndumps,fname,initc(ifile))
          call mfi(fname,ifile)
      enddo
      call setup_convect(3)
      if (nid /= 0) time=0
      time = glmax(time,1) ! Sync time across processors
#endif
      return
  endif

! use old reader (for ASCII + old binary support)
          
  if (param(67) < 1.0) then  ! zero only. should be abs.
      iffmat= .TRUE.  ! ascii
  else
      iffmat= .FALSE. ! binary
  endif

  do 6000 ifile=1,nfiles
      call sioflag(ndumps,fname,initc(ifile))
      ierr = 0
      if (nid == 0) then

          if (iffmat) then
              open (unit=91,file=fname,status='old',err=500)
          else
              len= ltrunc(fname,79)
          !           test for presence of file
              open (unit=91,file=hname &
              ,form='unformatted',status='old',err=500)
              close(unit=91)
              call byte_open(hname,ierr)
              if(ierr /= 0) goto 500
          ENDIF
          ifok = .TRUE. 
      endif

      500 continue
      call lbcast(ifok)
      if ( .NOT. ifok) goto 5000
               
      ndumps = 1
  
  !        Only NODE 0 reads from the disk.
  
      DO 1000 IDUMP=1,NDUMPS

          IF (NID == 0) THEN
          ! read header
              if (iffmat) then
                  ierr = 2
                  if(mod(param(67),1.0) == 0) then ! old header format
                      if(nelgt < 10000) then
                          read(91,91,err=10,end=10) &
                          neltr,nxr,nyr,nzr,rstime,istepr,(excoder1(i),i=1,30)
                          91 format(4i4,1x,g13.4,i5,1x,30a1)
                          ierr=0
                      else
                          read(91,92,err=10,end=10) &
                          neltr,nxr,nyr,nzr,rstime,istepr,(excoder1(i),i=1,30)
                          92 format(i10,3i4,1P1e18.9,i9,1x,30a1)
                          ierr=0
                      endif
                  else                          ! new head format
                      read(91,'(A132)',err=10,end=10) header
                      read(header,*) &
                      neltr,nxr,nyr,nzr,rstime,istepr,excoder
                      ierr=0
                  endif
              else
                  if(mod(param(67),1.0) == 0) then  ! old header format
                      call byte_read(hnami,20,ierr)
                      if(ierr /= 0) goto 10
                      icase = 2
                      if (nelgt < 10000) icase = 1
                      ipass = 0
                      93 continue  ! test each possible case  UGLY (7/31/07)
                      if(ipass < 2) then
                          ipass = ipass+1
                          if(icase == 1) then
                              read(hname,'(4i4,1x,g13.4,i5,1x,30a1)', &
                              err=94,end=94) &
                              neltr,nxr,nyr,nzr,rstime,istepr, &
                              (excoder1(i),i=1,30)
                              goto 95
                          else
                              read(hname,'(i10,3i4,1P1e18.9,i9,1x,30a1)', &
                              err=94,end=94) &
                              neltr,nxr,nyr,nzr,rstime,istepr, &
                              (excoder1(i),i=1,30)
                              goto 95
                          endif
                          94 icase = 3-icase  !  toggle: 2-->1  1-->2
                          goto 93
                      else
                          ierr=2
                          goto 10
                      endif
                      95 continue
                  else                         ! new head format
                      call byte_read(header,20,ierr)
                      if(ierr /= 0) goto 10
                      read(header,*) &
                      neltr,nxr,nyr,nzr,rstime,istepr,excoder
                  endif
                  call byte_read(bytetest,1,ierr)
              !                call byte_read2(bytetest,1,ierr)
                  if(ierr /= 0) goto 10
                  ifbytsw = if_byte_swap_test(bytetest,ierr)
                  if(ierr /= 0) goto 10
              endif
              mesg(1) = neltr
              mesg(2) = nxr
              mesg(3) = nyr
              mesg(4) = nzr
              write(6,*)  'Read mode: ', param(67)
              write(6,333)'neltr,nxr,nyr,nzr: ', neltr,nxr,nyr,nzr
              333 format(A,i9,3i4)
              call chcopy(mesg(5),excoder1,20)
              len  = 14*isize
          endif
          10 call err_chk(ierr,'Error reading restart header. Abort.$')

          IF (IDUMP == 1) THEN
              len  = 14*isize
              call bcast(mesg,len)
              neltr = mesg(1)
              nxr   = mesg(2)
              nyr   = mesg(3)
              nzr   = mesg(4)
              call   chcopy(excoder1,mesg(5),20)

              call lbcast(ifbytsw)

          !              Bounds checking on mapped data.
              IF (NXR > LXR) THEN
                  WRITE(6,20) NXR,NX1
                  20 FORMAT(//,2X, &
                  'ABORT:  Attempt to map from',I3, &
                  ' to N=',I3,'.',/,2X, &
                  'NEK5000 currently supports mapping from N+6 or less.' &
                  ,/,2X,'Increase N or LXR in IC.FOR.')
                  CALL EXITT
              ENDIF
          
          !              Figure out position of data in file "IFILE"
          
              NOUTS=0
              IPOSX=0
              IPOSY=0
              IPOSZ=0
              IPOSU=0
              IPOSV=0
              IPOSW=0
              IPOSP=0
              IPOST=0
              DO 40 I=1,NPSCAL
                  IPSPS(I)=0
              40 END DO

              IPS = 0
              NPS = 0
              DO 50 I=1, 30
                  IF (excoder1(i) == 'X') THEN
                      NOUTS=NOUTS + 1
                      IPOSX=NOUTS
                      NOUTS=NOUTS+1
                      IPOSY=NOUTS
                      IF (IF3D) THEN
                          NOUTS=NOUTS + 1
                          IPOSZ=NOUTS
                      ENDIF
                  ENDIF
                  IF (excoder1(i) == 'U') THEN
                      NOUTS=NOUTS + 1
                      IPOSU=NOUTS
                      NOUTS=NOUTS+1
                      IPOSV=NOUTS
                      IF (IF3D) THEN
                          NOUTS=NOUTS + 1
                          IPOSW=NOUTS
                      ENDIF
                  ENDIF
                  IF (excoder1(i) == 'P') THEN
                      NOUTS=NOUTS + 1
                      IPOSP=NOUTS
                  ENDIF
                  IF (excoder1(i) == 'T') THEN
                      NOUTS=NOUTS + 1
                      IPOST=NOUTS
                  ENDIF
                  IF(mod(param(67),1.0) == 0.0) THEN
                      i1 = i1_from_char(excoder1(i))
                      if (0 < i1 .AND. i1 < 10) then
                          if (i1 <= ldimt1) then
                              nouts=nouts + 1
                              ipsps(i1)=nouts
                          else
                              if (nid == 0) write(6,2) i1,i,excoder1(i)
                              2 format(2i4,a1,' PROBLEM W/ RESTART DATA')
                          endif
                      endif
                  ELSE
                      IF(excoder1(i) == 'S') THEN
                          READ(excoder1(i+1),'(I1)') NPS1
                          READ(excoder1(i+2),'(I1)') NPS0
                          NPS=10*NPS1 + NPS0
                          DO IS = 1, NPS
                              NOUTS=NOUTS + 1
                              IPSPS(IS)=NOUTS
                          ENDDO
                          GOTO 50
                      ENDIF
                  ENDIF
              50 END DO
               
              IF (NPS > (LDIMT-1)) THEN
                  IF (NID == 0) THEN
                      WRITE(*,'(A)') &
                      'ERROR: restart file has a NSPCAL > LDIMT'
                      WRITE(*,'(A,I2)') &
                      'Change LDIMT in SIZE'
                  ENDIF
                  CALL EXITT
              ENDIF

              lname=ltrunc(fname,132)
              if (nid == 0) write(6,61) (fname1(i),i=1,lname)
              if (nid == 0) write(6,62) &
              iposu,iposv,iposw,iposp,ipost,nps,nouts
              61 FORMAT(/,2X,'Restarting from file ',132A1)
              62 FORMAT(2X,'Columns for restart data U,V,W,P,T,S,N: ',7I4)

          !              Make sure the requested data is present in this file....
              if (iposx == 0) ifgetx= .FALSE. 
              if (iposy == 0) ifgetx= .FALSE. 
              if (iposz == 0) ifgetz= .FALSE. 
              if (iposu == 0) ifgetu= .FALSE. 
              if (iposv == 0) ifgetu= .FALSE. 
              if (iposw == 0) ifgetw= .FALSE. 
              if (iposp == 0) ifgetp= .FALSE. 
              if (ipost == 0) ifgett= .FALSE. 
              do 65 i=1,npscal
                  if (ipsps(i) == 0) ifgtps(i)= .FALSE. 
              65 END DO

          !              End of restart file header evaluation.

          ENDIF
      
      !           Read the error estimators
      !           not supported at the moment => just do dummy reading
      
          ifok = .FALSE. 
          IF(NID == 0)THEN
              if (iffmat) &
              READ(91,'(6G11.4)',END=15)(CDUMP,I=1,NELTR)
              ifok = .TRUE. 
          ENDIF
          15 continue
          call lbcast(ifok)
          if( .NOT. ifok) goto 1600
      
      !           Read the current dump, double buffer so that we can
      !           fit the data on a distributed memory machine,
      !           and so we won't have to read the restart file twice
      !           in case of an incomplete data file.
      
          NXYZR = NXR*NYR*NZR
      
      !           Read the data
      
          nelrr = min(nelgt,neltr) ! # of elements to _really_read
      ! why not just neltr?
          do 200 ieg=1,nelrr
              ifok = .FALSE. 
              IF (NID == 0) THEN
                  IF (MOD(IEG,10000) == 1) WRITE(6,*) 'Reading',IEG
                  IF (iffmat) THEN
                      READ(91,*,ERR=70,END=70) &
                      ((tdump(IXYZZ,II),II=1,NOUTS),IXYZZ=1,NXYZR)
                  ELSE
                      do ii=1,nouts
                          call byte_read(tdump(1,II),nxyzr,ierr)
                          if(ierr /= 0) then
                              write(6,*) "Error reading xyz restart data"
                              goto 70
                          endif
                      enddo
                  ENDIF
                  IFOK= .TRUE. 
              ENDIF
          
          !              Notify other processors that we've read the data OK.
          
              70 continue
              call lbcast(ifok)
              IF ( .NOT. IFOK) GOTO 1600
                write(*,*) "Oops: ifok" 
#if 0
          !              MAPDMP maps data from NXR to NX1
          !              (and sends data to the appropriate processor.)
          
          !              The buffer SDUMP is used so that if an incomplete dump
          !              file is found (e.g. due to UNIX io buffering!), then
          !              the previous read data stored in VX,VY,.., is not corrupted.
          
              IF (IFGETX) CALL MAPDMP &
              (SDUMP(1,1),TDUMP(1,IPOSX),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETX) CALL MAPDMP &
              (SDUMP(1,2),TDUMP(1,IPOSY),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETZ) CALL MAPDMP &
              (SDUMP(1,3),TDUMP(1,IPOSZ),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETU) CALL MAPDMP &
              (SDUMP(1,4),TDUMP(1,IPOSU),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETU) CALL MAPDMP &
              (SDUMP(1,5),TDUMP(1,IPOSV),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETW) CALL MAPDMP &
              (SDUMP(1,6),TDUMP(1,IPOSW),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETP) CALL MAPDMP &
              (SDUMP(1,7),TDUMP(1,IPOSP),IEG,NXR,NYR,NZR,ifbytsw)
              IF (IFGETT) CALL MAPDMP &
              (SDMP2(1,1),TDUMP(1,IPOST),IEG,NXR,NYR,NZR,ifbytsw)

          !              passive scalars
              do 100 ips=1,npscal
                  if (ifgtps(ips)) call mapdmp(sdmp2(1,ips+1) &
                  ,tdump(1,ipsps(ips)),ieg,nxr,nyr,nzr,IFBYTSW)
              100 END DO
#endif               
          200 END DO
      
      !           Successfully read a complete field, store it.
      
          nerr = 0              ! Count number of elements rec'd by nid
          do ieg=1,nelrr
              mid = gllnid(ieg)
              if (mid == nid) nerr = nerr+1
          enddo

          nxyz2=nx2*ny2*nz2
          nxyz1=nx1*ny1*nz1
          ntott=nerr*nxyz1
          ntotv=nerr*nxyz1   ! Problem for differing Vel. and Temp. counts!
      ! for now we read nelt dataset

          if (ifmhd .AND. ifile == 2) then
            write(*,*) "Oops: ifmhd"
#if 0
              if (ifgetu) call copy(bx,sdump(1,4),ntott)
              if (ifgetu) call copy(by,sdump(1,5),ntott)
              if (ifgetw) call copy(bz,sdump(1,6),ntott)
              if (ifgetp) then
                  if (nid == 0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                      call copy( pm,sdump(1,7),ntotv)
                  else
                      do iel=1,nelv
                          iiel = (iel-1)*nxyz1+1
                          call map12 (pm(1,1,1,iel),sdump(iiel,7),iel)
                      enddo
                  endif
              endif
              if (ifaxis .AND. ifgett) &
              call copy(t(1,1,1,1,2),sdmp2(1,1),ntott)
#endif
          elseif (ifpert .AND. ifile >= 2) then
            write(*,*) "Oops: ifpert"
#if 0
              j=ifile-1  ! pointer to perturbation field
              if (ifgetu) call copy(vxp(1,j),sdump(1,4),ntotv)
              if (ifgetu) call copy(vyp(1,j),sdump(1,5),ntotv)
              if (ifgetw) call copy(vzp(1,j),sdump(1,6),ntotv)
              if (ifgetp) then
                  if (nid == 0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                      call copy(prp(1,j),sdump(1,7),ntotv)
                  else
                      do ie=1,nelv
                          ie1 = (ie-1)*nxyz1+1
                          ie2 = (ie-1)*nxyz2+1
                          call map12 (prp(ie2,j),sdump(ie1,7),ie)
                      enddo
                  endif
              endif
              if (ifgett) call copy(tp(1,1,j),sdmp2(1,1),ntott)
          !              passive scalars
              do ips=1,NPSCAL
                  if (ifgtps(ips)) &
                  call copy(tp(1,ips+1,j),sdmp2(1,ips+1),ntott)
              enddo
#endif
          else  ! Std. Case
              allocate(sdump(LXYZT,7))
              if (ifgetx) call copy(xm1,sdump(1,1),ntott)
              if (ifgetx) call copy(ym1,sdump(1,2),ntott)
              if (ifgetz) call copy(zm1,sdump(1,3),ntott)
              if (ifgetu) call copy(vx ,sdump(1,4),ntotv)
              if (ifgetu) call copy(vy ,sdump(1,5),ntotv)
              if (ifgetw) call copy(vz ,sdump(1,6),ntotv)
!max              if (ifgetp) call copy(pm1,sdump(1,7),ntotv)
              if (ifgett) call copy(t,sdmp2(1,1),ntott)
          !              passive scalars
              allocate(sdmp2(lxyzt,ldimt))
              do i=1,NPSCAL
                  if (ifgtps(i)) &
                  call copy(t(1,1,1,1,i+1),sdmp2(1,i+1),ntott)
              enddo

!max              if (ifaxis) call axis_interp_ic(pm1)      ! Interpolate to axi mesh
!max              if (ifgetp) call map_pm1_to_pr(pm1,ifile) ! Interpolate pressure

              if (ifgtim) time=rstime
          endif

      1000 END DO
      GOTO 1600
       
      1600 CONTINUE
  
      IF (IDUMP == 1 .AND. NID == 0) THEN
          write(6,1700) fname
          write(6,1701) ieg,ixyzz
          write(6,1702) &
          ((tdump(jxyz,ii),ii=1,nouts),jxyz=ixyzz-1,ixyzz)
          1700 FORMAT(5X,'WARNING:  No data read in for file ',A132)
          1701 FORMAT(5X,'Failed on  element',I4,',  point',I5,'.')
          1702 FORMAT(5X,'Last read dump:',/,5G15.7)
          write(6,*) nid,'call exitt 1702a',idump
          call exitt
      ELSEIF (IDUMP == 1) THEN
          write(6,*) nid,'call exitt 1702b',idump
          call exitt
      ELSE
          IDUMP=IDUMP-1
          IF (NID == 0) WRITE(6,1800) IDUMP
          1800 FORMAT(2X,'Successfully read data from dump number',I3,'.')
      ENDIF
      if (iffmat) then
          if (nid == 0) close(unit=91)
      else
          ierr = 0
          if (nid == 0) call byte_close(ierr)
          call err_chk(ierr,'Error closing restart file in restart$')
      endif
      GOTO 6000

  !        Can't open file...
      5000 CONTINUE
      if (nid == 0) write(6,5001) fname
      5001 FORMAT(2X,'   *******   ERROR   *******    ' &
      ,/,2X,'   *******   ERROR   *******    ' &
      ,/,2X,'   Could not open restart file:' &
      ,/,A132 &
      ,//,2X,'Quitting in routine RESTART.')
      CLOSE(UNIT=91)
      call exitt
      IF (NID == 0) WRITE(6,5001) HNAME
      call exitt
  
  
  !     End of IFILE loop
  6000 END DO

  return
end subroutine restart_driver

!-----------------------------------------------------------------------
!> \brief Set IO flags according to Restart Options File, RSOPTS
subroutine sioflag(ndumps,fname,rsopts)
  use kinds, only : DP
  use size_m, only : ldimt, nfield
  use input, only : if3d, ifadvc, ifflow, ifheat
  use restart, only : ifgetx, ifgetu, ifgett, ifgetp, ifgetz, ifgetw
  use restart, only : ifgtps, ifgtim
  use string, only : indx1, ltrunc, indx_cut, ifgtrl, ljust, capit, csplit
  use tstep, only : time
  implicit none

  character(132) :: rsopts,fname
  character(2) ::  s2

!   Scratch variables..
  logical :: ifdeft,ifanyc
  CHARACTER(132) :: RSOPT     ,LINE
  CHARACTER(1) ::  RSOPT1(132),LINE1(132)
  EQUIVALENCE (RSOPT1,RSOPT)
  EQUIVALENCE (LINE1,LINE)

  integer :: len, len1, len4
  integer :: i, ndumps, ito, it1, it8, ita, itb, ixo, ivo, ipo
  real(DP) :: ttime, tdumps

!   Parse filename

!      CSPLIT splits S1 into two parts, delimited by S2.
!      S1 returns with 2nd part of S1.  CSPLIT returns 1st part.

  rsopt=rsopts
  call ljust(rsopt)
  call csplit(fname,rsopt,' ',1)
!   check fname for user supplied extension.
  if (indx1(fname,'.',1) == 0) then
      len=ltrunc(fname,132)
      len1=len+1
      len4=len+4
      fname(len1:len4)='.fld'
  endif

!   Parse restart options

!   set default flags

  ifgetx= .FALSE. 
  ifgetz= .FALSE. 
  ifgetu= .FALSE. 
  ifgetw= .FALSE. 
  ifgetp= .FALSE. 
  ifgett= .FALSE. 
  do 100 i=1,ldimt-1
      ifgtps(i)= .FALSE. 
  100 END DO
  ifgtim= .TRUE. 
  ndumps=0

!   Check for default case - just a filename given, no i/o options specified

  ifdeft= .TRUE. 

!   Parse file for i/o options and/or dump number

  CALL CAPIT(RSOPT,132)

  IF (LTRUNC(RSOPT,132) /= 0) THEN
  
  !        Check for explicit specification of restart TIME.
  
      ITO=INDX_CUT(RSOPT,'TIME',4)
      IFGTIM= .TRUE. 
      IF (ITO /= 0) THEN
      !           user has specified the time explicitly.
          IT1=INDX_CUT(RSOPT,'=',1)
          IT8=132-IT1
          CALL BLANK(LINE,132)
          CALL CHCOPY(LINE,RSOPT1(IT1),IT8)
          IF (IFGTRL(TTIME,LINE)) THEN
              IFGTIM= .FALSE. 
              TIME=TTIME
          ENDIF
      !           remove the user specified time from the RS options line.
          ITA=132-ITO+1
          CALL BLANK(RSOPT1(ITO),ITA)
          CALL LJUST(LINE)
          IT1=INDX1(LINE,' ',1)
          ITB=132-IT1+1
          CALL CHCOPY(RSOPT1(ITO),LINE1(IT1),ITB)
      ENDIF

  !        Parse field specifications.

      IXO=INDX_CUT(RSOPT,'X',1)
      IF (IXO /= 0) THEN
          ifdeft= .FALSE. 
          IFGETX= .TRUE. 
          IF (IF3D) IFGETZ= .TRUE. 
      ENDIF

      IVO=INDX_CUT(RSOPT,'U',1)
      IF (IVO /= 0) THEN
          ifdeft= .FALSE. 
          IFGETU= .TRUE. 
          IF (IF3D) IFGETW= .TRUE. 
      ENDIF

      IPO=INDX_CUT(RSOPT,'P',1)
      IF (IPO /= 0) THEN
          ifdeft= .FALSE. 
          IFGETP= .TRUE. 
      ENDIF

      ITO=INDX_CUT(RSOPT,'T',1)
      IF (ITO /= 0) THEN
          ifdeft= .FALSE. 
          IFGETT= .TRUE. 
      ENDIF

      do 300 i=1,ldimt-1
          write (s2,301) i
          ipo=indx_cut(rsopt,s2,2)
          if (ipo /= 0) then
              ifdeft= .FALSE. 
              ifgtps(i)= .TRUE. 
          endif
      300 END DO
      301 format('P',i1)

  !        Get number of dumps from remainder of user supplied line.
      if (ifgtrl(tdumps,rsopt)) ndumps=int(tdumps)
  endif

!   If no fields were explicitly specified, assume getting all fields.
  if (ifdeft) then
      IFGETX= .TRUE. 
      IF (IF3D) IFGETZ= .TRUE. 
      IFANYC= .FALSE. 
      DO 400 I=1,NFIELD
          IF (IFADVC(I)) IFANYC= .TRUE. 
      400 END DO
      IF (IFFLOW .OR. IFANYC) THEN
          IFGETU= .TRUE. 
          IF (IF3D) IFGETW= .TRUE. 
      ENDIF
      IF (IFFLOW) IFGETP= .TRUE. 
      IF (IFHEAT) IFGETT= .TRUE. 
      do 410 i=1,ldimt-1
          ifgtps(i)= .TRUE. 
      410 END DO
  ENDIF

  return
end subroutine sioflag

!----------------------------------------------------------------------
subroutine mapdmp(sdump,tdump,ieg,nxr,nyr,nzr,if_byte_sw)
  use kinds, only : DP, r4
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, nid, lelt
  use parallel, only : np, gllnid, nullpid, gllel
  implicit none

  integer, parameter :: LXYZ1=LX1*LY1*LZ1
  integer, PARAMETER :: LXR=LX1+6
  integer, PARAMETER :: LYR=LY1+6
  integer, PARAMETER :: LZR=LZ1+6
  integer, PARAMETER :: LXYZR=LXR*LYR*LZR

  REAL(DP), allocatable :: SDUMP(:,:)
  REAL(r4)   :: TDUMP(LXYZR)
  integer :: ieg, nxr, nyr, nzr
  logical :: if_byte_sw

  integer :: nxyz, nxyr, ierr, jnid, mtype, len, le1, ie
  real(DP) :: dummy


  NXYZ=NX1*NY1*NZ1
  NXYR=NXR*NYR*NZR
  ierr=0

!   Serial processor code:

  IF (NP == 1) THEN
      allocate(sdump(lxyz1,lelt))

      IF (if_byte_sw) call byte_reverse(TDUMP,NXYR,ierr)
      if(ierr /= 0) call exitti("Error in mapdmp")
      IF (NXR == NX1 .AND. NYR == NY1 .AND. NZR == NZ1) THEN
          CALL COPY4r(SDUMP(1,IEG),TDUMP,NXYZ)
      ELSE
      !           do the map    (assumes that NX=NY=NZ, or NX=NY, NZ=1)
          call mapab4r(sdump(1,ieg),tdump,nxr,1)
      ENDIF

  ELSE
  
  !     Parallel code - send data to appropriate processor and map.
  
      JNID=GLLNID(IEG)
      MTYPE=3333+IEG
      LEN=4*NXYR
      LE1=4
      IF (NID == 0 .AND. JNID /= 0) THEN
      !           hand-shake
          CALL CSEND(MTYPE,TDUMP,LE1,JNID,NULLPID)
          CALL CRECV(MTYPE,dummy,LE1)
          CALL CSEND(MTYPE,TDUMP,LEN,JNID,NULLPID)
      ELSEIF (NID /= 0 .AND. JNID == NID) THEN
      !           Receive data from node 0
          CALL CRECV(MTYPE,dummy,LE1)
          CALL CSEND(MTYPE,TDUMP,LE1,0,NULLPID)
          CALL CRECV(MTYPE,TDUMP,LEN)
      ENDIF
  
  !        If the data is targeted for this processor, then map
  !        to appropriate element.
  
      IF (JNID == NID) THEN
          allocate(sdump(lxyz1,lelt))
          IE=GLLEL(IEG)
          IF (if_byte_sw) call byte_reverse(TDUMP,NXYR,ierr)
          IF (NXR == NX1 .AND. NYR == NY1 .AND. NZR == NZ1) THEN
              CALL COPY4r(SDUMP(1,IE),TDUMP,NXYZ)
          ELSE
              call mapab4r(sdump(1,ie),tdump,nxr,1)
          ENDIF
      ENDIF
      call err_chk(ierr,'Error using byte_reverse in mapdmp.$')
  
  !        End of parallel distribution/map routine.
  
  ENDIF
  return
end subroutine mapdmp

!---------------------------------------------------------------
!> Interpolate Y(NXR,NYR,NZR,NEL) to X(NX1,NY1,NZ1,NEL)
!! (assumes that NXR=NYR=NZR, or NXR=NYR, NZR=1)
!! Input:  real*4,  Output:  default precision
!---------------------------------------------------------------
subroutine mapab4R(x,y,nxr,nel)
  use kinds, only : DP, r4
  use size_m, only : nid, ndim
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1
  use wz_m, only : zgm1
  implicit none

  integer :: nxr, nel
  REAL(DP) :: X(NX1,NY1,NZ1,NEL)
  REAL(r4) :: Y(NXR,NXR,NXR,NEL)

  integer, parameter :: LXR=LX1+6
  integer, parameter :: LYR=LY1+6
  integer, parameter :: LZR=LZ1+6
  integer, parameter :: LXYZR=LXR*LYR*LZR

  real(DP) ::   xa(lxyzr)      ,xb(lx1,ly1,lzr) ,xc(lxyzr) &
  , zgmr(lxr)      ,wgtr(lxr)
  real(DP), allocatable :: ires(:,:), itres(:,:)

  integer :: nzr, nyzr, nxy1, nxyzr
  integer :: ie, iz, izoff
  integer, save :: NOLD = 0
  
  allocate(ires(lxr,lxr), itres(lxr,lxr))

  NZR = NXR
  IF(NZ1 == 1) NZR=1
  NYZR = NXR*NZR
  NXY1 = NX1*NY1
  nxyzr = nxr*nxr*nzr

  IF (NXR /= NOLD) THEN
      NOLD=NXR
      CALL ZWGLL   (ZGMR,WGTR,NXR)
      CALL IGLLM   (IRES,ITRES,ZGMR,ZGM1,NXR,NX1,NXR,NX1)
      IF (NID == 0) WRITE(6,10) NXR,NX1
      10 FORMAT(2X,'Mapping restart data from Nold=',I2 &
      ,' to Nnew=',I2,'.')
  ENDIF

  DO 1000 IE=1,NEL
      call copy4r(xc,y(1,1,1,ie),nxyzr)
      CALL MXM (IRES,NX1,xc,NXR,XA,NYZR)
      DO 100 IZ=1,NZR
          IZOFF = 1 + (IZ-1)*NX1*NXR
          CALL MXM (XA(IZOFF),NX1,ITRES,NXR,XB(1,1,IZ),NY1)
      100 END DO
      IF (NDIM == 3) THEN
          CALL MXM (XB,NXY1,ITRES,NZR,X(1,1,1,IE),NZ1)
      ELSE
          CALL COPY(X(1,1,1,IE),XB,NXY1)
      ENDIF
  1000 END DO

  return
end subroutine mapab4R

!------------------------------------------------------------------
!> \brief  User specified fortran function (=0 if not specified)
!------------------------------------------------------------------
subroutine nekuic
  use size_m, only : nx1, ny1, nz1
  use input, only : ifmodel, ifkeps, ifldmhd
  use nekuse, only : turbk, turbe, ux, uy, uz, temp
  use parallel, only : lglel
  use soln, only : vx, vy, vz, bx, by, bz, t, vxp, vyp, vzp, tp, jp
  use tstep, only : ifield, nelfld
  use turbo, only : ifldk, iflde
  implicit none

  integer :: nel, ijke, i, j, k, iel, ieg

  NEL   = NELFLD(IFIELD)

  IF (IFMODEL .AND. IFKEPS .AND. IFIELD == IFLDK) THEN
  
      DO IEL=1,NEL
          ieg = lglel(iel)
          DO K=1,NZ1
              DO J=1,NY1
                  DO I=1,NX1
                      CALL NEKASGN (I,J,K,IEL)
                      CALL USERIC  (I,J,K,IEG)
                      T(I,J,K,IEL,IFIELD-1) = TURBK
                  enddo
              enddo
          enddo
      END DO
  
  ELSEIF (IFMODEL .AND. IFKEPS .AND. IFIELD == IFLDE) THEN
  
      DO IEL=1,NEL
          ieg = lglel(iel)
          DO K=1,NZ1
              DO J=1,NY1
                  DO I=1,NX1
                      CALL NEKASGN (I,J,K,IEL)
                      CALL USERIC  (I,J,K,IEG)
                      T(I,J,K,IEL,IFIELD-1) = TURBE
                  enddo
              enddo
          enddo
      END DO
  
  ELSE
  
      DO IEL=1,NEL
          ieg = lglel(iel)
          DO K=1,NZ1
              DO J=1,NY1
                  DO I=1,NX1
                      CALL NEKASGN (I,J,K,IEL)
                      CALL USERIC  (I,J,K,IEG)
                      if (jp == 0) then
                          IF (IFIELD == 1) THEN
                              VX(I,J,K,IEL) = UX
                              VY(I,J,K,IEL) = UY
                              VZ(I,J,K,IEL) = UZ
                          ELSEIF (IFIELD == ifldmhd) THEN
                              BX(I,J,K,IEL) = UX
                              BY(I,J,K,IEL) = UY
                              BZ(I,J,K,IEL) = UZ
                          ELSE
                              T(I,J,K,IEL,IFIELD-1) = TEMP
                          ENDIF
                      else
                          ijke = i+nx1*((j-1)+ny1*((k-1) + nz1*(iel-1)))
                          IF (IFIELD == 1) THEN
                              VXP(IJKE,jp) = UX
                              VYP(IJKE,jp) = UY
                              VZP(IJKE,jp) = UZ
                          ELSE
                              TP(IJKE,IFIELD-1,jp) = TEMP
                          ENDIF
                      endif
                  enddo
              enddo
          enddo
      END DO
  
  ENDIF

  return
end subroutine nekuic

!-----------------------------------------------------------------------
logical function if_byte_swap_test(bytetest,ierr)
  use kinds, only : DP, r4
  use size_m
  implicit none
   
  real(r4) :: bytetest
  integer :: ierr

  real(r4) :: test2
  real(r4), save :: test_pattern
  real(DP) :: eps, etest
   
  test_pattern = 6.54321_r4
  eps          = 0.00020
  etest        = abs(test_pattern-bytetest)
  if_byte_swap_test = .TRUE. 
  if (etest <= eps) if_byte_swap_test = .FALSE. 
   
  test2 = bytetest
  call byte_reverse(test2,1,ierr)
  if (nid == 0) &
  write(6,*) 'byte swap:',if_byte_swap_test,bytetest,test2
  return
end function if_byte_swap_test

!-----------------------------------------------------------------------
!> \brief Generate geometry data
subroutine geom_reset(icall)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, lx3
  use size_m, only : nx1, ny1, nz1, nelt, nid
  implicit none

  integer :: icall

  real(DP), allocatable :: XM3 (:,:,:,:), YM3(:,:,:,:), ZM3(:,:,:,:)

  integer :: ntot 

  if(nid == 0) write(6,*) 'regenerate geometry data',icall

  ntot = nx1*ny1*nz1*nelt

  if (lx3 == lx1) then
      CALL GEOM1 ()!XM1,YM1,ZM1)
  else
      allocate(xm3(lx1,ly1,lz1,lelt), ym3(lx1,ly1,lz1,lelt), zm3(lx1,ly1,lz1,lelt))
#if 0
      call map13_all(xm3,xm1)
      call map13_all(ym3,ym1)
      if (if3d) call map13_all(zm3,zm1)
      CALL GEOM1 (XM3,YM3,ZM3)
#endif
      deallocate(xm3,ym3,zm3)
  endif

  CALL GEOM2
  CALL UPDMSYS (1)
  CALL VOLUME
  CALL SETINVM
  CALL SETDEF
  CALL SFASTAX

  if(nid == 0) then
      write(6,*) 'done :: regenerate geometry data',icall
      write(6,*) ' '
  endif

  return
end subroutine geom_reset

!-----------------------------------------------------------------------
!> \brief Take direct stiffness avg of u
subroutine dsavg(u)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, nelt
  use size_m, only : lx1, ly1, lz1, lelt
  use input, only : ifflow
  use soln, ONLY : vmult, tmult
  use tstep, only : ifield

  implicit none
  real(DP) :: u(lx1,ly1,lz1,lelt)

  integer :: ifieldo, ntot


  ifieldo = ifield
  if (ifflow) then
      ifield = 1
      ntot = nx1*ny1*nz1*nelv
      call dssum(u)
      u = u * vmult
  else
      ifield = 2
      ntot = nx1*ny1*nz1*nelt
      call dssum(u)
      u = u * tmult(:,:,:,:,1)
  endif
  ifield = ifieldo

  return
end subroutine dsavg

!-----------------------------------------------------------------------
!> \brief open  blah000.fldnn
subroutine mbyte_open(hname,fid,ierr) 
  use size_m, only : nid
  use iso_c_binding, only : c_null_char
  use string, only : ltrunc, indx1
  use tstep, only : istep
  implicit none
   
  integer :: fid
  character(132) :: hname

  character(8) ::  fmt,s8
  character(8), save :: eight = "????????"

  character(132) :: fname

  integer :: len, ipass, k, i1, ierr

  len = ltrunc(hname,132)
  call chcopy (fname,hname,len)
  fname(len+1:len+1) = c_null_char

  do ipass=1,2      ! 2nd pass, in case 1 file/directory
      do k=8,1,-1
          i1 = indx1(fname,eight,k)
          if (i1 /= 0) then ! found k??? string
              write(fmt,1) k,k
              1 format('(i',i1,'.',i1,')')
              write(s8,fmt) fid
              call chcopy(fname(i1:),s8,k)
              goto 10
          endif
      enddo
      10 continue
  enddo
          
#ifdef MPIIO
  call byte_open_mpi(fname,ifh_mbyte,ierr)
  if(nid == 0) write(6,6) istep,(fname(k:k),k=1,len)
  6 format(1i8,' OPEN: ',132a1)
#else
  call byte_open(fname,ierr)
  write(6,6) nid,istep,(fname(k:k),k=1,len)
  6 format(2i8,' OPEN: ',132a1)
#endif

  return
end subroutine mbyte_open
!-----------------------------------------------------------------------
