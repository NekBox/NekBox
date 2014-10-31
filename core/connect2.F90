!-----------------------------------------------------------------------
!> \brief Read in data from preprocessor input file (.rea)
subroutine readat()
  use kinds, only : DP
  use size_m, only : nid, lelt, ndim, lx1, lx2, lx3
  use ctimer, only : dnekclock
  use input, only : reafle, ifmoab, h5mfle, numflu, numoth, ifheat
  use input, only : matids, matindx, imatie, ifchar, numsts
  use input, only : numbcs, ibcsts, bcf, bctyps
  use parallel, only : nelgv, nelgt, isize, np
  use string, only : cscan
  use zper, only : ifgtp
  implicit none
   
  logical :: ifre2
  character(132) :: tmp_string
  real(DP) :: etime_tmp
  integer :: idum(3*numsts+3)

  real(DP) :: edif, e1, e2, x
  integer :: i, iset, nelgs, maxrd, mread, iread

!   Test timer accuracy
  edif = 0.0
  do i = 1,10
      e1 = dnekclock()
      e2 = dnekclock()
      edif = edif + e2-e1
  enddo
  edif = edif/10.
  if(nid == 0) write(6,'(A,1pE15.7,A,/)') &
  ' timer accuracy: ', edif, ' sec'

  etime_tmp = dnekclock()

!   Open .rea file
  if(nid == 0) then
      write(6,*) 'read .rea file'
      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
  endif

!   Read parameters and logical flags
  CALL RDPARAM

!   Read Mesh Info
  if(nid == 0) then
      read(9,*)    ! xfac,yfac,xzero,yzero
      read(9,*)    ! dummy
      if (ifmoab) then
          read(9,*) h5mfle
      ! read fluid/solid material set ids
          read(9,*) numflu, numoth
          if (numflu+numoth > numsts) then
              write(6,'(A)') &
              'Number of fluid+other material sets too large.'
              write(6, '(A)') &
              'Need to increase NUMSTS in file INPUT.'
              call exitt
          else if (numoth > 0 .AND. .NOT. ifheat) then
              call exitt( &
              'Error: no. of other sets is non-zero but ifheat = false.')
          endif
          read(9,*) (matids(i), i = 1, numflu+numoth)
          do i = numflu+numoth+1, numsts
              matids(i) = -1
          enddo
          read(9,*) (matindx(i), i = 1, numflu+numoth)
          do i = numflu+numoth+1, numsts
              matindx(i) = -1
          enddo
          do i = 1, lelt
              imatie(i) = -1
          enddo
          read(9,*) numbcs
          if (numbcs > numsts) then
              write(6,'(A)') &
              'Number of BC sets too large.'
              write(6, '(A)') &
              'Need to increase NUMSTS in file INPUT.'
              call exitti
          endif
          do iset = 1, numbcs
              read(9,'(2I5,A3)') ibcsts(iset), bcf(iset), bctyps(iset)

          enddo
          nelgs = 0
          do iset = numbcs+1, numsts
              bcf(iset) = -1
              bctyps(iset) = 'E  '
              ibcsts(iset) = -1
          enddo
      else
          read(9,*)  nelgs,ndim,nelgv
          nelgt = abs(nelgs)
      endif
  endif
  call bcast(nelgs,ISIZE)
  call bcast(ndim ,ISIZE)
  call bcast(nelgv,ISIZE)
  call bcast(nelgt,ISIZE)
  call bcast(h5mfle,132)
  if (ifmoab) then
  ! pack into long int array and bcast as that
      if (nid == 0) then
          idum(1) = numflu
          idum(2) = numoth
          idum(3) = numbcs
          do iset = 1, numsts
              idum(3+iset) = matids(iset)
          enddo
          do iset = 1, numsts
              idum(3+numflu+numoth+iset) = ibcsts(iset)
          enddo
          do iset = 1, numsts
              idum(3+numflu+numoth+numbcs+iset) = matindx(iset)
          enddo
      endif
      call bcast(idum, ISIZE*(3+3*numsts))
      call bcast(bctyps, 3*numsts)
      call bcast(bcf, ISIZE*numsts)

      if (nid /= 0) then
          numflu = idum(1)
          numoth = idum(2)
          numbcs = idum(3)
          do iset = 1, numsts
              matids(iset) = idum(3+iset)
          enddo
          do iset = 1, numsts
              ibcsts(iset) = idum(3+numflu+numoth+iset)
          enddo
          do iset = 1, numsts
              matindx(iset) = idum(3+numflu+numoth+numbcs+iset)
          enddo
      endif
  endif
  ifre2 = .FALSE. 
  if(nelgs < 0) ifre2 = .TRUE.     ! use new .re2 reader

  ifgtp = .FALSE. 
  if (ndim < 0) ifgtp = .TRUE.     ! domain is a global tensor product

  if (ifmoab) then
#ifdef MOAB
      call nekMOAB_import
#endif
  else
!max        if (ifre2) call open_bin_file(ifbswap) ! rank0 will open and read
      if (nid == 0) then
          write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
          write(6,12) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
          12 format(1X,A,4I12,/,/)
      endif

      call chk_nel  ! make certain sufficient array sizes

      if ( .NOT. ifgtp) call mapelpr  ! read .map file, est. gllnid, etc.
      if (ifre2) then
        write(*,*) "Oops: ifre2"
!max            call bin_rd1(ifbswap) ! rank0 will read mesh data + distribute
      else

#if 1
          ! generate the mesh without reading
          call nekgsync()
          call genmesh
#else
          maxrd = 32            ! max # procs to read at once
          mread = (np-1)/maxrd+1   ! mod param
          iread = 0                ! mod param
          x     = 0
          do i=0,np-1,maxrd
              call nekgsync()
              if (mod(nid,mread) == iread) then
                  if (nid /= 0) then
                      open(UNIT=9,FILE=REAFLE,STATUS='OLD')
                      call cscan(tmp_string,'MESH DATA',9)
                      read(9,*) tmp_string
                  endif
                  if (ifgtp) then
!max                        call genbox
                  else
                      call rdmesh
                      call rdcurve !  Curved side data
                      call rdbdry  !  Boundary Conditions
                  endif
                  if (nid /= 0) close(unit=9)
              endif
              iread = iread + 1
          enddo
#endif
      endif
  endif

  if (nid == 0) then
    call cscan(tmp_string,'TAIL OPTS',9)
  endif

!   Read Restart options / Initial Conditions / Drive Force
  CALL RDICDF
!   Read materials property data
  CALL RDMATP
!   Read history data
  CALL RDHIST
!   Read output specs
  CALL RDOUT
!   Read objects
  CALL RDOBJ

  call nekgsync()

!   End of input data, close read file.
  IF(NID == 0) THEN
      CLOSE(UNIT=9)
      write(6,'(A,g13.5,A,/)')  ' done :: read .rea file ', &
      dnekclock()-etime_tmp,' sec'
  ENDIF

!   This is not an excellent place for this check, but will
!   suffice for now.   5/6/10
  if (ifchar .AND. (nelgv /= nelgt)) call exitti( &
  'ABORT: IFCHAR curr. not supported w/ conj. ht transfer$',nelgv)


  return
end subroutine readat

!-----------------------------------------------------------------------
!> \brief Read in parameters supplied by preprocessor and
!!   (eventually) echo check.
!!  .Broadcast run parameters to all processors
subroutine rdparam
  use size_m, only : nid, ndim, ldimt, ldimt1, npert
  use size_m, only : ldim, ly2, lx2, lz2, lgmres, lfdm, lbx1, lpx1  
  use size_m, only : lxd, lyd, lzd
  use size_m, only : lx1, lx1m, ly1, ly1m, lz1, lz1m
  use ctimer, only : ifsync
  use input, only : vnekton, param, npscal, cpfld
  use input, only : iftmsh, ifadvc, ifflow, ifheat, iftran, ifaxis, iflomach
  use input, only : ifusermv, ifcons, ifuservp, ifessr, ifmhd, ifvcoup
  use input, only : ifmgrid, ifpert, ifbase, ifschclob, ifexplvis, ifcyclic
  use input, only : ifmvbd, ifchar, if3d, ifsplit, ifldmhd
  use input, only : ifmoab, ifcoup, ifkeps, ifanls, ifmodel, ifstrs, ifaziv
  use parallel, only : isize, wdsize, csize
  use string, only : indx1, capit
  use zper, only : nelx, nely, nelz, ifzper, ifgfdm
  implicit none

  character(132) :: tmp_string(100)
  integer :: nparam, i, npscl1, npscl2, nskip, nlogic, ii, n_o, ktest

  VNEKTON = 3 ! dummy not really used anymore
        
  IF(NID == 0) THEN
      READ(9,*,ERR=400)
      READ(9,*,ERR=400)
      READ(9,*,ERR=400) NDIM
      READ(9,*,ERR=400) NPARAM
      DO 20 I=1,NPARAM
          READ(9,*,ERR=400) PARAM(I)
      20 END DO
  ENDIF
  call bcast(NDIM  ,ISIZE)
  call bcast(NPARAM,ISIZE)
  call bcast(PARAM ,200*WDSIZE)

  NPSCAL=INT(PARAM(23))
  NPSCL1=NPSCAL+1
  NPSCL2=NPSCAL+2

  IF (NPSCL1 > LDIMT) THEN
      if(nid == 0) then
          WRITE(6,21) LDIMT,NPSCL1
          21 FORMAT(//,2X,'Error: This NEKTON Solver has been compiled' &
          /,2X,'       for',I4,' passive scalars.  This run' &
          /,2X,'       requires that LDIMT be set to',I4,'.')
      endif
      call exitt
  ENDIF
     

!   Read in the passive scalar conduct and rhocp's:

!   fluid
!               .viscosity is PARAM(2)
!               .if it is negative, it indicates that Re has been input
!               .therefore, redefine PARAM(2) = -1.0/PARAM(2)

  if(param(2) < 0.0) param(2)  = -1.0/param(2)
  if(param(8) < 0.0) param(8)  = -1.0/param(8)
  if(param(29) < 0.0) param(29) = -1.0/param(29)

  CPFLD(1,1)=PARAM(2)
  CPFLD(1,2)=PARAM(1)
!   temperature
  CPFLD(2,1)=PARAM(8)
  CPFLD(2,2)=PARAM(7)
  CPFLD(2,3)=PARAM(9)

!   passive scalars

  IF(NID == 0) THEN
      READ(9,*,ERR=400) NSKIP
      IF (NSKIP > 0 .AND. NPSCAL > 0) THEN
          READ(9,*,ERR=400)(CPFLD(I,1),I=3,NPSCL2)
          IF(NPSCL2 < 9)READ(9,*)
          READ(9,*,ERR=400)(CPFLD(I,2),I=3,NPSCL2)
          IF(NPSCL2 < 9)READ(9,*)
          do i=3,npscl2
              if (cpfld(i,1) < 0) cpfld(i,1) = -1./cpfld(i,1)
              if (cpfld(i,2) < 0) cpfld(i,2) = -1./cpfld(i,2)
          enddo
      ELSE
          DO 25 I=1,NSKIP
              READ(9,*,ERR=500)
          25 END DO
      ENDIF
  ENDIF
  call bcast(cpfld,WDSIZE*LDIMT1*3)


!   Read logical equation type descriptors....

  IFTMSH(0) = .FALSE. 
  do i=1,NPSCL2
      IFTMSH(i) = .FALSE. 
      IFADVC(i) = .FALSE. 
  enddo
  IFFLOW    = .FALSE. 
  IFHEAT    = .FALSE. 
  IFTRAN    = .FALSE. 
  IFAXIS    = .FALSE. 
  IFAZIV    = .FALSE. 
  IFSTRS    = .FALSE. 
  IFLOMACH  = .FALSE. 
  IFMODEL   = .FALSE. 
  IFKEPS    = .FALSE. 
  IFMVBD    = .FALSE. 
  IFCHAR    = .FALSE. 
  IFANLS    = .FALSE. 
  IFMOAB    = .FALSE. 
  IFCOUP    = .FALSE. 
  IFVCOUP   = .FALSE. 
  IFMHD     = .FALSE. 
  IFESSR    = .FALSE. 
  IFTMSH(0) = .FALSE. 
  IFUSERVP  = .FALSE. 
  IFCONS    = .FALSE.    ! Use conservation form?
  IFUSERMV  = .FALSE. 
  IFCYCLIC  = .FALSE. 
  IFSYNC    = .FALSE. 
  IFEXPLVIS = .FALSE. 
  IFSCHCLOB = .FALSE. 
!   IFSPLIT   = .false.

  ifbase = .TRUE. 
  ifpert = .FALSE. 


  IF(NID == 0) READ(9,*,ERR=500) NLOGIC
  call bcast(NLOGIC,ISIZE)
  IF(NLOGIC > 100) THEN
      if(nid == 0) &
      write(6,*) 'ABORT: Too many logical switches', NLOGIC
      call exitt
  ENDIF

  if(nid == 0) READ(9,'(A132)',ERR=500) (tmp_string(i),i=1,NLOGIC)
  call bcast(tmp_string,100*132*CSIZE)

  do i = 1,NLOGIC
      call capit(tmp_string(i),132)
      if (indx1(tmp_string(i),'IFTMSH' ,6) > 0) then
          read(tmp_string(i),*,ERR=490) (IFTMSH(II),II=0,NPSCL2)
      elseif (indx1(tmp_string(i),'IFNAV'  ,5) > 0 .AND. &
          indx1(tmp_string(i),'IFADVC' ,6) > 0) then
          read(tmp_string(i),*,ERR=490) (IFADVC(II),II=1,NPSCL2)
      elseif (indx1(tmp_string(i),'IFADVC' ,6) > 0) then
          read(tmp_string(i),*,ERR=490) (IFADVC(II),II=1,NPSCL2)
      elseif (indx1(tmp_string(i),'IFFLOW' ,6) > 0) then
          read(tmp_string(i),*) IFFLOW
      elseif (indx1(tmp_string(i),'IFHEAT' ,6) > 0) then
          read(tmp_string(i),*) IFHEAT
      elseif (indx1(tmp_string(i),'IFTRAN' ,6) > 0) then
          read(tmp_string(i),*) IFTRAN
      elseif (indx1(tmp_string(i),'IFAXIS' ,6) > 0) then
          read(tmp_string(i),*) IFAXIS
      elseif (indx1(tmp_string(i),'IFAZIV' ,6) > 0) then
          read(tmp_string(i),*) IFAZIV
      elseif (indx1(tmp_string(i),'IFSTRS' ,6) > 0) then
          read(tmp_string(i),*) IFSTRS
      elseif (indx1(tmp_string(i),'IFLO'   ,4) > 0) then
          read(tmp_string(i),*) IFLOMACH
      elseif (indx1(tmp_string(i),'IFMGRID',7) > 0) then
      !             read(tmp_string(i),*) IFMGRID
      elseif (indx1(tmp_string(i),'IFKEPS' ,6) > 0) then
          read(tmp_string(i),*) IFKEPS
      elseif (indx1(tmp_string(i),'IFMODEL',7) > 0) then
          read(tmp_string(i),*) IFMODEL
      elseif (indx1(tmp_string(i),'IFMVBD' ,6) > 0) then
          read(tmp_string(i),*) IFMVBD
      elseif (indx1(tmp_string(i),'IFCHAR' ,6) > 0) then
          read(tmp_string(i),*) IFCHAR
      elseif (indx1(tmp_string(i),'IFANLS' ,6) > 0) then
          read(tmp_string(i),*) IFANLS
      elseif (indx1(tmp_string(i),'IFMOAB' ,6) > 0) then
          read(tmp_string(i),*) IFMOAB
      elseif (indx1(tmp_string(i),'IFCOUP' ,6) > 0) then
          read(tmp_string(i),*) IFCOUP
      elseif (indx1(tmp_string(i),'IFVCOUP' ,7) > 0) then
          read(tmp_string(i),*) IFVCOUP
      elseif (indx1(tmp_string(i),'IFMHD'  ,5) > 0) then
          read(tmp_string(i),*) IFMHD
      elseif (indx1(tmp_string(i),'IFCONS' ,6) > 0) then
          read(tmp_string(i),*) IFCONS
      elseif (indx1(tmp_string(i),'IFUSERVP',8) > 0) then
          read(tmp_string(i),*) IFUSERVP
      elseif (indx1(tmp_string(i),'IFUSERMV',8) > 0) then
          read(tmp_string(i),*) IFUSERMV
      elseif (indx1(tmp_string(i),'IFCYCLIC',8) > 0) then
          read(tmp_string(i),*) IFCYCLIC
      elseif (indx1(tmp_string(i),'IFPERT'  ,6) > 0) then
          read(tmp_string(i),*) IFPERT
      elseif (indx1(tmp_string(i),'IFBASE'  ,6) > 0) then
          read(tmp_string(i),*) IFBASE
      elseif (indx1(tmp_string(i),'IFSYNC'  ,6) > 0) then
          read(tmp_string(i),*) IFSYNC
      elseif (indx1(tmp_string(i),'IFEXPLVIS',9) > 0) then
          read(tmp_string(i),*) IFEXPLVIS
      elseif (indx1(tmp_string(i),'IFSCHCLOB',9) > 0) then
          read(tmp_string(i),*) IFSCHCLOB
      elseif (indx1(tmp_string(i),'IFSPLIT' ,7) > 0) then
      !              read(string,*) IFSPLIT
      else
          if(nid == 0) then
              write(6,'(1X,2A)') 'ABORT: Unknown logical flag', tmp_string
              write(6,'(30(A,/))') &
              ' Available logical flags:', &
              '   IFTMSH'   , &
              '   IFADVC'   , &
              '   IFFLOW'   , &
              '   IFHEAT'   , &
              '   IFTRAN'   , &
              '   IFAXIS'   , &
              '   IFCYCLIC' , &
              '   IFSTRS'   , &
              '   IFLOMACH' , &
              '   IFMGRID'  , &
              '   IFKEPS'   , &
              '   IFMVBD'   , &
              '   IFCHAR'   , &
              '   IFANLS'   , &
              '   IFUSERVP' , &
              '   IFUSERMV' , &
              '   IFSYNC'   , &
              '   IFCYCLIC' , &
              '   IFSPLIT'  , &
              '   IFEXPLVIS', &
              '   IFCONS'   , &
              '   IFMOAB'   , &
              '   IFCOUP'   , &
              '   IFVCOUP'
          endif
          call exitt
      endif
      490 continue
  enddo

  ifmgrid   = .FALSE. 
  if (ifsplit) ifmgrid   = .TRUE. 
  if (ifaxis ) ifmgrid   = .FALSE. 

  if (param(29) /= 0.) ifmhd  = .TRUE. 
  if (ifmhd)           ifessr = .TRUE. 
  if (ifmhd)           npscl1 = npscl1 + 1
  if (param(30) > 0)  ifuservp = .TRUE. 
  if (param(31) /= 0.) ifpert = .TRUE. 
  if (param(31) < 0.) ifbase = .FALSE.   ! don't time adv base flow
  npert = int(abs(param(31)))

  IF (NPSCL1 > LDIMT .AND. IFMHD) THEN
      if(nid == 0) then
          WRITE(6,22) LDIMT,NPSCL1
          22 FORMAT(/s,2X,'Error: This NEKTON Solver has been compiled' &
          /,2X,'       for',I4,' passive scalars.  A MHD run' &
          /,2X,'       requires that LDIMT be set to',I4,'.')
      endif
      call exitt
  ENDIF

  if (ifmvbd) then
      if (lx1 /= lx1m .OR. ly1 /= ly1m .OR. lz1 /= lz1m) &
      call exitti('Need lx1m=lx1 etc. in SIZE . $',lx1m)
  endif

  ifldmhd = npscal + 3
  if (ifmhd) then
      cpfld(ifldmhd,1) = param(29)  ! magnetic viscosity
      cpfld(ifldmhd,2) = param( 1)  ! magnetic rho same as for fluid
  endif

!   Set up default time dependent coefficients - NSTEPS,DT.

  if ( .NOT. iftran) then
      if (ifflow .AND. ifsplit) then
          iftran= .TRUE. 
      else
          param(11) = 1.0
          param(12) = 1.0
          param(19) = 0.0
      endif
  endif

!   Check here for global fast diagonalization method or z-homogeneity.
!   This is here because it influence the mesh read, which follows.
  nelx   = int(abs(param(116)))  ! check for global tensor-product structure
  nely   = int(abs(param(117)))
  nelz   = int(abs(param(118)))
  n_o    = 0

  if (n_o == 0) then
      ifzper= .FALSE. 
      ifgfdm= .FALSE. 
      if (nelz > 0) ifzper= .TRUE. 
      if (nelx > 0) ifgfdm= .TRUE. 
      if (nelx > 0) ifzper= .FALSE. 
  endif



!   Do some checks

  IF(NDIM /= LDIM)THEN
      IF(NID == 0) THEN
          WRITE(6,10) LDIM,NDIM
          10 FORMAT(//,2X,'ERROR: This NEKTON Solver has been compiled' &
          /,2X,'       for spatial dimension equal to',I2,'.' &
          /,2X,'       The data file has dimension',I2,'.')
      ENDIF
      call exitt
  ENDIF
  IF (NDIM == 3) IF3D= .TRUE. 
  IF (NDIM /= 3) IF3D= .FALSE. 

  if (if3d) then
      if (ly1 /= lx1 .OR. lz1 /= lx1) then
          if (nid == 0) write(6,13) lx1,ly1,lz1
          13 format('ERROR: lx1,ly1,lz1:',3i5,' must be equal for 3D')
          call exitt
      endif
      if (ly2 /= lx2 .OR. lz2 /= lx2) then
          if (nid == 0) write(6,14) lx2,ly2,lz2
          14 format('ERROR: lx2,ly2,lz2:',3i5,' must be equal for 3D')
          call exitt
      endif
  else
      if (ly1 /= lx1 .OR. lz1 /= 1) then
          if (nid == 0) write(6,12) lx1,ly1,lz1
          12 format('ERROR: ',3i5,' must have lx1=ly1; lz1=1, for 2D')
          call exitt
      endif
      if (ly2 /= lx2 .OR. lz2 /= 1) then
          if (nid == 0) write(6,11) lx2,ly2,lz2
          11 format('ERROR: ',3i5,' must have lx2=ly2; lz2=1, for 2D')
          call exitt
      endif
  endif

  if (lgmres < 5 .AND. param(42) == 0) then
      if(nid == 0) write(6,*) &
      'WARNING: lgmres might be too low!'
  endif

  if (ifsplit) then
      if (lx1 /= lx2) then
          if (nid == 0) write(6,43) lx1,lx2
          43 format('ERROR: lx1,lx2:',2i4,' must be equal for IFSPLIT=T')
          call exitt
      endif
  else
      if (lx2 < lx1-2) then
          if (nid == 0) write(6,44) lx1,lx2
          44 format('ERROR: lx1,lx2:',2i4,' lx2 must be lx-2 for IFSPLIT=F')
          call exitt
      endif
  endif

  if (ifmvbd .AND. ifsplit) then
      if(nid == 0) write(6,*) &
      'ABORT: Moving boundary in Pn-Pn is not supported'
      call exitt
  endif
  if (ifmoab .AND. .NOT. ifsplit) then
      if(nid == 0) write(6,*) &
      'ABORT: MOAB in Pn-Pn-2 is not supported'
      call exitt
  endif

  ktest = (lx1-lx1m) + (ly1-ly1m) + (lz1-lz1m)
  if (ifstrs .AND. ktest /= 0) then
      if(nid == 0) write(6,*) &
      'ABORT: Stress formulation requires lx1m=lx1, etc. in SIZE'
      call exitt
  endif

  if (ifgfdm .AND. ifsplit) call exitti &
  ('ERROR: FDM (p116>0) requires lx2=lx1-2 in SIZE$',lx2)

  if (ifgfdm .AND. lfdm == 0) call exitti &
  ('ERROR: FDM requires lfdm=1 in SIZE$',lfdm)

  if (ifsplit .AND. ifstrs) then
      if(nid == 0) write(6,*) &
      'ABORT: Stress formulation in Pn-Pn is not supported'
      call exitt
  endif

  if (ifsplit .AND. ifmhd) then
      if(nid == 0) write(6,*) &
      'ABORT: MHD in Pn-Pn is not supported'
      call exitt
  endif

  if (ifmhd .AND. lbx1 /= lx1) then
      if(nid == 0) write(6,*) &
      'ABORT: For MHD, need lbx1=lx1, etc.; Change SIZE '
      call exitt
  endif

  if (ifpert .AND. lpx1 /= lx1) then
      if(nid == 0) write(6,*) &
      'ABORT: For Lyapunov, need lpx1=lx1, etc.; Change SIZE '
  endif

  if (if3d) ifaxis = .FALSE. 

  if (iflomach .AND. .NOT. ifsplit) then
      if(nid == 0) write(6,*) &
      'ABORT: For lowMach, need lx2=lx1, etc.; Change SIZE '
      call exitt
  endif

  if (iflomach .AND. .NOT. ifheat) then
      if(nid == 0) write(6,*) &
      'ABORT For lowMach, need ifheat=true; Change IFHEAT'
      call exitt
  endif

!   if (ifsplit .and. param(55).ne.0) then
!      if(nid.eq.0) write(6,*)
!  $   'ABORT: Fixed mass flux not supported for Pn-Pn'
!      call exitt
!   endif


  if (ifmhd)           ifchar = .FALSE.   ! For now, at least.

!   set dealiasing handling
  if (param(99) < 0) then
      param(99) = -1       ! No  dealiasing
  else
      param(99) = 4        ! default
      if (ifaxis) param(99) = 3             ! For now, at least.
      if (ifmvbd) param(99) = 3             ! For now, at least.
  endif

  if (ifchar .AND. param(99) < 0) then
      if (nid == 0) write(6,*) &
      'ABORT: Characteristic scheme needs dealiasing!'
      call exitt
  endif

  if (param(99) > -1 .AND. (lxd < lx1 .OR. lyd < ly1 .OR. &
  lzd < lz1)) then
      if(nid == 0) write(6,*) &
      'ABORT: Dealiasing space too small; Check lxd,lyd,lzd in SIZE '
      call exitt
  endif

!   set I/O format handling
!   if (param(67).lt.0) then
!      param(67) = 0        ! ASCII
!   else ! elseif (param(67).ne.4) then
!      param(67) = 6        ! binary is default
!   endif

!   if (param(66).lt.0) then
!      param(66) = 0        ! ASCII
!   else ! elseif (param(66).ne.4) then
!      param(66) = 6        ! binary is default
!   endif

!   SET DEFAULT TO 6, ADJUSTED IN USR FILE ONLY
  param(66) = 6
  param(67) = 6

#ifndef MOAB
  if (ifmoab) then
      print *,"ABORT: ifmoab = .TRUE. in input but this ", &
      "version of nek not compiled with MOAB."
      call exitti
  endif
#endif

  return


!   Error handling:

  400 CONTINUE
  if(nid == 0) WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA' &
  ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
  call exitt

  500 CONTINUE
  if(nid == 0) WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA' &
  ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
  call exitt

  return
end subroutine rdparam

!-----------------------------------------------------------------------
!> \brief Generate local mesh elements.
!! 
!! Populate xc, yc, zc with element corner positions
!! Populate cbc, bc with element boundary conditions
!! \note Can't handle curves (zeros out)
subroutine genmesh
  use kinds, only : DP
  use size_m, only : ndim, nid, lelt
  use input, only : iffmtin, igroup, xc, yc, zc
  use input, only : curve, ccurve
  use input, only : bc, cbc
  use parallel, only : nelgt, gllnid, gllel, wdsize
  implicit none

  integer :: nsides, ieg, iel, lcbc, ldimt1
  integer :: shape_x(3)
  integer :: ix(3)
  real(DP) :: start_x(3)
  real(DP) :: end_x(3)
  real(DP) :: dx(3)
  real(DP) :: root(3)
  character(3) :: boundaries(6), tboundaries(6)

!   Read elemental mesh data, formatted
  iffmtin = .TRUE. 

  if (nid == 0) then
    read(9,*) start_x(1), end_x(1), shape_x(1)
    read(9,*) start_x(2), end_x(2), shape_x(2)
    read(9,*) start_x(3), end_x(3), shape_x(3)
    read(9,*) boundaries(1:6)
    read(9,*) tboundaries(1:6)
  endif

  call bcast(start_x,3*wdsize)  
  call bcast(end_x,  3*wdsize)  
  call bcast(shape_x,3*wdsize)  
  call bcast(boundaries,3*6)  
  call bcast(tboundaries,3*6)  

  dx = (end_x - start_x) / shape_x

  ldimt1 = 2
  curve = 0._dp
  CALL BLANK(CCURVE,12*LELT)
  LCBC=18*LELT*(LDIMT1 + 1)
  bc = 0._dp
  CALL BLANK(CBC,LCBC)

  NSIDES=NDIM*2
  DO IEG=1,NELGT
      IF (GLLNID(IEG) == NID) THEN
          IEL=GLLEL(IEG)

          ix(1) = mod(ieg - 1, shape_x(1))
          ix(2) = mod((ieg-1)/shape_x(1), shape_x(2))
          ix(3) = mod((ieg-1)/(shape_x(1)*shape_x(2)), shape_x(3))

          root = start_x + ix * dx

          igroup(iel) = 0
          XC(1,iel) = root(1)
          XC(2,iel) = root(1) + dx(1)
          XC(3,iel) = root(1) + dx(1)
          XC(4,iel) = root(1)
          XC(5,iel) = root(1)
          XC(6,iel) = root(1) + dx(1)
          XC(7,iel) = root(1) + dx(1)
          XC(8,iel) = root(1)

          YC(1,iel) = root(2)
          YC(2,iel) = root(2)
          YC(3,iel) = root(2) + dx(2)
          YC(4,iel) = root(2) + dx(2)
          YC(5,iel) = root(2)
          YC(6,iel) = root(2)
          YC(7,iel) = root(2) + dx(2)
          YC(8,iel) = root(2) + dx(2)

          ZC(1,iel) = root(3)
          ZC(2,iel) = root(3)
          ZC(3,iel) = root(3)
          ZC(4,iel) = root(3)
          ZC(5,iel) = root(3) + dx(3)
          ZC(6,iel) = root(3) + dx(3)
          ZC(7,iel) = root(3) + dx(3)
          ZC(8,iel) = root(3) + dx(3)

          CBC(:,IEL,:) = 'E'
          if (ix(2) == 0) then
            CBC(1,IEL,:) = boundaries(1)
            bc(1,1,iel,:) = ieg + (shape_x(2)-1)*shape_x(1)
          else
            bc(1,1,iel,:) = ieg - shape_x(1)
          endif

          if (ix(2) == shape_x(2) - 1) then
            CBC(3,IEL,:) = boundaries(3)
            bc(1,3,iel,:) = ieg - ix(2)*shape_x(1)
          else
            bc(1,3,iel,:) = ieg + shape_x(1)
          endif

          if (ix(1) == 0)  then
            CBC(4,IEL,:) = boundaries(4)
            bc(1,4,iel,:) = ieg + (shape_x(1) - 1)
          else
            bc(1,4,iel,:) = ieg - 1
          endif

          if (ix(1) == shape_x(1) - 1) then
            CBC(2,IEL,:) = boundaries(2)
            bc(1,2,iel,:) = ieg - ix(1)
          else
            bc(1,2,iel,:) = ieg +1
          endif

          if (ix(3) == 0) then
            CBC(5,IEL,1) = boundaries(5)
            CBC(5,IEL,2) = tboundaries(5)
            bc(1,5,iel,:) = ieg + (shape_x(3) - 1)*shape_x(2)*shape_x(1)
          else
            bc(1,5,iel,:) = ieg - shape_x(2)*shape_x(1)
          endif
          if (ix(3) == shape_x(3) - 1) then
            CBC(6,IEL,1) = boundaries(6)
            CBC(6,IEL,2) = tboundaries(6)
            bc(1,6,iel,:) = ieg - ix(3) * shape_x(2)*shape_x(1)
          else
            bc(1,6,iel,:) = ieg + shape_x(2)*shape_x(1)
          endif

          bc(2, 1, iel, :) = 3
          bc(2, 2, iel, :) = 4
          bc(2, 3, iel, :) = 1
          bc(2, 4, iel, :) = 2
          bc(2, 5, iel, :) = 6
          bc(2, 6, iel, :) = 5
       
      ENDIF
  END DO

  return
end subroutine genmesh

!-----------------------------------------------------------------------
!> \brief Read number of elements.
!! .Construct sequential element-processor partition according
!!  to number of elements and processors
!! .Selectively read mesh (defined by element vertices, and group numbers)
!!  on each processor
subroutine rdmesh
  use kinds, only : DP
  use size_m, only : ndim, nid
  use input, only : iffmtin, igroup, xc, yc, zc
  use parallel, only : nelgt, gllnid, gllel
  implicit none

  character(1) :: adum
  integer :: nsides, ieg, iel, ic

!   Read elemental mesh data, formatted
  iffmtin = .TRUE. 

  NSIDES=NDIM*2
  DO 40 IEG=1,NELGT
      IF (GLLNID(IEG) == NID) THEN
          IEL=GLLEL(IEG)

          igroup(iel) = 0
          read(9,30,err=31,end=600) igroup(iel)
          30 format(43x,i5)
      !           read(9,*,err=31,end=600) adum
          31 continue

      !           Read Corner data
          IF(NDIM == 2)THEN
              READ(9,*,ERR=500,END=600) (XC(IC,IEL),IC=1,4)
              READ(9,*,ERR=500,END=600) (YC(IC,IEL),IC=1,4)
              zc(:,iel) = 0._dp
          ELSE IF(NDIM == 3)THEN
              READ(9,*,ERR=500,END=600) (XC(IC,IEL),IC=1,4)
              READ(9,*,ERR=500,END=600) (YC(IC,IEL),IC=1,4)
              READ(9,*,ERR=500,END=600) (ZC(IC,IEL),IC=1,4)
              READ(9,*,ERR=500,END=600) (XC(IC,IEL),IC=5,8)
              READ(9,*,ERR=500,END=600) (YC(IC,IEL),IC=5,8)
              READ(9,*,ERR=500,END=600) (ZC(IC,IEL),IC=5,8)
          ENDIF
      ELSE
      !           Skip over this data for element NOT on this processor
          READ(9,41,ERR=500,END=600) ADUM
      !           Read Corner data
          IF(NDIM == 2)THEN
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
          ELSE IF(NDIM == 3)THEN
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
              READ(9,41,ERR=500,END=600) ADUM
          ENDIF
      ENDIF
  40 END DO
  41 FORMAT(A1)

!   End of mesh read.

  return

!   Error handling:

  if(nid == 0) WRITE(6,401)
  401 FORMAT(2X,'ERROR READING SCALE FACTORS, CHECK READ FILE' &
  ,/,2X,'ABORTING IN ROUTINE RDMESH.')
  call exitt

  500 CONTINUE
  if(nid == 0) WRITE(6,501) IEG
  501 FORMAT(2X,'ERROR READING MESH DATA NEAR ELEMENT',I12 &
  ,/,2X,'ABORTING IN ROUTINE RDMESH.')
  call exitt

  600 CONTINUE
  if(nid == 0) WRITE(6,601) IEG
  601 FORMAT(2X,'ERROR 2 READING MESH DATA NEAR ELEMENT',I12 &
  ,/,2X,'ABORTING IN ROUTINE RDMESH.')
  call exitt

  return
end subroutine rdmesh

!-----------------------------------------------------------------------
!> \brief .Read curve side data
!! .Disperse curve side data to all processors according
!!  to sequential partition scheme
subroutine rdcurve
  use kinds, only : DP
  use size_m, only : lelt, nid
  use input, only : iffmtin, curve, ccurve
  use parallel, only : nelgt, gllnid, gllel
  implicit none

  CHARACTER(1) :: ANS
  integer :: ncurve, icurve, iedg, ieg, iel
  real(DP) :: r1, r2, r3, r4, r5

  IF (IFFMTIN) THEN
  
  !     Read formatted curve side data
  
      READ(9,*)
      READ(9,*)NCURVE
      curve = 0._dp
      CALL BLANK(CCURVE,12*LELT)
      IF (NCURVE > 0) THEN
          DO 50 ICURVE=1,NCURVE
              IF (NELGT < 1000) THEN
                  READ(9,60,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
              ELSEIF (NELGT < 1000000) THEN
                  READ(9,61,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
              ELSE
                  READ(9,62,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
              ENDIF
              60 FORMAT(I3,I3 ,5G14.6,1X,A1)
              61 FORMAT(I2,I6 ,5G14.6,1X,A1)
              62 FORMAT(I2,I12,5G14.6,1X,A1)

              IF (GLLNID(IEG) == NID) THEN
                  IEL=GLLEL(IEG)
                  CURVE (1,IEDG,IEL)=R1
                  CURVE (2,IEDG,IEL)=R2
                  CURVE (3,IEDG,IEL)=R3
                  CURVE (4,IEDG,IEL)=R4
                  CURVE (5,IEDG,IEL)=R5
                  CCURVE(  IEDG,IEL)=ANS
              ENDIF
          50 END DO
      ENDIF
      return
  
  !     Error handling:
  
      500 CONTINUE
      if(nid == 0) WRITE(6,501)
      501 FORMAT(2X,'ERROR READING CURVE SIDE DATA' &
      ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
      return
  
  ELSE
  
  !     Read unformatted curve side data
  
      READ(8) NCURVE
      curve = 0._dp
      CALL BLANK(CCURVE,12*LELT)
      IF (NCURVE > 0) THEN
          DO 1050 ICURVE=1,NCURVE
              READ(8,ERR=1500,END=1500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
              IF (GLLNID(IEG) == NID) THEN
                  IEL=GLLEL(IEG)
                  CURVE (1,IEDG,IEL)=R1
                  CURVE (2,IEDG,IEL)=R2
                  CURVE (3,IEDG,IEL)=R3
                  CURVE (4,IEDG,IEL)=R4
                  CURVE (5,IEDG,IEL)=R5
                  CCURVE(  IEDG,IEL)=ANS
              ENDIF
          1050 END DO
      ENDIF
      return
  
  !     Error handling:
  
      1500 CONTINUE
      if(nid == 0) WRITE(6,1501)
      1501 FORMAT(2X,'ERROR READING unformatted CURVE SIDE DATA' &
      ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
  
      return
  ENDIF
  end subroutine rdcurve

!-----------------------------------------------------------------------
!> \brief Read Boundary Conditions (and connectivity data).
!! .Disperse boundary condition data to all processors
!!  according to sequential partition scheme
subroutine rdbdry
  use kinds, only : DP
  use size_m, only : ndim, lelt, ldimt1, nid
  use input, only : cbc, bc, vnekton, npscal
  use input, only : ifheat, ifmhd, ifflow, iffmtin, iftmsh
  use parallel, only : nelgt, nelgv, gllnid, gllel
  use scratch, only : cbcs, bcs
  use string, only : capit, indx1
  implicit none

  CHARACTER(1) :: CBC1, chtemp
  character(3) :: CBC3, chtmp3
  EQUIVALENCE (CHTEMP,CHTMP3)
  character(132) :: tmp_string
  integer :: nfldt, nbcs, ibcs, nsides, lcbc, ibcnew, ifield, nel, nbcrea, ieg
  integer :: ii, iside, icbc1, id2, id1, iel, lrbc

!   Set up TEMPORARY value for NFIELD - NFLDT

  NFLDT = 1
  IF (IFHEAT) NFLDT=2+NPSCAL
  if (ifmhd ) nfldt=2+npscal+1
  NBCS      = NFLDT
  IBCS      = 2
  IF (IFFLOW) IBCS = 1
  NSIDES    = 2*NDIM

!   Read boundary conditions for all fields

  LCBC=18*LELT*(LDIMT1 + 1)
  LRBC=30*LELT*(LDIMT1 + 1)
  bc = 0._dp
  CALL BLANK(CBC,LCBC)

!-----------------------------------------------------------------
!  Formatted Reads
!-----------------------------------------------------------------

  IF (IFFMTIN) THEN
  
      READ(9,*,ERR=500,END=500)  !   ***** BOUNDARY CONDITIONS *****
      ibcnew = 1
      DO 100 IFIELD=ibcnew,NBCS  !     DO 100 IFIELD=IBCS,NBCS
          NEL=NELGT
          if ( .NOT. iftmsh(ifield)) nel = nelgv
      !       Fluid and/or thermal
          read(9,81) tmp_string        !  ***** FLUID   BOUNDARY CONDITIONS *****
          call capit(tmp_string,132)

      !       write(6,*) 'reading BC:',ifield,ibcs,nbcs
      !       write(6,81) string
      !       if1 = indx1(string,'NO ',3)
      !       write(6,*) if1,' if NO.  quit.',ifield,ibcs,nbcs
      !       write(6,*) ifield,iftmsh(ifield),nel,' iftmsh'
      !       call exitt


          if (indx1(tmp_string,'NO ',3) == 0) then ! we have acitve bc info
 
              nbcrea = -1 ! below should catch
              IF(VNEKTON <= 2.52) NBCREA = 3
              IF(VNEKTON >= 2.55) NBCREA = 5
          
              DO IEG=1,NEL
                  DO ISIDE=1,NSIDES
                      IF (GLLNID(IEG) == NID) THEN
                          IEL=GLLEL(IEG)
                          IF (NELGT < 1000) THEN
                              READ(9,50,ERR=500,END=500) &
                              CHTEMP, &
                              CBC(ISIDE,IEL,IFIELD),ID1,ID2, &
                              (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                          !                 write(6,50)
                          !    $            CHTEMP,
                          !    $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
                          !    $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                              50 FORMAT(A1,A3,2I3,5G14.6)
                          ELSEIF (NELGT < 100000) THEN
                              READ(9,51,ERR=500,END=500) &
                              CHTEMP, &
                              CBC(ISIDE,IEL,IFIELD),ID1,ID2, &
                              (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                              51 FORMAT(A1,A3,I5,I1,5G14.6)
                          ELSEIF (NELGT < 1000000) THEN
                              READ(9,52,ERR=500,END=500) &
                              CHTEMP, &
                              CBC(ISIDE,IEL,IFIELD),ID1, &
                              (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                              52 FORMAT(A1,A3,I6,5G14.6)
                          ELSE
                              READ(9,53,ERR=500,END=500) &
                              CHTEMP, &
                              CBC(ISIDE,IEL,IFIELD),ID1, &
                              (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                              53 FORMAT(A1,A3,I12,5G18.11)
                          ENDIF
                      !              Mesh B.C.'s in 1st column of 1st field
                          IF (CHTEMP /= ' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
                      !              check for fortran function as denoted by lower case bc's:
                          CBC3=CBC(ISIDE,IEL,IFIELD)
                          ICBC1=ICHAR(CBC3(1:1))
                      !              IF (ICBC1.GE.97.AND.ICBC1.LE.122) THEN
                      !                 IF(CBC3(3:3).NE.'i')NLINES=BC(1,ISIDE,IEL,IFIELD)
                      !                 IF(CBC3(3:3).EQ.'i')NLINES=BC(4,ISIDE,IEL,IFIELD)
                      !                 DO 60 I=1,NLINES
                      !  60             READ(9,*,ERR=500,END=500)
                      !              ENDIF
                      ELSE
                          READ(9,*,ERR=500,END=500)   cbc1  ! dummy read, pff 4/28/05
                      ENDIF
                  enddo
              END DO
          endif
          81 format(a132)
      100 END DO
  
  !     END OF BC READ
  
  !     Check for dummy line:  "NO THERMAL B.C.'S"
      IF (NFLDT == 1) READ(9,*,ERR=500,END=500)
  
      return
  
  !     Error handling:
  
      500 CONTINUE
      if(nid == 0) WRITE(6,501) IFIELD,IEG
      501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I12 &
      ,/,2X,'ABORTING IN ROUTINE RDBDRY.')
      call exitt
      return
  
  
  ELSE
  
  !-----------------------------------------------------------------
  !  UNformatted Reads
  !-----------------------------------------------------------------
  
  !     READ(8,ERR=500,END=500)
      DO IFIELD=IBCS,NBCS
          NEL=NELGT
      !        Fluid and/or thermal
          NBCREA = 5
      
          DO IEG=1,NEL
              DO ISIDE=1,NSIDES
                  IF (GLLNID(IEG) == NID) THEN
                      IEL=GLLEL(IEG)
                      READ(8,ERR=1500,END=1500) &
                      CHTMP3, &
                      CBC(ISIDE,IEL,IFIELD),ID1,ID2, &
                      (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                  
                  !              Mesh B.C.'s in 1st column of 1st field
                      IF (CHTEMP /= ' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
                  !              check for fortran function as denoted by lower case bc's:
                  ELSE
                      IEL=1
                      READ(8,ERR=1500,END=1500) CHTMP3, &
                      CBCS(ISIDE,IEL),ID1,ID2,(BCS(II,ISIDE,IEL),II=1,NBCREA)
                  !              check for fortran function as denoted by lower case bcs:
                  ENDIF
              enddo
          END DO
      END DO
  
  !     END OF BC READ
  
      return
  
  !     Error handling:
  
      1500 CONTINUE
      if(nid == 0) WRITE(6,1501) IFIELD,IEG
      1501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I12 &
      ,/,2X,'(unformatted) ABORTING IN ROUTINE RDBDRY.')
      call exitt
  ENDIF

  return
end subroutine rdbdry

!-----------------------------------------------------------------------
!> \brief Read Initial Conditions / Drive Force.
!! Broadcast ICFILE to all processors
subroutine rdicdf
  use size_m, only : nid
  use input, only : initc
  use parallel, only : csize
  use string, only : capit, indx1, ifgtil
  implicit none

  character(132) :: line
  integer :: ierr, nskip, i
  integer, external :: iglmax

  ierr = 0

  if (nid == 0) then   !  Read names of restart files

      call blank(initc,15*132)
      read (9,80,err=200,end=200) line
      call capit(line,132)
      if (indx1(line,'RESTART',7) /= 0) then
          if ( .NOT. ifgtil(nskip,line)) goto 200
      !          read(line,*,err=200,end=200) nskip
          do 50 i=1,nskip
              read(9,80,err=200,end=200) initc(i)
          50 END DO
          read(9,80,err=200,end=200) line
      endif
      80 format(a132)

      if ( .NOT. ifgtil(nskip,line)) goto 200

  !       Read initial conditions
      do 100 i=1,nskip
          read(9,80,err=200,end=200) line
      100 END DO

  !       Read drive force data
      read(9,*,err=200,end=200)
      read(9,*,err=200,end=200) nskip
      do 110 i=1,nskip
          read(9,80,err=200,end=200) line
      110 END DO
  endif

  ierr = iglmax(ierr,1)
  if (ierr == 0) then
      call bcast(initc,15*132*csize)
      return
  else
      goto 210
  endif

!   Error handling:

  200 ierr = 1
  ierr = iglmax(ierr,1)
        
  210 continue
  if (nid == 0) write(6,300)
  300 format(2x,'Error reading initial condition/drive force data' &
  ,/,2x,'aborting in routine rdicdf.')
  call exitti('rdicdf error$',ierr)

  return
end subroutine rdicdf

!-----------------------------------------------------------------------
!> \brief .Read materials property data
!!  .Disperse material properties to all processors according
!!   to sequential partition scheme
subroutine rdmatp
  use kinds, only : DP
  use size_m, only : ldimt1, nid
  use input, only : matype, cpgrp, ifvps
  use parallel, only : isize, wdsize
  implicit none

  CHARACTER(132) :: LINE
  integer :: nskip, npacks, iig, igrp, ifld, itype, iprop

  matype = 0
  cpgrp = 0._dp

!   Read material property data

  IF(NID == 0) THEN
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NSKIP
      READ(9,*,ERR=200,END=200) NPACKS
      DO IIG=1,NPACKS
          IFVPS= .TRUE. 
          READ(9,*)IGRP,IFLD,ITYPE
          MATYPE(IGRP,IFLD)=ITYPE
          DO IPROP=1,3
              IF(ITYPE == 1) READ(9,* ) CPGRP(IGRP,IFLD,IPROP)
              IF(ITYPE == 2) READ(9,80) LINE
              80 FORMAT(A132)
          enddo
      END DO
  ENDIF

  CALL BCAST(MATYPE,16*LDIMT1*ISIZE)
  CALL BCAST(CPGRP ,48*LDIMT1*WDSIZE)

  return

!   Error handling:

  200 CONTINUE
  if(nid == 0) WRITE(6,201)
  201 FORMAT(2X,'ERROR READING MATERIAL PROPERTIES DATA' &
  ,/,2X,'ABORTING IN ROUTINE RDMATP.')
  call exitt

  return
end subroutine rdmatp

!-----------------------------------------------------------------------
!> \brief .Read history data
!! .Broadcast to all processors
subroutine rdhist
  use size_m, only : lhis, nid, nx1, ny1, nz1
  use input, only : lochis, hcode, nhis
  use parallel, only : nelgt, isize, csize
  implicit none

  integer :: ierr, i, ii, i2

  CALL BLANK (HCODE ,11*lhis)
  lochis = 0

  ierr=0
  IF(NID == 0) THEN
  !       Read history data
      READ (9,*)
      READ (9,*,ERR=200,END=200) NHIS
      if (nhis > lhis) then
          write(6,*) nid,' Too many history pts. RESET LHIS.',nhis,lhis
          ierr=1
      endif
       
      if(ierr == 0) then
      !       HCODE(10) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
          if (nhis > 0) then
              do i=1,nhis
                  if (nelgt < 100000) then
                      read(9,130,err=200,end=200) &
                      (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
                      130 format(1x,11a1,1x,4i5)
                  else
                      read(9,131,err=200,end=200) &
                      (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
                      131 format(1x,11a1,1x,3i5,i10)
                  endif
              
              !           threshold lochis locations to allow easy specification of "NX,NY,NZ"
              !           pff 1/7/97
              
                  if (hcode(10,i) == 'H') then
                      lochis(1,i) = min(lochis(1,i),nx1)
                      lochis(2,i) = min(lochis(2,i),ny1)
                      lochis(3,i) = min(lochis(3,i),nz1)
                  
                  !              if lochis_k = -1, set it to nxk/2   pff 8/21/03
                  
                      if (lochis(1,i) == -1) lochis(1,i) = (nx1+1)/2
                      if (lochis(2,i) == -1) lochis(2,i) = (ny1+1)/2
                      if (lochis(3,i) == -1) lochis(3,i) = (nz1+1)/2
                  endif
              enddo
          endif
      endif
  ENDIF
  call err_chk(ierr,' Too many histroy pts. RESET LHIS$')

  call bcast(NHIS  ,ISIZE)
  call bcast(HCODE ,11*LHIS*CSIZE)
  call bcast(LOCHIS,4*LHIS*ISIZE)

  return

!   Error handling:

  200 CONTINUE
  if(nid == 0) WRITE(6,201)
  201 FORMAT(2X,'ERROR READING HISTORY DATA' &
  ,/,2X,'ABORTING IN ROUTINE RDHIST.')
  call exitt

  return
end subroutine rdhist

!-----------------------------------------------------------------------
!> \brief Read output specs, broadcast to all processors
subroutine rdout
  use size_m, only : nid, ldimt1
  use input, only : ifpsco, ifxyo, ifvo, ifpo, ifto, ifbo, ipsco
  use parallel, only : lsize, isize
  implicit none

  logical :: lbuf(5+ldimt1)
  integer :: iflag, nouts, k, i
  integer, external :: iglmax

  call lfalse(lbuf,5+ldimt1)
  iflag = 0                           ! Check for valid ipsco read

  IF(NID == 0) THEN                   ! Read output specs

      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOUTS
      READ(9,*,ERR=200,END=200) IFXYO
      READ(9,*,ERR=200,END=200) IFVO
      READ(9,*,ERR=200,END=200) IFPO
      READ(9,*,ERR=200,END=200) IFTO
      READ(9,*,ERR=200,END=200) IFBO   !  IFTGO

      lbuf(1) = IFXYO
      lbuf(2) = IFVO
      lbuf(3) = IFPO
      lbuf(4) = IFTO
      lbuf(5) = IFBO

      k = 5

      call lfalse(ifpsco,ldimt1)
      read(9,*,err=200,end=200) ipsco
      if (ipsco > 0) then
          if (ipsco > ldimt1) then    ! Invalid ifpsco read
              iflag = 1
          else
              do i=1,ipsco
                  read(9,*,err=200,end=200) ifpsco(i)
                  k = k+1
                  lbuf(k) = ifpsco(i)
              enddo
          endif
      endif

  endif


  iflag = iglmax(iflag,1)                       ! Check for valid ipsco read
  if (iflag > 0) call exitti                    & ! Invalid ifpsco read
  ('Error in rdout.  Increase ldimt1 in SIZE to$',ipsco)

  k = 5+ldimt1
  call bcast(lbuf ,LSIZE*k)
  call bcast(IPSCO,ISIZE  )

  ifxyo = lbuf(1)
  ifvo  = lbuf(2)
  ifpo  = lbuf(3)
  ifto  = lbuf(4)
  ifbo  = lbuf(5)

  k = 5
  do i=1,ipsco
      k = k+1
      ifpsco(i) = lbuf(k)
  enddo

  return


!   Error handling:

  200 CONTINUE
  WRITE(6,201)
  201 FORMAT(2X,'ERROR READING OUTPUT SPECIFICATION DATA' &
  ,/,2X,'ABORTING IN ROUTINE RDOUT.')
  call exitt

  return
  end subroutine rdout

!-----------------------------------------------------------------------
!> \brief Read objects, Broadcast to all processors
subroutine rdobj
  use size_m, only : nid, maxobj, maxmbr
  use input, only : nobj, nmember, object
  use parallel, only : isize
  implicit none

  integer :: ierr, iobj, member, k

!   Default if no data is read No Objects
  ierr=0
  IF(NID == 0) THEN
      NOBJ=0
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOBJ
       
      IF(NOBJ > MAXOBJ) ierr=1
       
      if(ierr == 0) then
          DO 10 IOBJ = 1,NOBJ
              READ(9,*,ERR=200,END=200) NMEMBER(IOBJ)
              IF(NMEMBER(IOBJ) > MAXMBR)THEN
                  PRINT*,'ERROR: Too many members in object ',IOBJ
                  ierr=2
              ENDIF
              if(ierr == 0) then
                  DO 5 MEMBER=1,NMEMBER(IOBJ)
                      READ(9,*,ERR=200,END=200) OBJECT(IOBJ,MEMBER,1), &
                      OBJECT(IOBJ,MEMBER,2)
                  5 END DO
              endif
          10 END DO
          write(6,*) nobj,' objects found' &
          ,(nmember(k),k=1,nobj)
      endif
  endif
  call err_chk(ierr,'ERROR, too many objects:$')
   
  call bcast(NOBJ   ,ISIZE)
  call bcast(NMEMBER,MAXOBJ*ISIZE)
  call bcast(OBJECT ,MAXOBJ*MAXMBR*2*ISIZE)

   
  return

!   Error handling:  For old versions, default to no objects

  200 CONTINUE
  NOBJ=0
   
  return
end subroutine rdobj

!=====================================================================
!> \brief Verify that mesh and dssum are properly defined by performing
!!    a direct stiffness operation on the X,Y and Z coordinates.
!! Note that periodic faces are not checked here.
!! \todo clean up ta,tb,qmask
!=====================================================================
subroutine vrdsmsh()
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelt, ndim, nid
  use geom, only : xm1, ym1, zm1
  use input, only : ifheat, cbc, if3d
  use parallel, only : lglel
  use soln, only : tmult, vmult
  use tstep, only : ifield
  implicit none

  real(DP), allocatable :: TA(:,:,:,:),TB(:,:,:,:) &
  ,QMASK(:,:,:,:)
  real(DP) :: tmp(2)

  CHARACTER(3) :: CB
  integer :: ierr, nxyz1, ntot, nfaces, ie, ieg, ix, iy, iz, iface, iel
  real(DP) :: eps, xmx, xmn, ymx, ymn, zmx, zmn, xscmax, xscmin
  real(DP) :: scal1, scal2, scal3, xscale, yscmax, yscmin, zscale
  real(DP) :: yscale, zscmax, zscmin
  real(DP), external :: glmin, glmax
  integer :: iglsum

!    call  vrdsmshx  ! verify mesh topology

  allocate(TA(LX1,LY1,LZ1,LELT),TB(LX1,LY1,LZ1,LELT) &
  ,QMASK(LX1,LY1,LZ1,LELT))

  if(nid == 0) write(*,*) 'verify mesh topology'

  IERR      = 0
  EPS       = 1.0e-04
  EPS       = 1.0e-03
  IFIELD    = 1
  IF (IFHEAT) IFIELD = 2
  NXYZ1     = NX1*NY1*NZ1
  NTOT      = NX1*NY1*NZ1*NELT
  NFACES    = 2*NDIM

  xmx = glmax(xm1,ntot)
  xmn = glmin(xm1,ntot)
  ymx = glmax(ym1,ntot)
  ymn = glmin(ym1,ntot)
  zmx = glmax(zm1,ntot)
  zmn = glmin(zm1,ntot)
  if (nid == 0) write(6,*) xmn,xmx,' Xrange'
  if (nid == 0) write(6,*) ymn,ymx,' Yrange'
  if (nid == 0) write(6,*) zmn,zmx,' Zrange'
!   return

!   First check - use 1/Multiplicity

  IF (IFHEAT) THEN
      CALL COPY(TA,TMULT,NTOT)
  ELSE
      CALL COPY(TA,VMULT,NTOT)
  ENDIF

!   write(6,1)
!  $(nid,'tab4',lglel(ie),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
! 1 format(i3,a4,i3,16f5.2)

  CALL DSSUM(TA)

!   write(6,1)
!  $(nid,'taaf',lglel(ie),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)

  tb = 1._dp - ta
  DO IE=1,NELT
      ieg=lglel(ie)
      DO IZ=1,NZ1
          DO IY=1,NY1
              DO IX=1,NX1
                  IF (ABS(TB(IX,IY,IZ,IE)) > EPS ) THEN
                      WRITE(6,1005) IX,IY,IZ,IEG &
                      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE) &
                      ,TA(IX,IY,IZ,IE),eps
                  !           WRITE(7,1005) IX,IY,IZ,IEG
                  !    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
                  !    $      ,QMASK(IX,IY,IZ,IE)
                      1005 FORMAT(2X,'WARNING: DSSUM problem at:',/ &
                      ,1X,'I,J,K,IE:',3I5,i12,/ &
                      ,2X,'Near X =',3G16.8,', d:',2G16.8)
                      IERR=4
                  ENDIF
              enddo
          enddo
      enddo
  END DO

!   Set up QMASK quickly to annihilate checks on periodic bc's

  qmask = 1._dp
  DO IEL=1,NELT
      DO IFACE=1,NFACES
          CB =CBC(IFACE,IEL,IFIELD)
          IF (CB == 'P  ' .OR. cb == 'p  ') &
          CALL FACEV(QMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
      enddo
  END DO
  CALL DSOP(QMASK,'MUL')

!    xxmin = glmin(xm1,ntot)
!    yymin = glmin(ym1,ntot)
!    zzmin = glmin(zm1,ntot)
!    xxmax = glmax(xm1,ntot)
!    yymax = glmax(ym1,ntot)
!    zzmax = glmax(zm1,ntot)
!    if (nid.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
!  7 format('xyz minmx2:',6g13.5)




!   X-component

  CALL COPY(TA,XM1,NTOT)
  CALL COPY(TB,XM1,NTOT)
  CALL DSOP(TA,'MIN')
  CALL DSOP(TB,'MAX')
  ta = (ta - xm1) * qmask
  tb = (tb - xm1) * qmask
  DO IE=1,NELT
      XSCMAX = maxval(XM1(:,:,:,IE))
      XSCMIN = minval(XM1(:,:,:,IE))
      SCAL1=ABS(XSCMAX-XSCMIN)
      SCAL2=ABS(XSCMAX)
      SCAL3=ABS(XSCMIN)
      SCAL1=MAX(SCAL1,SCAL2)
      SCAL1=MAX(SCAL1,SCAL3)
      XSCALE = 1./SCAL1
      ieg=lglel(ie)
      DO IZ=1,NZ1
          DO IY=1,NY1
              DO IX=1,NX1
                  if (abs(ta(ix,iy,iz,ie)*xscale) > eps .OR. &
                  abs(tb(ix,iy,iz,ie)*xscale) > eps ) then
                      write(6,1105) ix,iy,iz,ieg &
                      ,xm1(ix,iy,iz,ie),ym1(ix,iy,iz,ie),zm1(ix,iy,iz,ie) &
                      ,tb(ix,iy,iz,ie),ta(ix,iy,iz,ie),XSCALE
                      1105 format(1x,'WARNING1 Element mesh mismatch at:',/ &
                      ,1x,'i,j,k,ie:',3i5,I12,/ &
                      ,1X,'Near X =',3G16.8,', d:',3G16.8)
                      ierr=1
                  endif
              enddo
          enddo
      enddo
  enddo

!   Y-component

  CALL COPY(TA,YM1,NTOT)
  CALL COPY(TB,YM1,NTOT)
  CALL DSOP(TA,'MIN')
  CALL DSOP(TB,'MAX')
  ta = (ta - ym1) * qmask
  tb = (tb - ym1) * qmask
  DO IE=1,NELT
      YSCMAX = MAXVAL(YM1(:,:,:,IE))
      YSCMIN = MINVAL(YM1(:,:,:,IE))
      SCAL1=ABS(YSCMAX-YSCMIN)
      SCAL2=ABS(YSCMAX)
      SCAL3=ABS(YSCMIN)
      SCAL1=MAX(SCAL1,SCAL2)
      SCAL1=MAX(SCAL1,SCAL3)
      YSCALE = 1./SCAL1
      ieg=lglel(ie)
      DO IZ=1,NZ1
          DO IY=1,NY1
              DO IX=1,NX1
                  IF (ABS(TA(IX,IY,IZ,IE)*YSCALE) > EPS .OR. &
                  ABS(TB(IX,IY,IZ,IE)*YSCALE) > EPS ) THEN
                      WRITE(6,1205) IX,IY,IZ,IEG &
                      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE) &
                      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),yscale
                      1205 FORMAT(1X,'WARNING2 Element mesh mismatch at:',/ &
                      ,1X,'I,J,K,IE:',3I5,i12,/ &
                      ,1X,'Near Y =',3G16.8,', d:',3G16.8)
                      IERR=2
                  ENDIF
              enddo
          enddo
      enddo
  enddo

!   Z-component

  IF (IF3D) THEN
      CALL COPY(TA,ZM1,NTOT)
      CALL COPY(TB,ZM1,NTOT)
      CALL DSOP(TA,'MIN')
      CALL DSOP(TB,'MAX')
      ta = (ta - zm1) * qmask
      tb = (tb - zm1) * qmask
      DO IE=1,NELT
          ZSCMAX = MAXVAL(ZM1(:,:,:,IE))
          ZSCMIN = MINVAL(ZM1(:,:,:,IE))
          SCAL1=ABS(ZSCMAX-ZSCMIN)
          SCAL2=ABS(ZSCMAX)
          SCAL3=ABS(ZSCMIN)
          SCAL1=MAX(SCAL1,SCAL2)
          SCAL1=MAX(SCAL1,SCAL3)
          ZSCALE = 1./SCAL1
          ieg=lglel(ie)
          DO IZ=1,NZ1
              DO IY=1,NY1
                  DO IX=1,NX1
                      IF (ABS(TA(IX,IY,IZ,IE)*ZSCALE) > EPS .OR. &
                      ABS(TB(IX,IY,IZ,IE)*ZSCALE) > EPS ) THEN
                          WRITE(6,1305) IX,IY,IZ,IEG &
                          ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE) &
                          ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),zscale
                          1305 FORMAT(1X,'WARNING3 Element mesh mismatch at:',/ &
                          ,1X,'I,J,K,IE:',3I5,i12,/ &
                          ,1X,'Near Z =',3G16.8,', d:',3G16.8)
                          IERR=3
                      ENDIF
                  enddo
              enddo
          enddo
      enddo
  ENDIF
   
  ierr = iglsum(ierr,1)
  IF (IERR > 0) THEN
      if(nid == 0) WRITE(6,1400)
      1400 FORMAT &
      (' Mesh consistency check failed.  EXITING in VRDSMSH.')
      call exitt
  ENDIF
   
  tmp(1)=ierr
  CALL GOP(tmp,tmp(2),'M  ',1)
  IF (tmp(1) >= 4.0) THEN
      WRITE(6,1400) &
      (' Mesh consistency check failed.  EXITING in VRDSMSH.')
      call exitt
  ENDIF
   
  if(nid == 0) then
      write(6,*) 'done :: verify mesh topology'
      write(6,*) ' '
  endif

  return
end subroutine vrdsmsh

!-----------------------------------------------------------------------
subroutine chk_nel
  use size_m, only : lelt, lelv, lelg, nid, nelt
  use parallel, only : np, nelgt, nelgv, nelgt_max
  implicit none

  integer :: neltmx, nelvmx, lelt_needed
  integer, external :: iglmax

  neltmx=np*lelt
  nelvmx=np*lelv

  neltmx=min(neltmx,lelg)
  nelvmx=min(nelvmx,lelg)

  nelgt = iglmax(nelgt,1)
  nelgv = iglmax(nelgv,1)

!   write(6,*) nid,' inside chk_nel',nelgt,neltmx,nelvmx

  if (nelgt > neltmx .OR. nelgv > nelvmx) then
      if (nid == 0) then
          lelt_needed = nelgt/np
          if (mod(nelgt,np) /= 0) lelt_needed = lelt_needed + 1
          write(6,12) lelt,lelg,lelt_needed,np,nelgt
          12 format(//,2X,'ABORT: Problem size too large!' &
          ,/,2X &
          ,/,2X,'This solver has been compiled for:' &
          ,/,2X,'   number of elements/proc  (lelt):',i12 &
          ,/,2X,'   total number of elements (lelg):',i12 &
          ,/,2X &
          ,/,2X,'Recompile with the following SIZE  parameters:' &
          ,/,2X,'   lelt >= ',i12,'  for np = ',i12 &
          ,/,2X,'   lelg >= ',i12,/)
      !           write(6,*)'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
      !           write(6,*)'help:',lelt,lelv,lelgv
      endif
      call exitt
  endif

  if(nelgt > nelgt_max) then
      if(nid == 0) write(6,*) &
      'ABORT: Total number of elements too large!', &
      '       nel_max = ', nelgt_max
      call exitt
  endif

  if (nelt > lelt) then
      write(6,'(A,3I12)') 'ABORT: nelt>lelt!', nid, nelt, lelt
      call exitt
  endif

  return
end subroutine chk_nel

