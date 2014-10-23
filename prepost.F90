!-----------------------------------------------------------------------
!> \brief Store results for later postprocessing.
!!  Recent updates:
!!  p65 now indicates the number of parallel i/o files; iff p66 >= 6
subroutine prepost(ifdoin,prefin)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nid
  use ctimer, only : icalld, nprep, etime1, dnekclock, tprep
  use input, only : schfle, ifschclob, ifpsco
  use tstep, only : iostep, timeio, istep, nsteps, lastep, time, ntdump
  implicit none

  logical :: ifdoin
  character(3) :: prefin

  real(DP) :: tdmp(4)

  character(3) :: prefix

  logical, save :: ifdoit = .FALSE.

  real(DP), allocatable ::  pm1(:,:,:,:)

  logical :: ifhis

  integer, save :: maxstep = 999999999
  integer :: ierr, iiidmp, idummy
  real(DP) :: timdump

  if (iostep < 0 .OR. timeio < 0) return

  icalld=icalld+1
  nprep=icalld

#ifndef NOTIMER
  etime1=dnekclock()
#endif

!   Trigger history output only if prefix = 'his'   pff 8/18/05
  ifhis  = .FALSE. 
  prefix = prefin
  if (prefin == 'his') ifhis  = .TRUE. 
  if (prefix == 'his') prefix = '   '

  if(icalld == 1) then
      ierr = 0
      if (nid == 0) then
          write(6,*) 'schfile:',schfle
          if (ifschclob) then
              open(unit=26,file=schfle,err=44,form='formatted')
          else
              open(unit=26,file=schfle,err=44,form='formatted', &
              status='new')
          endif
          goto 45
          44 ierr = 1
      endif
      45 continue
      call err_chk(ierr,'.sch file already exists. Use IFSCHCLOB=F to &
     & disable this check BUT BEWARE!!!!!!$')
  endif

  allocate(pm1(lx1,ly1,lz1,lelv))
  call prepost_map(0, pm1) ! map pr and axisymm. arrays

  if(istep >= nsteps) lastep=1

  timdump=0
  if(timeio /= 0.0)then
      if(time >= (ntdump + 1) * timeio) then
          timdump=1.
          ntdump=ntdump+1
      endif
  endif

  if (istep > 0 .AND. iostep > 0) then
      if(mod(istep,iostep) == 0) ifdoit= .TRUE. 
  endif


! check for io request in file 'ioinfo'
  iiidmp=0
  if (nid == 0 .AND. (mod(istep,10) == 0 .OR. istep < 200)) then
      open(unit=87,file='ioinfo',status='old',err=88)
      read(87,*,end=87,err=87) idummy
      if (iiidmp == 0) iiidmp=idummy
      if (idummy /= 0) then  ! overwrite last i/o request
          rewind(87)
          write(87,86)
          86 format(' 0')
      endif
      87 continue
      close(unit=87)
      88 continue
      if (iiidmp /= 0) write(6,*) 'Output:',iiidmp
  endif

  tdmp(1)=iiidmp
  call gop(tdmp,tdmp(3),'+  ',1)
  iiidmp= int(tdmp(1))
  if (iiidmp < 0) maxstep=abs(iiidmp)
  if (istep >= maxstep .OR. iiidmp == -2) lastep=1
  if (iiidmp == -2) return
  if (iiidmp < 0) iiidmp = 0

  if (ifdoin) ifdoit= .TRUE. 
  if (iiidmp /= 0 .OR. lastep == 1 .OR. timdump == 1.) ifdoit= .TRUE. 


  if (ifdoit .AND. nid == 0)write(6,*)'call outfld: ifpsco:',ifpsco(1)
  if (ifdoit) call outfld(prefix, pm1)

  call outhis(ifhis, pm1)

  call prepost_map(1, pm1) ! map back axisymm. arrays

  if (lastep == 1 .AND. nid == 0) close(unit=26)

#ifndef NOTIMER
  tprep=tprep+dnekclock()-etime1
#endif

  ifdoit= .FALSE. 
  return

end subroutine prepost

!-----------------------------------------------------------------------
!> \brief Store results for later postprocessing
subroutine prepost_map(isave, pm1) ! isave=0-->fwd, isave=1-->bkwd
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, ly2, lz2
  use size_m, only : nx1, ny1, nz1, nelv, nx2, ny2, nz2
  use input, only : ifaxis, ifsplit
  use ixyz, only : ixm21, iytm21, iztm21
  use soln, only : pr
  use tstep, only : if_full_pres
  implicit none

  integer, intent(in) :: isave
  real(DP), intent(out) ::  pm1    (lx1,ly1,lz1,lelv)

!   Work arrays and temporary arrays
  real(DP) :: pa(lx1,ly2,lz2),pb(lx1,ly1,lz2)
  integer :: e
  integer :: ntot1, nyz2, nxy1, nxyz, nxyz2, iz

  if (isave == 0) then ! map to GLL grid

      if (ifaxis) then
        write(*,*) "Oops: ifaxis"
#if 0
          ntotm1 = nx1*ny1*nelt
          call copy (yax,ym1,ntotm1)
          do 5 e=1,nelt
              if (ifrzer(e)) then
                  call mxm  (ym1(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                  call copy (ym1(1,1,1,e),pb,nx1*ny1)
              endif
          5 END DO
          if (ifflow) then
              ntotm1 = nx1*ny1*nelv
              ntotm2 = nx2*ny2*nelv
              call copy (vxax,vx,ntotm1)
              call copy (vyax,vy,ntotm1)
              call copy (prax,pr,ntotm2)
              do 10 e=1,nelv
                  if (ifrzer(e)) then
                      call mxm  (vx(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                      call copy (vx(1,1,1,e),pb,nx1*ny1)
                      call mxm  (vy(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                      call copy (vy(1,1,1,e),pb,nx1*ny1)
                      call mxm  (pr(1,1,1,e),nx2,iatjl2,ny2,pb,ny2)
                      call copy (pr(1,1,1,e),pb,nx2*ny2)
                  endif
              10 END DO
          endif
          if (ifheat) then
              ntotm1 = nx1*ny1*nelt
              do 15 ifldt=1,npscal+1
                  call copy (tax(1,1,1,ifldt),t(1,1,1,1,ifldt),ntotm1)
              15 END DO
              do 30 e=1,nelt
                  if (ifrzer(e)) then
                      do 25 ifldt=1,npscal+1
                          call mxm  (t(1,1,1,e,ifldt),nx1,iatjl1,ny1, &
                          pb,ny1)
                          call copy (t(1,1,1,e,ifldt),pb,nx1*ny1)
                      25 END DO
                  endif
              30 END DO
          endif
#endif
      endif
  !        Map the pressure onto the velocity mesh
  
      ntot1 = nx1*ny1*nz1*nelv
      nyz2  = ny2*nz2
      nxy1  = nx1*ny1
      nxyz  = nx1*ny1*nz1
      nxyz2 = nx2*ny2*nz2
  
      if (ifsplit) then
          call copy(pm1,pr,ntot1)
      elseif (if_full_pres) then
          pm1 = 0._dp
          do e=1,nelv
              call copy(pm1(1,1,1,e),pr(1,1,1,e),nxyz2)
          enddo
      else
          do 1000 e=1,nelv
              call mxm (ixm21,nx1,pr(1,1,1,e),nx2,pa(1,1,1),nyz2)
              do 100 iz=1,nz2
                  call mxm (pa(1,1,iz),nx1,iytm21,ny2,pb(1,1,iz),ny1)
              100 END DO
              call mxm (pb(1,1,1),nxy1,iztm21,nz2,pm1(1,1,1,e),nz1)
          1000 END DO
      endif

  else       ! map back
      if (ifaxis) then
        write(*,*) "Oops: ifaxis"
#if 0
          ntot1 = nx1*ny1*nelt
          call copy (ym1,yax,ntot1)
          if (ifflow) then
              ntot1 = nx1*ny1*nelv
              ntot2 = nx2*ny2*nelv
              call copy (vx,vxax,ntot1)
              call copy (vy,vyax,ntot1)
              call copy (pr,prax,ntot2)
          endif
          if (ifheat) then
              ntot1 = nx1*ny1*nelt
              do 3000 ifldt=1,npscal+1
                  call copy (t(1,1,1,1,ifldt),tax(1,1,1,ifldt),ntot1)
              3000 END DO
          endif
#endif
      endif

  endif
  return
end subroutine prepost_map

!-----------------------------------------------------------------------
!> \brief output .fld file
subroutine outfld(prefix, pm1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nid
  use input, only : param
  use tstep, only : istep, time
  implicit none

!   Work arrays and temporary arrays

  character(3), intent(in) ::    prefix

  real(DP), intent(in) ::  pm1    (lx1,ly1,lz1,lelv)

  real(DP) :: p66

  if(nid == 0) then
      WRITE(6,1001) istep,time
      1001 FORMAT(/,i9,1pe12.4,' Write checkpoint:')
  endif
  call nekgsync()

  p66 = abs(param(66))
  if (p66 == 6) then
      call mfo_outfld(prefix, pm1)
      call nekgsync                ! avoid race condition w/ outfld
      return
  endif

  write(*,*) "Oops: p66 /= 6"
#if 0
  ifxyo_s = ifxyo              ! Save ifxyo

  iprefix = i_find_prefix(prefix,99)

  ierr = 0
  if (nid == 0) then

  !       Open new file for each dump on /cfs
      nopen(iprefix)=nopen(iprefix)+1

      if (prefix == '   ' .AND. nopen(iprefix) == 1) ifxyo = .TRUE. ! 1st file

      if (prefix == 'rst' .AND. max_rst > 0) &
      nopen(iprefix) = mod1(nopen(iprefix),max_rst) ! restart

      call file2(nopen(iprefix),prefix)
      if (p66 < 1.0) then
          open(unit=24,file=fldfle,form='formatted',status='unknown')
      else
          fldfilei = 0
          len = ltrunc   (fldfle,131)
          call chcopy    (fldfile2,fldfle,len)
          call byte_open (fldfile2,ierr)
      !          write header as character string
          call blank(fhdfle,132)
      endif
  endif
  call bcast(ifxyo,lsize)
  if(p66 >= 1.0) &
  call err_chk(ierr,'Error opening file in outfld. Abort. $')

!   Figure out what goes in EXCODE
  CALL BLANK(EXCODE,30)
  NDUMPS=NDUMPS+1
  i=1
  if (mod(p66,1.0) == 0.0) then !old header format
      IF(IFXYO) then
          EXCODE(1)='X'
          EXCODE(2)=' '
          EXCODE(3)='Y'
          EXCODE(4)=' '
          i = 5
          IF(IF3D) THEN
              EXCODE(i)  ='Z'
              EXCODE(i+1)=' '
              i = i + 2
          ENDIF
      ENDIF
      IF(IFVO) then
          EXCODE(i)  ='U'
          EXCODE(i+1)=' '
          i = i + 2
      ENDIF
      IF(IFPO) THEN
          EXCODE(i)='P'
          EXCODE(i+1)=' '
          i = i + 2
      ENDIF
      IF(IFTO) THEN
          EXCODE(i)='T '
          EXCODE(i+1)=' '
          i = i + 1
      ENDIF
      do iip=1,ldimt1
          if (ifpsco(iip)) then
              write(excode(iip+I)  ,'(i1)') iip
              write(excode(iip+I+1),'(a1)') ' '
              i = i + 1
          endif
      enddo
  else
  ! ew header format
      IF (IFXYO) THEN
          EXCODE(i)='X'
          i = i + 1
      ENDIF
      IF (IFVO) THEN
          EXCODE(i)='U'
          i = i + 1
      ENDIF
      IF (IFPO) THEN
          EXCODE(i)='P'
          i = i + 1
      ENDIF
      IF (IFTO) THEN
          EXCODE(i)='T'
          i = i + 1
      ENDIF
      IF (LDIMT > 1) THEN
          NPSCALO = 0
          do k = 1,ldimt-1
              if(ifpsco(k)) NPSCALO = NPSCALO + 1
          enddo
          IF (NPSCALO > 0) THEN
              EXCODE(i) = 'S'
              WRITE(EXCODE(i+1),'(I1)') NPSCALO/10
              WRITE(EXCODE(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
          ENDIF
      ENDIF
  endif
       

!   Dump header
  ierr = 0
  if (nid == 0) call dump_header(excode,p66,ierr)
  call err_chk(ierr,'Error dumping header in outfld. Abort. $')

  call get_id(id)

  nxyz  = nx1*ny1*nz1

  ierr = 0
  do ieg=1,nelgt

      jnid = gllnid(ieg)
      ie   = gllel (ieg)

      if (nid == 0) then
          if (jnid == 0) then
              call fill_tmp(tdump,id,ie)
          else
              mtype=2000+ieg
              len=4*id*nxyz
              dum1=0.
              call csend (mtype,dum1,wdsize,jnid,nullpid)
              call crecv (mtype,tdump,len)
          endif
          if(ierr == 0) call out_tmp(id,p66,ierr)
      elseif (nid == jnid) then
          call fill_tmp(tdump,id,ie)
          dum1=0.

          mtype=2000+ieg
          len=4*id*nxyz
          call crecv (mtype,dum1,wdsize)
          call csend (mtype,tdump,len,node0,nullpid)
      endif
  enddo
  call err_chk(ierr,'Error writing file in outfld. Abort. $')

  ifxyo = ifxyo_s           ! restore ifxyo

  if (nid == 0) call close_fld(p66,ierr)
  call err_chk(ierr,'Error closing file in outfld. Abort. $')
#endif

  return
end subroutine outfld

!-----------------------------------------------------------------------
!> \brief output time history info.
subroutine outhis(ifhis, pm1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv
  use size_m, only : nx1, ny1, nz1, nelv, nid
  use geom, only : xm1, ym1, zm1
  use input, only : param, nhis, hcode, lochis, if3d, qinteg
  use parallel, only : gllnid, gllel, np, wdsize, node0, nullpid
  use soln, only : jp, vx, vy, vz, t
  use tstep, only : istep, dt, time
  implicit none

  real(DP) :: pm1    (lx1,ly1,lz1,lelv)

  real(DP) :: hdump(25)
  real(DP) :: xpart(10),ypart(10),zpart(10)
  logical :: ifhis

  integer, save :: icalld = 0
  integer :: iohis, nvar, ih, ii, ihisps, mtype, len, iobj, isk, iq
  integer :: ipart, i, iel, k, j, ip, kp, ielp, ix, iy, iz, ieg, jnid, ie
  real(DP) :: rmin, x, y, z, r, one
  real(DP), external :: glmax

  iohis=1
  if (param(52) >= 1) iohis = int(param(52))
  if (mod(istep,iohis) == 0 .AND. ifhis) then
      if (nhis > 0) then
          IPART=0
          DO 2100 I=1,NHIS
              IF(HCODE(10,I) == 'P')then       ! Do particle paths
                  if (ipart <= 10) ipart=ipart+1
                  if (istep == 0) then          ! Particle has original coordinates
                  !               Restarts?
                      XPART(IPART)= &
                      XM1(LOCHIS(1,I),LOCHIS(2,I),LOCHIS(3,I),LOCHIS(4,I))
                      YPART(IPART)= &
                      YM1(LOCHIS(1,I),LOCHIS(2,I),LOCHIS(3,I),LOCHIS(4,I))
                      ZPART(IPART)= &
                      ZM1(LOCHIS(1,I),LOCHIS(2,I),LOCHIS(3,I),LOCHIS(4,I))
                  ELSE
                  !               Kludge: Find Closest point
                      RMIN=1.0E7
                      ip = 1; jp = 1; kp = 1; ielp = 1
                      DO IEL=1,NELV
                          DO K=1,NZ1
                              DO J=1,NY1
                                  DO II=1,NX1
                                      X = XM1(II,J,K,IEL)
                                      Y = YM1(II,J,K,IEL)
                                      Z = ZM1(II,J,K,IEL)
                                      R=SQRT( (X-XPART(IPART))**2 + (Y-YPART(IPART))**2 &
                                      +       (Z-ZPART(IPART))**2 )
                                      IF(R < RMIN) then
                                          RMIN=R
                                          IP=II
                                          JP=J
                                          KP=K
                                          IELP=IEL
                                      ENDIF
                                  enddo
                              enddo
                          enddo
                      END DO
                      XPART(IPART) = XPART(IPART) + DT * VX(IP,JP,KP,IELP)
                      YPART(IPART) = YPART(IPART) + DT * VY(IP,JP,KP,IELP)
                      ZPART(IPART) = ZPART(IPART) + DT * VZ(IP,JP,KP,IELP)
                  ENDIF
              !            Dump particle data for history point first
              !            Particle data is Time, x,y,z.
                  WRITE(26,'(4G14.6,A10)')TIME,XPART(IPART),YPART(IPART) &
                  ,ZPART(IPART),'  Particle'
              ENDIF
          !         Figure out which fields to dump
              NVAR=0
              IF(HCODE(10,I) == 'H')then
              !           Do histories
              
              
                  IX =LOCHIS(1,I)
                  IY =LOCHIS(2,I)
                  IZ =LOCHIS(3,I)
                  IEG=LOCHIS(4,I)
                  JNID=GLLNID(IEG)
                  IE  =GLLEL(IEG)
              
              !------------------------------------------------------------------------
              !           On first call, write out XYZ location of history points
              
                  if (icalld == 0) then
                      one = glmax(one,1)           ! Force synch.  pff 04/16/04
                      IF (NID == JNID) then
                          IF (NP > 1 .AND. .NOT. IF3D) &
                          WRITE(6,22) NID,I,IX,IY,ie,IEG &
                          ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE)
                          IF (NP > 1 .AND. IF3D) &
                          WRITE(6,23) NID,I,IX,IY,IZ,ie,IEG,XM1(IX,IY,IZ,IE) &
                          ,YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
                          IF (NP == 1 .AND. .NOT. IF3D) &
                          WRITE(6,32) I,IX,IY,ie,IEG &
                          ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE)
                          IF (NP == 1 .AND. IF3D) &
                          WRITE(6,33) I,IX,IY,IZ,ie,IEG,XM1(IX,IY,IZ,IE) &
                          ,YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
                          22 FORMAT(i6,' History point:',I3,' at (',2(I2,','), &
                          &          2(I4,','),'); X,Y,Z = (',G12.4,',',G12.4,',',G12.4,').')
                          23 FORMAT(i6,' History point:',I3,' at (',3(I2,','), &
                          &          2(I4,','),'); X,Y,Z = (',G12.4,',',G12.4,',',G12.4,').')
                          32 FORMAT(2X,' History point:',I3,' at (',2(I2,','), &
                          &          2(I4,','),'); X,Y,Z = (',G12.4,',',G12.4,',',G12.4,').')
                          33 FORMAT(2X,' History point:',I3,' at (',3(I2,','), &
                          &          2(I4,','),'); X,Y,Z = (',G12.4,',',G12.4,',',G12.4,').')
                      ENDIF
                  ENDIF
              !------------------------------------------------------------------------
              
                  IF(HCODE(1,I) == 'U')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=VX(IX,IY,IZ,IE)
                  elseif(HCODE(1,I) == 'X')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=XM1(IX,IY,IZ,IE)
                  ENDIF
                  IF(HCODE(2,I) == 'V')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=VY(IX,IY,IZ,IE)
                  elseif(HCODE(2,I) == 'Y')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=YM1(IX,IY,IZ,IE)
                  ENDIF
                  IF(HCODE(3,I) == 'W')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=VZ(IX,IY,IZ,IE)
                  elseif(HCODE(3,I) == 'Z')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=ZM1(IX,IY,IZ,IE)
                  ENDIF
                  IF(HCODE(4,I) == 'P')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=PM1(IX,IY,IZ,IE)
                  ENDIF
                  IF(HCODE(5,I) == 'T')then
                      NVAR=NVAR+1
                      HDUMP(NVAR)=T (IX,IY,IZ,IE,1)
                  ENDIF
                  IF(HCODE(6,I) /= ' ' .AND. HCODE(6,I) /= '0') then
                      READ(HCODE(6,I),'(I1)',ERR=13)IHISPS
                      13 CONTINUE
                  !              Passive scalar data here
                      NVAR=NVAR+1
                      HDUMP(NVAR)=T (IX,IY,IZ,IE,IHISPS+1)
                  ENDIF
              
              !--------------------------------------------------------------
              !           Dump out the NVAR values for this history point
              !--------------------------------------------------------------
                  MTYPE=2200+I
                  LEN=WDSIZE*NVAR
              
              !           If point on this processor, send data to node 0
                  IF (NVAR > 0 .AND. NID /= 0 .AND. JNID == NID) &
                  call csend (mtype,hdump,len,node0,nullpid)
              
              !           If processor 0, recv data (unless point resides on node0).
                  IF (NVAR > 0 .AND. NID == 0 .AND. JNID /= NID) &
                  call crecv (mtype,hdump,len)
              
                  IF (NVAR > 0 .AND. NID == 0) &
                  WRITE(26,'(1p6e16.8)')TIME,(HDUMP(II),II=1,NVAR)
              
              !--------------------------------------------------------------
              !         End of history points
              !--------------------------------------------------------------
              ENDIF
          2100 END DO
      !        Now do Integrated quantities (Drag, Lift, Flux, etc.)
      !        Find out which groups are to be dumped
!max            IF (IFINTQ) CALL INTGLQ
          DO 2200 IH=1,NHIS
              IF(HCODE(10,IH) == 'I') then
                  IOBJ=LOCHIS(1,IH)
                  ISK=0
                  DO 2205 IQ=1,3
                      IF (HCODE(IQ,IH) /= ' ') ISK=ISK + 1
                  2205 END DO
                  DO 2207 IQ=5,6
                      IF (HCODE(IQ,IH) /= ' ') ISK=ISK + 1
                  2207 END DO
                  IF (NID == 0) &
                  WRITE(26,'(1p6e16.8)')TIME,(QINTEG(II,IOBJ),II=1,ISK)
              ENDIF
          2200 END DO
      ENDIF

  endif

  icalld = icalld+1

  return
end subroutine outhis

!=======================================================================
subroutine copyX4(a,b,n)
  use kinds, only : DP, r4
  implicit none
  integer, intent(in) :: n
  REAL(r4), intent(out) :: A(n)
  REAL(DP), intent(in) :: B(n)
  integer :: i
  DO 100 I = 1, N
      A(I) = real(B(I), kind=r4)
  100 END DO
  return
end subroutine copyX4

!=======================================================================
subroutine copy4r(a,b,n)
  use kinds, only : DP, r4
  implicit none
  real(DP) :: a(1)
  real(r4) :: b(1)
  integer :: n, i
  do i = 1, n
      a(i) = b(i)
  enddo
  return
end subroutine copy4r

!=======================================================================
integer function i_find_prefix(prefix,imax)
  implicit none

  character(3) :: prefix
  integer :: imax

  character(3), save :: prefixes(99)
  data prefixes /99*'...'/
  integer, save :: nprefix = 0
  integer :: i

!   Scan existing list of prefixes for a match to "prefix"
  do i=1,nprefix
      if (prefix == prefixes(i)) then
          i_find_prefix = i
          return
      endif
  enddo

!   If we're here, we didn't find a match.. bump list and return

  nprefix                = nprefix + 1
  prefixes(nprefix)      = prefix
  i_find_prefix          = nprefix

!   Array bounds check on prefix list

  if (nprefix > 99 .OR. nprefix > imax) then
      write(6,*) 'Hey! nprefix too big! ABORT in i_find_prefix' &
      ,nprefix,imax
      call exitt
  endif

  return
end function i_find_prefix

!-----------------------------------------------------------------------
!> \brief mult-file output
subroutine mfo_outfld(prefix, pm1) 
  use kinds, only : DP, i8
  use ctimer, only : dnekclock_sync
  use size_m, only : lx1, ly1, lz1, lelv, ldimt, lxo
  use size_m, only : nx1, ny1, nz1, nelt, nid, ndim
  use restart, only : nxo, nyo, nzo, nrg, iheadersize, pid0, nelb, wdsizo
  use restart, only : ifh_mbyte, nfileo
  use geom, only : xm1, ym1, zm1
  use input, only : ifxyo, ifxyo_, ifreguo, if3d, ifvo, ifpo, ifto, ifpsco
  use parallel, only : isize, nelgt, lsize
  use soln, only : vx, vy, vz, t
  use tstep, only : istep, time
  implicit none

  character(3), intent(in) :: prefix

  real(DP), intent(in) :: pm1 (lx1,ly1,lz1,lelv)  ! mapped pressure

  integer(i8) :: offs0,offs,stride,strideB,nxyzo8
  logical :: ifxyo_s

  real(DP) :: tiostart, dnbyte, tio
  real(DP), external :: glsum
  integer :: nout, ierr, ioflds, k

  tiostart=dnekclock_sync()

  ifxyo_s = ifxyo
  ifxyo_  = ifxyo
  nout = nelt
  nxo  = nx1
  nyo  = ny1
  nzo  = nz1
  if (ifreguo) then ! dump on regular (uniform) mesh
      if (nrg > lxo) then
          if (nid == 0) write(6,*) &
          'WARNING: nrg too large, reset to lxo!'
          nrg = lxo
      endif
      nxo  = nrg
      nyo  = nrg
      nzo  = 1
      if(if3d) nzo = nrg
  endif
  offs0 = iHeaderSize + 4 + isize*nelgt

  ierr=0
  if (nid == pid0) then
      call mfo_open_files(prefix,ierr)         ! open files on i/o nodes
  endif
  call err_chk(ierr,'Error opening file in mfo_open_files. $')
  call bcast(ifxyo_,lsize)
  ifxyo = ifxyo_
  call mfo_write_hdr                     ! create element mapping +

!   call exitti('this is wdsizo A:$',wdsizo)
! write hdr
  nxyzo8  = nxo*nyo*nzo
  strideB = nelB * nxyzo8*wdsizo
  stride  = nelgt* nxyzo8*wdsizo

  ioflds = 0
! dump all fields based on the t-mesh to avoid different
! topologies in the post-processor
  if (ifxyo) then
      offs = offs0 + ndim*strideB
      call byte_set_view(offs,ifh_mbyte)
      if (ifreguo) then
        write(*,*) "Oops: ifreguo"
#if 0
          call map2reg(ur1,nrg,xm1,nout)
          call map2reg(ur2,nrg,ym1,nout)
          if (if3d) call map2reg(ur3,nrg,zm1,nout)
          call mfo_outv(ur1,ur2,ur3,nout,nxo,nyo,nzo)
#endif
      else
          call mfo_outv(xm1,ym1,zm1,nout,nxo,nyo,nzo)
      endif
      ioflds = ioflds + ndim
  endif
  if (ifvo ) then
      offs = offs0 + ioflds*stride + ndim*strideB
      call byte_set_view(offs,ifh_mbyte)
      if (ifreguo) then
        write(*,*) "Oops: ifreguo"
#if 0
          call map2reg(ur1,nrg,vx,nout)
          call map2reg(ur2,nrg,vy,nout)
          if (if3d) call map2reg(ur3,nrg,vz,nout)
          call mfo_outv(ur1,ur2,ur3,nout,nxo,nyo,nzo)
#endif
      else
          call mfo_outv(vx,vy,vz,nout,nxo,nyo,nzo)  ! B-field handled thru outpost
      endif
      ioflds = ioflds + ndim
  endif
  if (ifpo ) then
      offs = offs0 + ioflds*stride + strideB
      call byte_set_view(offs,ifh_mbyte)
      if (ifreguo) then
        write(*,*) "Oops: ifreguo"
#if 0
          call map2reg(ur1,nrg,pm1,nout)
          call mfo_outs(ur1,nout,nxo,nyo,nzo)
#endif
      else
          call mfo_outs(pm1,nout,nxo,nyo,nzo)
      endif
      ioflds = ioflds + 1
  endif
  if (ifto ) then
      offs = offs0 + ioflds*stride + strideB
      call byte_set_view(offs,ifh_mbyte)
      if (ifreguo) then
        write(*,*) "Oops: ifreguo"
#if 0
          call map2reg(ur1,nrg,t,nout)
          call mfo_outs(ur1,nout,nxo,nyo,nzo)
#endif
      else
          call mfo_outs(t,nout,nxo,nyo,nzo)
      endif
      ioflds = ioflds + 1
  endif
  do k=1,ldimt-1
      if(ifpsco(k)) then
          offs = offs0 + ioflds*stride + strideB
          call byte_set_view(offs,ifh_mbyte)
          if (ifreguo) then
        write(*,*) "Oops: ifreguo"
#if 0
              call map2reg(ur1,nrg,t(1,1,1,1,k+1),nout)
              call mfo_outs(ur1,nout,nxo,nyo,nzo)
#endif
          else
              call mfo_outs(t(1,1,1,1,k+1),nout,nxo,nyo,nzo)
          endif
          ioflds = ioflds + 1
      endif
  enddo
  dnbyte = 1.*ioflds*nout*wdsizo*nxo*nyo*nzo

  if (if3d) then
      offs0   = offs0 + ioflds*stride
      strideB = nelB *2*4   ! min/max single precision
      stride  = nelgt*2*4
      ioflds  = 0
  ! add meta data to the end of the file
      if (ifxyo) then
          offs = offs0 + ndim*strideB
          call byte_set_view(offs,ifh_mbyte)
          call mfo_mdatav(xm1,ym1,zm1,nout)
          ioflds = ioflds + ndim
      endif
      if (ifvo ) then
          offs = offs0 + ioflds*stride + ndim*strideB
          call byte_set_view(offs,ifh_mbyte)
          call mfo_mdatav(vx,vy,vz,nout)
          ioflds = ioflds + ndim
      endif
      if (ifpo ) then
          offs = offs0 + ioflds*stride + strideB
          call byte_set_view(offs,ifh_mbyte)
          call mfo_mdatas(pm1,nout)
          ioflds = ioflds + 1
      endif
      if (ifto ) then
          offs = offs0 + ioflds*stride + strideB
          call byte_set_view(offs,ifh_mbyte)
          call mfo_mdatas(t,nout)
          ioflds = ioflds + 1
      endif
      do k=1,ldimt-1
          offs = offs0 + ioflds*stride + strideB
          call byte_set_view(offs,ifh_mbyte)
          if(ifpsco(k)) call mfo_mdatas(t(1,1,1,1,k+1),nout)
          ioflds = ioflds + 1
      enddo
      dnbyte = dnbyte + 2.*ioflds*nout*wdsizo
  endif

  ierr = 0
  if (nid == pid0) then
#ifdef MPIIO 
    call byte_close_mpi(ifh_mbyte,ierr)
#else 
    call byte_close(ierr)
#endif
  endif
  call err_chk(ierr,'Error closing file in mfo_outfld. Abort. $')

  tio = dnekclock_sync()-tiostart
  dnbyte = glsum(dnbyte,1)
  dnbyte = dnbyte + iHeaderSize + 4. + isize*nelgt
  dnbyte = dnbyte/1024/1024
  if(nid == 0) write(6,7) istep,time,dnbyte,dnbyte/tio, &
  nfileo
  7 format(/,i9,1pe12.4,' done :: Write checkpoint',/, &
  &        30X,'file size = ',3pG12.2,'MB',/, &
  &        30X,'avg data-throughput = ',0pf7.1,'MB/s',/, &
  &        30X,'io-nodes = ',i5,/)

  ifxyo = ifxyo_s ! restore old value

  return
end subroutine mfo_outfld

!-----------------------------------------------------------------------
subroutine io_init ! determine which nodes will output
  use kinds, only : DP
  use size_m, only : nid, lxo, nelt
  use input, only : param, ifreguo
  use parallel, only : np, wdsize
  use restart, only : ifdiro, nfileo, nproc_o
  use restart, only : fid0, pid0, pid1, wdsizo, nrg, nelB, pid00
  implicit none

  integer :: nn
  integer, external :: igl_running_sum
  real(DP), external :: glmin

  ifdiro = .FALSE. 

#ifdef MPIIO

#ifdef MPIIO_NOCOL
  nfileo  = abs(param(65))
  if(nfileo == 0) nfileo = 1
  if(np < nfileo) nfileo=np
  nproc_o = np / nfileo              !  # processors pointing to pid0
  fid0    = nid/nproc_o              !  file id
  pid0    = nproc_o*fid0             !  my parent i/o node
  pid1    = min(np-1,pid0+nproc_o-1) !  range of sending procs
  fid0    = 0
#else
  nfileo  = np
  nproc_o = 1
  fid0    = 0
  pid0    = nid
  pid1    = 0
#endif

#else
  if(param(65) < 0) ifdiro = .TRUE. !  p65 < 0 --> multi subdirectories
  nfileo  = int(abs(param(65)))
  if(nfileo == 0) nfileo = 1
  if(np < nfileo) nfileo=np
  nproc_o = np / nfileo              !  # processors pointing to pid0
  fid0    = nid/nproc_o              !  file id
  pid0    = nproc_o*fid0             !  my parent i/o node
  pid1    = min(np-1,pid0+nproc_o-1) !  range of sending procs
#endif

  call nek_comm_io(nfileo)

  wdsizo = 4                             ! every proc needs this
  if (param(63) > 0) wdsizo = 8         ! 64-bit .fld file
  if (wdsizo > wdsize) then
      if(nid == 0) write(6,*) 'ABORT: wdsizo > wdsize!'
      call exitt
  endif

  ifreguo = .FALSE.   ! by default we dump the data based on the GLL mesh
  nrg = lxo

! how many elements are present up to rank nid
  nn = nelt
  nelB = igl_running_sum(nn)
  nelB = nelB - nelt
       
  pid00 = int(glmin(pid0,1))

  return
end subroutine io_init

!-----------------------------------------------------------------------
subroutine mfo_open_files(prefix,ierr) ! open files
  use kinds, only : DP
  use input, only : ifreguo, ifxyo_, session
  use restart, only : max_rst, nfileo, ifdiro, fid0
  use string, only : ltrunc
  implicit none

  character(3) :: prefix
  integer :: ierr

  character(132) ::  fname
  character(1) ::   fnam1(132)
  equivalence  (fnam1,fname)

  character(6), save ::  six = "??????"
  character(6) :: str

  character(1), save :: slash = '/', dot = '.'

  integer, save :: nopen(99,2)
  data    nopen  / 198*0 /

  integer :: iprefix, nfld, k, len, ndigit
  integer, external :: i_find_prefix, mod1
  real(DP) :: rfileo

  call blank(fname,132)      !  zero out for byte_open()

  iprefix        = i_find_prefix(prefix,99)
  if (ifreguo) then
      nopen(iprefix,2) = nopen(iprefix,2)+1
      nfld             = nopen(iprefix,2)
  else
      nopen(iprefix,1) = nopen(iprefix,1)+1
      nfld             = nopen(iprefix,1)
  endif

  ! check for full-restart request
  if (prefix == 'rst' .AND. max_rst > 0) nfld = mod1(nfld,max_rst)

  call restart_nfld( nfld, prefix ) ! Check for Restart option.
  if (prefix == '   ' .AND. nfld == 1) ifxyo_ = .TRUE. ! 1st file

#ifdef MPIIO
  rfileo = 1
#else
  rfileo = nfileo
#endif
  ndigit = int(log10(rfileo) + 1)
       
  k = 1
  if (ifdiro) then                                  !  Add directory
      call chcopy(fnam1(1),'A',1)
      call chcopy(fnam1(2),six,ndigit)  ! put ???? in string
      k = 2 + ndigit
      call chcopy(fnam1(k),slash,1)
      k = k+1
  endif

  if (prefix(1:1) /= ' ' .AND. prefix(2:2) /= ' ' .AND.     & !  Add prefix
  prefix(3:3) /= ' ') then
  call chcopy(fnam1(k),prefix,3)
  k = k+3
  endif

  len=ltrunc(session,132)                           !  Add SESSION
  call chcopy(fnam1(k),session,len)
  k = k+len
       
  if (ifreguo) then
      len=4
      call chcopy(fnam1(k),'_reg',len)
      k = k+len
  endif

  call chcopy(fnam1(k),six,ndigit)                  !  Add file-id holder
  k = k + ndigit

  call chcopy(fnam1(k  ),dot,1)                     !  Add .f appendix
  call chcopy(fnam1(k+1),'f',1)
  k = k + 2

  write(str,4) nfld                                 !  Add nfld number
  4 format(i5.5)
  call chcopy(fnam1(k),str,5)
  k = k + 5

  call mbyte_open(fname,fid0,ierr)                  !  Open blah000.fnnnn
!    write(6,*) nid,fid0,' FILE:',fname
   
  return
end subroutine mfo_open_files

!-----------------------------------------------------------------------
!> \brief Check for Restart option and return proper nfld value.
!!   Also, convenient spot to explain restart strategy.
!!   The approach is as follows:
!!       Prefix rs4 would indicate 4 files in the restart cycle.
!!       This would be normal usage for velocity only, with
!!       checkpoints taking place in synch with standard io.
!!       The resultant restart sequence might look like:
!!       blah.fld09           Step 0
!!       rs4blah.fld01             1
!!       rs4blah.fld02             2
!!       which implies that fld09 would be used as the i.c.
!!       in the restart, rs4blah.fld01 would overwrite the
!!       solution at Step 1, and rs4blah.fld02 would overwrite
!!       Step 2.   Net result is that Steps 0-2 of the restart
!!       session have solutions identical to those computed in
!!       the prior run.   (It's important that both runs use
!!       the same dt in this case.)
!!       Another equally possible restart sequence would be:
!!       blah.fld10           Step 0
!!       rs4blah.fld03             1
!!       rs4blah.fld04             2
!!       Why the 3 & 4 ?   If one were to use only 1 & 2, there
!!       is a risk that the system crashes while writing, say,
!!       rs4blah.fld01, in which case the restart is compromised --
!!       very frustrating at the end of a run that has been queued
!!       for a week.  By providing a toggled sequence in pairs such as
!!       (1,2),   (3,4),  (1,2), ...
!!       ensures that one always has at least one complete restart
!!       sequence.   In the example above, the following files would
!!       be written, in order:
!!       :
!!       :
!!       blah.fld09
!!       rs4blah.fld01
!!       rs4blah.fld02
!!       blah.fld10
!!       rs4blah.fld03
!!       rs4blah.fld04
!!       blah.fld11
!!       rs4blah.fld01       (overwriting existing rs4blah.fld01)
!!       rs4blah.fld02       (    "           "        "  .fld02)
!!       blah.fld12
!!       rs4blah.fld03       (   etc.  )
!!       rs4blah.fld04
!!       :
!!       :
!!       Other strategies are possible, according to taste.
!!       Here is a data-intensive one:
!!       MHD + double-precision restart, but single-precision std files
!!       In this case, single-precision files are kept as the running
!!       file sequence (i.e., for later post-processing) but dbl-prec.
!!       is required for restart.  A total of 12 temporary restart files
!!       must be saved:  (3 for velocity, 3 for B-field) x 2 for redundancy.
!!       This is expressed, using hexadecimal notation (123456789abc...),
!!       as prefix='rsc'.
subroutine restart_nfld( nfld, prefix )
  use string, only : ltrunc, indx1
  implicit none
  character(3) :: prefix
  character(16), save :: kst = '0123456789abcdef' 
  character(1) ::  ks1(0:15),kin
  equivalence (ks1,kst)

  integer :: kfld, nfld, nfln
  integer, external :: mod1

  if (indx1(prefix,'rs',2) == 1) then
      read(prefix,3) kin
      3 format(2x,a1)
      do kfld=1,15
          if (ks1(kfld) == kin) goto 10
      enddo
      10 if (kfld == 16) kfld=4 ! std. default
      nfln = mod1(nfld,kfld) ! Restart A (1,2) and B (3,4)
      write(6,*) nfln,nfld,kfld,' kfld'
      nfld = nfln
  endif

  return
end subroutine restart_nfld

!-----------------------------------------------------------------------
subroutine outpost(v1,v2,v3,vp,vt,name3)
  use kinds, only : DP
  use input, only : ifto
  implicit none

  real(DP) :: v1(*),v2(*),v3(*),vp(*),vt(*)
  character(3) :: name3

  integer :: itmp
  itmp=0
  if (ifto) itmp=1
  call outpost2(v1,v2,v3,vp,vt,itmp,name3)

  return
end subroutine outpost

!-----------------------------------------------------------------------
subroutine outpost2(v1,v2,v3,vp,vt,nfldt,name3)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, lelv, ldimt
  use size_m, only : lx2, ly2, lz2
  use size_m, only : nx1, ny1, nz1, nx2, ny2, nz2, nelv, nelt
  use input, only : ifto, ifpsco
  use soln, only : vx, vy, vz, pr, t
  implicit none

  integer, parameter :: ltot1=lx1*ly1*lz1*lelt
  integer, parameter :: ltot2=lx2*ly2*lz2*lelv
  real(DP), allocatable :: w1(:), w2(:), w3(:), wp(:), wt(:,:)
  real(DP) :: v1(1),v2(1),v3(1),vp(1),vt(ltot1,1)
  character(3) :: name3
  logical :: if_save(ldimt)

  integer :: ntot1, ntot1t, ntot2, nfldt, i

  allocate(w1(ltot1),w2(ltot1),w3(ltot1),wp(ltot2),wt(ltot1,ldimt))

  ntot1  = nx1*ny1*nz1*nelv
  ntot1t = nx1*ny1*nz1*nelt
  ntot2  = nx2*ny2*nz2*nelv

  if(nfldt > ldimt) then
      write(6,*) 'ABORT: outpost data too large (nfldt>ldimt)!'
      call exitt
  endif

! store solution
  call copy(w1,vx,ntot1)
  call copy(w2,vy,ntot1)
  call copy(w3,vz,ntot1)
  call copy(wp,pr,ntot2)
  do i = 1,nfldt
      call copy(wt(1,i),t(1,1,1,1,i),ntot1t)
  enddo

! swap with data to dump
  call copy(vx,v1,ntot1)
  call copy(vy,v2,ntot1)
  call copy(vz,v3,ntot1)
  call copy(pr,vp,ntot2)
  do i = 1,nfldt
      call copy(t(1,1,1,1,i),vt(1,i),ntot1t)
  enddo

! dump data
  if_save(1) = ifto
  ifto = .FALSE. 
  if(nfldt > 0) ifto = .TRUE. 
  do i = 1,ldimt-1
      if_save(i+1) = ifpsco(i)
      ifpsco(i) = .FALSE. 
      if(i+1 <= nfldt) ifpsco(i) = .TRUE. 
  enddo

  call prepost( .TRUE. ,name3)

  ifto = if_save(1)
  do i = 1,ldimt-1
      ifpsco(i) = if_save(i+1)
  enddo

! restore solution data
  call copy(vx,w1,ntot1)
  call copy(vy,w2,ntot1)
  call copy(vz,w3,ntot1)
  call copy(pr,wp,ntot2)
  do i = 1,nfldt
      call copy(t(1,1,1,1,i),wt(1,i),ntot1t)
  enddo

  return
end subroutine outpost2

!-----------------------------------------------------------------------
subroutine mfo_mdatav(u,v,w,nel)
  use kinds, only : DP, r4
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1, ndim, nelt, nid, lelt
  use input, only : if3d
  use restart, only : pid0, pid1
  implicit none

  real(DP), intent(in) :: u(lx1*ly1*lz1,*),v(lx1*ly1*lz1,*),w(lx1*ly1*lz1,*)
  integer, intent(in) :: nel

  real(r4) :: buffer(1+6*lelt)

  integer :: e, inelp, mtype, k, idum, nout, j, ierr, leo, len, nxyz, n

  call nekgsync() ! clear outstanding message queues.

  nxyz = nx1*ny1*nz1
  n    = 2*ndim
  len  = 4 + 4*(n*lelt)   ! recv buffer size
  leo  = 4 + 4*(n*nelt)
  ierr = 0

! Am I an I/O node?
  if (nid == pid0) then
      j = 1
      do e=1,nel
          buffer(j+0) = real(minval(u(:,e)), kind=r4)
          buffer(j+1) = real(maxval(u(:,e)), kind=r4)
          buffer(j+2) = real(minval(v(:,e)), kind=r4)
          buffer(j+3) = real(maxval(v(:,e)), kind=r4)
          j = j + 4
          if(if3d) then
              buffer(j+0) = real(minval(w(:,e)), kind=r4)
              buffer(j+1) = real(maxval(w(:,e)), kind=r4)
              j = j + 2
          endif
      enddo

  ! write out my data
      nout = n*nel
      if(ierr == 0) then
#ifdef MPIIO 
      call byte_write_mpi(buffer,nout,-1,ifh_mbyte,ierr)
#else 
      call byte_write(buffer,nout,ierr)
#endif
      endif

  ! write out the data of my childs
      idum  = 1
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)           ! handshake
          call crecv(mtype,buffer,len)
          inelp = int(buffer(1))
          nout  = n*inelp
          if(ierr == 0) then
#ifdef MPIIO
          call byte_write_mpi(buffer(2),nout,-1,ifh_mbyte,ierr)
#else
          call byte_write(buffer(2),nout,ierr)
#endif
          endif
      enddo
  else
      j = 1
      buffer(j) = nel
      j = j + 1
      do e=1,nel
          buffer(j+0) = real(minval(u(:,e)), kind=r4)
          buffer(j+1) = real(maxval(u(:,e)), kind=r4)
          buffer(j+2) = real(minval(v(:,e)), kind=r4)
          buffer(j+3) = real(maxval(v(:,e)), kind=r4)
          j = j + 4
          if(n == 6) then
              buffer(j+0) = real(minval(w(:,e)), kind=r4)
              buffer(j+1) = real(maxval(w(:,e)), kind=r4)
              j = j + 2
          endif
      enddo

  ! send my data to my pararent I/O node
      mtype = nid
      call crecv(mtype,idum,4)                ! hand-shake
      call csend(mtype,buffer,leo,pid0,0)     ! u4 :=: u8
  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_mdatav. $')

  return
end subroutine mfo_mdatav

!-----------------------------------------------------------------------
subroutine mfo_mdatas(u,nel)
  use kinds, only : DP, r4
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, nelt, nid
  use restart, only : pid0, pid1
  implicit none

  real(DP), intent(in) :: u(lx1*ly1*lz1,*)
  integer, intent(in) :: nel

  real(r4) :: buffer(1+2*lelt)

  integer :: e, inelp, mtype, k, idum, nout, j, ierr, leo, len, n, nxyz

  call nekgsync() ! clear outstanding message queues.

  nxyz = nx1*ny1*nz1
  n    = 2
  len  = 4 + 4*(n*lelt)    ! recv buffer size
  leo  = 4 + 4*(n*nelt)
  ierr = 0

! Am I an I/O node?
  if (nid == pid0) then
      j = 1
      do e=1,nel
          buffer(j+0) = real(minval(u(:,e)), kind=r4)
          buffer(j+1) = real(maxval(u(:,e)), kind=r4)
          j = j + 2
      enddo

  ! write out my data
      nout = n*nel
      if(ierr == 0) then
#ifdef MPIIO
      call byte_write_mpi(buffer,nout,-1,ifh_mbyte,ierr)
#else
      call byte_write(buffer,nout,ierr)
#endif
      endif
  ! write out the data of my childs
      idum  = 1
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)           ! handshake
          call crecv(mtype,buffer,len)
          inelp = int(buffer(1))
          nout  = n*inelp
          if(ierr == 0) then
#ifdef MPIIO
          call byte_write_mpi(buffer(2),nout,-1,ifh_mbyte,ierr)
#else
          call byte_write(buffer(2),nout,ierr)
#endif
          endif
      enddo
  else
      j = 1
      buffer(j) = nel
      j = j + 1
      do e=1,nel
          buffer(j+0) = real(minval(u(:,e)), kind=r4)
          buffer(j+1) = real(maxval(u(:,e)), kind=r4)
          j = j + 2
      enddo

  ! send my data to my pararent I/O node
      mtype = nid
      call crecv(mtype,idum,4)                ! hand-shake
      call csend(mtype,buffer,leo,pid0,0)     ! u4 :=: u8
  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_mdatas. $')

  return
end subroutine mfo_mdatas

!-----------------------------------------------------------------------
!> \brief output a scalar field
subroutine mfo_outs(u,nel,mx,my,mz)
  use kinds, only : DP, r4
  use size_m, only : nid, lelt, lxo
  use restart, only : wdsizo, pid0, pid1
  implicit none

  integer, intent(in) :: nel, mx, my, mz
  real(DP), intent(in) :: u(mx,my,mz,1)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, len, leo, ntot, idum, ierr, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  len  = 8 + 8*(lelt*nxyz)  ! recv buffer size
  leo  = 8 + wdsizo*(nel*nxyz)
  ntot = nxyz*nel

  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*2*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*1*lelt))
  endif

  if (nid == pid0) then

      if (wdsizo == 4) then             ! 32-bit output
          call copyx4 (u4,u,ntot)
      else
          call copy   (u8,u,ntot)
      endif
      if(wdsizo == 4 .and. ierr == 0) then
      nout = wdsizo/4 * ntot
#ifdef MPIIO 
      call byte_write_mpi(u4,nout,-1,ifh_mbyte,ierr)
#else 
      call byte_write(u4,nout,ierr)          ! u4 :=: u8
#endif
      elseif(ierr == 0) then
      nout = wdsizo/4 * ntot
#ifdef MPIIO 
      call byte_write_mpi(u8,nout,-1,ifh_mbyte,ierr)
#else 
      call byte_write(u8,nout,ierr)          ! u4 :=: u8
#endif
      endif

  ! write out the data of my childs
      idum  = 1
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)       ! handshake
          if (wdsizo == 4 .AND. ierr == 0) then
            call crecv(mtype,u4,len)
            nout  = wdsizo/4 * nxyz * int(u4(1))
#ifdef MPIIO
            call byte_write_mpi(u4(3),nout,-1,ifh_mbyte,ierr)
#else
            call byte_write(u4(3),nout,ierr)
#endif
          elseif(ierr == 0) then
            call crecv(mtype,u8,len)
            nout  = wdsizo/4 * nxyz * int(u8(1))
#ifdef MPIIO
            call byte_write_mpi(u8(2),nout,-1,ifh_mbyte,ierr)
#else
            call byte_write(u8(2),nout,ierr)
#endif
          endif
      enddo

  else

      if (wdsizo == 4) then             ! 32-bit output
        u4(1)= nel
        call copyx4 (u4(3),u,ntot)
        mtype = nid
        call crecv(mtype,idum,4)            ! hand-shake
        call csend(mtype,u4,leo,pid0,0)     ! u4 :=: u8
      else
        u8(1)= nel
        call copy   (u8(2),u,ntot)
        mtype = nid
        call crecv(mtype,idum,4)            ! hand-shake
        call csend(mtype,u8,leo,pid0,0)     ! u4 :=: u8
      endif

  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_outs. $')

  return
end subroutine mfo_outs

!-----------------------------------------------------------------------
!> \brief output a vector field
subroutine mfo_outv(u,v,w,nel,mx,my,mz) 
  use kinds, only : DP, r4
  use size_m, only : nid, ndim, lxo, lelt
  use input, only : if3d
  use restart, only : wdsizo, pid0, pid1
  implicit none
 
  integer, intent(in) :: mx, my, mz
  real(DP), intent(in) :: u(mx*my*mz,*),v(mx*my*mz,*),w(mx*my*mz,*)

  real(r4), allocatable :: u4(:)
  real(DP), allocatable :: u8(:)

  integer :: nxyz, len, leo, nel, idum, ierr
  integer :: j, iel, nout, k, mtype

  call nekgsync() ! clear outstanding message queues.
  if(mx > lxo .OR. my > lxo .OR. mz > lxo) then
      if(nid == 0) write(6,*) 'ABORT: lxo too small'
      call exitt
  endif

  nxyz = mx*my*mz
  len  = 8 + 8*(lelt*nxyz*ndim)   ! recv buffer size (u4)
  leo  = 8 + wdsizo*(nel*nxyz*ndim)
  idum = 1
  ierr = 0

  if (wdsizo == 4) then
    allocate(u4(2+lxo*lxo*lxo*6*lelt))
  else
    allocate(u8(1+lxo*lxo*lxo*3*lelt))
  endif
  
  if (nid == pid0) then
      j = 0
      if (wdsizo == 4) then             ! 32-bit output
          do iel = 1,nel
              call copyx4   (u4(j+1),u(1,iel),nxyz)
              j = j + nxyz
              call copyx4   (u4(j+1),v(1,iel),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copyx4 (u4(j+1),w(1,iel),nxyz)
                  j = j + nxyz
              endif
          enddo
      else
          do iel = 1,nel
              call copy     (u8(j+1),u(1,iel),nxyz)
              j = j + nxyz
              call copy     (u8(j+1),v(1,iel),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copy   (u8(j+1),w(1,iel),nxyz)
                  j = j + nxyz
              endif
          enddo
      endif
      nout = wdsizo/4 * ndim*nel * nxyz
      if (wdsizo == 4 .and. ierr == 0) then
#ifdef MPIIO 
      call byte_write_mpi(u4,nout,-1,ifh_mbyte,ierr)
#else 
      call byte_write(u4,nout,ierr)          ! u4 :=: u8
#endif
      elseif (ierr == 0) then
#ifdef MPIIO 
      call byte_write_mpi(u8,nout,-1,ifh_mbyte,ierr)
#else 
      call byte_write(u8,nout,ierr)          ! u4 :=: u8
#endif
      endif
  ! write out the data of my childs
      do k=pid0+1,pid1
          mtype = k
          call csend(mtype,idum,4,k,0)           ! handshake

          if (wdsizo == 4 .AND. ierr == 0) then
              call crecv(mtype,u4,len)
              nout  = wdsizo/4 * ndim*nxyz * int(u4(1))
#ifdef MPIIO
              call byte_write_mpi(u4(3),nout,-1,ifh_mbyte,ierr)
#else
              call byte_write(u4(3),nout,ierr)
#endif
          elseif(ierr == 0) then
              call crecv(mtype,u8,len)
              nout  = wdsizo/4 * ndim*nxyz * int(u8(1))
#ifdef MPIIO
              call byte_write_mpi(u8(2),nout,-1,ifh_mbyte,ierr)
#else
              call byte_write(u8(2),nout,ierr)
#endif
          endif
      enddo
  else

      if (wdsizo == 4) then             ! 32-bit output
          u4(1) = nel
          j = 2
          do iel = 1,nel
              call copyx4   (u4(j+1),u(1,iel),nxyz)
              j = j + nxyz
              call copyx4   (u4(j+1),v(1,iel),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copyx4 (u4(j+1),w(1,iel),nxyz)
                  j = j + nxyz
              endif
          enddo
        mtype = nid
        call crecv(mtype,idum,4)            ! hand-shake
        call csend(mtype,u4,leo,pid0,0)     ! u4 :=: u8
      else
          u8(1) = nel
          j = 1
          do iel = 1,nel
              call copy     (u8(j+1),u(1,iel),nxyz)
              j = j + nxyz
              call copy     (u8(j+1),v(1,iel),nxyz)
              j = j + nxyz
              if(if3d) then
                  call copy   (u8(j+1),w(1,iel),nxyz)
                  j = j + nxyz
              endif
          enddo
        mtype = nid
        call crecv(mtype,idum,4)            ! hand-shake
        call csend(mtype,u8,leo,pid0,0)     ! u4 :=: u8
      endif
  endif

  call err_chk(ierr,'Error writing data to .f00 in mfo_outv. $')
  return
end subroutine mfo_outv

!-----------------------------------------------------------------------
!> \brief write hdr, byte key, els.
subroutine mfo_write_hdr          
  use kinds, only : r4
  use size_m, only : nid, nelt, lelt, ldimt
  use input, only : ifxyo, ifvo, ifpo, ifto, ifpsco
  use parallel, only : nelgt, lglel
  use restart, only : nfileo, pid0, pid1, rdcode1, wdsizo, nxo, nyo, nzo
  use restart, only : fid0, iheadersize
  use tstep, only : istep, time
  implicit none

  real(r4) :: test_pattern
  real(r4), allocatable :: padding(:)
  integer :: lglist(0:lelt)

  character(132) :: hdr

  integer :: idum, nfileoo, nelo, j, mtype, inelp, ierr, i, npscalo, k
  integer :: ibsw_out, len
  integer :: pad_size

  call nekgsync()
  idum = 1

#ifdef MPIIO
  nfileoo = 1   ! all data into one file
  nelo = nelgt
#else
  nfileoo = nfileo
  if(nid == pid0) then                ! how many elements to dump
      nelo = nelt
      do j = pid0+1,pid1
          mtype = j
          call csend(mtype,idum,4,j,0)   ! handshake
          call crecv(mtype,inelp,4)
          nelo = nelo + inelp
      enddo
  else
      mtype = nid
      call crecv(mtype,idum,4)          ! hand-shake
      call csend(mtype,nelt,4,pid0,0)   ! u4 :=: u8
  endif
#endif

  ierr = 0
  if(nid == pid0) then

      call blank(hdr,132)              ! write header
      call blank(rdcode1,10)
      i = 1
      IF (IFXYO) THEN
          rdcode1(i)='X'
          i = i + 1
      ENDIF
      IF (IFVO) THEN
          rdcode1(i)='U'
          i = i + 1
      ENDIF
      IF (IFPO) THEN
          rdcode1(i)='P'
          i = i + 1
      ENDIF
      IF (IFTO) THEN
          rdcode1(i)='T'
          i = i + 1
      ENDIF
      IF (LDIMT > 1) THEN
          NPSCALO = 0
          do k = 1,ldimt-1
              if(ifpsco(k)) NPSCALO = NPSCALO + 1
          enddo
          IF (NPSCALO > 0) THEN
              rdcode1(i) = 'S'
              WRITE(rdcode1(i+1),'(I1)') NPSCALO/10
              WRITE(rdcode1(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
          ENDIF
      ENDIF
       
      write(hdr,1) wdsizo,nxo,nyo,nzo,nelo,nelgt,time,istep,fid0,nfileoo &
      ,   (rdcode1(i),i=1,10)        ! 74+20=94
      1 format('#max',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13, &
      &        1x,i9,1x,i6,1x,i6,1x,10a)

  ! if we want to switch the bytes for output
  ! switch it again because the hdr is in ASCII
      call get_bytesw_write(ibsw_out)
  !      if (ibsw_out.ne.0) call set_bytesw_write(ibsw_out)
      if (ibsw_out /= 0) call set_bytesw_write(0)

      test_pattern = 6.54321_r4           ! write test pattern for byte swap

      pad_size = (8 * (2**20) - (iHeaderSize + 4) ) / 4
      allocate(padding(pad_size)); padding = 0.
#ifdef MPIIO
  ! only rank0 (pid00) will write hdr + test_pattern
      call byte_write_mpi(hdr,iHeaderSize/4,pid00,ifh_mbyte,ierr)
      call byte_write_mpi(test_pattern,1,pid00,ifh_mbyte,ierr)
      call byte_write_mpi(padding, pad_size, pid00, ifh_mbyte, ierr)
#else
      call byte_write(hdr,iHeaderSize/4,ierr)
      call byte_write(test_pattern,1,ierr)
      ! pad up to 8MB
      call byte_write(padding, pad_size, ierr)
#endif
      deallocate(padding)

  endif

  call err_chk(ierr,'Error writing header in mfo_write_hdr. $')

! write global element numbering for this group
  if(nid == pid0) then
#ifdef MPIIO
      ioff = iHeaderSize + 4 + 4*pad_size + nelB*isize
      call byte_set_view (ioff,ifh_mbyte)
      call byte_write_mpi(lglel,nelt,-1,ifh_mbyte,ierr)
#else
      call byte_write(lglel,nelt,ierr)
#endif
      pad_size = -nelt
      do j = pid0+1,pid1
          mtype = j
          call csend(mtype,idum,4,j,0)   ! handshake
          len = 4*(lelt+1)
          call crecv(mtype,lglist,len)
          if(ierr == 0) then
#ifdef MPIIO
              call byte_write_mpi(lglist(1),lglist(0),-1,ifh_mbyte,ierr)
#else
              call byte_write(lglist(1),lglist(0),ierr)
#endif
          pad_size = pad_size - lglist(0)
          endif
      enddo

    ! pad up to 8MB
    do while (pad_size < 0) 
      pad_size = pad_size + (8 * (2**20)) / 4
    enddo
    allocate(padding(pad_size)); padding = 0.
#ifdef MPIIO
    call byte_write_mpi(padding, pad_size, pid00, ifh_mbyte, ierr)
#else
    call byte_write(padding, pad_size, ierr)
#endif
    deallocate(padding)

  else
      mtype = nid
      call crecv(mtype,idum,4)          ! hand-shake
              
      lglist(0) = nelt
      lglist(1:nelt) = lglel(1:nelt)

      len = 4*(nelt+1)
      call csend(mtype,lglist,len,pid0,0)
  endif

 

  call err_chk(ierr,'Error writing global nums in mfo_write_hdr$')
  return
end subroutine mfo_write_hdr

!-----------------------------------------------------------------------
