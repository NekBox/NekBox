!-----------------------------------------------------------------------
!> \brief Subroutine to initialize logical flags
SUBROUTINE SETLOG()
  use size_m, only : ndim, nelv, nelt, nfield, nid, lxd
  use ctimer, only : ifsync
  use geom, only : ifvcor, ifgeom, ifsurt, ifwcno, ifeppm, ifqinp, ifmelt
  use input, only : param, ifadvc, iftmsh, ifneknek, ifmoab, ifkeps, ifmodel
  use input, only : ifnonl, cbc, ifheat, ifmvbd, ifflow, ifnav, ifstrs, ifprint
  use input, only : ifaxis, ifcyclic, ifchar, ifusermv, ifuservp, iflomach
  use input, only : ifsplit, iftran, nobj, lochis, ifintq, hcode, nhis
  use tstep, only : ifield, nelfld
  use turbo, only : ifswall, ifcwuz
  implicit none

  CHARACTER(3) CB
  LOGICAL :: IFALGN,IFNORX,IFNORY,IFNORZ

  integer :: nface, nmxv, nmxt, iel, ifc, iq, ih, iobj

  NFACE  = 2*NDIM
  NMXV   = NFACE*NELV
  NMXT   = NFACE*NELT

  IFPRINT = .TRUE. 
  IFVCOR  = .TRUE. 
  IFGEOM  = .FALSE. 
  IFINTQ  = .FALSE. 
  IFSURT  = .FALSE. 
  IFWCNO  = .FALSE. 
  IFSWALL = .FALSE. 
  DO IFIELD=1,NFIELD
      IFNONL(IFIELD) = .FALSE. 
  END DO

  CALL LFALSE (IFEPPM,NMXV)
  CALL LFALSE (IFQINP,NMXV)

!max    IF (IFMODEL) CALL SETSHL

  IF (IFMVBD) THEN
      IFGEOM = .TRUE. 
      IF ( IFFLOW .AND. .NOT. IFNAV )       IFWCNO          = .TRUE. 
      IF ( IFMELT .AND. .NOT. IFFLOW )      IFWCNO          = .TRUE. 
  ENDIF

  IF (IFFLOW) THEN

  ! k         call check_cyclic  ! fow now; set in .rea file

      IFIELD = 1
      DO IEL=1,NELV
          DO IFC=1,NFACE
              CB = CBC(IFC,IEL,IFIELD)
              CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
              IF ( .NOT. IFSTRS ) CALL CHKCBC  (CB,IEL,IFC,IFALGN)
              IF  (CB == 'O  ' .OR. CB == 'o  ' .OR. &
              CB == 'ON ' .OR. CB == 'on ' .OR. &
              CB == 'S  ' .OR. CB == 's  ' .OR. &
              CB == 'SL ' .OR. CB == 'sl ' .OR. &
              CB == 'MM ' .OR. CB == 'mm ' .OR. &
              CB == 'MS ' .OR. CB == 'ms ')  THEN
                  IFVCOR          = .FALSE. 
                  IFEPPM(IFC,IEL) = .TRUE. 
              ENDIF
              IF  (CB == 'VL ' .OR. CB == 'vl ' .OR. &
              CB == 'WSL' .OR. CB == 'wsl' .OR. &
              CB == 'SL ' .OR. CB == 'sl ' .OR. &
              CB == 'SHL' .OR. CB == 'shl' .OR. &
              CB == 'MM ' .OR. CB == 'mm ' .OR. &
              CB == 'MS ' .OR. CB == 'ms ' .OR. &
              CB == 'O  ' .OR. CB == 'o  ' .OR. &
              CB == 'ON ' .OR. CB == 'on ')  THEN
                  IFQINP(IFC,IEL) = .TRUE. 
              ENDIF
              IF  (CB == 'MS ' .OR. CB == 'ms ' .OR. &
              CB == 'MM ' .OR. CB == 'mm ' .OR. &
              CB == 'MSI' .OR. CB == 'msi' ) THEN
                  IFSURT          = .TRUE. 
              ENDIF
              IF  (CB == 'WS ' .OR. CB == 'ws ' .OR. &
              CB == 'WSL' .OR. CB == 'wsl') THEN
                  IFSWALL         = .TRUE. 
                  IFCWUZ          = .TRUE. 
              ENDIF
          END DO
      enddo
  ENDIF

  IF (IFHEAT) THEN
  
      DO IFIELD=2,NFIELD
          DO IEL=1,NELFLD(IFIELD)
              DO IFC=1,NFACE
                  CB=CBC(IFC,IEL,IFIELD)
                  IF  (CB == 'r  ' .OR. CB == 'R  ') THEN
                      IFNONL(IFIELD)  = .TRUE. 
                  ENDIF
              enddo
          enddo
      END DO
  
  ENDIF

!max    if (ifmhd) call set_ifbcor

  IF (NHIS > 0) THEN
      IQ = 0
      DO IH=1, NHIS
          IF ( HCODE(10,IH) == 'I' ) THEN
              IFINTQ = .TRUE. 
              IOBJ   = LOCHIS(1,IH)
              IQ     = IQ + 1
              IF (IOBJ > NOBJ .OR. IOBJ < 0)  THEN
                  WRITE (6,*) &
                  'ERROR : Undefined Object for integral',IQ
                  call exitt
              ENDIF
          ENDIF
      END DO
  ENDIF

!   Establish global consistency of LOGICALS amongst all processors.

  CALL GLLOG(IFVCOR , .FALSE. )
  CALL GLLOG(IFSURT , .TRUE. )
  CALL GLLOG(IFSWALL, .TRUE. )
  CALL GLLOG(IFCWUZ , .TRUE. )
  CALL GLLOG(IFWCNO , .TRUE. )
  DO IFIELD=2,NFIELD
      CALL GLLOG(IFNONL(IFIELD), .TRUE. )
  END DO

  IF (NID == 0) THEN
      WRITE (6,*) 'IFTRAN   =',IFTRAN
      WRITE (6,*) 'IFFLOW   =',IFFLOW
      WRITE (6,*) 'IFHEAT   =',IFHEAT
      WRITE (6,*) 'IFSPLIT  =',IFSPLIT
      WRITE (6,*) 'IFLOMACH =',IFLOMACH
      WRITE (6,*) 'IFUSERVP =',IFUSERVP
      WRITE (6,*) 'IFUSERMV =',IFUSERMV
      WRITE (6,*) 'IFSTRS   =',IFSTRS
      WRITE (6,*) 'IFCHAR   =',IFCHAR
      WRITE (6,*) 'IFCYCLIC =',IFCYCLIC
      WRITE (6,*) 'IFAXIS   =',IFAXIS
      WRITE (6,*) 'IFMVBD   =',IFMVBD
      WRITE (6,*) 'IFMELT   =',IFMELT
      WRITE (6,*) 'IFMODEL  =',IFMODEL
      WRITE (6,*) 'IFKEPS   =',IFKEPS
      WRITE (6,*) 'IFMOAB   =',IFMOAB
      WRITE (6,*) 'IFNEKNEK =',IFNEKNEK
      WRITE (6,*) 'IFSYNC   =',IFSYNC
      WRITE (6,*) '  '
      WRITE (6,*) 'IFVCOR   =',IFVCOR
      WRITE (6,*) 'IFINTQ   =',IFINTQ
      WRITE (6,*) 'IFCWUZ   =',IFCWUZ
      WRITE (6,*) 'IFSWALL  =',IFSWALL
      WRITE (6,*) 'IFGEOM   =',IFGEOM
      WRITE (6,*) 'IFSURT   =',IFSURT
      WRITE (6,*) 'IFWCNO   =',IFWCNO
      DO IFIELD=1,NFIELD
          WRITE (6,*) '  '
          WRITE (6,*) 'IFTMSH for field',IFIELD,'   = ',IFTMSH(IFIELD)
          WRITE (6,*) 'IFADVC for field',IFIELD,'   = ',IFADVC(IFIELD)
          WRITE (6,*) 'IFNONL for field',IFIELD,'   = ',IFNONL(IFIELD)
      END DO
      WRITE (6,*) '  '
      if (param(99) > 0) write(6,*) 'Dealiasing enabled, lxd=', lxd
  ENDIF

  RETURN
END SUBROUTINE SETLOG

!-----------------------------------------------------------------------
!> \brief Check direction of normal of an element face for
!!  alignment with the X, Y, or Z axis.
SUBROUTINE CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
  use kinds, only : DP
  use size_m, only : ndim, nx1, ny1
  use geom, only : unx, uny, unz
  implicit none

  LOGICAL :: IFALGN,IFNORX,IFNORY,IFNORZ
  integer :: ifc, iel

  integer :: ix, iy, ncpf
  real(DP) :: sumx, sumy, sumz, tolnor

  SUMX    = 0.0
  SUMY    = 0.0
  SUMZ    = 0.0
  TOLNOR  = 1.0e-3
  IFALGN  = .FALSE. 
  IFNORX  = .FALSE. 
  IFNORY  = .FALSE. 
  IFNORZ  = .FALSE. 

  IF (NDIM == 2) THEN
  
      NCPF = NX1
      DO IX=1,NX1
          SUMX = SUMX + ABS( ABS(UNX(IX,1,IFC,IEL)) - 1.0 )
          SUMY = SUMY + ABS( ABS(UNY(IX,1,IFC,IEL)) - 1.0 )
      END DO
      SUMX = SUMX / NCPF
      SUMY = SUMY / NCPF
      IF ( SUMX < TOLNOR ) THEN
          IFNORX  = .TRUE. 
          IFALGN = .TRUE. 
      ENDIF
      IF ( SUMY < TOLNOR ) THEN
          IFNORY  = .TRUE. 
          IFALGN = .TRUE. 
      ENDIF
  
  ELSE
  
      NCPF = NX1*NX1
      DO IX=1,NX1
          DO IY=1,NY1
              SUMX = SUMX + ABS( ABS(UNX(IX,IY,IFC,IEL)) - 1.0 )
              SUMY = SUMY + ABS( ABS(UNY(IX,IY,IFC,IEL)) - 1.0 )
              SUMZ = SUMZ + ABS( ABS(UNZ(IX,IY,IFC,IEL)) - 1.0 )
          enddo
      END DO
      SUMX = SUMX / NCPF
      SUMY = SUMY / NCPF
      SUMZ = SUMZ / NCPF
      IF ( SUMX < TOLNOR ) THEN
          IFNORX  = .TRUE. 
          IFALGN = .TRUE. 
      ENDIF
      IF ( SUMY < TOLNOR ) THEN
          IFNORY  = .TRUE. 
          IFALGN = .TRUE. 
      ENDIF
      IF ( SUMZ < TOLNOR ) THEN
          IFNORZ  = .TRUE. 
          IFALGN = .TRUE. 
      ENDIF
  
  ENDIF

  RETURN
END SUBROUTINE CHKNORD

!-----------------------------------------------------------------------
SUBROUTINE CHKAXCB()
  use size_m, only : ndim, nelv
  use input, only : cbc
  implicit none

  CHARACTER(3) CB
  integer :: ifld, nface, iel, ifc

  IFLD  = 1
  NFACE = 2*NDIM

  DO IEL=1,NELV
      DO IFC=1,NFACE
          CB = CBC(IFC,IEL,IFLD)
          IF  (CB == 'A  ' .AND. IFC /= 1)  GOTO 9000
      enddo
  END DO

  RETURN

  9000 WRITE (6,*) ' Element face on the axis of symmetry must be FACE 1'
  WRITE (6,*) ' Element',IEL,'   face',IFC,'  is on the axis.'
  call exitt

END SUBROUTINE CHKAXCB
!-----------------------------------------------------------------------
!> \brief Check for illegal boundary conditions
SUBROUTINE CHKCBC (CB,IEL,IFC,IFALGN)
  implicit none
  CHARACTER(3) CB
  LOGICAL :: IFALGN
  integer :: iel, ifc

!   Laplacian formulation only

  IF  (CB == 'SH ' .OR.  CB == 'sh ' .OR. &
  CB == 'SHL' .OR.  CB == 'shl' .OR. &
  CB == 'S  ' .OR.  CB == 's  ' .OR. &
  CB == 'SL ' .OR.  CB == 'sl ' .OR. &
  CB == 'MM ' .OR.  CB == 'mm ' .OR. &
  CB == 'MS ' .OR.  CB == 'ms ' .OR. &
  CB == 'MSI' .OR.  CB == 'msi'    )                GOTO 9001
  IF ( .NOT. IFALGN .AND. &
  (CB == 'ON ' .OR.  CB == 'on ' .OR. CB == 'SYM') ) GOTO 9010
  RETURN

  9001 WRITE (6,*) ' Illegal traction boundary conditions detected for'
  GOTO 9999
  9010 WRITE (6,*) ' Mixed B.C. on a side nonaligned with either the X,Y, &
  &or Z axis detected for'
  9999 WRITE (6,*) ' Element',IEL,'   side',IFC,'.'
  WRITE (6,*) ' Selected option only allowed for STRESS FORMULATION'
  WRITE (6,*) ' Execution terminates'
  call exitt
END SUBROUTINE CHKCBC
!-----------------------------------------------------------------------
!> \brief Zero out masks corresponding to Dirichlet boundary points.
SUBROUTINE BCMASK
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, ndim, nfield
  use input, only : ifmvbd, ifflow, cbc, ifstrs, ifheat, ipscal, ifmhd, ifaziv
  use input, only : ifaxis
  use soln, only : v1mask, v2mask, v3mask, pmask, omask, tmask
  use tstep, only : ifield, nelfld
  implicit none

  character(3) :: cb
  character(1) :: cb1(3)
  equivalence (cb1,cb)

  logical :: ifalgn,ifnorx,ifnory,ifnorz

  integer :: nfaces, nxyz, nel, ntot, iel, iface

  NFACES=2*NDIM
  NXYZ  =NX1*NY1*NZ1


!   Masks for moving mesh

  IF (IFMVBD) THEN
#if 0
      IFIELD = 0
      CALL STSMASK (W1MASK,W2MASK,W3MASK)
      do e=1,nelv
          do f=1,nfaces
              if (cbc(f,e,1) == 'msi' .OR. cbc(f,e,1) == 'msi') then
                  call facev(w1mask,e,f,0.0,nx1,ny1,nz1)
                  call facev(w2mask,e,f,0.0,nx1,ny1,nz1)
                  call facev(w3mask,e,f,0.0,nx1,ny1,nz1)
              endif
          enddo
      enddo
#endif
  ENDIF

!   Masks for flow variables

  IF (IFFLOW) THEN
      IFIELD = 1
      NEL    = NELFLD(IFIELD)
      NTOT   = NXYZ*NEL
  
  !        Pressure mask
  
      pmask = 1._dp
#if 0
      DO IEL=1,NELV
          DO IFACE=1,NFACES
              CB=CBC(IFACE,IEL,IFIELD)
              IF (CB == 'O  ' .OR. CB == 'ON ') &
              CALL FACEV(PMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
          enddo
      END DO
  
  !        Zero out mask at Neumann-Dirichlet interfaces
  
      CALL DSOP(PMASK,'MUL',NX1,NY1,NZ1)
#endif
  
  !        Velocity masks
  
  !        write(6,*) 'MASK ifstrs',ifstrs,ifield
  !        call exitt
      IF (IFSTRS) THEN
!max            CALL STSMASK (V1MASK,V2MASK,V3MASK)
      ELSE
      
          v1mask = 1._dp
          v2mask = 1._dp
          v3mask = 1._dp
!max          if (ifaxis) CALL RONE( OMASK,NTOT)
      
          DO IEL=1,NELV
              DO IFACE=1,NFACES
                  CB =CBC(IFACE,IEL,IFIELD)
                  CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFACE,IEL)
              
              !            All-Dirichlet boundary conditions
              
                  IF (CB == 'v  ' .OR. CB == 'V  ' .OR. CB == 'vl ' .OR. &
                  CB == 'VL ' .OR. CB == 'W  ') THEN
                      CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      cycle
                  ENDIF
              
              !        Mixed-Dirichlet-Neumann boundary conditions
              
                  IF (CB == 'SYM') THEN
                      IF ( .NOT. IFALGN .OR. IFNORX ) &
                      CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      IF ( IFNORY ) &
                      CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      IF ( IFNORZ ) &
                      CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      cycle
                  ENDIF
                  IF (CB == 'ON ') THEN
                      IF ( IFNORY .OR. IFNORZ ) &
                      CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      IF ( .NOT. IFALGN .OR. IFNORX .OR. IFNORZ ) &
                      CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      IF ( .NOT. IFALGN .OR. IFNORX .OR. IFNORY ) &
                      CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      cycle
                  ENDIF
                  IF (CB == 'A  ') THEN
                      CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                      CALL FACEV ( OMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
                  ENDIF
              enddo
          END DO

          if (ifaxis) CALL DSOP  ( OMASK,'MUL',NX1,NY1,NZ1)
          call opdsop(v1mask,v2mask,v3mask,'MUL') ! no rotation for mul


      ENDIF
  
  ENDIF

!   Masks for passive scalars +
!   k and e if k-e turbulence modem:
!   k = nfield-1
!   e = nfield

  IF (IFHEAT) THEN
  
      DO IFIELD=2,NFIELD
          IPSCAL = IFIELD-1
          NEL    = NELFLD(IFIELD)
          NTOT   = NXYZ*NEL
          tmask(:,:,:,:,ipscal) = 1._dp
          DO IEL=1,NEL
              DO IFACE=1,NFACES
                  CB =CBC(IFACE,IEL,IFIELD)
              
              !           Assign mask values.
              
                  IF  (CB == 'T  ' .OR. CB == 't  ' .OR. &
                  (CB == 'A  ' .AND. IFAZIV)    .OR. &
                  CB == 'MCI' .OR. CB == 'MLI' .OR. &
                  CB == 'KD ' .OR. CB == 'kd ' .OR. &
                  CB == 'ED ' .OR. CB == 'ed ' .OR. &
                  CB == 'KW ' .OR. CB == 'KWS' .OR. CB == 'EWS') &
                  CALL FACEV (TMASK(1,1,1,1,IPSCAL), &
                  IEL,IFACE,0.0,NX1,NY1,NZ1)
              enddo
          END DO
          CALL DSOP (TMASK(1,1,1,1,IPSCAL),'MUL',NX1,NY1,NZ1)
      END DO
  
  ENDIF

!   Masks for B-field

  if (ifmhd) then
    write(*,*) "Oops: ifmhd"
#if 0
      ifield = ifldmhd
      nel    = nelfld(ifield)
      ntot   = nxyz*nel
  
  !        B-field pressure mask
  
      call rone(bpmask,ntot)
      do iel=1,nelv
          do iface=1,nfaces
              cb=cbc(iface,iel,ifield)
              if (cb == 'O  ' .OR. cb == 'ON ') &
              call facev(bpmask,iel,iface,0.0,nx1,ny1,nz1)
          enddo
      enddo
  
  !        Zero out mask at Neumann-Dirichlet interfaces
  
      call dsop(bpmask,'MUL')
  
  !        B-field masks
  
      if (ifstrs) then
!max            call stsmask (b1mask,b2mask,b3mask)
      else
      
          call rone(b1mask,ntot)
          call rone(b2mask,ntot)
          call rone(b3mask,ntot)
      
          do iel=1,nelv
              do iface=1,nfaces
                  cb =cbc(iface,iel,ifield)
                  call chknord (ifalgn,ifnorx,ifnory,ifnorz,iface,iel)
              
                  if (cb == 'v  ' .OR. cb == 'V  ' .OR. cb == 'vl ' .OR. &
                  cb == 'VL ' .OR. cb == 'W  ') then
                  
                  !               All-Dirichlet boundary conditions
                  
                      call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                      call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                      call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                  
                  elseif (cb == 'SYM') then
                  
                  !               Mixed-Dirichlet-Neumann boundary conditions
                  
                      if ( .NOT. ifalgn .OR. ifnorx ) &
                      call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( ifnory ) &
                      call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( ifnorz ) &
                      call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                  
                  elseif (cb == 'ON ') then
                  
                  !               Mixed-Dirichlet-Neumann boundary conditions
                  
                      if ( ifnory .OR. ifnorz ) &
                      call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( .NOT. ifalgn .OR. ifnorx .OR. ifnorz ) &
                      call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( .NOT. ifalgn .OR. ifnorx .OR. ifnory ) &
                      call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                  
                  elseif (cb == 'A  ') then
                  
                  !               axisymmetric centerline
                  
                      call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                  
                  else
                  
                      if ( cb1(1) == 'd' ) &
                      call facev (b1mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( cb1(2) == 'd' ) &
                      call facev (b2mask,iel,iface,0.0,nx1,ny1,nz1)
                      if ( cb1(3) == 'd' .AND. if3d ) &
                      call facev (b3mask,iel,iface,0.0,nx1,ny1,nz1)
                  
                  endif
              enddo
          enddo
      
          call dsop(b1mask,'MUL')
          call dsop(b2mask,'MUL')
          if (ndim == 3) call dsop(b3mask,'MUL')
      endif
#endif
  endif

  RETURN
END SUBROUTINE BCMASK
!-----------------------------------------------------------------------
!> \brief Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
!! Use IFIELD as a guide to which boundary conditions are to be applied.
!!
!! \attention Most of the utility of this routine has been removed to avoid
!! temporary allocations.  This modification is not valid for the  boundaries
!!  - v
!!  - vl
!!  - ws
!!  - wsl
!!  - mv
!!  - mvn
!!  - d
!!  - on
SUBROUTINE BCDIRVC(V1,V2,V3,mask1,mask2,mask3)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, ndim, lelv
  use ctimer, only : icalld, tusbc, nusbc, etime1, dnekclock
  use input, only : if3d, cbc, bc, ifstrs
  use tstep, only : ifield, nelfld
  implicit none

  REAL(DP) :: V1(NX1,NY1,NZ1,LELV),V2(NX1,NY1,NZ1,LELV) &
  ,V3(NX1,NY1,NZ1,LELV)
  real(DP) :: mask1(nx1,ny1,nz1,lelv),mask2(nx1,ny1,nz1,lelv) &
  ,mask3(nx1,ny1,nz1,lelv)

  character(3) cb
  character(1) :: cb1(3)
  equivalence (cb1,cb)

  logical :: ifonbc

  integer :: nfaces, nxyz, nel, ntot, isweep, ie, iface
  real(DP) :: bc1, bc2, bc3

  ifonbc = .FALSE. 

#ifndef NOTIMER
  if (icalld == 0) then
      tusbc=0.0
      nusbc=0
      icalld=icalld+1
  endif
  nusbc=nusbc+1
  etime1=dnekclock()
#endif


  NFACES=2*NDIM
  NXYZ  =NX1*NY1*NZ1
  NEL   =NELFLD(IFIELD)
  NTOT  =NXYZ*NEL

#if 0
  allocate(TMP1(nx1,ny1,nz1,nelfld(ifield)) &
  , TMP2(nx1,ny1,nz1,nelfld(ifield)) &
  , TMP3(nx1,ny1,nz1,nelfld(ifield)) )

  tmp1 = 0._dp; tmp2 = 0._dp
  IF (IF3D) tmp3 = 0._dp
#endif

!   Velocity boundary conditions

!   write(6,*) 'BCDIRV: ifield',ifield
  DO ISWEEP=1,2
      DO IE=1,NEL
          DO IFACE=1,NFACES
              CB  = CBC(IFACE,IE,IFIELD)
              BC1 = BC(1,IFACE,IE,IFIELD)
              BC2 = BC(2,IFACE,IE,IFIELD)
              BC3 = BC(3,IFACE,IE,IFIELD)

              IF (CB == 'V  ' .OR. CB == 'VL '  .OR. &
              CB == 'WS ' .OR. CB == 'WSL') THEN
                  write(*,*) "Oops: CB = v, vl, ws, wsl"
#if 0
                  CALL FACEV (TMP1,IE,IFACE,BC1,NX1,NY1,NZ1)
                  CALL FACEV (TMP2,IE,IFACE,BC2,NX1,NY1,NZ1)
                  IF (IF3D) CALL FACEV (TMP3,IE,IFACE,BC3,NX1,NY1,NZ1)
                  IF ( IFQINP(IFACE,IE) ) &
                  CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                  TMP3(1,1,1,IE),IE,IFACE)
#endif
              ENDIF

              IF (CB == 'v  ' .OR. CB == 'vl ' .OR. &
              CB == 'ws ' .OR. CB == 'wsl' .OR. &
              CB == 'mv ' .OR. CB == 'mvn' .OR. &
              cb1(1) == 'd' .OR. cb1(2) == 'd' .OR. cb1(3) == 'd') then
                  write(*,*) "Oops: CB = something bad"
#if 0
                  call faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie), &
                  tmp3(1,1,1,ie),ie,iface,nx1,ny1,nz1)

                  IF ( IFQINP(IFACE,IE) ) &
                  CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                  TMP3(1,1,1,IE),IE,IFACE)
#endif
              ENDIF

              IF (CB == 'ON ' .OR. CB == 'on ') then   ! 5/21/01 pff
                  write(*,*) "Oops: CB = ON"
#if 0
                  ifonbc = .TRUE. 
                  CALL FACEIV ('v  ',TMP1(1,1,1,IE),TMP2(1,1,1,IE), &
                  TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
#endif
              ENDIF
          enddo
      END DO
#if 0
      DO IE=1,NEL
          DO IFACE=1,NFACES
              IF (CBC(IFACE,IE,IFIELD) == 'W  ') THEN
                  CALL FACEV (TMP1,IE,IFACE,0.0,NX1,NY1,NZ1)
                  CALL FACEV (TMP2,IE,IFACE,0.0,NX1,NY1,NZ1)
                  IF (IF3D) CALL FACEV (TMP3,IE,IFACE,0.0,NX1,NY1,NZ1)
              ENDIF
          enddo
      END DO
 
  !        Take care of Neumann-Dirichlet shared edges...
  
      if (isweep == 1) then
          call opdsop(tmp1,tmp2,tmp3,'MXA')
      else
          call opdsop(tmp1,tmp2,tmp3,'MNA')
      endif
#endif
  END DO
!   Copy temporary array to velocity arrays.

  IF ( .NOT. IFSTRS ) THEN
      v1 = v1 * mask1
      v2 = v2 * mask2
      IF (IF3D) v3 = v3 * mask3 
      if (ifonbc) then
        write(*,*) "Oops: ifonbc"
#if 0
          call antimsk1(tmp1,mask1,ntot)
          call antimsk1(tmp2,mask2,ntot)
          if (if3d) call antimsk1(tmp3,mask3,ntot)
#endif
      endif
  ELSE
    write(*,*) "Oops: ifstrs"
#if 0
      IF (IFMODEL) THEN
          CALL COPY (TMQ1,TMP1,NTOT)
          CALL COPY (TMQ2,TMP2,NTOT)
          IF (NDIM == 3) CALL COPY (TMQ3,TMP3,NTOT)
          CALL AMASK (TMP1,TMP2,TMP3,TMQ1,TMQ2,TMQ3,NELV)
      ENDIF
      CALL RMASK (V1,V2,V3,NELV)
#endif
  ENDIF

#if 0
  v1 = v1 + tmp1; v2 = v2 + tmp2 
  IF (IF3D) v3 = v3 + tmp3
#endif


#ifndef NOTIMER
  tusbc=tusbc+(dnekclock()-etime1)
#endif

  RETURN
END SUBROUTINE BCDIRVC
!-----------------------------------------------------------------------
!> \brief Apply Dirichlet boundary conditions to surface of scalar, S.
!! Use IFIELD as a guide to which boundary conditions are to be applied.
SUBROUTINE BCDIRSC(S)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, ndim, nfield
  use ctimer, only : icalld, tusbc, nusbc, etime1, dnekclock
  use input, only : cbc, bc
  use soln, only : tmask
  use tstep, only : ifield, nelfld
  implicit none

  real(DP) :: S(LX1,LY1,LZ1,LELT)
  real(DP), allocatable :: tmp(:,:,:,:) 

  CHARACTER(3) CB

  integer :: ifld, nfaces, nxyz, nel, ntot, isweep, ie, iface, nfldt
  real(DP) :: BC1, BC2, BC3, BC4, BCK, BCE

#ifndef NOTIMER
  if (icalld == 0) then
      tusbc=0.0
      nusbc=0
      icalld=icalld+1
  endif
  nusbc=nusbc+1
  etime1=dnekclock()
#endif

  IFLD   = 1
  NFACES = 2*NDIM
  NXYZ   = NX1*NY1*NZ1
  NEL    = NELFLD(IFIELD)
  NTOT   = NXYZ*NEL
  NFLDT  = NFIELD - 1

  allocate(tmp(nx1,ny1,nz1,nelfld(ifield)))
  tmp = 0._dp

!     Temperature boundary condition

  DO ISWEEP=1,2
    
#if 0
      IF (IFMODEL .AND. IFKEPS .AND. IFIELD >= NFLDT) &
      CALL TURBWBC (TMP,TMA,SMU)
#endif
    
      DO IE=1,NEL
          DO IFACE=1,NFACES
              CB=CBC(IFACE,IE,IFIELD)
              BC1=BC(1,IFACE,IE,IFIELD)
              BC2=BC(2,IFACE,IE,IFIELD)
              BC3=BC(3,IFACE,IE,IFIELD)
              BC4=BC(4,IFACE,IE,IFIELD)
              BCK=BC(4,IFACE,IE,IFLD)
              BCE=BC(5,IFACE,IE,IFLD)
              IF (CB == 'T  ') CALL FACEV (TMP,IE,IFACE,BC1,NX1,NY1,NZ1)
              IF (CB == 'MCI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
              IF (CB == 'MLI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
              IF (CB == 'KD ') CALL FACEV (TMP,IE,IFACE,BCK,NX1,NY1,NZ1)
              IF (CB == 'ED ') CALL FACEV (TMP,IE,IFACE,BCE,NX1,NY1,NZ1)
              IF (CB == 't  ' .OR. CB == 'kd ' .OR. CB == 'ed ') then
                write(*,*) "Oops: CB = t, kd, or ed"
!max              CALL FACEIS (CB,TMP(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
              ENDIF
          END DO
      enddo
  
  !        Take care of Neumann-Dirichlet shared edges...
  
      IF (ISWEEP == 1) CALL DSOP(TMP,'MXA',NX1,NY1,NZ1)
      IF (ISWEEP == 2) CALL DSOP(TMP,'MNA',NX1,NY1,NZ1)
  END DO

!     Copy temporary array to temperature array.

  s = s * tmask(:,:,:,:,ifield-1) + tmp

#ifndef NOTIMER
  tusbc=tusbc+(dnekclock()-etime1)
#endif

  RETURN
END SUBROUTINE BCDIRSC

!-----------------------------------------------------------------------
!> \brief Apply Neumann boundary conditions to surface of scalar, S.
!!  Use IFIELD as a guide to which boundary conditions are to be applied.
!!  If ITYPE = 1, then S is returned as the rhs contribution to the
!!                volumetric flux.
!!  If ITYPE =-1, then S is returned as the lhs contribution to the
!!                diagonal of A.
SUBROUTINE BCNEUSC(S,ITYPE)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt, nx1, ny1, nz1, ndim
  use ctimer, only : icalld, tusbc, nusbc, etime1, dnekclock
  use nekuse, only : hc, tinf, hrad, flux
  use geom, only : area
  use input, only : cbc, bc
  use geom, only : bm1
  use parallel, only : lglel
  use soln, only : t
  use tstep, only : ifield, nelfld
  implicit none

  real(DP), intent(out) :: S(LX1,LY1,LZ1,LELT)
  integer :: itype

  CHARACTER(3) CB

  real(DP) :: ts
  integer :: nfaces, nxyz, nel, ntot
  integer :: ie, iface, ieg, ia, ix, iy, iz
  integer :: kx1, kx2, ky1, ky2, kz1, kz2

#ifndef NOTIMER
  if (icalld == 0) then
      tusbc=0.0
      nusbc=0
      icalld=icalld+1
  endif
  nusbc=nusbc+1
  etime1=dnekclock()
#endif

  NFACES=2*NDIM
  NXYZ  =NX1*NY1*NZ1
  NEL   =NELFLD(IFIELD)
  NTOT  =NXYZ*NEL
  s = 0._dp

  IF (ITYPE == -1) THEN
  
  !        Compute diagonal contributions to accomodate Robin boundary conditions
  
      DO IE=1,NEL
          DO IFACE=1,NFACES
              ieg=lglel(ie)
              CB =CBC(IFACE,IE,IFIELD)
              IF (CB == 'C  ' .OR. CB == 'c  ' .OR. &
              CB == 'R  ' .OR. CB == 'r  ') THEN
              
                  IF (CB == 'C  ') HC   = BC(2,IFACE,IE,IFIELD)
                  IF (CB == 'R  ') THEN
                      TINF = BC(1,IFACE,IE,IFIELD)
                      HRAD = BC(2,IFACE,IE,IFIELD)
                  ENDIF
                  IA=0
              
              ! IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
              
                  CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
                  DO IZ=KZ1,KZ2
                      DO IY=KY1,KY2
                          DO IX=KX1,KX2
                              IA = IA + 1
                              TS = T(IX,IY,IZ,IE,IFIELD-1)
                              IF (CB == 'c  ' .OR. CB == 'r  ') THEN
                                  CALL NEKASGN (IX,IY,IZ,IE)
                                  CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                              ENDIF
                              IF (CB == 'r  ' .OR. CB == 'R  ') &
                              HC = HRAD * (TINF**2 + TS**2) * (TINF + TS)
                              S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) + &
                              HC*AREA(IA,1,IFACE,IE)/BM1(IX,IY,IZ,IE)
                          enddo
                      enddo
                  END DO
              ENDIF
          enddo
      END DO
  ENDIF
  IF (ITYPE == 1) THEN
  
  !        Add passive scalar fluxes to rhs
  
      DO IE=1,NEL
          DO IFACE=1,NFACES
              ieg=lglel(ie)
              CB =CBC(IFACE,IE,IFIELD)
              IF (CB == 'F  ' .OR. CB == 'f  ' .OR. &
              CB == 'C  ' .OR. CB == 'c  ' .OR. &
              CB == 'R  ' .OR. CB == 'r  ' ) THEN
              
                  IF (CB == 'F  ') FLUX=BC(1,IFACE,IE,IFIELD)
                  IF (CB == 'C  ') FLUX=BC(1,IFACE,IE,IFIELD) &
                  *BC(2,IFACE,IE,IFIELD)
                  IF (CB == 'R  ') THEN
                      TINF=BC(1,IFACE,IE,IFIELD)
                      HRAD=BC(2,IFACE,IE,IFIELD)
                  ENDIF
              
              !              Add local weighted flux values to rhs, S.
              
              ! IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
                  IA=0
                  CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
                  DO IZ=KZ1,KZ2
                      DO IY=KY1,KY2
                          DO IX=KX1,KX2
                              IA = IA + 1
                              TS = T(IX,IY,IZ,IE,IFIELD-1)
                              IF (CB == 'f  ') THEN
                                  CALL NEKASGN (IX,IY,IZ,IE)
                                  CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                              ENDIF
                              IF (CB == 'c  ') THEN
                                  CALL NEKASGN (IX,IY,IZ,IE)
                                  CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                                  FLUX = TINF*HC
                              ENDIF
                              IF (CB == 'r  ') THEN
                                  CALL NEKASGN (IX,IY,IZ,IE)
                                  CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                              ENDIF
                              IF (CB == 'R  ' .OR. CB == 'r  ') &
                              FLUX = HRAD*(TINF**2 + TS**2)*(TINF + TS) * TINF
                          
                          !                 Add computed fluxes to boundary surfaces:
                          
                              S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) &
                              + FLUX*AREA(IA,1,IFACE,IE)
                          enddo
                      enddo
                  enddo
              ENDIF
          enddo
      END DO
  ENDIF

#ifndef NOTIMER
    tusbc=tusbc+(dnekclock()-etime1)
#endif

  RETURN
END SUBROUTINE BCNEUSC
!-----------------------------------------------------------------------
!> \brief Assign NEKTON variables for definition (by user) of
!!   boundary conditions at collocation point (IX,IY,IZ)
!!   of element IEL.
!!     X             X-coordinate
!!     Y             Y-coordinate
!!     Z             Z-coordinate
!!     UX            X-velocity
!!     UY            Y-velocity
!!     UZ            Z-velocity
!!     TEMP          Temperature
!!     PS1           Passive scalar No. 1
!!     PS2           Passive scalar No. 2
!!      .             .
!!      .             .
!!     PS9           Passive scalar No. 9
!!     SI2           Strainrate invariant II
!!     SI3           Strainrate invariant III
!!   Variables to be defined by user for imposition of
!!   boundary conditions :
!!     SH1           Shear component No. 1
!!     SH2           Shear component No. 2
!!     TRX           X-traction
!!     TRY           Y-traction
!!     TRZ           Z-traction
!!     SIGMA         Surface-tension coefficient
!!     FLUX          Flux
!!     HC            Convection heat transfer coefficient
!!     HRAD          Radiation  heat transfer coefficient
!!     TINF          Temperature at infinity
SUBROUTINE NEKASGN (IX,IY,IZ,IEL)
  use size_m
  use geom
  use input
  use nekuse
  use soln
  use tstep
  implicit none

  integer, intent(in) :: ix, iy, iz, iel

  CHARACTER(3) CB

  integer :: ips

  X     = XM1(IX,IY,IZ,IEL)
  Y     = YM1(IX,IY,IZ,IEL)
  Z     = ZM1(IX,IY,IZ,IEL)
  R     = X**2+Y**2
  IF (R > 0.0) R=SQRT(R)
  IF (X /= 0.0 .OR. Y /= 0.0) THETA = ATAN2(Y,X)

  UX    = VX(IX,IY,IZ,IEL)
  UY    = VY(IX,IY,IZ,IEL)
  UZ    = VZ(IX,IY,IZ,IEL)
  TEMP  = T(IX,IY,IZ,IEL,1)
  DO 100 IPS=1,NPSCAL
      PS(IPS) = T(IX,IY,IZ,IEL,IPS+1)
  100 END DO
  UDIFF = VDIFF (IX,IY,IZ,IEL,IFIELD)
  UTRANS= VTRANS(IX,IY,IZ,IEL,IFIELD)

  cbu   = cb

  RETURN
END SUBROUTINE NEKASGN
!-----------------------------------------------------------------------
SUBROUTINE LFALSE (IFA,N)
  implicit none
  LOGICAL :: IFA(*)
  integer :: n, i
  DO 100 I=1,N
      IFA(I)= .FALSE. 
  100 END DO
  RETURN
END SUBROUTINE LFALSE
!-----------------------------------------------------------------------
SUBROUTINE UNITVEC (X,Y,Z,N)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  real(DP), intent(inout) :: X(n),Y(n),Z(n)
  real(DP) :: xlngth
  integer :: i
  DO I=1,N
      XLNGTH = SQRT( X(I)**2 + Y(I)**2 + Z(I)**2 )
      IF (XLNGTH /= 0.0) THEN
          X(I) = X(I)/XLNGTH
          Y(I) = Y(I)/XLNGTH
          Z(I) = Z(I)/XLNGTH
      ENDIF
  END DO
  RETURN
END SUBROUTINE UNITVEC
!-----------------------------------------------------------------------
