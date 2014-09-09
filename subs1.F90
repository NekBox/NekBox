!-----------------------------------------------------------------------
!> \brief Set the new time step. All cases covered.
subroutine setdt
  use kinds, only : DP
  use size_m, only : nid
  use input, only : param, ifflow, ifneknek, ifprint
  use soln, only : vx, vy, vz
  use tstep, only : dt, dtinit, time, fintim, lastep, courno, ctarg
  implicit none

  real(DP), save :: umax = 0._dp
  REAL(DP), save :: DTOLD = 0._dp
  REAL(DP), save :: DTOpf = 0._dp
  logical, save :: iffxdt = .FALSE.

  real(DP) :: dtcfl, dtfs, dtmin, dtmax
  real(DP), external :: uglmin

  if (param(12) < 0 .OR. iffxdt) then
      iffxdt    = .TRUE. 
      param(12) = abs(param(12))
      dt        = param(12)
      dtopf     = dt
      call compute_cfl(umax,vx,vy,vz,1.0)
      goto 200
  else IF (PARAM(84) /= 0.0) THEN
      if (dtold == 0.0) then
          dt   =param(84)
          dtold=param(84)
          dtopf=param(84)
          return
      else
          dtold=dt
          dtopf=dt
          dt=dtopf*param(85)
          dt=min(dt,param(12))
      endif
  endif

!   Find DT=DTCFL based on CFL-condition (if applicable)

  CALL SETDTC(umax)
  DTCFL = DT

!   Find DTFS based on surface tension (if applicable)

!  CALL SETDTFS (DTFS)

!   Select appropriate DT
  dtfs = 0.
  IF ((DT == 0.) .AND. (DTFS > 0.)) THEN
      DT = DTFS
  ELSEIF ((DT > 0.) .AND. (DTFS > 0.)) THEN
      DT = MIN(DT,DTFS)
  ELSEIF ((DT == 0.) .AND. (DTFS == 0.)) THEN
      DT = 0.
      IF (IFFLOW .AND. NID == 0 .AND. IFPRINT) THEN
          WRITE (6,*) 'WARNING: CFL-condition & surface tension'
          WRITE (6,*) '         are not applicable'
      endif
  ELSEIF ((DT > 0.) .AND. (DTFS == 0.)) THEN
      DT = DT
  ELSE
      DT = 0.
      IF (NID == 0) WRITE (6,*) 'WARNING: DT<0 or DTFS<0'
      IF (NID == 0) WRITE (6,*) '         Reset DT      '
  endif

!   Check DT against user-specified input, DTINIT=PARAM(12).

  IF ((DT > 0.) .AND. (DTINIT > 0.)) THEN
      DT = MIN(DT,DTINIT)
  ELSEIF ((DT == 0.) .AND. (DTINIT > 0.)) THEN
      DT = DTINIT
  ELSEIF ((DT > 0.) .AND. (DTINIT == 0.)) THEN
      DT = DT
  ELSEIF ( .NOT. iffxdt) THEN
      DT = 0.001
      IF(NID == 0)WRITE (6,*) 'WARNING: Set DT=0.001 (arbitrarily)'
  endif

!   Check if final time (user specified) has been reached.

  200 IF (TIME+DT >= FINTIM .AND. FINTIM /= 0.0) THEN
  !        Last step
      LASTEP = 1
      DT = FINTIM-TIME
      IF (NID == 0) WRITE (6,*) 'Final time step = ',DT
  endif

  COURNO = DT*UMAX
  IF (NID == 0 .AND. IFPRINT .AND. DT /= DTOLD) &
  WRITE (6,100) DT,DTCFL,DTFS,DTINIT
  100 FORMAT(5X,'DT/DTCFL/DTFS/DTINIT',4E12.3)

!   Put limits on how much DT can change.

  IF (DTOLD /= 0.0) THEN
      DTMIN=0.8*DTOLD
      DTMAX=1.2*DTOLD
      DT = MIN(DTMAX,DT)
      DT = MAX(DTMIN,DT)
  endif
  DTOLD=DT

!    IF (PARAM(84).NE.0.0) THEN
!          dt=dtopf*param(85)
!          dt=min(dt,param(12))
!    endif

  if (iffxdt) dt=dtopf
  COURNO = DT*UMAX

! synchronize time step for multiple sessions
  if (ifneknek) dt=uglmin(dt,1)

  if (iffxdt .AND. abs(courno) > 10.*abs(ctarg)) then
      if (nid == 0) write(6,*) 'CFL, Ctarg!',courno,ctarg
!max      call emerxit
  endif

  return
end subroutine setdt


!----------------------------------------------------------------------
!> \brief Check convergence for non-linear passisve scalar solver.
!!  Relevant for solving heat transport problems with radiation b.c.
!----------------------------------------------------------------------
subroutine cvgnlps (ifconv)
  use kinds, only : DP
  use input, only : ifnonl
  use tstep, only : ifield, tnrmh1, tolnl
  implicit none

  LOGICAL ::  IFCONV
  real(DP) :: tnorm1, tnorm2, eps

  IF (IFNONL(IFIELD)) THEN
      IFCONV = .FALSE. 
  ELSE
      IFCONV = .TRUE. 
      return
  endif

  TNORM1 = TNRMH1(IFIELD-1)
  CALL UNORM
  TNORM2 = TNRMH1(IFIELD-1)
  EPS = ABS((TNORM2-TNORM1)/TNORM2)
  IF (EPS < TOLNL) IFCONV = .TRUE. 

  return
end subroutine cvgnlps

!---------------------------------------------------------------------
!> \brief Norm calculation.
!---------------------------------------------------------------------
subroutine unorm
  use kinds, only : DP
  use soln, only : vx, vy, vz, t
  use tstep, only : ifield, imesh, tnrml8, tnrmsm, tnrmh1, time, dt, tnrml2
  use tstep, only : vnrmh1, vnrmsm, vnrml2, vnrml8, istep, dtinvm, vmean, tmean
  implicit none

  real(DP) :: tden, arg

  IF (IFIELD == 1) THEN
  
  !        Compute norms of the velocity.
  !        Compute time mean (L2) of the inverse of the time step.
  !        Compute L2 in time, H1 in space of the velocity.
  
      CALL NORMVC (VNRMH1,VNRMSM,VNRML2,VNRML8,VX,VY,VZ)
      IF (ISTEP == 0) return
      IF (ISTEP == 1) THEN
          DTINVM = 1./DT
          VMEAN  = VNRML8
      ELSE
          tden   = time
          if (time <= 0) tden = abs(time)+1.e-9
          arg    = ((TIME-DT)*DTINVM**2+1./DT)/tden
          if (arg > 0) DTINVM = SQRT(arg)
          arg    = ((TIME-DT)*VMEAN**2+DT*VNRMH1**2)/tden
          if (arg > 0) VMEAN  = SQRT(arg)
      endif
  ELSE
  
  !     Compute norms of a passive scalar
  
      CALL NORMSC (TNRMH1(IFIELD-1),TNRMSM(IFIELD-1), &
      TNRML2(IFIELD-1),TNRML8(IFIELD-1), &
      T(1,1,1,1,IFIELD-1),IMESH)
      TMEAN(IFIELD-1) = 0.
  endif

  return
end subroutine unorm

!--------------------------------------------------------------
!> \brief Compute new timestep based on CFL-condition
!--------------------------------------------------------------
subroutine setdtc(umax)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv
  use size_m, only : nx1, ny1, nz1, nelv, ndim, nid
  use geom, only : ifwcno, xm1, ym1, zm1
  use input, only : param, iftran, ifflow, ifnav, ifheat, ifadvc
  use input, only : ipscal, npscal
  use mass, only : binvm1
  use soln, only : vx, vy, vz, v1mask, v2mask, v3mask, bfx, bfy, bfz
  use tstep, only : lastep, dt, ifield, courno, ctarg, avtran, dtinit
  implicit none

  real(DP) :: umax
  real(DP) ::  u(lx1,ly1,lz1,lelv), v(lx1,ly1,lz1,lelv), w(lx1,ly1,lz1,lelv)

  REAL, save :: VCOUR

  INTEGER, save :: IFIRST = 0
  integer :: irst, iconv, ntot, ntotl, ntotd, i
  real(DP) :: dtold, cold, cmax, cmin, fmax, density, amax, dxchar, vold, cpred
  real(DP) :: a, b, c, discr, dtlow, dthi
  real(DP), external :: glmax, glmin

!   Steady state => all done
  IF ( .NOT. IFTRAN) THEN
      IFIRST=1
      LASTEP=1
      return
  endif

  irst = param(46)
  if (irst > 0) ifirst=1

!   First time around

  IF (IFIRST == 0) THEN
      DT     = DTINIT
      IF (IFFLOW) THEN
          IFIELD = 1
          CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      endif
  endif
  IFIRST=IFIRST+1

  DTOLD = DT

!   Convection ?

!   Don't enforce Courant condition if there is no convection.


  ICONV=0
  IF (IFFLOW .AND. IFNAV) ICONV=1
  IF (IFWCNO)             ICONV=1
  IF (IFHEAT) THEN
      DO 10 IPSCAL=0,NPSCAL
          IF (IFADVC(IPSCAL+2)) ICONV=1
      10 END DO
  endif
  IF (ICONV == 0) THEN
      DT=0.
      return
  endif


!   Find Courant and Umax


  NTOT   = NX1*NY1*NZ1*NELV
  NTOTL  = LX1*LY1*LZ1*LELV
  NTOTD  = NTOTL*NDIM
  COLD   = COURNO
  CMAX   = 1.2*CTARG
  CMIN   = 0.8*CTARG
  CALL CUMAX (VX,VY,VZ,u,v,w,UMAX)

!   Zero DT

  IF (DT == 0.0) THEN
  
      IF (UMAX /= 0.0) THEN
          DT = CTARG/UMAX
          VCOUR = UMAX
      ELSEIF (IFFLOW) THEN
      
      !           We'll use the body force to predict max velocity
      
          CALL SETPROP
          IFIELD = 1
      
          CALL MAKEUF
          CALL OPDSSUM (BFX,BFY,BFZ)
          CALL OPCOLV  (BFX,BFY,BFZ,BINVM1)
          FMAX=0.0
          CALL RZERO (U,NTOTD)
          DO 600 I=1,NTOT
              U(I,1,1,1) = ABS(BFX(I,1,1,1))
              V(I,1,1,1) = ABS(BFY(I,1,1,1))
              W(I,1,1,1) = ABS(BFZ(I,1,1,1))
          600 END DO
          FMAX    = GLMAX (U,NTOTD)
          DENSITY = AVTRAN(1)
          AMAX    = FMAX/DENSITY
          DXCHAR  = SQRT( (XM1(1,1,1,1)-XM1(2,1,1,1))**2 + &
          (YM1(1,1,1,1)-YM1(2,1,1,1))**2 + &
          (ZM1(1,1,1,1)-ZM1(2,1,1,1))**2 )
          DXCHAR  = GLMIN (dxchar,1)
          IF (AMAX /= 0.) THEN
              DT = SQRT(CTARG*DXCHAR/AMAX)
          ELSE
              IF (NID == 0) &
              WRITE (6,*) 'CFL: Zero velocity and body force'
              DT = 0.0
              return
          endif
      ELSEIF (IFWCNO) THEN
          IF (NID == 0) &
          WRITE (6,*) ' Stefan problem with no fluid flow'
          DT = 0.0
          return
      endif
  
  ELSEIF ((DT > 0.0) .AND. (UMAX /= 0.0)) THEN
  
  
  !     Nonzero DT & nonzero velocity
  
  
      COURNO = DT*UMAX
      VOLD   = VCOUR
      VCOUR  = UMAX
      IF (IFIRST == 1) THEN
          COLD = COURNO
          VOLD = VCOUR
      endif
      CPRED  = 2.*COURNO-COLD
  
  !     Change DT if it is too big or if it is too small
  
  !     if (nid.eq.0)
  !    $write(6,917) dt,umax,vold,vcour,cpred,cmax,courno,cmin
  ! 917 format(' dt',4f9.5,4f10.6)
      IF(COURNO > CMAX .OR. CPRED > CMAX .OR. COURNO < CMIN) THEN
      
          A=(VCOUR-VOLD)/DT
          B=VCOUR
      !           -C IS Target Courant number
          C=-CTARG
          DISCR=B**2-4*A*C
          DTOLD=DT
          IF(DISCR <= 0.0)THEN
              if (nid == 0) &
              PRINT*,'Problem calculating new DT Discriminant=',discr
              DT=DT*(CTARG/COURNO)
          !               IF(DT.GT.DTOLD) DT=DTOLD
          ELSE IF(ABS((VCOUR-VOLD)/VCOUR) < 0.001)THEN
          !              Easy: same v as before (LINEARIZED)
              DT=DT*(CTARG/COURNO)
          !     if (nid.eq.0)
          !    $write(6,918) dt,dthi,dtlow,discr,a,b,c
          ! 918 format(' d2',4f9.5,4f10.6)
          ELSE
              DTLOW=(-B+SQRT(DISCR) )/(2.0*A)
              DTHI =(-B-SQRT(DISCR) )/(2.0*A)
              IF(DTHI > 0.0 .AND. DTLOW > 0.0)THEN
                  DT = MIN (DTHI,DTLOW)
              !     if (nid.eq.0)
              !    $write(6,919) dt,dthi,dtlow,discr,a,b,c
              ! 919 format(' d3',4f9.5,4f10.6)
              ELSE IF(DTHI <= 0.0 .AND. DTLOW <= 0.0)THEN
              !                 PRINT*,'DTLOW,DTHI',DTLOW,DTHI
              !                 PRINT*,'WARNING: Abnormal DT from CFL-condition'
              !                 PRINT*,'         Keep going'
                  DT=DT*(CTARG/COURNO)
              ELSE
              !                 Normal case; 1 positive root, one negative root
                  DT = MAX (DTHI,DTLOW)
              !     if (nid.eq.0)
              !    $write(6,929) dt,dthi,dtlow,discr,a,b,c
              ! 929 format(' d4',4f9.5,4f10.6)
              endif
          endif
      !           We'll increase gradually-- make it the geometric mean between
      !     if (nid.eq.0)
      !    $write(6,939) dt,dtold
      ! 939 format(' d5',4f9.5,4f10.6)
          IF (DTOLD/DT < 0.2) DT = DTOLD*5
      endif
  
  endif

  return
end subroutine setdtc

subroutine cumax (v1,v2,v3,u,v,w,umax)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv, ndim
  use geom, only : rxm1, rym1, sxm1, sym1, rzm1, szm1, txm1, tym1, tzm1
  use input, only : ifaxis
  use wz_m, only : zgm1
  implicit none

  real(DP) :: V1(LX1,LY1,LZ1,1), V2(LX1,LY1,LZ1,1), V3(LX1,LY1,LZ1,1)
  real(DP) ::  u    (lx1,ly1,lz1,lelv) &
  ,             v    (lx1,ly1,lz1,lelv) &
  ,             w    (lx1,ly1,lz1,lelv)    
  real(DP) :: umax

  real(DP) ::  xrm1 (lx1,ly1,lz1,lelv) &
  ,            xsm1 (lx1,ly1,lz1,lelv) &
  ,            xtm1 (lx1,ly1,lz1,lelv) &
  ,            yrm1 (lx1,ly1,lz1,lelv) &
  ,            ysm1 (lx1,ly1,lz1,lelv) &
  ,            ytm1 (lx1,ly1,lz1,lelv)

  real(DP) :: zrm1 (lx1,ly1,lz1,lelv) &
  ,             zsm1 (lx1,ly1,lz1,lelv) &
  ,             ztm1 (lx1,ly1,lz1,lelv)


  real(DP) :: x(lx1,ly1,lz1,lelv), r(lx1,ly1,lz1,lelv)

  real(DP) :: drst(lx1), drsti(lx1)

  real(DP) :: U3(3)
  INTEGER, save :: ICALLD = 0

  integer :: ntot, ntotl, ntotd
  integer :: i, ie, ix, iy, iz
  real(DP), external :: vlmax, glmax


  NTOT  = NX1*NY1*NZ1*NELV
  NTOTL = LX1*LY1*LZ1*LELV
  NTOTD = NTOTL*NDIM

!   Compute isoparametric partials.

  CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1, &
  IFAXIS)

!   Compute maximum U/DX

  IF (ICALLD == 0) THEN
      ICALLD=1
      DRST (1)=ABS(ZGM1(2,1)-ZGM1(1,1))
      DRSTI(1)=1.0/DRST(1)
      DO 400 I=2,NX1-1
          DRST (I)=ABS(ZGM1(I+1,1)-ZGM1(I-1,1))/2.0
          DRSTI(I)=1.0/DRST(I)
      400 END DO
      DRST (NX1)=DRST(1)
      DRSTI(NX1)=1.0/DRST(NX1)
  endif

!   Zero out scratch arrays U,V,W for ALL declared elements...

  CALL RZERO3 (U,V,W,NTOTL)

  IF (NDIM == 2) THEN

      CALL VDOT2  (U,V1  ,V2  ,RXM1,RYM1,NTOT)
      CALL VDOT2  (R,RXM1,RYM1,RXM1,RYM1,NTOT)
      CALL VDOT2  (X,XRM1,YRM1,XRM1,YRM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(U,R,NTOT)
  
      CALL VDOT2  (V,V1  ,V2  ,SXM1,SYM1,NTOT)
      CALL VDOT2  (R,SXM1,SYM1,SXM1,SYM1,NTOT)
      CALL VDOT2  (X,XSM1,YSM1,XSM1,YSM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(V,R,NTOT)
  
  ELSE
  
      CALL VDOT3  (U,V1  ,V2  ,V3  ,RXM1,RYM1,RZM1,NTOT)
      CALL VDOT3  (R,RXM1,RYM1,RZM1,RXM1,RYM1,RZM1,NTOT)
      CALL VDOT3  (X,XRM1,YRM1,ZRM1,XRM1,YRM1,ZRM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(U,R,NTOT)
  
      CALL VDOT3  (V,V1  ,V2  ,V3  ,SXM1,SYM1,SZM1,NTOT)
      CALL VDOT3  (R,SXM1,SYM1,SZM1,SXM1,SYM1,SZM1,NTOT)
      CALL VDOT3  (X,XSM1,YSM1,ZSM1,XSM1,YSM1,ZSM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(V,R,NTOT)
  
      CALL VDOT3  (W,V1  ,V2  ,V3  ,TXM1,TYM1,TZM1,NTOT)
      CALL VDOT3  (R,TXM1,TYM1,TZM1,TXM1,TYM1,TZM1,NTOT)
      CALL VDOT3  (X,XTM1,YTM1,ZTM1,XTM1,YTM1,ZTM1,NTOT)
      CALL COL2   (R,X,NTOT)
      CALL VSQRT  (R,NTOT)
      CALL INVCOL2(W,R,NTOT)
  
  endif

  DO IE=1,NELV
      DO IX=1,NX1
          DO IY=1,NY1
              DO IZ=1,NZ1
                  U(IX,IY,IZ,IE)=ABS( U(IX,IY,IZ,IE)*DRSTI(IX) )
                  V(IX,IY,IZ,IE)=ABS( V(IX,IY,IZ,IE)*DRSTI(IY) )
                  W(IX,IY,IZ,IE)=ABS( W(IX,IY,IZ,IE)*DRSTI(IZ) )
              enddo
          enddo
      enddo
  END DO

  U3(1)   = VLMAX(U,NTOT)
  U3(2)   = VLMAX(V,NTOT)
  U3(3)   = VLMAX(W,NTOT)
  UMAX    = GLMAX(U3,3)

  return
end subroutine cumax

!> \brief Collocate B with A on the surface IFACE1 of element IE.
!!  A is a (NX,NY,NZ) data structure
!!  B is a (NX,NY,IFACE) data structure
!!  IFACE1 is in the preprocessor notation
!!  IFACE  is the dssum notation.
!!  5 Jan 1989 15:12:22      PFF
subroutine faccl2(a,b,iface1)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nx1, ny1, nz1
  use topol, only : eface1, skpdat
  implicit none

  real(DP) :: A(LX1,LY1,LZ1),B(LX1,LY1)
  real(DP) :: af(lx1*ly1*lz1), bf(lx1*ly1)
  integer :: iface1

  integer :: j1, j2, i, jskip2, jf2, js2, jskip1, jf1, js1, iface

!   Set up counters

  CALL DSSET(NX1,NY1,NZ1)
  IFACE  = EFACE1(IFACE1)
  JS1    = SKPDAT(1,IFACE)
  JF1    = SKPDAT(2,IFACE)
  JSKIP1 = SKPDAT(3,IFACE)
  JS2    = SKPDAT(4,IFACE)
  JF2    = SKPDAT(5,IFACE)
  JSKIP2 = SKPDAT(6,IFACE)

  af = reshape(a, (/lx1*ly1*lz1/))
  bf = reshape(b, (/lx1*ly1/))
  I = 0
  DO J2=JS2,JF2,JSKIP2
      DO J1=JS1,JF1,JSKIP1
          I = I+1
          af(J1+lx1*(J2-1)) = Af(J1+lx1*(J2-1))*Bf(I)
      enddo
  END DO
  a = reshape(af, (/lx1, ly1, lz1/))

  return
end subroutine faccl2

!-----------------------------------------------------------------------
!> \brief Set the variable property arrays H1 and H2
!! in the Helmholtz equation.
!! (associated with variable IFIELD)
!! INTLOC =      integration type
subroutine sethlm (h1,h2,intloc)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use input, only : iftran, param
  use soln, only : vdiff, vtrans
  use tstep, only : ifield, nelfld, bd, dt
  implicit none

  real(DP) :: h1(1),h2(1)
  integer :: intloc

  integer :: nel, ntot1, i
  real(DP) :: dtbd

  nel   = nelfld(ifield)
  ntot1 = nx1*ny1*nz1*nel

  if (iftran) then
      dtbd = bd(1)/dt
      call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
      if (intloc == 0) then
          call rzero (h2,ntot1)
      else
          if (ifield == 1 .OR. param(107) == 0) then

              call cmult2 (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

          else   ! unsteady reaction-diffusion type equation

              do i=1,ntot1
                  h2(i) = dtbd*vtrans(i,1,1,1,ifield) + param(107)
              enddo

          endif

      endif

  !        if (ifield.eq.1 .and. ifanls) then   ! this should be replaced
  !           const = 2.                        ! with a correct stress
  !           call cmult (h1,const,ntot1)       ! formulation
  !        endif

  ELSE
      CALL COPY  (H1,VDIFF (1,1,1,1,IFIELD),NTOT1)
      CALL RZERO (H2,NTOT1)
      if (param(107) /= 0) then
          write(6,*) 'SPECIAL SETHLM!!',param(107)
      !           call cfill (h2,param(107),ntot1)
          call copy  (h2,vtrans(1,1,1,1,ifield),ntot1)
      endif
  endif

  return
end subroutine sethlm

!-----------------------------------------------------------------------
!> \brief Set material properties
!!  Material type: 0 for default  (PARAM and PCOND/PRHOCP)
!!                 1 for constant props;
!!                 2 for fortran function;
!-----------------------------------------------------------------------
subroutine vprops
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use input, only : ifuservp, ifvarp, matype, igroup, iflomach
  use input, only : cpfld, cpgrp
  use soln, only : vdiff, vtrans
  use tstep, only : ifield, nelfld, istep
  implicit none

  integer :: nxyz1, nel, ntot1, iel, igrp, itype, itest
  real(DP) :: cdiff, ctrans
  integer, external :: iglmax

  NXYZ1 = NX1*NY1*NZ1
  NEL   = NELFLD(IFIELD)
  NTOT1 = NXYZ1*NEL

  IF (ISTEP == 0) THEN
  
  !        First time around, set defaults
  
      ifvarp(ifield) = .FALSE. 
      if (iflomach) ifvarp(ifield) = .TRUE. 

      if ( .NOT. ifvarp(ifield)) then ! check all groups
          do iel=1,nel
              igrp  = igroup(iel)
              itype = matype(igrp,ifield)
              if(itype /= 0) ifvarp(ifield) = .TRUE. 
          enddo
      endif

      itest = 0                        ! test against all processors
      if (ifvarp(ifield)) itest = 1
      itest = iglmax(itest,1)
      if (itest > 0) ifvarp(ifield) = .TRUE. 

  endif

!   Fill up property arrays every time step

!   First, check for turbulence models

#if 0
  IF (IFMODEL .AND. IFKEPS) THEN
      CALL TURBFLD (IFKFLD,IFEFLD)
      IF (IFKFLD)           CALL TPROPK
      IF (IFEFLD)           CALL TPROPE
      IF (IFKFLD .OR. IFEFLD) return
  endif
#endif

!...  No turbulence models, OR current field is not k or e.

  DO IEL=1,NEL
  
      IGRP=IGROUP(IEL)

      if (ifuservp) then
#if 0
      
      !           User specified fortran function   (pff 2/13/01)
          CALL NEKUVP (IEL)
          DIFMIN = VLMIN(VDIFF(1,1,1,IEL,IFIELD),NXYZ1)
          IF (DIFMIN <= 0.0) THEN
              WRITE (6,100) DIFMIN,IFIELD,IGRP
              CALL EXITT
          endif
#endif        
      ELSE IF(MATYPE(IGRP,IFIELD) == 1)THEN
      
      !           Constant property within groups of elements
      
          CDIFF  = CPGRP(IGRP,IFIELD,1)
          CTRANS = CPGRP(IGRP,IFIELD,2)
          CALL CFILL(VDIFF (1,1,1,IEL,IFIELD),CDIFF,NXYZ1)
          CALL CFILL(VTRANS(1,1,1,IEL,IFIELD),CTRANS,NXYZ1)
          IF (CDIFF <= 0.0) THEN
              WRITE(6,100) CDIFF,IFIELD,IGRP
              100 FORMAT(2X,'ERROR:  Non-positive diffusivity (' &
              ,G12.3,') specified for field',I2,', group',I2 &
              ,' element',I4,'.' &
              ,/,'ABORTING in VPROPS',//)
              CALL EXITT
          endif
      
      ELSE IF(MATYPE(IGRP,IFIELD) == 2)THEN
        write(*,*) "Oops: matype" 
#if 0
      !           User specified fortran function
      
          CALL NEKUVP (IEL)
      
          DIFMIN = VLMIN(VDIFF(1,1,1,IEL,IFIELD),NXYZ1)
          IF (DIFMIN <= 0.0) THEN
              WRITE (6,100) DIFMIN,IFIELD,IGRP
              CALL EXITT
          endif
#endif        
      ELSE IF(MATYPE(IGRP,IFIELD) == 0)THEN
      
      !           Default constant property
      
          CDIFF  = CPFLD(IFIELD,1)
          CTRANS = CPFLD(IFIELD,2)
      !           write(6,*) 'vdiff:',ifield,cdiff,ctrans
          CALL CFILL(VDIFF (1,1,1,IEL,IFIELD),CDIFF,NXYZ1)
          CALL CFILL(VTRANS(1,1,1,IEL,IFIELD),CTRANS,NXYZ1)
          IF (CDIFF <= 0.0) THEN
              WRITE(6,200) CDIFF,IFIELD
              200 FORMAT(2X,'ERROR:  Non-positive diffusivity (' &
              ,G12.3,') specified for field',I2,'.',/ &
              ,'ABORTING in VPROPS',//)
              CALL EXITT
          endif
      endif
  
  END DO

!     Turbulence models --- sum eddy viscosity/diffusivity
#if 0
  IF (IFMODEL .AND. (IFIELD == 1 .OR. IFIELD == 2)) &
  CALL TVISCOS
#endif

  return
end subroutine vprops

!-----------------------------------------------------------------------
!> Set ifsolv = .FALSE.
subroutine setsolv
  use mesh, only : ifsolv
  implicit none
    
  IFSOLV = .FALSE. 
  return
end subroutine setsolv

!-----------------------------------------------------------------------
subroutine flush_io
  return
end subroutine flush_io

!-----------------------------------------------------------------------
!!     IOP = 0
!!     Extract vector (A1,A2,A3) from (B1,B2,B3) on face IFACE1.
!!     IOP = 1
!!     Extract vector (B1,B2,B3) from (A1,A2,A3) on face IFACE1.
!!     A1, A2, A3 have the (NX,NY,NFACE) data structure
!!     B1, B2, B3 have the (NX,NY,NZ)    data structure
!!     IFACE1 is in the preprocessor notation
!!     IFACE  is the dssum notation.
SUBROUTINE FACEXV (A1,A2,A3,B1,B2,B3,IFACE1,IOP)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx1, ly1, lz1
  use topol, only : eface1, skpdat
  implicit none

  real(DP) :: A1(LX1,LY1),    A2(LX1,LY1),    A3(LX1,LY1)
  real(DP) :: B1(LX1,LY1,LZ1),B2(LX1,LY1,LZ1),B3(LX1,LY1,LZ1)
  real(DP) :: A1f(lx1*ly1), A2f(lx1*ly1), A3f(lx1*ly1)
  real(DP) :: B1f(lx1*ly1*lz1), b2f(lx1*ly1*lz1), b3f(lx1*ly1*lz1)
  integer :: iface1, iop

  integer :: j1, j2, i, jskip2, jf2, js2, jskip1, jf1, js1, iface

  CALL DSSET(NX1,NY1,NZ1)
  IFACE  = EFACE1(IFACE1)
  JS1    = SKPDAT(1,IFACE)
  JF1    = SKPDAT(2,IFACE)
  JSKIP1 = SKPDAT(3,IFACE)
  JS2    = SKPDAT(4,IFACE)
  JF2    = SKPDAT(5,IFACE)
  JSKIP2 = SKPDAT(6,IFACE)
  I = 0

  IF (IOP == 0) THEN
      b1f = reshape(b1, (/lx1*ly1*lz1/))
      b2f = reshape(b2, (/lx1*ly1*lz1/))
      b3f = reshape(b3, (/lx1*ly1*lz1/))
      DO J2=JS2,JF2,JSKIP2
          DO J1=JS1,JF1,JSKIP1
              I = I+1
              A1f(I) = B1f(J1+lx1*(J2-1))
              A2f(I) = B2f(J1+lx1*(J2-1))
              A3f(I) = B3f(J1+lx1*(J2-1))
          enddo
      END DO
      A1 = reshape(A1f, (/lx1, ly1/))
      A2 = reshape(A2f, (/lx1, ly1/))
      A3 = reshape(A3f, (/lx1, ly1/))
  ELSE
      a1f = reshape(a1, (/lx1*ly1/))
      a2f = reshape(a2, (/lx1*ly1/))
      a3f = reshape(a3, (/lx1*ly1/))
      DO J2=JS2,JF2,JSKIP2
          DO J1=JS1,JF1,JSKIP1
              I = I+1
              B1f(J1+lx1*(J2-1)) = A1f(I)
              B2f(J1+lx1*(J2-1)) = A2f(I)
              B3f(J1+lx1*(J2-1)) = A3f(I)
          enddo
      END DO
      B1 = reshape(B1f, (/lx1, ly1, lz1/))
      B2 = reshape(B2f, (/lx1, ly1, lz1/))
      B3 = reshape(B3f, (/lx1, ly1, lz1/))
  ENDIF

  RETURN
END SUBROUTINE FACEXV

SUBROUTINE CMULT2 (A,B,CONST,N)
  use kinds, only : DP
  implicit none
  real(DP) :: A(1),B(1), const
  integer :: n, i
  DO I=1,N
      A(I)=B(I)*CONST
  END DO
  RETURN
END SUBROUTINE CMULT2

SUBROUTINE UPDMSYS (IFLD)
  use size_m
  use geom, only : iflmsf
  implicit none

  integer, intent(in) :: ifld

  IF ( .NOT. IFLMSF(IFLD)) RETURN
#if 0
  NEL  = NELFLD(IFLD)
  CALL SETHMSK (HVMASK,HFMASK,IFLD,NEL)
  CALL SETCSYS (HVMASK,HFMASK,NEL)
#endif

  RETURN
END SUBROUTINE UPDMSYS

SUBROUTINE SETCDOF
  use size_m, only : ndim, nelt
  use input, only : cbc, cdof
  implicit none 

  integer :: nface, ifc, iel
  NFACE = 2*NDIM

  DO IEL=1,NELT
      DO IFC=1,NFACE
          CDOF(IFC,IEL)=CBC(IFC,IEL,0)(1:1)
      enddo
  END DO

  RETURN
END SUBROUTINE SETCDOF
