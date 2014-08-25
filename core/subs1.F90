!-----------------------------------------------------------------------
    subroutine setdt

!     Set the new time step. All cases covered.

    use size_m
    use input
    use soln
    use tstep
    common /cprint/ ifprint
    logical ::         ifprint
    common /udxmax/ umax
    REAL ::     DTOLD
    SAVE     DTOLD
    DATA     DTOLD /0.0/
    REAL ::     DTOpf
    SAVE     DTOpf
    DATA     DTOpf /0.0/
    logical :: iffxdt
    save    iffxdt
    data    iffxdt / .FALSE. /

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

!     Find DT=DTCFL based on CFL-condition (if applicable)

    CALL SETDTC
    DTCFL = DT

!     Find DTFS based on surface tension (if applicable)

!    CALL SETDTFS (DTFS)

!     Select appropriate DT

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

!     Check DT against user-specified input, DTINIT=PARAM(12).

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

!     Check if final time (user specified) has been reached.

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

!     Put limits on how much DT can change.

    IF (DTOLD /= 0.0) THEN
        DTMIN=0.8*DTOLD
        DTMAX=1.2*DTOLD
        DT = MIN(DTMAX,DT)
        DT = MAX(DTMIN,DT)
    endif
    DTOLD=DT

!      IF (PARAM(84).NE.0.0) THEN
!            dt=dtopf*param(85)
!            dt=min(dt,param(12))
!      endif

    if (iffxdt) dt=dtopf
    COURNO = DT*UMAX

! synchronize time step for multiple sessions
    if (ifneknek) dt=uglmin(dt,1)

    if (iffxdt .AND. abs(courno) > 10.*abs(ctarg)) then
        if (nid == 0) write(6,*) 'CFL, Ctarg!',courno,ctarg
        call emerxit
    endif


    return
    end subroutine setdt

!--------------------------------------------------------

    subroutine cvgnlps (ifconv)
!----------------------------------------------------------------------

!     Check convergence for non-linear passisve scalar solver.
!     Relevant for solving heat transport problems with radiation b.c.

!----------------------------------------------------------------------
    use size_m
    use input
    use tstep
    LOGICAL ::  IFCONV

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

    subroutine unorm
!---------------------------------------------------------------------

!     Norm calculation.

!---------------------------------------------------------------------
    use size_m
    use soln
    use tstep

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

    subroutine setdtc
!--------------------------------------------------------------

!     Compute new timestep based on CFL-condition

!--------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use mvgeom
    use soln
    use tstep

    common /ctmp1/ u(lx1,ly1,lz1,lelv) &
    ,             v(lx1,ly1,lz1,lelv) &
    ,             w(lx1,ly1,lz1,lelv)
    common /ctmp0/ x(lx1,ly1,lz1,lelv) &
    ,             r(lx1,ly1,lz1,lelv)
    common /udxmax/ umax


    REAL :: VCOUR
    SAVE VCOUR

    INTEGER :: IFIRST
    SAVE    IFIRST
    DATA    IFIRST/0/


!     Steady state => all done


    IF ( .NOT. IFTRAN) THEN
        IFIRST=1
        LASTEP=1
        return
    endif

    irst = param(46)
    if (irst > 0) ifirst=1

!     First time around

    IF (IFIRST == 0) THEN
        DT     = DTINIT
        IF (IFFLOW) THEN
            IFIELD = 1
            CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
        endif
    endif
    IFIRST=IFIRST+1

    DTOLD = DT

!     Convection ?

!     Don't enforce Courant condition if there is no convection.


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


!     Find Courant and Umax


    NTOT   = NX1*NY1*NZ1*NELV
    NTOTL  = LX1*LY1*LZ1*LELV
    NTOTD  = NTOTL*NDIM
    COLD   = COURNO
    CMAX   = 1.2*CTARG
    CMIN   = 0.8*CTARG
    CALL CUMAX (VX,VY,VZ,UMAX)

!     Zero DT

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

    subroutine cumax (v1,v2,v3,umax)

    use size_m
    use geom
    use input
    use wz_m

    common /scrns/ xrm1 (lx1,ly1,lz1,lelv) &
    ,             xsm1 (lx1,ly1,lz1,lelv) &
    ,             xtm1 (lx1,ly1,lz1,lelv) &
    ,             yrm1 (lx1,ly1,lz1,lelv) &
    ,             ysm1 (lx1,ly1,lz1,lelv) &
    ,             ytm1 (lx1,ly1,lz1,lelv)
    common /scrmg/ zrm1 (lx1,ly1,lz1,lelv) &
    ,             zsm1 (lx1,ly1,lz1,lelv) &
    ,             ztm1 (lx1,ly1,lz1,lelv)
    common /ctmp1/ u    (lx1,ly1,lz1,lelv) &
    ,             v    (lx1,ly1,lz1,lelv) &
    ,             w    (lx1,ly1,lz1,lelv)
    common /ctmp0/ x    (lx1,ly1,lz1,lelv) &
    ,             r    (lx1,ly1,lz1,lelv)
    common /delrst/ drst(lx1),drsti(lx1)

    DIMENSION V1(LX1,LY1,LZ1,1) &
    , V2(LX1,LY1,LZ1,1) &
    , V3(LX1,LY1,LZ1,1)
    DIMENSION U3(3)
    INTEGER :: ICALLD
    SAVE    ICALLD
    DATA    ICALLD /0/

    NTOT  = NX1*NY1*NZ1*NELV
    NTOTL = LX1*LY1*LZ1*LELV
    NTOTD = NTOTL*NDIM

!     Compute isoparametric partials.

    CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1, &
    IFAXIS)

!     Compute maximum U/DX

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

!     Zero out scratch arrays U,V,W for ALL declared elements...

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

    DO 500 IE=1,NELV
        DO 500 IX=1,NX1
            DO 500 IY=1,NY1
                DO 500 IZ=1,NZ1
                    U(IX,IY,IZ,IE)=ABS( U(IX,IY,IZ,IE)*DRSTI(IX) )
                    V(IX,IY,IZ,IE)=ABS( V(IX,IY,IZ,IE)*DRSTI(IY) )
                    W(IX,IY,IZ,IE)=ABS( W(IX,IY,IZ,IE)*DRSTI(IZ) )
    500 END DO

    U3(1)   = VLMAX(U,NTOT)
    U3(2)   = VLMAX(V,NTOT)
    U3(3)   = VLMAX(W,NTOT)
    UMAX    = GLMAX(U3,3)

    return
    end subroutine cumax

    subroutine cdxmin2 (dtst,rhosig,iel,ifc,ifaxis)

    use size_m
    use dxyz
    use geom
    common /delrst/ drst(lx1),drsti(lx1)
    common /ctmp0/  xfm1(lx1),yfm1(lx1),t1xf(lx1),t1yf(lx1)
    DIMENSION DTST(LX1,1)
    LOGICAL :: IFAXIS

    DELTA = 1.E-9
    X     = 1.+DELTA
    Y     = 1.
    DIFF  = ABS(X-Y)
    IF (DIFF == 0.) EPS = 1.E-6
    IF (DIFF > 0.) EPS = 1.E-13

    CALL FACEC2 (XFM1,YFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),IFC)

    IF (IFC == 1 .OR. IFC == 3) THEN
        CALL MXM (DXM1,NX1,XFM1,NX1,T1XF,1)
        CALL MXM (DXM1,NX1,YFM1,NX1,T1YF,1)
    ELSE
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
        CALL MXM (DYM1,NY1,XFM1,NY1,T1XF,1)
        CALL MXM (DYM1,NY1,YFM1,NY1,T1YF,1)
    endif

    IF (IFAXIS) THEN
        DO 100 IX=1,NX1
            IF (YFM1(IX) < EPS) THEN
                DTST(IX,1) = 1.e+10
            ELSE
                XJ = SQRT( T1XF(IX)**2 + T1YF(IX)**2 )*DRST(IX)
                DTST(IX,1) = RHOSIG * SQRT( XJ**3 ) * YFM1(IX)
            endif
        100 END DO
    ELSE
        DO 200 IX=1,NX1
            XJ = SQRT( T1XF(IX)**2 + T1YF(IX)**2 )*DRST(IX)
            DTST(IX,1) = RHOSIG * SQRT( XJ**3 )
        200 END DO
    endif

    return
    end subroutine cdxmin2
    subroutine cdxmin3 (dtst,rhosig,iel,ifc)

    use size_m
    use dxyz
    use geom
    common /delrst/ drst(lx1),drsti(lx1)
    common /ctmp0/  xfm1(lx1,ly1),yfm1(lx1,ly1),zfm1(lx1,ly1)
    common /ctmp1/  drm1(lx1,lx1),drtm1(lx1,ly1) &
    ,  dsm1(lx1,lx1),dstm1(lx1,ly1)
    common /scrmg/  xrm1(lx1,ly1),yrm1(lx1,ly1),zrm1(lx1,ly1) &
    ,  xsm1(lx1,ly1),ysm1(lx1,ly1),zsm1(lx1,ly1)
    dimension dtst(lx1,ly1)

    call facexv (xfm1,yfm1,zfm1,xm1(1,1,1,iel),ym1(1,1,1,iel), &
    zm1(1,1,1,iel),ifc,0)
    call setdrs (drm1,drtm1,dsm1,dstm1,ifc)

    CALL MXM (DRM1,NX1, XFM1,NX1,XRM1,NY1)
    CALL MXM (DRM1,NX1, YFM1,NX1,YRM1,NY1)
    CALL MXM (DRM1,NX1, ZFM1,NX1,ZRM1,NY1)
    CALL MXM (XFM1,NX1,DSTM1,NY1,XSM1,NY1)
    CALL MXM (YFM1,NX1,DSTM1,NY1,YSM1,NY1)
    CALL MXM (ZFM1,NX1,DSTM1,NY1,ZSM1,NY1)

    DO 100 IX=1,NX1
        DO 100 IY=1,NY1
            DELR = XRM1(IX,IY)**2 + YRM1(IX,IY)**2 + ZRM1(IX,IY)**2
            DELS = XSM1(IX,IY)**2 + YSM1(IX,IY)**2 + ZSM1(IX,IY)**2
            DELR = SQRT( DELR )*DRST(IX)
            DELS = SQRT( DELS )*DRST(IY)
            XJ   = MIN( DELR,DELS )
            DTST(IX,IY) = RHOSIG * SQRT( XJ**3 )
    100 END DO

    return
    end subroutine cdxmin3

    FUNCTION FACDOT(A,B,IFACE1)

!     Take the dot product of A and B on the surface IFACE1 of element IE.

!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.
!         5 Jan 1989 15:12:22      PFF

    use size_m
    use topol
    DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    SUM=0.0
    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            SUM = SUM + A(J1,J2,1)*B(I,1)
    100 END DO

    FACDOT = SUM

    return
    END FUNCTION FACDOT

    subroutine fcaver(xaver,a,iel,iface1)
!------------------------------------------------------------------------

!     Compute the average of A over the face IFACE1 in element IEL.

!         A is a (NX,NY,NZ) data structure
!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.
!------------------------------------------------------------------------
    use size_m
    use geom
    use topol
    REAL :: A(LX1,LY1,LZ1,1)

    FCAREA = 0.
    XAVER  = 0.

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            FCAREA = FCAREA+AREA(I,1,IFACE1,IEL)
            XAVER  = XAVER +AREA(I,1,IFACE1,IEL)*A(J1,J2,1,IEL)
    100 END DO

    XAVER = XAVER/FCAREA
    return
    end subroutine fcaver
    subroutine faccl2(a,b,iface1)

!     Collocate B with A on the surface IFACE1 of element IE.

!         A is a (NX,NY,NZ) data structure
!         B is a (NX,NY,IFACE) data structure
!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.
!         5 Jan 1989 15:12:22      PFF

    use size_m
    use topol
    DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A(J1,J2,1) = A(J1,J2,1)*B(I,1)
    100 END DO

    return
    end subroutine faccl2

    subroutine faccl3(a,b,c,iface1)

!     Collocate B with A on the surface IFACE1 of element IE.

!         A is a (NX,NY,NZ) data structure
!         B is a (NX,NY,IFACE) data structure
!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.
!         5 Jan 1989 15:12:22      PFF

    use size_m
    use topol
    DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A(J1,J2,1) = B(J1,J2,1)*C(I,1)
    100 END DO

    return
    end subroutine faccl3
    subroutine faddcl3(a,b,c,iface1)

!     Collocate B with C and add to A on the surface IFACE1 of element IE.

!         A is a (NX,NY,NZ) data structure
!         B is a (NX,NY,NZ) data structure
!         C is a (NX,NY,IFACE) data structure
!         IFACE1 is in the preprocessor notation
!         IFACE  is the dssum notation.
!         29 Jan 1990 18:00 PST   PFF

    use size_m
    use topol
    DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)

!     Set up counters

    CALL DSSET(NX1,NY1,NZ1)
    IFACE  = EFACE1(IFACE1)
    JS1    = SKPDAT(1,IFACE)
    JF1    = SKPDAT(2,IFACE)
    JSKIP1 = SKPDAT(3,IFACE)
    JS2    = SKPDAT(4,IFACE)
    JF2    = SKPDAT(5,IFACE)
    JSKIP2 = SKPDAT(6,IFACE)

    I = 0
    DO 100 J2=JS2,JF2,JSKIP2
        DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A(J1,J2,1) = A(J1,J2,1) + B(J1,J2,1)*C(I,1)
    100 END DO

    return
    end subroutine faddcl3
!-----------------------------------------------------------------------
    subroutine sethlm (h1,h2,intloc)
     
!     Set the variable property arrays H1 and H2
!     in the Helmholtz equation.
!     (associated with variable IFIELD)
!     INTLOC =      integration type

    use size_m
    use input
    use soln
    use tstep

    real :: h1(1),h2(1)

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

    subroutine vprops
!-----------------------------------------------------------------------

!     Set material properties

!     Material type: 0 for default  (PARAM and PCOND/PRHOCP)
!                    1 for constant props;
!                    2 for fortran function;

!-----------------------------------------------------------------------
    use size_m
    use input
    use soln
    use tstep
    LOGICAL ::  IFKFLD,IFEFLD

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

!     Fill up property arrays every time step

!     First, check for turbulence models

#if 0
    IF (IFMODEL .AND. IFKEPS) THEN
        CALL TURBFLD (IFKFLD,IFEFLD)
        IF (IFKFLD)           CALL TPROPK
        IF (IFEFLD)           CALL TPROPE
        IF (IFKFLD .OR. IFEFLD) return
    endif
#endif

!...  No turbulence models, OR current field is not k or e.

    DO 1000 IEL=1,NEL
    
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
    
    1000 END DO

!     Turbulence models --- sum eddy viscosity/diffusivity
#if 0
    IF (IFMODEL .AND. (IFIELD == 1 .OR. IFIELD == 2)) &
    CALL TVISCOS
#endif

    return
    end subroutine vprops

    subroutine diagnos
    return
    end subroutine diagnos

    subroutine setsolv
    use size_m
    common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV
    IFSOLV = .FALSE. 
    return
    end subroutine setsolv

    subroutine chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binv, &
    vol,tol,nel)
!-------------------------------------------------------------------

!     Check that the tolerances are not too small for the CG-solver.
!     Important when calling the CG-solver (Gauss-Lobatto mesh) with
!     zero Neumann b.c.

!-------------------------------------------------------------------
    use size_m
    use eigen
    use input
    use mass
    common /cprint/ ifprint
    logical ::         ifprint
    common /ctmp0/ wa (lx1,ly1,lz1,lelt)

    dimension r1    (lx1,ly1,lz1,1) &
    , r2    (lx1,ly1,lz1,1) &
    , r3    (lx1,ly1,lz1,1) &
    , rmask1(lx1,ly1,lz1,1) &
    , rmask2(lx1,ly1,lz1,1) &
    , rmask3(lx1,ly1,lz1,1) &
    , rmult (lx1,ly1,lz1,1) &
    , binv  (lx1,ly1,lz1,1)

    NTOT1 = NX1*NY1*NZ1*NEL

    IF (EIGAA /= 0.0) THEN
        ACONDNO = EIGGA/EIGAA
    ELSE
        ACONDNO = 10.0
    endif

!     Check Single or double precision

    DELTA = 1.0E-9
    X     = 1.0 + DELTA
    Y     = 1.0
    DIFF  = ABS(X - Y)
    IF (DIFF == 0.0) EPS = 1.0E-6
    IF (DIFF > 0.0) EPS = 1.0E-13

    CALL OPDOT (WA,R1,R2,R3,R1,R2,R3,NTOT1)
    RINIT = GLSC3(WA,BINV,RMULT,NTOT1)
    RINIT = SQRT (RINIT/VOL)
    RMIN  = EPS*RINIT

    IF (TOL < RMIN) THEN
        TOLOLD = TOL
        TOL = RMIN
        IF (NID == 0 .AND. IFPRINT) &
        WRITE(6,*)'New CG1(stress)-tolerance (RINIT*epsm) = ',TOL,TOLOLD
    endif

    IF (NDIM == 2) THEN
        CALL ADD3 (WA,RMASK1,RMASK2,NTOT1)
    ELSE
        CALL ADD4 (WA,RMASK1,RMASK2,RMASK3,NTOT1)
    endif
    BCNEU1 = GLSC2 (WA,RMULT,NTOT1)
    BCNEU2 = (NDIM) * GLSUM(RMULT,NTOT1)
    BCTEST = ABS(BCNEU1 - BCNEU2)
    IF (BCTEST < 0.1) THEN
        IF (NDIM == 2) THEN
            CALL ADD3 (WA,R1,R2,NTOT1)
        ELSE
            CALL ADD4 (WA,R1,R2,R3,NTOT1)
        endif
        OTR    = GLSC2(WA,RMULT,NTOT1) / ( NDIM )
        TOLMIN = ABS(OTR) * ACONDNO
        IF (TOL < TOLMIN) THEN
            TOLOLD = TOL
            TOL = TOLMIN
            IF (NID == 0) &
            WRITE (6,*) 'New CG1(stress)-tolerance (OTR) = ',TOL,TOLOLD
        endif
    endif

    return
    end subroutine chktcgs
!-----------------------------------------------------------------------
    subroutine stnrate (u1,u2,u3,nel,matmod)

!     Compute strainrates

!     CAUTION : Stresses and strainrates share the same scratch commons

    use size_m
    use geom
    use input
    use tstep

    common /ctmp0/ exz(lx1*ly1*lz1*lelt) &
    , eyz(lx1*ly1*lz1*lelt)
    common /ctmp1/ exx(lx1*ly1*lz1*lelt) &
    , exy(lx1*ly1*lz1*lelt) &
    , eyy(lx1*ly1*lz1*lelt) &
    , ezz(lx1*ly1*lz1*lelt)

    dimension u1(lx1,ly1,lz1,1) &
    , u2(lx1,ly1,lz1,1) &
    , u3(lx1,ly1,lz1,1)

    NTOT1 = NX1*NY1*NZ1*NEL

    CALL RZERO3 (EXX,EYY,EZZ,NTOT1)
    CALL RZERO3 (EXY,EXZ,EYZ,NTOT1)

    CALL UXYZ  (U1,EXX,EXY,EXZ,NEL)
    CALL UXYZ  (U2,EXY,EYY,EYZ,NEL)
    IF (NDIM == 3) CALL UXYZ   (U3,EXZ,EYZ,EZZ,NEL)

    CALL INVCOL2 (EXX,JACM1,NTOT1)
    CALL INVCOL2 (EXY,JACM1,NTOT1)
    CALL INVCOL2 (EYY,JACM1,NTOT1)
     
    IF (IFAXIS) CALL AXIEZZ (U2,EYY,EZZ,NEL)

    IF (NDIM == 3) THEN
        CALL INVCOL2 (EXZ,JACM1,NTOT1)
        CALL INVCOL2 (EYZ,JACM1,NTOT1)
        CALL INVCOL2 (EZZ,JACM1,NTOT1)
    endif

    return
    end subroutine stnrate
!-----------------------------------------------------------------------
    subroutine stress (h1,h2,nel,matmod,ifaxis)

!     MATMOD.GE.0        Fluid material models
!     MATMOD.LT.0        Solid material models

!     CAUTION : Stresses and strainrates share the same scratch commons

    use size_m
    common /ctmp1/ txx(lx1,ly1,lz1,lelt) &
    , txy(lx1,ly1,lz1,lelt) &
    , tyy(lx1,ly1,lz1,lelt) &
    , tzz(lx1,ly1,lz1,lelt)
    common /ctmp0/ txz(lx1,ly1,lz1,lelt) &
    , tyz(lx1,ly1,lz1,lelt)
    common /scrsf/ t11(lx1,ly1,lz1,lelt) &
    , t22(lx1,ly1,lz1,lelt) &
    , t33(lx1,ly1,lz1,lelt) &
    , hii(lx1,ly1,lz1,lelt)

    DIMENSION H1(LX1,LY1,LZ1,1),H2(LX1,LY1,LZ1,1)
    LOGICAL :: IFAXIS

    NTOT1 = NX1*NY1*NZ1*NEL

    IF (MATMOD == 0) THEN

    !        Newtonian fluids

        CONST = 2.0
        CALL CMULT2 (HII,H1,CONST,NTOT1)
        CALL COL2   (TXX,HII,NTOT1)
        CALL COL2   (TXY,H1 ,NTOT1)
        CALL COL2   (TYY,HII,NTOT1)
        IF (IFAXIS .OR. NDIM == 3) CALL COL2 (TZZ,HII,NTOT1)
        IF (NDIM == 3) THEN
            CALL COL2 (TXZ,H1 ,NTOT1)
            CALL COL2 (TYZ,H1 ,NTOT1)
        endif
    
    ELSEIF (MATMOD == -1) THEN
    
    !        Elastic solids
    
        CONST = 2.0
        CALL ADD3S   (HII,H1,H2,CONST,NTOT1)
        CALL COPY    (T11,TXX,NTOT1)
        CALL COPY    (T22,TYY,NTOT1)
        CALL COL3    (TXX,HII,T11,NTOT1)
        CALL ADDCOL3 (TXX,H1 ,T22,NTOT1)
        CALL COL3    (TYY,H1 ,T11,NTOT1)
        CALL ADDCOL3 (TYY,HII,T22,NTOT1)
        CALL COL2    (TXY,H2     ,NTOT1)
        IF (IFAXIS .OR. NDIM == 3) THEN
            CALL COPY (T33,TZZ,NTOT1)
            CALL COL3    (TZZ,H1 ,T11,NTOT1)
            CALL ADDCOL3 (TZZ,H1 ,T22,NTOT1)
            CALL ADDCOL3 (TZZ,HII,T33,NTOT1)
            CALL ADDCOL3 (TXX,H1 ,T33,NTOT1)
            CALL ADDCOL3 (TYY,H1 ,T33,NTOT1)
        endif
        IF (NDIM == 3) THEN
            CALL COL2 (TXZ,H2     ,NTOT1)
            CALL COL2 (TYZ,H2     ,NTOT1)
        endif
    
    endif

    return
    end subroutine stress
!-----------------------------------------------------------------------
    subroutine aijuj (au1,au2,au3,nel,ifaxis)

    use size_m
    common /ctmp1/ txx(lx1,ly1,lz1,lelt) &
    , txy(lx1,ly1,lz1,lelt) &
    , tyy(lx1,ly1,lz1,lelt) &
    , tzz(lx1,ly1,lz1,lelt)
    common /ctmp0/ txz(lx1,ly1,lz1,lelt) &
    , tyz(lx1,ly1,lz1,lelt)

    DIMENSION AU1(LX1,LY1,LZ1,1) &
    , AU2(LX1,LY1,LZ1,1) &
    , AU3(LX1,LY1,LZ1,1)
    LOGICAL :: IFAXIS

    CALL TTXYZ (AU1,TXX,TXY,TXZ,NEL)
    CALL TTXYZ (AU2,TXY,TYY,TYZ,NEL)
    IF (IFAXIS)    CALL AXITZZ (AU2,TZZ,NEL)
    IF (NDIM == 3) CALL TTXYZ  (AU3,TXZ,TYZ,TZZ,NEL)

    return
    end subroutine aijuj
!-----------------------------------------------------------------------
    subroutine uxyz (u,ex,ey,ez,nel)

    use size_m
    use geom
    common /scrsf/ ur(lx1,ly1,lz1,lelt) &
    , us(lx1,ly1,lz1,lelt) &
    , ut(lx1,ly1,lz1,lelt)

    dimension u (lx1,ly1,lz1,1) &
    , ex(lx1,ly1,lz1,1) &
    , ey(lx1,ly1,lz1,1) &
    , ez(lx1,ly1,lz1,1)

    NTOT1 = NX1*NY1*NZ1*NEL

    CALL URST (U,UR,US,UT,NEL)

    CALL ADDCOL3 (EX,RXM1,UR,NTOT1)
    CALL ADDCOL3 (EX,SXM1,US,NTOT1)
    CALL ADDCOL3 (EY,RYM1,UR,NTOT1)
    CALL ADDCOL3 (EY,SYM1,US,NTOT1)

    IF (NDIM == 3) THEN
        CALL ADDCOL3 (EZ,RZM1,UR,NTOT1)
        CALL ADDCOL3 (EZ,SZM1,US,NTOT1)
        CALL ADDCOL3 (EZ,TZM1,UT,NTOT1)
        CALL ADDCOL3 (EX,TXM1,UT,NTOT1)
        CALL ADDCOL3 (EY,TYM1,UT,NTOT1)
    endif

    return
    end subroutine uxyz
    subroutine urst (u,ur,us,ut,nel)

    use size_m
    use geom
    use input

    DIMENSION U (LX1,LY1,LZ1,1) &
    , UR(LX1,LY1,LZ1,1) &
    , US(LX1,LY1,LZ1,1) &
    , UT(LX1,LY1,LZ1,1)

    DO 100 IEL=1,NEL
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
        CALL DDRST (U (1,1,1,IEL),UR(1,1,1,IEL), &
        US(1,1,1,IEL),UT(1,1,1,IEL))
    100 END DO

    return
    end subroutine urst
    subroutine ddrst (u,ur,us,ut)

    use size_m
    use dxyz

    DIMENSION U (LX1,LY1,LZ1) &
    , UR(LX1,LY1,LZ1) &
    , US(LX1,LY1,LZ1) &
    , UT(LX1,LY1,LZ1)

    NXY1 = NX1*NY1
    NYZ1 = NY1*NZ1

    CALL MXM (DXM1,NX1,U,NX1,UR,NYZ1)
    IF (NDIM == 2) THEN
        CALL MXM (U,NX1,DYTM1,NY1,US,NY1)
    ELSE
        DO 10 IZ=1,NZ1
            CALL MXM (U(1,1,IZ),NX1,DYTM1,NY1,US(1,1,IZ),NY1)
        10 END DO
        CALL MXM (U,NXY1,DZTM1,NZ1,UT,NZ1)
    endif

    return
    end subroutine ddrst
    subroutine axiezz (u2,eyy,ezz,nel)

    use size_m
    use geom

    DIMENSION U2 (LX1,LY1,LZ1,1) &
    , EYY(LX1,LY1,LZ1,1) &
    , EZZ(LX1,LY1,LZ1,1)

    NXYZ1  = NX1*NY1*NZ1

    DO 100 IEL=1,NEL
        IF ( IFRZER(IEL) ) THEN
            DO 200 IX=1,NX1
                EZZ(IX, 1,1,IEL) = EYY(IX,1,1,IEL)
                DO 200 IY=2,NY1
                    EZZ(IX,IY,1,IEL) = U2(IX,IY,1,IEL) / YM1(IX,IY,1,IEL)
            200 END DO
        ELSE
            CALL INVCOL3 (EZZ(1,1,1,IEL),U2(1,1,1,IEL),YM1(1,1,1,IEL), &
            NXYZ1)
        endif
    100 END DO

    return
    end subroutine axiezz
!-----------------------------------------------------------------------
    subroutine flush_io
    return
    end subroutine flush_io
!-----------------------------------------------------------------------
    subroutine axsf_e_3d(au,av,aw,u,v,w,h1,h2,ur,e)

!                                         du_i
!     Compute the gradient tensor G_ij := ----  ,  for element e
!                                         du_j

    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    real :: au(1),av(1),aw(1),u(1),v(1),w(1),h1(1),h2(1)

    parameter (l=lx1*ly1*lz1)
    real :: ur(l,ndim,ndim)

    integer :: e,p

    p = nx1-1      ! Polynomial degree
    n = nx1*ny1*nz1


    call local_grad3(ur(1,1,1),ur(1,2,1),ur(1,3,1),u,p,1,dxm1,dxtm1)
    call local_grad3(ur(1,1,2),ur(1,2,2),ur(1,3,2),v,p,1,dxm1,dxtm1)
    call local_grad3(ur(1,1,3),ur(1,2,3),ur(1,3,3),w,p,1,dxm1,dxtm1)

    do i=1,n

        u1 = ur(i,1,1)*rxm1(i,1,1,e) + ur(i,2,1)*sxm1(i,1,1,e) &
        + ur(i,3,1)*txm1(i,1,1,e)
        u2 = ur(i,1,1)*rym1(i,1,1,e) + ur(i,2,1)*sym1(i,1,1,e) &
        + ur(i,3,1)*tym1(i,1,1,e)
        u3 = ur(i,1,1)*rzm1(i,1,1,e) + ur(i,2,1)*szm1(i,1,1,e) &
        + ur(i,3,1)*tzm1(i,1,1,e)

        v1 = ur(i,1,2)*rxm1(i,1,1,e) + ur(i,2,2)*sxm1(i,1,1,e) &
        + ur(i,3,2)*txm1(i,1,1,e)
        v2 = ur(i,1,2)*rym1(i,1,1,e) + ur(i,2,2)*sym1(i,1,1,e) &
        + ur(i,3,2)*tym1(i,1,1,e)
        v3 = ur(i,1,2)*rzm1(i,1,1,e) + ur(i,2,2)*szm1(i,1,1,e) &
        + ur(i,3,2)*tzm1(i,1,1,e)

        w1 = ur(i,1,3)*rxm1(i,1,1,e) + ur(i,2,3)*sxm1(i,1,1,e) &
        + ur(i,3,3)*txm1(i,1,1,e)
        w2 = ur(i,1,3)*rym1(i,1,1,e) + ur(i,2,3)*sym1(i,1,1,e) &
        + ur(i,3,3)*tym1(i,1,1,e)
        w3 = ur(i,1,3)*rzm1(i,1,1,e) + ur(i,2,3)*szm1(i,1,1,e) &
        + ur(i,3,3)*tzm1(i,1,1,e)

        dj  = h1(i)*w3m1(i,1,1)*jacmi(i,e)
        s11 = dj*(u1 + u1)! S_ij
        s12 = dj*(u2 + v1)
        s13 = dj*(u3 + w1)
        s21 = dj*(v1 + u2)
        s22 = dj*(v2 + v2)
        s23 = dj*(v3 + w2)
        s31 = dj*(w1 + u3)
        s32 = dj*(w2 + v3)
        s33 = dj*(w3 + w3)

    !        Sum_j : (r_p/x_j) h1 J S_ij

        ur(i,1,1)=rxm1(i,1,1,e)*s11+rym1(i,1,1,e)*s12+rzm1(i,1,1,e)*s13
        ur(i,2,1)=sxm1(i,1,1,e)*s11+sym1(i,1,1,e)*s12+szm1(i,1,1,e)*s13
        ur(i,3,1)=txm1(i,1,1,e)*s11+tym1(i,1,1,e)*s12+tzm1(i,1,1,e)*s13

        ur(i,1,2)=rxm1(i,1,1,e)*s21+rym1(i,1,1,e)*s22+rzm1(i,1,1,e)*s23
        ur(i,2,2)=sxm1(i,1,1,e)*s21+sym1(i,1,1,e)*s22+szm1(i,1,1,e)*s23
        ur(i,3,2)=txm1(i,1,1,e)*s21+tym1(i,1,1,e)*s22+tzm1(i,1,1,e)*s23

        ur(i,1,3)=rxm1(i,1,1,e)*s31+rym1(i,1,1,e)*s32+rzm1(i,1,1,e)*s33
        ur(i,2,3)=sxm1(i,1,1,e)*s31+sym1(i,1,1,e)*s32+szm1(i,1,1,e)*s33
        ur(i,3,3)=txm1(i,1,1,e)*s31+tym1(i,1,1,e)*s32+tzm1(i,1,1,e)*s33

    enddo
    call local_grad3_t &
    (au,ur(1,1,1),ur(1,2,1),ur(1,3,1),p,1,dxm1,dxtm1,av)
    call local_grad3_t &
    (av,ur(1,1,2),ur(1,2,2),ur(1,3,2),p,1,dxm1,dxtm1,ur)
    call local_grad3_t &
    (aw,ur(1,1,3),ur(1,2,3),ur(1,3,3),p,1,dxm1,dxtm1,ur)

    do i=1,n
        au(i)=au(i) + h2(i)*bm1(i,1,1,e)*u(i)
        av(i)=av(i) + h2(i)*bm1(i,1,1,e)*v(i)
        aw(i)=aw(i) + h2(i)*bm1(i,1,1,e)*w(i)
    enddo

    return
    end subroutine axsf_e_3d
!-----------------------------------------------------------------------
    subroutine axsf_fast(au,av,aw,u,v,w,h1,h2,ifld)
    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    parameter (l=lx1*ly1*lz1)
    real :: au(l,1),av(l,1),aw(l,1),u(l,1),v(l,1),w(l,1),h1(l,1),h2(l,1)

    common /btmp0/ ur(l,ldim,ldim)

    integer :: e

    nel = nelfld(ifld)

    if (if3d) then
        do e=1,nel
            call axsf_e_3d(au(1,e),av(1,e),aw(1,e),u(1,e),v(1,e),w(1,e) &
            ,h1(1,e),h2(1,e),ur,e)
        enddo
    else
#if 0
        do e=1,nel
            call axsf_e_2d(au(1,e),av(1,e),u(1,e),v(1,e) &
            ,h1(1,e),h2(1,e),ur,e)
        enddo
#endif
    endif

    return
    end subroutine axsf_fast
!-----------------------------------------------------------------------
    subroutine ttxyz (ff,tx,ty,tz,nel)

    use size_m
    use dxyz
    use geom
    use input
    use tstep
    use wz_m

    DIMENSION TX(LX1,LY1,LZ1,1) &
    , TY(LX1,LY1,LZ1,1) &
    , TZ(LX1,LY1,LZ1,1) &
    , FF(LX1*LY1*LZ1,1)

    common /scrsf/ fr(lx1*ly1*lz1,lelt) &
    , fs(lx1*ly1*lz1,lelt) &
    , ft(lx1*ly1*lz1,lelt)
    real ::           wa(lx1,ly1,lz1,lelt)
    equivalence   (wa,ft)
    real :: ys(lx1)

    NXYZ1 = NX1*NY1*NZ1
    NTOT1 = NXYZ1*NEL

    CALL COL3    (FR,RXM1,TX,NTOT1)
    CALL ADDCOL3 (FR,RYM1,TY,NTOT1)
    CALL COL3    (FS,SXM1,TX,NTOT1)
    CALL ADDCOL3 (FS,SYM1,TY,NTOT1)

    IF (NDIM == 3) THEN
        CALL ADDCOL3 (FR,RZM1,TZ,NTOT1)
        CALL ADDCOL3 (FS,SZM1,TZ,NTOT1)
        CALL COL3    (FT,TXM1,TX,NTOT1)
        CALL ADDCOL3 (FT,TYM1,TY,NTOT1)
        CALL ADDCOL3 (FT,TZM1,TZ,NTOT1)
    endif

    IF (IFAXIS) THEN
        DO 100 IEL=1,NEL
            IF ( IFRZER(IEL) ) THEN
                CALL MXM (YM1(1,1,1,IEL),NX1,DATM1,NY1,YS,1)
                DO 120 IX=1,NX1
                    IY = 1
                    WA(IX,IY,1,IEL)=YS(IX)*W2AM1(IX,IY)
                    DO 120 IY=2,NY1
                        DNR = 1.0 + ZAM1(IY)
                        WA(IX,IY,1,IEL)=YM1(IX,IY,1,IEL)*W2AM1(IX,IY)/DNR
                120 END DO
            ELSE
                CALL COL3 (WA(1,1,1,IEL),YM1(1,1,1,IEL),W2CM1,NXYZ1)
            endif
        100 END DO
        CALL COL2 (FR,WA,NTOT1)
        CALL COL2 (FS,WA,NTOT1)
    else
        do 180 iel=1,nel
            call col2(fr(1,iel),w3m1,nxyz1)
            call col2(fs(1,iel),w3m1,nxyz1)
            call col2(ft(1,iel),w3m1,nxyz1)
        180 END DO
    endif


    DO 200 IEL=1,NEL
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
        CALL TTRST (FF(1,IEL),FR(1,IEL),FS(1,IEL), &
        FT(1,IEL),FR(1,IEL)) ! FR work array
    200 END DO

    return
    end subroutine ttxyz
!-----------------------------------------------------------------------
    subroutine ttrst (ff,fr,fs,ft,ta)

    use size_m
    use dxyz
    use tstep

    DIMENSION FF(LX1,LY1,LZ1) &
    , FR(LX1,LY1,LZ1) &
    , FS(LX1,LY1,LZ1) &
    , FT(LX1,LY1,LZ1) &
    , TA(LX1,LY1,LZ1)

    NXY1  = NX1*NY1
    NYZ1  = NY1*NZ1
    NXYZ1 = NXY1*NZ1

    CALL MXM (DXTM1,NX1,FR,NX1,FF,NYZ1)
    IF (NDIM == 2) THEN
        CALL MXM  (FS,NX1,DYM1,NY1,TA,NY1)
        CALL ADD2 (FF,TA,NXYZ1)
    ELSE
        DO 10 IZ=1,NZ1
            CALL MXM  (FS(1,1,IZ),NX1,DYM1,NY1,TA(1,1,IZ),NY1)
        10 END DO
        CALL ADD2 (FF,TA,NXYZ1)
        CALL MXM  (FT,NXY1,DZM1,NZ1,TA,NZ1)
        CALL ADD2 (FF,TA,NXYZ1)
    endif

    return
    end subroutine ttrst
!-----------------------------------------------------------------------
    subroutine axitzz (vfy,tzz,nel)

    use size_m
    use dxyz
    use geom
    use wz_m
    common /ctmp0/ phi(lx1,ly1)

    DIMENSION VFY(LX1,LY1,LZ1,1) &
    , TZZ(LX1,LY1,LZ1,1)

    NXYZ1 = NX1*NY1*NZ1

    DO 100 IEL=1,NEL
        CALL SETAXW1 ( IFRZER(IEL) )
        CALL COL4 (PHI,TZZ(1,1,1,IEL),JACM1(1,1,1,IEL),W3M1,NXYZ1)
        IF ( IFRZER(IEL) ) THEN
            DO 120 IX=1,NX1
                DO 120 IY=2,NY1
                    DNR = PHI(IX,IY)/( 1.0 + ZAM1(IY) )
                    DDS = WXM1(IX) * WAM1(1) * DATM1(IY,1) * &
                    JACM1(IX,1,1,IEL) * TZZ(IX,1,1,IEL)
                    VFY(IX,IY,1,IEL)=VFY(IX,IY,1,IEL) + DNR + DDS
            120 END DO
        ELSE
            CALL ADD2 (VFY(1,1,1,IEL),PHI,NXYZ1)
        endif
    100 END DO

    return
    end subroutine axitzz
!-----------------------------------------------------------------------
    subroutine setaxdy (ifaxdy)

    use size_m
    use dxyz

    LOGICAL :: IFAXDY

    IF (IFAXDY) THEN
        CALL COPY (DYM1 ,DAM1 ,NY1*NY1)
        CALL COPY (DYTM1,DATM1,NY1*NY1)
    ELSE
        CALL COPY (DYM1 ,DCM1 ,NY1*NY1)
        CALL COPY (DYTM1,DCTM1,NY1*NY1)
    endif

    return
    end subroutine setaxdy
!-----------------------------------------------------------------------
