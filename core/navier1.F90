!-----------------------------------------------------------------------
    subroutine ctolspl (tolspl,respr)

!     Compute the pressure tolerance

    use size_m
    use mass
    use tstep
    REAL ::           RESPR (LX2,LY2,LZ2,LELV)
    COMMON /SCRMG/ WORK  (LX1,LY1,LZ1,LELV)

    NTOT1 = NX1*NY1*NZ1*NELV
    CALL INVCOL3 (WORK,RESPR,BM1,NTOT1)
    CALL COL2    (WORK,RESPR,NTOT1)
    RINIT  = SQRT (GLSUM (WORK,NTOT1)/VOLVM1)
    IF (TOLPDF > 0.) THEN
        TOLSPL = TOLPDF
        TOLMIN = TOLPDF
    ELSE
        TOLSPL = TOLPS/DT
        TOLMIN = RINIT*PRELAX
    ENDIF
    IF (TOLSPL < TOLMIN) THEN
        TOLOLD = TOLSPL
        TOLSPL = TOLMIN
        IF (NID == 0) &
        WRITE (6,*) 'Relax the pressure tolerance ',TOLSPL,TOLOLD
    ENDIF
    return
    end subroutine ctolspl
!------------------------------------------------------------------------
    subroutine ortho (respr)

!     Orthogonalize the residual in the pressure solver with respect
!     to (1,1,...,1)T  (only if all Dirichlet b.c.).

    use size_m
    use geom
    use input
    use parallel
    use soln
    use tstep
    real :: respr (lx2,ly2,lz2,lelv)
    integer*8 :: ntotg,nxyz2

    nxyz2 = nx2*ny2*nz2
    ntot  = nxyz2*nelv
    ntotg = nxyz2*nelgv

    if (ifield == 1) then
        if (ifvcor) then
            rlam  = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)
        endif
    elseif (ifield == ifldmhd) then
        if (ifbcor) then
            rlam = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)
        endif
    else
        call exitti('ortho: unaccounted ifield = $',ifield)
    endif

    return
    end subroutine ortho
!------------------------------------------------------------------------
    subroutine cdabdtp (ap,wp,h1,h2,h2inv,intype)

!     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
!     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

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
    REAL ::           AP    (LX2,LY2,LZ2,1)
    REAL ::           WP    (LX2,LY2,LZ2,1)
    REAL ::           H1    (LX1,LY1,LZ1,1)
    REAL ::           H2    (LX1,LY1,LZ1,1)
    REAL ::           H2INV (LX1,LY1,LZ1,1)

    COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV) &
    ,             TA2 (LX1,LY1,LZ1,LELV) &
    ,             TA3 (LX1,LY1,LZ1,LELV) &
    ,             TB1 (LX1,LY1,LZ1,LELV) &
    ,             TB2 (LX1,LY1,LZ1,LELV) &
    ,             TB3 (LX1,LY1,LZ1,LELV)

    CALL OPGRADT (TA1,TA2,TA3,WP)
    IF ((INTYPE == 0) .OR. (INTYPE == -1)) THEN
        TOLHIN=TOLHS
        CALL OPHINV (TB1,TB2,TB3,TA1,TA2,TA3,H1,H2,TOLHIN,NMXH)
    ELSE
        if (ifanls) then
            dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
            CALL OPBINV1(TB1,TB2,TB3,TA1,TA2,TA3,dtbdi)
        else
            CALL OPBINV (TB1,TB2,TB3,TA1,TA2,TA3,H2INV)
        endif
    ENDIF
    CALL OPDIV  (AP,TB1,TB2,TB3)

    return
    end subroutine cdabdtp

!-----------------------------------------------------------------------
    subroutine opgrad (out1,out2,out3,inp)
!---------------------------------------------------------------------

!     Compute OUTi = Di*INP, i=1,2,3.
!     the gradient of the scalar field INP.
!     Note: OUTi is defined on the pressure mesh !!!

!---------------------------------------------------------------------
    use size_m
    use geom
    use input

    REAL :: OUT1 (LX2,LY2,LZ2,1)
    REAL :: OUT2 (LX2,LY2,LZ2,1)
    REAL :: OUT3 (LX2,LY2,LZ2,1)
    REAL :: INP  (LX1,LY1,LZ1,1)

    iflg = 0

    if (ifsplit .AND. .NOT. ifaxis) then
        call wgradm1(out1,out2,out3,inp,nelv) ! weak grad on FLUID mesh
        return
    endif

    NTOT2 = NX2*NY2*NZ2*NELV
    CALL MULTD (OUT1,INP,RXM2,SXM2,TXM2,1,iflg)
    CALL MULTD (OUT2,INP,RYM2,SYM2,TYM2,2,iflg)
    IF (NDIM == 3) &
    CALL MULTD (OUT3,INP,RZM2,SZM2,TZM2,3,iflg)

    return
    end subroutine opgrad
!-----------------------------------------------------------------------
    subroutine cdtp (dtx,x,rm2,sm2,tm2,isd)
!-------------------------------------------------------------

!     Compute DT*X (entire field)

!-------------------------------------------------------------
    use ctimer
    use size_m
    use dxyz
    use esolv
    use geom
    use input
    use ixyz
    use mass
    use wz_m

    real :: dtx  (lx1*ly1*lz1,lelv)
    real :: x    (lx2*ly2*lz2,lelv)
    real :: rm2  (lx2*ly2*lz2,lelv)
    real :: sm2  (lx2*ly2*lz2,lelv)
    real :: tm2  (lx2*ly2*lz2,lelv)

    common /ctmp1/ wx  (lx1*ly1*lz1) &
    ,             ta1 (lx1*ly1*lz1) &
    ,             ta2 (lx1*ly1*lz1) &
    ,             ta3 (lx1*ly1,lz1)

    REAL ::           DUAX(LX1)

    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV

    integer :: e

#ifndef NOTIMER
    if (icalld == 0) tcdtp=0.0
    icalld=icalld+1
    ncdtp=icalld
    etime1=dnekclock()
#endif

    nxyz1 = nx1*ny1*nz1
    nxyz2 = nx2*ny2*nz2
    nyz2  = ny2*nz2
    nxy1  = nx1*ny1

    n1    = nx1*ny1
    n2    = nx1*ny2

    do e=1,nelv

    !       Use the appropriate derivative- and interpolation operator in
    !       the y-direction (= radial direction if axisymmetric).
        if (ifaxis) then
            ny12   = ny1*ny2
            if (ifrzer(e)) then
                call copy (iym12,iam12,ny12)
                call copy (dym12,dam12,ny12)
                call copy (w3m2,w2am2,nxyz2)
            else
                call copy (iym12,icm12,ny12)
                call copy (dym12,dcm12,ny12)
                call copy (w3m2,w2cm2,nxyz2)
            endif
        endif
    
    !      Collocate with weights
    
        if(ifsplit) then
            call col3 (wx,bm1(1,1,1,e),x(1,e),nxyz1)
            call invcol2(wx,jacm1(1,1,1,e),nxyz1)
        else
            if ( .NOT. ifaxis) call col3 (wx,w3m2,x(1,e),nxyz2)
        
            if (ifaxis) then
                if (ifrzer(e)) then
                    call col3    (wx,x(1,e),bm2(1,1,1,e),nxyz2)
                    call invcol2 (wx,jacm2(1,1,1,e),nxyz2)
                else
                    call col3    (wx,w3m2,x(1,e),nxyz2)
                    call col2    (wx,ym2(1,1,1,e),nxyz2)
                endif
            endif
        endif
    
        if (ndim == 2) then
            write(*,*) "Whoops! cdtp"
#if 0
            if ( .NOT. ifdfrm(e) .AND. ifalgn(e)) then
            
                if (      ifrsxy(e) .AND. isd == 1  .OR. &
                 .NOT. ifrsxy(e) .AND. isd == 2) then
                
                    call col3 (ta1,wx,rm2(1,e),nxyz2)
                    call mxm  (dxtm12,nx1,ta1,nx2,ta2,nyz2)
                    call mxm  (ta2,nx1,iym12,ny2,dtx(1,e),ny1)
                else
                    call col3 (ta1,wx,sm2(1,e),nxyz2)
                    call mxm  (ixtm12,nx1,ta1,nx2,ta2,nyz2)
                    call mxm  (ta2,nx1,dym12,ny2,dtx(1,e),ny1)
                endif
            else
                call col3 (ta1,wx,rm2(1,e),nxyz2)
                call mxm  (dxtm12,nx1,ta1,nx2,ta2,nyz2)
                call mxm  (ta2,nx1,iym12,ny2,dtx(1,e),ny1)

                call col3 (ta1,wx,sm2(1,e),nxyz2)
                call mxm  (ixtm12,nx1,ta1,nx2,ta2,nyz2)
                call mxm  (ta2,nx1,dym12,ny2,ta1,ny1)

                call add2 (dtx(1,e),ta1,nxyz1)
            endif
#endif
        else
            if (ifsplit) then
                call col3 (ta1,wx,rm2(1,e),nxyz2)
                call mxm  (dxtm12,nx1,ta1,nx2,dtx(1,e),nyz2)
                call col3 (ta1,wx,sm2(1,e),nxyz2)
                i1 = 1
                i2 = 1
                do iz=1,nz2
                    call mxm  (ta1(i2),nx1,dym12,ny2,ta2(i1),ny1)
                    i1 = i1 + n1
                    i2 = i2 + n2
                enddo
                call add2 (dtx(1,e),ta2,nxyz1)
                call col3 (ta1,wx,tm2(1,e),nxyz2)
                call mxm  (ta1,nxy1,dzm12,nz2,ta2,nz1)
                call add2 (dtx(1,e),ta2,nxyz1)
            else
                write(*,*) "Whoops! cdtp"
#if 0
                call col3 (ta1,wx,rm2(1,e),nxyz2)
                call mxm  (dxtm12,nx1,ta1,nx2,ta2,nyz2)
                i1 = 1
                i2 = 1
                do iz=1,nz2
                    call mxm  (ta2(i2),nx1,iym12,ny2,ta1(i1),ny1)
                    i1 = i1 + n1
                    i2 = i2 + n2
                enddo
                call mxm  (ta1,nxy1,izm12,nz2,dtx(1,e),nz1)

                call col3 (ta1,wx,sm2(1,e),nxyz2)
                call mxm  (ixtm12,nx1,ta1,nx2,ta2,nyz2)
                i1 = 1
                i2 = 1
                do iz=1,nz2
                    call mxm  (ta2(i2),nx1,dym12,ny2,ta1(i1),ny1)
                    i1 = i1 + n1
                    i2 = i2 + n2
                enddo
                call mxm  (ta1,nxy1,izm12,nz2,ta2,nz1)
                call add2 (dtx(1,e),ta2,nxyz1)

                call col3 (ta1,wx,tm2(1,e),nxyz2)
                call mxm  (ixtm12,nx1,ta1,nx2,ta2,nyz2)
                i1 = 1
                i2 = 1
                do iz=1,nz2
                    call mxm  (ta2(i2),nx1,iym12,ny2,ta1(i1),ny1)
                    i1 = i1 + n1
                    i2 = i2 + n2
                enddo
                call mxm  (ta1,nxy1,dzm12,nz2,ta2,nz1)
                call add2 (dtx(1,e),ta2,nxyz1)
#endif
            endif

        endif
    
    !     If axisymmetric, add an extra diagonal term in the radial
    !     direction (only if solving the momentum equations and ISD=2)
    !     NOTE: NZ1=NZ2=1
    
    
        if(ifsplit) then

            if (ifaxis .AND. (isd == 4)) then
                call copy    (ta1,x(1,e),nxyz1)
                if (ifrzer(e)) THEN
                    call rzero(ta1, nx1)
                    call mxm  (x  (1,e),nx1,datm1,ny1,duax,1)
                    call copy (ta1,duax,nx1)
                endif
                call col2    (ta1,baxm1(1,1,1,e),nxyz1)
                call add2    (dtx(1,e),ta1,nxyz1)
            endif

        else

            if (ifaxis .AND. (isd == 2)) then
                call col3    (ta1,x(1,e),bm2(1,1,1,e),nxyz2)
                call invcol2 (ta1,ym2(1,1,1,e),nxyz2)
                call mxm     (ixtm12,nx1,ta1,nx2,ta2,ny2)
                call mxm     (ta2,nx1,iym12,ny2,ta1,ny1)
                call add2    (dtx(1,e),ta1,nxyz1)
            endif

        endif

    enddo

#ifndef NOTIMER
    tcdtp=tcdtp+(dnekclock()-etime1)
#endif
    return
    end subroutine cdtp

    subroutine multd (dx,x,rm2,sm2,tm2,isd,iflg)
!---------------------------------------------------------------------

!     Compute D*X
!     X    : input variable, defined on M1
!     DX   : output variable, defined on M2 (note: D is rectangular)
!     RM2 : RXM2, RYM2 or RZM2
!     SM2 : SXM2, SYM2 or SZM2
!     TM2 : TXM2, TYM2 or TZM2
!     ISD : spatial direction (x=1,y=2,z=3)
!     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)

!---------------------------------------------------------------------
    use ctimer
    use size_m
    use dxyz
    use esolv
    use geom
    use input
    use ixyz
    use mass
    use wz_m

    real ::           dx   (lx2*ly2*lz2,lelv)
    real ::           x    (lx1*ly1*lz1,lelv)
    real ::           rm2  (lx2*ly2*lz2,lelv)
    real ::           sm2  (lx2*ly2*lz2,lelv)
    real ::           tm2  (lx2*ly2*lz2,lelv)

    common /ctmp1/ ta1 (lx1*ly1*lz1) &
    ,             ta2 (lx1*ly1*lz1) &
    ,             ta3 (lx1*ly1*lz1)

    real ::           duax(lx1)

    common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
    logical :: ifdfrm, iffast, ifh2, ifsolv

    integer :: e

#ifndef NOTIMER
    if (icalld == 0) tmltd=0.0
    icalld=icalld+1
    nmltd=icalld
    etime1=dnekclock()
#endif

    nyz1  = ny1*nz1
    nxy2  = nx2*ny2
    nxyz1 = nx1*ny1*nz1
    nxyz2 = nx2*ny2*nz2

    n1    = nx2*ny1
    n2    = nx2*ny2

    do e=1,nelv

    !        Use the appropriate derivative- and interpolation operator in
    !        the y-direction (= radial direction if axisymmetric).
        if (ifaxis) then
            ny12   = ny1*ny2
            if (ifrzer(e)) then
                call copy (iytm12,iatm12,ny12)
                call copy (dytm12,datm12,ny12)
                call copy (w3m2,w2am2,nxyz2)
            else
                call copy (iytm12,ictm12,ny12)
                call copy (dytm12,dctm12,ny12)
                call copy (w3m2,w2cm2,nxyz2)
            endif
        endif

        if (ndim == 2) then
            write(*,*) "Whoops! multd"
#if 0
            if ( .NOT. ifdfrm(e) .AND. ifalgn(e)) then
            
                if (      ifrsxy(e) .AND. isd == 1  .OR. &
                 .NOT. ifrsxy(e) .AND. isd == 2) then
                    call mxm     (dxm12,nx2,x(1,e),nx1,ta1,nyz1)
                    call mxm     (ta1,nx2,iytm12,ny1,dx(1,e),ny2)
                    call col2    (dx(1,e),rm2(1,e),nxyz2)
                else
                    call mxm     (ixm12,nx2,x(1,e),nx1,ta1,nyz1)
                    call mxm     (ta1,nx2,dytm12,ny1,dx(1,e),ny2)
                    call col2    (dx(1,e),sm2(1,e),nxyz2)
                endif
            else
                call mxm     (dxm12,nx2,x(1,e),nx1,ta1,nyz1)
                call mxm     (ta1,nx2,iytm12,ny1,dx(1,e),ny2)
                call col2    (dx(1,e),rm2(1,e),nxyz2)
                call mxm     (ixm12,nx2,x(1,e),nx1,ta1,nyz1)
                call mxm     (ta1,nx2,dytm12,ny1,ta3,ny2)
                call addcol3 (dx(1,e),ta3,sm2(1,e),nxyz2)
            endif
#endif
        else  ! 3D

        !           if (ifsplit) then
        
        !             call mxm  (dxm12,nx2,x(1,e),nx1,dx(1,e),nyz1)
        !             call col2 (dx(1,e),rm2(1,e),nxyz2)
        !             i1=1
        !             i2=1
        !             do iz=1,nz1
        !                call mxm (x(1,e),nx2,dytm12,ny1,ta1(i2),ny2)
        !                i1=i1+n1
        !                i2=i2+n2
        !             enddo
        !             call addcol3 (dx(1,e),ta1,sm2(1,e),nxyz2)
        !             call mxm (x(1,e),nxy2,dztm12,nz1,ta1,nz2)
        !             call addcol3 (dx(1,e),ta1,tm2(1,e),nxyz2)

        !           else ! PN - PN-2

            call mxm (dxm12,nx2,x(1,e),nx1,ta1,nyz1)
            i1=1
            i2=1
            do iz=1,nz1
                call mxm (ta1(i1),nx2,iytm12,ny1,ta2(i2),ny2)
                i1=i1+n1
                i2=i2+n2
            enddo
            call mxm  (ta2,nxy2,iztm12,nz1,dx(1,e),nz2)
            call col2 (dx(1,e),rm2(1,e),nxyz2)

            call mxm  (ixm12,nx2,x(1,e),nx1,ta3,nyz1) ! reuse ta3 below
            i1=1
            i2=1
            do iz=1,nz1
                call mxm (ta3(i1),nx2,dytm12,ny1,ta2(i2),ny2)
                i1=i1+n1
                i2=i2+n2
            enddo
            call mxm     (ta2,nxy2,iztm12,nz1,ta1,nz2)
            call addcol3 (dx(1,e),ta1,sm2(1,e),nxyz2)

        !            call mxm (ixm12,nx2,x(1,e),nx1,ta1,nyz1) ! reuse ta3 from above
            i1=1
            i2=1
            do iz=1,nz1
                call mxm (ta3(i1),nx2,iytm12,ny1,ta2(i2),ny2)
                i1=i1+n1
                i2=i2+n2
            enddo
            call mxm (ta2,nxy2,dztm12,nz1,ta3,nz2)
            call addcol3 (dx(1,e),ta3,tm2(1,e),nxyz2)
        !           endif
        endif
    
    !        Collocate with the weights on the pressure mesh


        if(ifsplit) then
            call col2   (dx(1,e),bm1(1,1,1,e),nxyz1)
            call invcol2(dx(1,e),jacm1(1,1,1,e),nxyz1)
        else
            write(*,*) "Whoops! multd"
#if 0
            if ( .NOT. ifaxis) call col2 (dx(1,e),w3m2,nxyz2)
            if (ifaxis) then
                if (ifrzer(e)) then
                    call col2    (dx(1,e),bm2(1,1,1,e),nxyz2)
                    call invcol2 (dx(1,e),jacm2(1,1,1,e),nxyz2)
                else
                    call col2    (dx(1,e),w3m2,nxyz2)
                    call col2    (dx(1,e),ym2(1,1,1,e),nxyz2)
                endif
            endif
#endif
        endif

    !        If axisymmetric, add an extra diagonal term in the radial
    !        direction (ISD=2).
    !        NOTE: NZ1=NZ2=1

        if(ifsplit) then

            if (ifaxis .AND. (isd == 2) .AND. iflg == 1) then
                write(*,*) "Whoops! multd"
#if 0
                call copy    (ta3,x(1,e),nxyz1)
                if (ifrzer(e)) then
                    call rzero(ta3, nx1)
                    call mxm  (x(1,e),nx1,datm1,ny1,duax,1)
                    call copy (ta3,duax,nx1)
                endif
                call col2    (ta3,baxm1(1,1,1,e),nxyz1)
                call add2    (dx(1,e),ta3,nxyz2)
#endif
            endif

        else
            write(*,*) "Whoops! multd"
#if 0
            if (ifaxis .AND. (isd == 2)) then
                call mxm     (ixm12,nx2,x(1,e),nx1,ta1,ny1)
                call mxm     (ta1,nx2,iytm12,ny1,ta2,ny2)
                call col3    (ta3,bm2(1,1,1,e),ta2,nxyz2)
                call invcol2 (ta3,ym2(1,1,1,e),nxyz2)
                call add2    (dx(1,e),ta3,nxyz2)
            endif
#endif
        endif
    
    enddo

#ifndef NOTIMER
    tmltd=tmltd+(dnekclock()-etime1)
#endif
    return
    end subroutine multd

    subroutine ophinv (out1,out2,out3,inp1,inp2,inp3,h1,h2,tolh,nmxi)
!----------------------------------------------------------------------

!     OUT = (H1*A+H2*B)-1 * INP  (implicit)

!----------------------------------------------------------------------
    use size_m
    use input
    use soln
    use tstep
    REAL :: OUT1 (LX1,LY1,LZ1,1)
    REAL :: OUT2 (LX1,LY1,LZ1,1)
    REAL :: OUT3 (LX1,LY1,LZ1,1)
    REAL :: INP1 (LX1,LY1,LZ1,1)
    REAL :: INP2 (LX1,LY1,LZ1,1)
    REAL :: INP3 (LX1,LY1,LZ1,1)
    REAL :: H1   (LX1,LY1,LZ1,1)
    REAL :: H2   (LX1,LY1,LZ1,1)

    IMESH = 1

    if (ifstrs) then
#if 0
        MATMOD = 0
        if (ifield == ifldmhd) then
            CALL HMHZSF  ('NOMG',OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2, &
            B1MASK,B2MASK,B3MASK,VMULT, &
            TOLH,NMXI,MATMOD)
        else
            CALL HMHZSF  ('NOMG',OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2, &
            V1MASK,V2MASK,V3MASK,VMULT, &
            TOLH,NMXI,MATMOD)
        endif
#endif
    elseif (ifcyclic) then
#if 0
        matmod = 0
        if (ifield == ifldmhd) then
            call hmhzsf  ('bxyz',out1,out2,out3,inp1,inp2,inp3,h1,h2, &
            b1mask,b2mask,b3mask,vmult, &
            tolh,nmxi,matmod)
        else
            call hmhzsf  ('vxyz',out1,out2,out3,inp1,inp2,inp3,h1,h2, &
            v1mask,v2mask,v3mask,vmult, &
            tolh,nmxi,matmod)
        endif
#endif
    else
        if (ifield == ifldmhd) then
            CALL HMHOLTZ ('BX  ',OUT1,INP1,H1,H2,B1MASK,VMULT, &
            IMESH,TOLH,NMXI,1)
            CALL HMHOLTZ ('BY  ',OUT2,INP2,H1,H2,B2MASK,VMULT, &
            IMESH,TOLH,NMXI,2)
            IF (NDIM == 3) &
            CALL HMHOLTZ ('BZ  ',OUT3,INP3,H1,H2,B3MASK,VMULT, &
            IMESH,TOLH,NMXI,3)
        else
            CALL HMHOLTZ ('VELX',OUT1,INP1,H1,H2,V1MASK,VMULT, &
            IMESH,TOLH,NMXI,1)
            CALL HMHOLTZ ('VELY',OUT2,INP2,H1,H2,V2MASK,VMULT, &
            IMESH,TOLH,NMXI,2)
            IF (NDIM == 3) &
            CALL HMHOLTZ ('VELZ',OUT3,INP3,H1,H2,V3MASK,VMULT, &
            IMESH,TOLH,NMXI,3)
        endif
    ENDIF

    return
    end subroutine ophinv

    subroutine ophx (out1,out2,out3,inp1,inp2,inp3,h1,h2)
!----------------------------------------------------------------------

!     OUT = (H1*A+H2*B) * INP

!----------------------------------------------------------------------
    use size_m
    use input
    use soln
    REAL :: OUT1 (LX1,LY1,LZ1,1)
    REAL :: OUT2 (LX1,LY1,LZ1,1)
    REAL :: OUT3 (LX1,LY1,LZ1,1)
    REAL :: INP1 (LX1,LY1,LZ1,1)
    REAL :: INP2 (LX1,LY1,LZ1,1)
    REAL :: INP3 (LX1,LY1,LZ1,1)
    REAL :: H1   (LX1,LY1,LZ1,1)
    REAL :: H2   (LX1,LY1,LZ1,1)

    IMESH = 1

    IF (IFSTRS) THEN
#if 0
        MATMOD = 0
        CALL AXHMSF (OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2,MATMOD)
#endif
    ELSE
        CALL AXHELM (OUT1,INP1,H1,H2,IMESH,1)
        CALL AXHELM (OUT2,INP2,H1,H2,IMESH,2)
        IF (NDIM == 3) &
        CALL AXHELM (OUT3,INP3,H1,H2,IMESH,3)
    ENDIF

    return
    end subroutine ophx
!-----------------------------------------------------------------------
    subroutine opbinv (out1,out2,out3,inp1,inp2,inp3,h2inv)
!--------------------------------------------------------------------

!     Compute OUT = (H2*B)-1 * INP   (explicit)

!--------------------------------------------------------------------
    use size_m
    use input
    use mass
    use opctr
    use soln

    REAL :: OUT1  (1)
    REAL :: OUT2  (1)
    REAL :: OUT3  (1)
    REAL :: INP1  (1)
    REAL :: INP2  (1)
    REAL :: INP3  (1)
    REAL :: H2INV (1)



    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opbinv'
    endif

    call opmask  (inp1,inp2,inp3)
    call opdssum (inp1,inp2,inp3)

    NTOT=NX1*NY1*NZ1*NELV

    isbcnt = ntot*(1+ndim)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)

    call invcol3 (out1,bm1,h2inv,ntot)  ! this is expensive and should
    call dssum   (out1,nx1,ny1,nz1)     ! be changed (pff, 3/18/09)
    if (if3d) then
        do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
            out3(i)=inp3(i)*tmp
        enddo
    else
        do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
        enddo
    endif

    return
    end subroutine opbinv
!-----------------------------------------------------------------------
    subroutine opbinv1(out1,out2,out3,inp1,inp2,inp3,SCALE)
!--------------------------------------------------------------------

!     Compute OUT = (B)-1 * INP   (explicit)

!--------------------------------------------------------------------
    use size_m
    use input
    use mass
    use opctr
    use soln

    REAL :: OUT1  (1)
    REAL :: OUT2  (1)
    REAL :: OUT3  (1)
    REAL :: INP1  (1)
    REAL :: INP2  (1)
    REAL :: INP3  (1)



    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opbnv1'
    endif

    CALL OPMASK  (INP1,INP2,INP3)
    CALL OPDSSUM (INP1,INP2,INP3)

    NTOT=NX1*NY1*NZ1*NELV

    isbcnt = ntot*(1+ndim)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)

    IF (IF3D) THEN
        DO 100 I=1,NTOT
            TMP    =BINVM1(I,1,1,1)*scale
            OUT1(I)=INP1(I)*TMP
            OUT2(I)=INP2(I)*TMP
            OUT3(I)=INP3(I)*TMP
        100 END DO
    ELSE
        DO 200 I=1,NTOT
            TMP    =BINVM1(I,1,1,1)*scale
            OUT1(I)=INP1(I)*TMP
            OUT2(I)=INP2(I)*TMP
        200 END DO
    ENDIF

    return
    end subroutine opbinv1
!-----------------------------------------------------------------------
    subroutine uzprec (rpcg,rcg,h1m1,h2m1,intype,wp)
!--------------------------------------------------------------------

!     Uzawa preconditioner

!--------------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use parallel
    use tstep

    REAL ::           RPCG (LX2,LY2,LZ2,LELV)
    REAL ::           RCG  (LX2,LY2,LZ2,LELV)
    REAL ::           WP   (LX2,LY2,LZ2,LELV)
    REAL ::           H1M1 (LX1,LY1,LZ1,LELV)
    REAL ::           H2M1 (LX1,LY1,LZ1,LELV)
    COMMON /SCRCH/ H1M2 (LX2,LY2,LZ2,LELV) &
    ,             H2M2 (LX2,LY2,LZ2,LELV)

    integer :: kstep
    save    kstep
    data    kstep/-1/

    integer*8 :: ntotg,nxyz2


    if (solver_type == 'pdm') then
        write(*,*) "Oops, gfdm"
#if 0
        call gfdm_pres_solv(rpcg,rcg,h1m2,h2m2, .TRUE. ,0.0)
        return
#endif
    endif

    NTOT2 = NX2*NY2*NZ2*NELV
    if (istep /= kstep .AND. .NOT. ifanls) then
        kstep=istep
        DO 100 IE=1,NELV
            CALL MAP12 (H1M2(1,1,1,IE),H1M1(1,1,1,IE),IE)
            CALL MAP12 (H2M2(1,1,1,IE),H2M1(1,1,1,IE),IE)
        100 END DO
    endif

    IF (INTYPE == 0) THEN
        CALL COL3        (WP,RCG,H1M2,NTOT2)
        CALL COL3        (RPCG,WP,BM2INV,NTOT2)
    ELSEIF (INTYPE == -1) THEN
        CALL EPREC       (WP,RCG)
        CALL COL2        (WP,H2M2,NTOT2)
        CALL COL3        (RPCG,RCG,BM2INV,NTOT2)
        CALL COL2        (RPCG,H1M2,NTOT2)
        CALL ADD2        (RPCG,WP,NTOT2)
    ELSEIF (INTYPE == 1) THEN
        if (ifanls) then
            CALL EPREC2      (RPCG,RCG)
            DTBD = BD(1)/DT
            CALL cmult       (RPCG,DTBD,ntot2)
        else
            CALL EPREC2      (RPCG,RCG)
        !           CALL COL2        (RPCG,H2M2,NTOT2)
        endif
    ELSE
        CALL COPY        (RPCG,RCG,NTOT2)
    ENDIF

    call ortho (rpcg)

    return
    end subroutine uzprec

    subroutine eprec (z2,r2)
!----------------------------------------------------------------

!     Precondition the explicit pressure operator (E) with
!     a Neumann type (H1) Laplace operator: JT*A*J.
!     Invert A by conjugate gradient iteration or multigrid.

!     NOTE: SCRNS is used.

!----------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use parallel
    use soln
    use tstep
    REAL ::           Z2   (LX2,LY2,LZ2,LELV)
    REAL ::           R2   (LX2,LY2,LZ2,LELV)
    COMMON /SCRNS/ MASK (LX1,LY1,LZ1,LELV) &
    ,R1   (LX1,LY1,LZ1,LELV) &
    ,X1   (LX1,LY1,LZ1,LELV) &
    ,W2   (LX2,LY2,LZ2,LELV) &
    ,H1   (LX1,LY1,LZ1,LELV) &
    ,H2   (LX1,LY1,LZ1,LELV)
    REAL ::    MASK
    COMMON /CPRINT/ IFPRINT, IFHZPC
    LOGICAL ::         IFPRINT, IFHZPC
    integer*8 :: ntotg,nxyz
     
    nxyz   = nx1*ny1*nz1
    ntotg  = nxyz*nelgv
    ntot1  = nxyz*nelv
    ntot2  = nx2*ny2*nz2*nelv
    nfaces = 2*ndim

!     Set the tolerance for the preconditioner

    CALL COL3 (W2,R2,BM2INV,NTOT2)
    RINIT  = SQRT(GLSC2(W2,R2,NTOT2)/VOLVM2)
    EPS    = 0.02
    TOL    = EPS*RINIT
!     if (param(91).gt.0) tol=param(91)*rinit
!     if (param(91).lt.0) tol=-param(91)

    DO 100 IEL=1,NELV
        CALL MAP21E (R1(1,1,1,IEL),R2(1,1,1,IEL),IEL)
    100 END DO

    if (ifvcor) then
        otr1 = glsum (r1,ntot1)
        call rone  (x1,ntot1)
        c2   = -otr1/ntotg
        call add2s2 (r1,x1,c2,ntot1)
    
        OTR1 = GLSUM (R1,NTOT1)
        TOLMIN = 10.*ABS(OTR1)
        IF (TOL < TOLMIN) THEN
            if (nid == 0) &
            write(6,*) 'Resetting tol in EPREC:(old,new)',tol,tolmin
            TOL = TOLMIN
        ENDIF
    ENDIF

    CALL RONE    (H1,NTOT1)
    CALL RZERO   (H2,NTOT1)
    IFHZPC = .TRUE. 
    CALL HMHOLTZ ('PREC',X1,R1,H1,H2,PMASK,VMULT,IMESH,TOL,NMXH,1)
    IFHZPC = .FALSE. 

    DO 200 IEL=1,NELV
        CALL MAP12 (Z2(1,1,1,IEL),X1(1,1,1,IEL),IEL)
    200 END DO

    return
    end subroutine eprec

    subroutine convprn (iconv,rbnorm,rrpt,res,z,tol)
!-----------------------------------------------------------------
!                                               T
!     Convergence test for the pressure step;  r z

!-----------------------------------------------------------------
    use size_m
    use mass
    REAL ::           RES  (1)
    REAL ::           Z    (1)
    REAL ::           wrk1(2),wrk2(2)

    ntot2   = nx2*ny2*nz2*nelv
    wrk1(1) = vlsc21 (res,bm2inv,ntot2)  !  res*bm2inv*res
    wrk1(2) = vlsc2  (res,z     ,ntot2)  !  res*z
    call gop(wrk1,wrk2,'+  ',2)
    rbnorm  = sqrt(wrk1(1)/volvm2)
    rrpt    = wrk1(2)

!     CALL CONVPR (RCG,tolpss,ICONV,RNORM)
!     RRP1 = GLSC2 (RPCG,RCG,NTOT2)
!     RBNORM = SQRT (GLSC2 (BM2,TB,NTOT2)/VOLVM2)

    ICONV  = 0
    IF (RBNORM < TOL) ICONV=1
    return
    end subroutine convprn

    subroutine chktcg2 (tol,res,iconv)
!-------------------------------------------------------------------

!     Check that the tolerances are not too small for the CG-solver.
!     Important when calling the CG-solver (Gauss  mesh) with
!     all Dirichlet velocity b.c. (zero Neumann for the pressure).

!-------------------------------------------------------------------
    use size_m
    use geom
    use input
    use mass
    use tstep
    REAL ::           RES (LX2,LY2,LZ2,LELV)
    COMMON /CTMP0/ TA  (LX2,LY2,LZ2,LELV) &
    ,             TB  (LX2,LY2,LZ2,LELV)
    COMMON /CPRINT/ IFPRINT
    LOGICAL ::         IFPRINT

    ICONV = 0

!     Single or double precision???

    DELTA = 1.E-9
    X     = 1.+DELTA
    Y     = 1.
    DIFF  = ABS(X-Y)
    IF (DIFF == 0.) EPS = 1.E-5
    IF (DIFF > 0.) EPS = 1.E-10

!     Relaxed pressure iteration; maximum decrease in the residual (ER)

    IF (PRELAX /= 0.) EPS = PRELAX

    NTOT2 = NX2*NY2*NZ2*NELV
    CALL COL3     (TA,RES,BM2INV,NTOT2)
    CALL COL3     (TB,TA,TA,NTOT2)
    CALL COL2     (TB,BM2,NTOT2)
    RINIT = SQRT( GLSUM (TB,NTOT2)/VOLVM2 )
    IF (RINIT < TOL) THEN
        ICONV = 1
        return
    ENDIF
    IF (TOLPDF > 0.) THEN
        RMIN = TOLPDF
    ELSE
        RMIN  = EPS*RINIT
    ENDIF
    IF (TOL < RMIN) THEN
        TOLOLD = TOL
        TOL = RMIN
        IF (NID == 0 .AND. IFPRINT) WRITE (6,*) &
        'New CG2-tolerance (RINIT*10-5/10-10) = ',TOL,TOLOLD
    ENDIF
    IF (IFVCOR) THEN
        OTR = GLSUM (RES,NTOT2)
        TOLMIN = ABS(OTR)*100.
        IF (TOL < TOLMIN) THEN
            TOLOLD = TOL
            TOL = TOLMIN
            IF (NID == 0 .AND. IFPRINT) &
            WRITE (6,*) 'New CG2-tolerance (OTR) = ',TOLMIN,TOLOLD
        ENDIF
    ENDIF
    return
    end subroutine chktcg2

    subroutine dudxyz (du,u,rm1,sm1,tm1,jm1,imsh,isd)
!--------------------------------------------------------------

!     DU   - dU/dx or dU/dy or dU/dz
!     U    - a field variable defined on mesh 1
!     RM1  - dr/dx or dr/dy or dr/dz
!     SM1  - ds/dx or ds/dy or ds/dz
!     TM1  - dt/dx or dt/dy or dt/dz
!     JM1  - the Jacobian
!     IMESH - topology: velocity (1) or temperature (2) mesh

!--------------------------------------------------------------
    use size_m
    use dxyz
    use geom
    use input
    use tstep

    REAL ::  DU  (LX1,LY1,LZ1,1)
    REAL ::  U   (LX1,LY1,LZ1,1)
    REAL ::  RM1 (LX1,LY1,LZ1,1)
    REAL ::  SM1 (LX1,LY1,LZ1,1)
    REAL ::  TM1 (LX1,LY1,LZ1,1)
    REAL ::  JM1 (LX1,LY1,LZ1,1)

    COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
    LOGICAL :: IFDFRM, IFFAST, IFH2, IFSOLV

    REAL ::  DRST(LX1,LY1,LZ1)

    IF (imsh == 1) NEL = NELV
    IF (imsh == 2) NEL = NELT
    NXY1  = NX1*NY1
    NYZ1  = NY1*NZ1
    NXYZ1 = NX1*NY1*NZ1
    NTOT  = NXYZ1*NEL

    DO 1000 IEL=1,NEL
    
        IF (IFAXIS) CALL SETAXDY (IFRZER(IEL) )
    
        IF (NDIM == 2) THEN
            CALL MXM     (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
            CALL COL2    (DU(1,1,1,IEL),RM1(1,1,1,IEL),NXYZ1)
            CALL MXM     (U(1,1,1,IEL),NX1,DYTM1,NY1,DRST,NY1)
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,SM1(1,1,1,IEL),NXYZ1)
        ELSE
            CALL MXM   (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
            CALL COL2  (DU(1,1,1,IEL),RM1(1,1,1,IEL),NXYZ1)
            DO 20 IZ=1,NZ1
                CALL MXM  (U(1,1,IZ,IEL),NX1,DYTM1,NY1,DRST(1,1,IZ),NY1)
            20 END DO
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,SM1(1,1,1,IEL),NXYZ1)
            CALL MXM     (U(1,1,1,IEL),NXY1,DZTM1,NZ1,DRST,NZ1)
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,TM1(1,1,1,IEL),NXYZ1)
        ENDIF
    
    1000 END DO

!     CALL INVCOL2 (DU,JM1,NTOT)
    CALL COL2 (DU,JACMI,NTOT)

    return
    end subroutine dudxyz

!-----------------------------------------------------------------------
    subroutine makef
!---------------------------------------------------------------------

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)
!                      (2) driving force due to natural convection
!                      (3) convection term

!     !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
!              current time step is completed.

!----------------------------------------------------------------------
    use size_m
    use input
    use mass
    use soln
    use tstep

    call makeuf
!    if (ifnatc)                               call natconv
!    if (ifexplvis .AND. ifsplit)                call explstrs
    if (ifnav .AND. ( .NOT. ifchar))              call advab
!    if (ifmvbd)                               call admeshv
    if (iftran)                               call makeabf
    if ((iftran .AND. .NOT. ifchar) .OR. &
    (iftran .AND. .NOT. ifnav .AND. ifchar))   call makebdf
    if (ifnav .AND. ifchar .AND. ( .NOT. ifmvbd))   call advchar
#if 0
    if (ifmodel)                              call twallsh
#endif

    return
    end subroutine makef

    subroutine makeuf
!---------------------------------------------------------------------

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

!----------------------------------------------------------------------
    use size_m
    use mass
    use soln
    use tstep

    TIME = TIME-DT
    CALL NEKUF   (BFX,BFY,BFZ)
    CALL OPCOLV (BFX,BFY,BFZ,BM1)
    TIME = TIME+DT

    return
    end subroutine makeuf

    subroutine nekuf (f1,f2,f3)
    use size_m
    use nekuse
    use parallel
    REAL :: F1 (LX1,LY1,LZ1,LELV)
    REAL :: F2 (LX1,LY1,LZ1,LELV)
    REAL :: F3 (LX1,LY1,LZ1,LELV)
    CALL OPRZERO (F1,F2,F3)
    DO 100 IEL=1,NELV
        ielg = lglel(iel)
        DO 100 K=1,NZ1
            DO 100 J=1,NY1
                DO 100 I=1,NX1
                    CALL NEKASGN (I,J,K,IEL)
                    CALL USERF   (I,J,K,IELG)
                    F1(I,J,K,IEL) = FFX
                    F2(I,J,K,IEL) = FFY
                    F3(I,J,K,IEL) = FFZ
    100 END DO
    return
    end subroutine nekuf

    subroutine advab
!---------------------------------------------------------------

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

!---------------------------------------------------------------
    use size_m
    use mass
    use soln
    use tstep

    COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV) &
    ,             TA2 (LX1,LY1,LZ1,LELV) &
    ,             TA3 (LX1,LY1,LZ1,LELV)

    NTOT1 = NX1*NY1*NZ1*NELV
    CALL CONVOP  (TA1,VX)
    CALL CONVOP  (TA2,VY)
    CALL SUBCOL3 (BFX,BM1,TA1,NTOT1)
    CALL SUBCOL3 (BFY,BM1,TA2,NTOT1)
    IF (NDIM == 2) THEN
        CALL RZERO (TA3,NTOT1)
    ELSE
        CALL CONVOP  (TA3,VZ)
        CALL SUBCOL3 (BFZ,BM1,TA3,NTOT1)
    ENDIF


    return
    end subroutine advab
!-----------------------------------------------------------------------
    subroutine makebdf

!     Add contributions to F from lagged BD terms.

    use size_m
    use geom
    use input
    use mass
    use soln
    use tstep

    COMMON /SCRNS/ TA1(LX1,LY1,LZ1,LELV) &
    ,             TA2(LX1,LY1,LZ1,LELV) &
    ,             TA3(LX1,LY1,LZ1,LELV) &
    ,             TB1(LX1,LY1,LZ1,LELV) &
    ,             TB2(LX1,LY1,LZ1,LELV) &
    ,             TB3(LX1,LY1,LZ1,LELV) &
    ,             H2 (LX1,LY1,LZ1,LELV)

    NTOT1 = NX1*NY1*NZ1*NELV
    CONST = 1./DT
    CALL CMULT2(H2,vtrans(1,1,1,1,ifield),CONST,NTOT1)
    CALL OPCOLV3c (TB1,TB2,TB3,VX,VY,VZ,BM1,bd(2))

    DO 100 ILAG=2,NBD
        IF (IFGEOM) THEN
            CALL OPCOLV3c(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1), &
            VYLAG (1,1,1,1,ILAG-1), &
            VZLAG (1,1,1,1,ILAG-1), &
            BM1LAG(1,1,1,1,ILAG-1),bd(ilag+1))
        ELSE
            CALL OPCOLV3c(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1), &
            VYLAG (1,1,1,1,ILAG-1), &
            VZLAG (1,1,1,1,ILAG-1), &
            BM1                   ,bd(ilag+1))
        ENDIF
        CALL OPADD2  (TB1,TB2,TB3,TA1,TA2,TA3)
    100 END DO
    CALL OPADD2col (BFX,BFY,BFZ,TB1,TB2,TB3,h2)

    return
    end subroutine makebdf
!-----------------------------------------------------------------------
    subroutine makeabf
!-----------------------------------------------------------------------

!     Sum up contributions to kth order extrapolation scheme.
!     NOTE: rho^{n+1} should multiply all the Sum_q{beta_q} term
!           if rho is not constant!


!-----------------------------------------------------------------------
    use size_m
    use soln
    use tstep

    COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV) &
    ,             TA2 (LX1,LY1,LZ1,LELV) &
    ,             TA3 (LX1,LY1,LZ1,LELV)

    NTOT1 = NX1*NY1*NZ1*NELV

    AB0 = AB(1)
    AB1 = AB(2)
    AB2 = AB(3)
    CALL ADD3S2 (TA1,ABX1,ABX2,AB1,AB2,NTOT1)
    CALL ADD3S2 (TA2,ABY1,ABY2,AB1,AB2,NTOT1)
    CALL COPY   (ABX2,ABX1,NTOT1)
    CALL COPY   (ABY2,ABY1,NTOT1)
    CALL COPY   (ABX1,BFX,NTOT1)
    CALL COPY   (ABY1,BFY,NTOT1)
    CALL ADD2S1 (BFX,TA1,AB0,NTOT1)
    CALL ADD2S1 (BFY,TA2,AB0,NTOT1)
    CALL COL2   (BFX,VTRANS,NTOT1)          ! multiply by density
    CALL COL2   (BFY,VTRANS,NTOT1)
    IF (NDIM == 3) THEN
        CALL ADD3S2 (TA3,ABZ1,ABZ2,AB1,AB2,NTOT1)
        CALL COPY   (ABZ2,ABZ1,NTOT1)
        CALL COPY   (ABZ1,BFZ,NTOT1)
        CALL ADD2S1 (BFZ,TA3,AB0,NTOT1)
        CALL COL2   (BFZ,VTRANS,NTOT1)
    ENDIF

    return
    end subroutine makeabf

!-----------------------------------------------------------------------
    subroutine setab3 (ab0,ab1,ab2)

!     Set coefficients for 3rd order Adams-Bashforth scheme
!     (variable time step).

    use size_m
    use tstep

    IF (ISTEP <= 2) THEN
        AB0 = 1.
        AB1 = 0.
        AB2 = 0.
    ELSE
        DT0 = DTLAG(1)
        DT1 = DTLAG(2)
        DT2 = DTLAG(3)
        AB2  = (DT0/DT2)*((DT0/3.+DT1/2.)/(DT1+DT2))
        AB1  = -(DT0/DT1)*(0.5+(DT0/3.+DT1/2.)/DT2)
        AB0  = 1.-AB1-AB2
    ENDIF
    return
    end subroutine setab3

    subroutine setabbd (ab,dtlag,nab,nbd)
!-----------------------------------------------------------------------

!     Compute Adams-Bashforth coefficients (order NAB, less or equal to 3)

!     NBD .EQ. 1
!     Standard Adams-Bashforth coefficients

!     NBD .GT. 1
!     Modified Adams-Bashforth coefficients to be used in con-
!     junction with Backward Differentiation schemes (order NBD)

!-----------------------------------------------------------------------
    REAL :: AB(NAB),DTLAG(NAB)

    DT0 = DTLAG(1)
    DT1 = DTLAG(2)
    DT2 = DTLAG(3)

    IF ( NAB == 1 ) THEN
    
        AB(1) = 1.0
    
    ELSEIF ( NAB == 2 ) THEN
    
        DTA =  DT0/DT1
    
        IF ( NBD == 1 ) THEN
        
            AB(2) = -0.5*DTA
            AB(1) =  1.0 - AB(2)
        
        ELSEIF ( NBD == 2 ) THEN
        
            AB(2) = -DTA
            AB(1) =  1.0 - AB(2)
        
        ENDIF
    
    ELSEIF ( NAB == 3 ) THEN
    
        DTS =  DT1 + DT2
        DTA =  DT0 / DT1
        DTB =  DT1 / DT2
        DTC =  DT0 / DT2
        DTD =  DTS / DT1
        DTE =  DT0 / DTS
    
        IF ( NBD == 1 ) THEN
        
            AB(3) =  DTE*( 0.5*DTB + DTC/3. )
            AB(2) = -0.5*DTA - AB(3)*DTD
            AB(1) =  1.0 - AB(2) - AB(3)
        
        ELSEIF ( NBD == 2 ) THEN
        
            AB(3) =  2./3.*DTC*(1./DTD + DTE)
            AB(2) = -DTA - AB(3)*DTD
            AB(1) =  1.0 - AB(2) - AB(3)
        
        ELSEIF ( NBD == 3 ) THEN
        
            AB(3) =  DTE*(DTB + DTC)
            AB(2) = -DTA*(1.0 + DTB + DTC)
            AB(1) =  1.0 - AB(2) - AB(3)
        
        ENDIF
    
    ENDIF

    return
    end subroutine setabbd

    subroutine setbd (bd,dtbd,nbd)
!-----------------------------------------------------------------------

!     Compute bacward-differentiation coefficients of order NBD

!-----------------------------------------------------------------------
    PARAMETER (NDIM = 10)
    REAL :: BDMAT(NDIM,NDIM),BDRHS(NDIM)
    INTEGER :: IR(NDIM),IC(NDIM)
    REAL :: BD(1),DTBD(1)

    CALL RZERO (BD,10)
    IF (NBD == 1) THEN
        BD(1) = 1.
        BDF   = 1.
    ELSEIF (NBD >= 2) THEN
        NSYS = NBD+1
        CALL BDSYS (BDMAT,BDRHS,DTBD,NBD,NDIM)
        CALL LU    (BDMAT,NSYS,NDIM,IR,IC)
        CALL SOLVE (BDRHS,BDMAT,1,NSYS,NDIM,IR,IC)
        DO 30 I=1,NBD
            BD(I) = BDRHS(I)
        30 END DO
        BDF = BDRHS(NBD+1)
    ENDIF

!     Normalize

    DO 100 IBD=NBD,1,-1
        BD(IBD+1) = BD(IBD)
    100 END DO
    BD(1) = 1.
    DO 200 IBD=1,NBD+1
        BD(IBD) = BD(IBD)/BDF
    200 END DO
!     write(6,1) (bd(k),k=1,nbd+1)
!   1 format('bd:',1p8e13.5)

    return
    end subroutine setbd

    subroutine bdsys (a,b,dt,nbd,ndim)
    REAL :: A(NDIM,9),B(9),DT(9)
    CALL RZERO (A,NDIM**2)
    N = NBD+1
    DO 10 J=1,NBD
        A(1,J) = 1.
    10 END DO
    A(1,NBD+1) = 0.
    B(1) = 1.
    DO 20 J=1,NBD
        SUMDT = 0.
        DO 25 K=1,J
            SUMDT = SUMDT+DT(K)
        25 END DO
        A(2,J) = SUMDT
    20 END DO
    A(2,NBD+1) = -DT(1)
    B(2) = 0.
    DO 40 I=3,NBD+1
        DO 30 J=1,NBD
            SUMDT = 0.
            DO 35 K=1,J
                SUMDT = SUMDT+DT(K)
            35 END DO
            A(I,J) = SUMDT**(I-1)
        30 END DO
        A(I,NBD+1) = 0.
        B(I) = 0.
    40 END DO
    return
    end subroutine bdsys
    subroutine tauinit (tau,ilag)
!-------------------------------------------------------------------

!     Set initial time for subintegration

!-------------------------------------------------------------------
    use size_m
    use tstep
    TAU   = 0.
    DO 10 I=NBD,ILAG+1,-1
        TAU = TAU+DTLAG(I)
    10 END DO
    return
    end subroutine tauinit

    subroutine velconv (vxn,vyn,vzn,tau)
!--------------------------------------------------------------------

!     Compute convecting velocity field (linearization)

!--------------------------------------------------------------------
    use size_m
    use soln
    use tstep
    REAL :: VXN (LX1,LY1,LZ1,LELV)
    REAL :: VYN (LX1,LY1,LZ1,LELV)
    REAL :: VZN (LX1,LY1,LZ1,LELV)
    CALL VELCHAR (VX,VXN,VXLAG,NBD,TAU,DTLAG)
    CALL VELCHAR (VY,VYN,VYLAG,NBD,TAU,DTLAG)
    IF (NDIM == 3) &
    CALL VELCHAR (VZ,VZN,VZLAG,NBD,TAU,DTLAG)
    return
    end subroutine velconv

    subroutine frkconv (y,x,mask)
!--------------------------------------------------------------------

!     Evaluate right-hand-side for Runge-Kutta scheme in the case of
!     pure convection.

!--------------------------------------------------------------------
    use size_m
    use mass
    use tstep
    REAL :: Y    (LX1,LY1,LZ1,1)
    REAL :: X    (LX1,LY1,LZ1,1)
    REAL :: MASK (LX1,LY1,LZ1,1)

    IF (IMESH == 1) NEL=NELV
    IF (IMESH == 2) NEL=NELT
    NTOT1 = NX1*NY1*NZ1*NEL
    CALL CONVOP (Y,X)
    CALL COL2   (Y,BM1,NTOT1)
    CALL DSSUM  (Y,NX1,NY1,NZ1)
    IF (IMESH == 1) CALL COL2 (Y,BINVM1,NTOT1)
    IF (IMESH == 2) CALL COL2 (Y,BINTM1,NTOT1)
    CALL COL2   (Y,MASK,NTOT1)

    return
    end subroutine frkconv

    subroutine velchar (vel,vn,vlag,nbd,tau,dtbd)
!-----------------------------------------------------------------------

!     Compute linearized velocity field.

!-----------------------------------------------------------------------
    use size_m
    REAL :: VEL  (LX1,LY1,LZ1,LELV)
    REAL :: VN   (LX1,LY1,LZ1,LELV)
    REAL :: VLAG (LX1,LY1,LZ1,LELV,9)
    REAL :: DTBD (NBD)

    NTOT1 = NX1*NY1*NZ1*NELV
    IF (NBD == 1) THEN
        CALL COPY (VEL,VN,NTOT1)
        return
    ELSEIF (NBD == 2) THEN
        C1 = TAU/DTBD(2)
        C2 = 1.-C1
        CALL ADD3S2 (VEL,VN,VLAG,C1,C2,NTOT1)
    ELSEIF (NBD == 3) THEN
        F1 = TAU**2-DTBD(3)*TAU
        F2 = TAU**2-(DTBD(2)+DTBD(3))*TAU
        F3 = DTBD(2)*DTBD(3)
        F4 = DTBD(2)*(DTBD(2)+DTBD(3))
        R1 = F2/F3
        R2 = F1/F4
        C1 = R2
        C2 = -R1
        C3 = 1+R1-R2
        CALL ADD3S2 (VEL,VLAG(1,1,1,1,1),VLAG(1,1,1,1,2),C2,C3,NTOT1)
        CALL ADD2S2 (VEL,VN,C1,NTOT1)
    ELSE
        WRITE (6,*) 'Need higher order expansion in VELCHAR'
        call exitt
    ENDIF

    return
    end subroutine velchar

    subroutine lagvel
!-----------------------------------------------------------------------

!     Keep old velocity field(s)

!-----------------------------------------------------------------------
    use size_m
    use input
    use soln
    use tstep

    NTOT1 = NX1*NY1*NZ1*NELV

!      DO 100 ILAG=NBDINP-1,2,-1
    DO 100 ILAG=3-1,2,-1
        CALL COPY (VXLAG (1,1,1,1,ILAG),VXLAG (1,1,1,1,ILAG-1),NTOT1)
        CALL COPY (VYLAG (1,1,1,1,ILAG),VYLAG (1,1,1,1,ILAG-1),NTOT1)
        IF (NDIM == 3) &
        CALL COPY (VZLAG (1,1,1,1,ILAG),VZLAG (1,1,1,1,ILAG-1),NTOT1)
    100 END DO

    CALL OPCOPY (VXLAG,VYLAG,VZLAG,VX,VY,VZ)

    return
    end subroutine lagvel

    subroutine setordbd
!----------------------------------------------------------------------

!     Set up parameters for backward differentiation scheme.

!----------------------------------------------------------------------
    use size_m
    use input
    use tstep

!     IF (IFSPLIT .OR. NBDINP.EQ.0) THEN     undid hardwire, 3/6/92 pff
    IF ( NBDINP < 1) THEN
        NBD = 1
    ELSE
        IF ((ISTEP == 0) .OR. (ISTEP == 1))        NBD = 1
        IF ((ISTEP > 1) .AND. (ISTEP <= NBDINP))  NBD = ISTEP
        IF (ISTEP > NBDINP)                     NBD = NBDINP
    ENDIF

    return
    end subroutine setordbd

!-----------------------------------------------------------------------
!---------------------------------------------------------------
!> \brief Compute error norms of a (scalar) field variable X
!! defined on mesh 1 or mesh 2.
!! The error norms are normalized with respect to the volume
!! (except for Linf).
!---------------------------------------------------------------
subroutine normsc (h1,semi,l2,linf,x,imesh)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelt
  use size_m, only : nx1, ny1, nz1, nelv, nelt, ndim
  use mass, only : bm1, voltm1, volvm1
  implicit none

  real(DP) :: h1, semi, l2, linf
  real(DP) ::           X  (LX1,LY1,LZ1,1)
  integer :: imesh

  real(DP), allocatable, dimension(:,:,:,:) :: y, ta1, ta2
  REAL(DP) :: LENGTH, vol
  integer :: nel, nxyz1, ntot1
  real(DP), external :: glamax, glsum

  allocate(Y(LX1,LY1,LZ1,LELT), TA1(LX1,LY1,LZ1,LELT), TA2(LX1,LY1,LZ1,LELT))

  IF (IMESH == 1) THEN
      NEL = NELV
      VOL = VOLVM1
  ELSEIF (IMESH == 2) THEN
      NEL = NELT
      VOL = VOLTM1
  ENDIF
  LENGTH = VOL**(1./(NDIM))
  NXYZ1  = NX1*NY1*NZ1
  NTOT1  = NXYZ1*NEL

  H1     = 0.
  SEMI   = 0.
  L2     = 0.
  LINF   = 0.

  LINF = GLAMAX (X,NTOT1)

  CALL COL3   (TA1,X,X,NTOT1)
  CALL COL2   (TA1,BM1,NTOT1)
  L2   = GLSUM  (TA1,NTOT1)
  IF (L2 < 0.0) L2 = 0.

  CALL RONE   (TA1,NTOT1)
  CALL RZERO  (TA2,NTOT1)
  CALL AXHELM (Y,X,TA1,TA2,IMESH,1)
  CALL COL3   (TA1,Y,X,NTOT1)
  SEMI = GLSUM  (TA1,NTOT1)
  IF (SEMI < 0.0) SEMI = 0.

  H1   = SQRT((SEMI*LENGTH**2+L2)/VOL)
  SEMI = SQRT(SEMI/VOL)
  L2   = SQRT(L2/VOL)
  IF (H1 < 0.) H1 = 0.

  return
end subroutine normsc

    subroutine normvc (h1,semi,l2,linf,x1,x2,x3)
!---------------------------------------------------------------

!     Compute error norms of a (vector) field variable (X1,X2,X3)
!     defined on mesh 1 (velocity mesh only !)
!     The error norms are normalized with respect to the volume
!     (except for Linf).

!---------------------------------------------------------------
    use size_m
    use mass

    REAL ::           X1 (LX1,LY1,LZ1,1)
    REAL ::           X2 (LX1,LY1,LZ1,1)
    REAL ::           X3 (LX1,LY1,LZ1,1)
    COMMON /SCRMG/ Y1 (LX1,LY1,LZ1,LELT) &
    ,Y2 (LX1,LY1,LZ1,LELT) &
    ,Y3 (LX1,LY1,LZ1,LELT) &
    ,TA1(LX1,LY1,LZ1,LELT)
    COMMON /SCRCH/ TA2(LX1,LY1,LZ1,LELT)
    REAL :: H1,SEMI,L2,LINF
    REAL :: LENGTH

    IMESH  = 1
    NEL    = NELV
    VOL    = VOLVM1
    LENGTH = VOL**(1./(NDIM))
    NXYZ1  = NX1*NY1*NZ1
    NTOT1  = NXYZ1*NEL

    H1     = 0.
    SEMI   = 0.
    L2     = 0.
    LINF   = 0.

    CALL COL3 (TA1,X1,X1,NTOT1)
    CALL COL3 (TA2,X2,X2,NTOT1)
    CALL ADD2 (TA1,TA2,NTOT1)
    IF (NDIM == 3) THEN
        CALL COL3 (TA2,X3,X3,NTOT1)
        CALL ADD2 (TA1,TA2,NTOT1)
    ENDIF
    LINF = GLAMAX (TA1,NTOT1)
    LINF = SQRT( LINF )

    CALL COL3 (TA1,X1,X1,NTOT1)
    CALL COL3 (TA2,X2,X2,NTOT1)
    CALL ADD2 (TA1,TA2,NTOT1)
    IF (NDIM == 3) THEN
        CALL COL3 (TA2,X3,X3,NTOT1)
        CALL ADD2 (TA1,TA2,NTOT1)
    ENDIF
    CALL COL2 (TA1,BM1,NTOT1)
    L2 = GLSUM  (TA1,NTOT1)
    IF (L2 < 0.0) L2 = 0.

    CALL RONE  (TA1,NTOT1)
    CALL RZERO (TA2,NTOT1)
    CALL OPHX  (Y1,Y2,Y3,X1,X2,X3,TA1,TA2)
    CALL COL3  (TA1,Y1,X1,NTOT1)
    CALL COL3  (TA2,Y2,X2,NTOT1)
    CALL ADD2  (TA1,TA2,NTOT1)
    IF (NDIM == 3) THEN
        CALL COL3  (TA2,Y3,X3,NTOT1)
        CALL ADD2  (TA1,TA2,NTOT1)
    ENDIF
    SEMI = GLSUM (TA1,NTOT1)
    IF (SEMI < 0.0) SEMI = 0.

    H1   = SQRT((SEMI*LENGTH**2+L2)/VOL)
    SEMI = SQRT(SEMI/VOL)
    L2   = SQRT(L2  /VOL)
    IF (H1 < 0.) H1 = 0.

    return
    end subroutine normvc

    subroutine convsp (ifstsp)
    LOGICAL :: IFSTSP
    IFSTSP = .FALSE. 
    return
    end subroutine convsp

    subroutine opmask (res1,res2,res3)
!----------------------------------------------------------------------

!     Mask the residual arrays.

!----------------------------------------------------------------------
    use size_m
    use input
    use soln
    use tstep
    REAL :: RES1(1),RES2(1),RES3(1)

    NTOT1 = NX1*NY1*NZ1*NELV

!     sv=glsum(v3mask,ntot1)
!     sb=glsum(b3mask,ntot1)
!     write(6,*) istep,' ifld:',ifield,intype,sv,sb
    IF (IFSTRS) THEN
!max        CALL RMASK (RES1,RES2,RES3,NELV)
    ELSE
        if (ifield == ifldmhd) then
            CALL COL2 (RES1,B1MASK,NTOT1)
            CALL COL2 (RES2,B2MASK,NTOT1)
            IF (NDIM == 3) &
            CALL COL2 (RES3,B3MASK,NTOT1)
        else
            CALL COL2 (RES1,V1MASK,NTOT1)
            CALL COL2 (RES2,V2MASK,NTOT1)
            IF (NDIM == 3) &
            CALL COL2 (RES3,V3MASK,NTOT1)
        endif
    ENDIF

    return
    end subroutine opmask

    subroutine opadd2 (a1,a2,a3,b1,b2,b3)
    use size_m
    REAL :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL ADD2(A1,B1,NTOT1)
    CALL ADD2(A2,B2,NTOT1)
    IF(NDIM == 3)CALL ADD2(A3,B3,NTOT1)
    return
    end subroutine opadd2

    subroutine opsub2 (a1,a2,a3,b1,b2,b3)
    use size_m
    REAL :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL SUB2(A1,B1,NTOT1)
    CALL SUB2(A2,B2,NTOT1)
    IF(NDIM == 3)CALL SUB2(A3,B3,NTOT1)
    return
    end subroutine opsub2

    subroutine opcolv3(a1,a2,a3,b1,b2,b3,c)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1)
    REAL :: B1(1),B2(1),B3(1)
    REAL :: C (1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opcolv'
    endif

    isbcnt = ntot1*ndim
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            A1(I)=B1(I)*C(I)
            A2(I)=B2(I)*C(I)
            A3(I)=B3(I)*C(I)
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            A1(I)=B1(I)*C(I)
            A2(I)=B2(I)*C(I)
        200 END DO
    ENDIF
    return
    end subroutine opcolv3

    subroutine opcolv (a1,a2,a3,c)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1),C(1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opcolv'
    endif

    isbcnt = ntot1*ndim
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            A1(I)=A1(I)*C(I)
            A2(I)=A2(I)*C(I)
            A3(I)=A3(I)*C(I)
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            A1(I)=A1(I)*C(I)
            A2(I)=A2(I)*C(I)
        200 END DO
    ENDIF
    return
    end subroutine opcolv

    subroutine opcol2 (a1,a2,a3,b1,b2,b3)
    use size_m
    REAL :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL COL2(A1,B1,NTOT1)
    CALL COL2(A2,B2,NTOT1)
    IF(NDIM == 3)CALL COL2(A3,B3,NTOT1)
    return
    end subroutine opcol2

    subroutine opchsgn (a,b,c)
    use size_m
    REAL :: A(1),B(1),C(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL CHSIGN(A,NTOT1)
    CALL CHSIGN(B,NTOT1)
    IF(NDIM == 3)CALL CHSIGN(C,NTOT1)
    return
    end subroutine opchsgn

    subroutine opcopy (a1,a2,a3,b1,b2,b3)
    use size_m
    REAL :: A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL COPY(A1,B1,NTOT1)
    CALL COPY(A2,B2,NTOT1)
    IF(NDIM == 3)CALL COPY(A3,B3,NTOT1)
    return
    end subroutine opcopy

!-----------------------------------------------------------------------
    subroutine rotate_cyc(r1,r2,r3,idir)

    use size_m
    use geom
    use input
    use parallel
    use tstep

    real :: r1(lx1,ly1,lz1,1) &
    , r2(lx1,ly1,lz1,1) &
    , r3(lx1,ly1,lz1,1)

    integer :: e,f
    logical :: ifxy
     
     
!     (1) Face n-t transformation


    nface = 2*ndim
    do e=1,nelfld(ifield)
        do f=1,nface

            if(cbc(f,e,ifield) == 'P  ' .OR. cbc(f,e,ifield) == 'p  ')then
                write(*,*) "Oops: cyclic"
!                call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,f)
                if (idir == 1) then
                    k=0
                    do j2=js2,jf2,jskip2
                        do j1=js1,jf1,jskip1
                            k=k+1

                            dotprod = unx(k,1,f,e)*ym1(j1,j2,1,e) &
                            -uny(k,1,f,e)*xm1(j1,j2,1,e)
                            ifxy = .FALSE. 
                            if (abs(unz(k,1,f,e)) < 0.0001) ifxy = .TRUE. 

                            cost =  unx(k,1,f,e)
                            sint =  uny(k,1,f,e)
                            rnor = ( r1(j1,j2,1,e)*cost + r2(j1,j2,1,e)*sint )
                            rtn1 = (-r1(j1,j2,1,e)*sint + r2(j1,j2,1,e)*cost )

                            if (ifxy .AND. dotprod >= 0.0) then
                                r1(j1,j2,1,e) = rnor
                                r2(j1,j2,1,e) = rtn1
                            elseif (ifxy) then
                                r1(j1,j2,1,e) =-rnor
                                r2(j1,j2,1,e) =-rtn1
                            endif
                        enddo
                    enddo

                else    ! reverse rotate

                    k=0
                    do j2=js2,jf2,jskip2
                        do j1=js1,jf1,jskip1
                            k=k+1

                            dotprod = unx(k,1,f,e)*ym1(j1,j2,1,e) &
                            -uny(k,1,f,e)*xm1(j1,j2,1,e)
                            ifxy = .FALSE. 
                            if (abs(unz(k,1,f,e)) < 0.0001) ifxy = .TRUE. 

                            cost =  unx(k,1,f,e)
                            sint =  uny(k,1,f,e)
                            rnor = ( r1(j1,j2,1,e)*cost - r2(j1,j2,1,e)*sint )
                            rtn1 = ( r1(j1,j2,1,e)*sint + r2(j1,j2,1,e)*cost )

                            if(ifxy .AND. dotprod >= 0.0) then
                                r1(j1,j2,1,e) = rnor
                                r2(j1,j2,1,e) = rtn1
                            elseif (ifxy) then
                                r1(j1,j2,1,e) =-rnor
                                r2(j1,j2,1,e) =-rtn1
                            endif
                        enddo
                    enddo
                endif

            endif

        enddo
    enddo

    return
    end subroutine rotate_cyc
!-----------------------------------------------------------------------
    subroutine opdssum (a,b,c)! NOTE: opdssum works on FLUID/MHD arrays only!

    use size_m
    use geom
    use input
    use parallel
    use tstep

    real :: a(1),b(1),c(1)

    if (ifcyclic) then
        call rotate_cyc  (a,b,c,1)
        call vec_dssum   (a,b,c,nx1,ny1,nz1)
        call rotate_cyc  (a,b,c,0)
    else
        call vec_dssum   (a,b,c,nx1,ny1,nz1)
    endif

    return
    end subroutine opdssum
!-----------------------------------------------------------------------
    subroutine opdsop (a,b,c,op)! opdsop works on FLUID/MHD arrays only!

    use size_m
    use geom
    use geom
    use input
    use parallel
    use tstep

    real :: a(1),b(1),c(1)
    character(3) :: op

    if (ifcyclic) then

        if (op == '*  ' .OR. op == 'mul' .OR. op == 'MUL') then
            call vec_dsop    (a,b,c,nx1,ny1,nz1,op)
        else
            call rotate_cyc  (a,b,c,1)
            call vec_dsop    (a,b,c,nx1,ny1,nz1,op)
            call rotate_cyc  (a,b,c,0)
        endif

    else

        call vec_dsop    (a,b,c,nx1,ny1,nz1,op)

    endif

    return
    end subroutine opdsop
!-----------------------------------------------------------------------
    subroutine oprzero (a,b,c)
    use size_m
    REAL :: A(1),B(1),C(1)
    NTOT1=NX1*NY1*NZ1*NELV
    CALL RZERO(A,NTOT1)
    CALL RZERO(B,NTOT1)
    IF(NDIM == 3) CALL RZERO(C,NTOT1)
    return
    end subroutine oprzero
!-----------------------------------------------------------------------
    subroutine opcolv2c(a1,a2,a3,b1,b2,c)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1)
    REAL :: B1(1),B2(1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opcv2c'
    endif

    isbcnt = ntot1*(ndim+2)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            tmp = c*b1(i)*b2(i)
            A1(I)=A1(I)*tmp
            A2(I)=A2(I)*tmp
            A3(I)=A3(I)*tmp
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            tmp = c*b1(i)*b2(i)
            A1(I)=A1(I)*tmp
            A2(I)=A2(I)*tmp
        200 END DO
    ENDIF
    return
    end subroutine opcolv2c
!-----------------------------------------------------------------------
    subroutine opcolv2(a1,a2,a3,b1,b2)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1)
    REAL :: B1(1),B2(1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opclv2'
    endif

    isbcnt = ntot1*(ndim+1)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            tmp = b1(i)*b2(i)
            A1(I)=A1(I)*tmp
            A2(I)=A2(I)*tmp
            A3(I)=A3(I)*tmp
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            tmp = b1(i)*b2(i)
            A1(I)=A1(I)*tmp
            A2(I)=A2(I)*tmp
        200 END DO
    ENDIF
    return
    end subroutine opcolv2
!-----------------------------------------------------------------------
    subroutine opadd2col(a1,a2,a3,b1,b2,b3,c)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1)
    REAL :: B1(1),B2(1),B3(1),C(1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opa2cl'
    endif

    isbcnt = ntot1*(ndim*2)
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            A1(I)=A1(I)+b1(i)*c(i)
            A2(I)=A2(I)+b2(i)*c(i)
            A3(I)=A3(I)+b3(i)*c(i)
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            A1(I)=A1(I)+b1(i)*c(i)
            A2(I)=A2(I)+b2(i)*c(i)
        200 END DO
    ENDIF
    return
    end subroutine opadd2col
!-----------------------------------------------------------------------
    subroutine opcolv3c(a1,a2,a3,b1,b2,b3,c,d)
    use size_m
    use opctr
    REAL :: A1(1),A2(1),A3(1)
    REAL :: B1(1),B2(1),B3(1)
    REAL :: C (1)

    NTOT1=NX1*NY1*NZ1*NELV

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'opcv3c'
    endif

    isbcnt = ntot1*ndim*2
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

    IF (NDIM == 3) THEN
        DO 100 I=1,NTOT1
            A1(I)=B1(I)*C(I)*d
            A2(I)=B2(I)*C(I)*d
            A3(I)=B3(I)*C(I)*d
        100 END DO
    ELSE
        DO 200 I=1,NTOT1
            A1(I)=B1(I)*C(I)*d
            A2(I)=B2(I)*C(I)*d
        200 END DO
    ENDIF
    return
    end subroutine opcolv3c
!-----------------------------------------------------------------------

    subroutine uzawa (rcg,h1,h2,h2inv,intype,iter)
!-----------------------------------------------------------------------

!     Solve the pressure equation by (nested) preconditioned
!     conjugate gradient iteration.
!     INTYPE =  0  (steady)
!     INTYPE =  1  (explicit)
!     INTYPE = -1  (implicit)

!-----------------------------------------------------------------------
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
    COMMON  /CTOLPR/ DIVEX
    COMMON  /CPRINT/ IFPRINT
    LOGICAL ::          IFPRINT
    REAL ::             RCG  (LX2,LY2,LZ2,LELV)
    REAL ::             H1   (LX1,LY1,LZ1,LELV)
    REAL ::             H2   (LX1,LY1,LZ1,LELV)
    REAL ::             H2INV(LX1,LY1,LZ1,LELV)
    COMMON /SCRUZ/   WP   (LX2,LY2,LZ2,LELV) &
    ,               XCG  (LX2,LY2,LZ2,LELV) &
    ,               PCG  (LX2,LY2,LZ2,LELV) &
    ,               RPCG (LX2,LY2,LZ2,LELV)
     
    real*8 :: etime1,dnekclock
    integer*8 :: ntotg,nxyz2


    etime1 = dnekclock()
    DIVEX = 0.
    ITER  = 0

    CALL CHKTCG2 (TOLPS,RCG,ICONV)
    if (param(21) > 0 .AND. tolps > abs(param(21))) &
    TOLPS = abs(param(21))

!      IF (ICONV.EQ.1) THEN
!         IF (NID.EQ.0) WRITE(6,9999) ITER,DIVEX,TOLPS
!         return
!      ENDIF

    nxyz2 = nx2*ny2*nz2
    ntot2 = nxyz2*nelv
    ntotg = nxyz2*nelgv

    CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
    RRP1 = GLSC2 (RPCG,RCG,NTOT2)
    CALL COPY    (PCG,RPCG,NTOT2)
    CALL RZERO   (XCG,NTOT2)
    if (rrp1 == 0) return
    BETA = 0.
    div0=0.

    tolpss = tolps
    DO 1000 ITER=1,NMXP
    
    !        CALL CONVPR  (RCG,tolpss,ICONV,RNORM)
        call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)

        if (iter == 1)      div0   = rnorm
        if (param(21) < 0) tolpss = abs(param(21))*div0

        ratio = rnorm/div0
        IF (IFPRINT .AND. NID == 0) &
        WRITE (6,66) iter,tolpss,rnorm,div0,ratio,istep
        66 format(i5,1p4e12.5,i8,' Divergence')
    
        IF (ICONV == 1 .AND. iter > 1) GOTO 9000
    !        IF (ICONV.EQ.1.and.(iter.gt.1.or.istep.le.2)) GOTO 9000
    !        IF (ICONV.EQ.1) GOTO 9000
    !        if (ratio.le.1.e-5) goto 9000


        IF (ITER /= 1) THEN
            BETA = RRP1/RRP2
            CALL ADD2S1 (PCG,RPCG,BETA,NTOT2)
        ENDIF

        CALL CDABDTP  (WP,PCG,H1,H2,H2INV,INTYPE)
        PAP   = GLSC2 (PCG,WP,NTOT2)

        IF (PAP /= 0.) THEN
            ALPHA = RRP1/PAP
        ELSE
            pcgmx = glamax(pcg,ntot2)
            wp_mx = glamax(wp ,ntot2)
            ntot1 = nx1*ny1*nz1*nelv
            h1_mx = glamax(h1 ,ntot1)
            h2_mx = glamax(h2 ,ntot1)
            if (nid == 0) write(6,*) 'ERROR: pap=0 in uzawa.' &
            ,iter,pcgmx,wp_mx,h1_mx,h2_mx
            call exitt
        ENDIF
        CALL ADD2S2 (XCG,PCG,ALPHA,NTOT2)
        CALL ADD2S2 (RCG,WP,-ALPHA,NTOT2)

        if (iter == -1) then
            call convprn (iconv,rnrm1,rrpx,rcg,rpcg,tolpss)
            if (iconv == 1) then
                rnorm = rnrm1
                ratio = rnrm1/div0
                if (nid == 0) &
                write (6,66) iter,tolpss,rnrm1,div0,ratio,istep
                goto 9000
            endif
        endif

        call ortho(rcg)

        RRP2 = RRP1
        CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
    !        RRP1 = GLSC2 (RPCG,RCG,NTOT2)

    1000 END DO
    if (nid == 0) WRITE (6,3001) ITER,RNORM,tolpss
    if (istep > 20) CALL EMERXIT
    3001 FORMAT(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
    9000 CONTINUE

    divex = rnorm
    iter  = iter-1

    if (iter > 0) call copy (rcg,xcg,ntot2)
    call ortho(rcg)

    etime1 = dnekclock()-etime1
    IF (NID == 0) WRITE(6,9999) ISTEP,ITER,DIVEX,tolpss,div0,etime1
    9999 FORMAT(I10,' U-Press std. : ',I6,1p4E13.4)
    19999 FORMAT(I10,' U-Press 1.e-5: ',I6,1p4E13.4)


    return
    end subroutine uzawa
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine setmap(n1,nd)

    use size_m
    use dealias

    parameter(lx=80)
    real :: z1(lx),zd(lx),w(lx)

    if (n1 > lx .OR. nd > lx) then
        write(6,*)'ERROR: increase lx in setmap to max:',n1,nd
        call exitt
    endif

    call zwgll(z1,w,n1)
    call zwgll(zd,w,nd)
    call igllm(im1d,im1dt,z1,zd,n1,nd,n1,nd)
    call igllm(imd1,imd1t,zd,z1,nd,n1,nd,n1)

    return
    end subroutine setmap
!-----------------------------------------------------------------------
    subroutine transpose(a,lda,b,ldb)
    real :: a(lda,1),b(ldb,1)

    do j=1,ldb
        do i=1,lda
            a(i,j) = b(j,i)
        enddo
    enddo
    return
    end subroutine transpose
!-----------------------------------------------------------------------
    subroutine convop(conv,fi)

!     Compute the convective term CONV for a passive scalar field FI
!     using the skew-symmetric formulation.
!     The field variable FI is defined on mesh M1 (GLL) and
!     the velocity field is assumed given.

!     IMPORTANT NOTE: Use the scratch-arrays carefully!!!!!

!     The common-block SCRNS is used in CONV1 and CONV2.
!     The common-blocks CTMP0 and CTMP1 are also used as scratch-arrays
!     since there is no direct stiffness summation or Helmholtz-solves.

    use ctimer
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

!     Use the common blocks CTMP0 and CTMP1 as work space.

    COMMON /SCRCH/  CMASK1 (LX1,LY1,LZ1,LELV) &
    ,              CMASK2 (LX1,LY1,LZ1,LELV)
    COMMON /CTMP1/  MFI    (LX1,LY1,LZ1,LELV) &
    ,              DMFI   (LX1,LY1,LZ1,LELV) &
    ,              MDMFI  (LX1,LY1,LZ1,LELV)
    REAL ::   MFI,DMFI,MDMFI

!     Arrays in parameter list

    REAL ::    CONV (LX1,LY1,LZ1,1)
    REAL ::    FI   (LX1,LY1,LZ1,1)

#ifndef NOTIMER
    if (icalld == 0) tadvc=0.0
    icalld=icalld+1
    nadvc=icalld
    etime1=dnekclock()
#endif

    NXYZ1 = NX1*NY1*NZ1
    NTOT1 = NX1*NY1*NZ1*NELV
    NTOTZ = NX1*NY1*NZ1*nelfld(ifield)

    CALL RZERO  (CONV,NTOTZ)

    if (param(86) /= 0.0) then  ! skew-symmetric form
!max        call convopo(conv,fi)
        goto 100
    endif

!     write(6,*) istep,param(99),' CONVOP',ifpert
!     ip99 = param(99)
!     if (istep.gt.5) call exitti(' CONVOP dbg: $',ip99)

    if (param(99) == 2 .OR. param(99) == 3) then
!max        call conv1d(conv,fi)  !    use dealiased form
    elseif (param(99) == 4) then
        if (ifpert) then
            call convect_new (conv,fi, .FALSE. ,vx,vy,vz, .FALSE. )
        else
            call convect_new (conv,fi, .FALSE. ,vxd,vyd,vzd, .TRUE. )
        endif
        call invcol2     (conv,bm1,ntot1)  ! local mass inverse
    elseif (param(99) == 5) then
!max        call convect_cons(conv,fi, .FALSE. ,vx,vy,vz, .FALSE. )
        call invcol2     (conv,bm1,ntot1)  ! local mass inverse
    else
!max        call conv1 (conv,fi)  !    use the convective form
    endif

    100 continue

#ifndef NOTIMER
    tadvc=tadvc+(dnekclock()-etime1)
#endif

    return
    end subroutine convop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine opdiv(outfld,inpx,inpy,inpz)
!---------------------------------------------------------------------

!     Compute OUTFLD = SUMi Di*INPi,
!     the divergence of the vector field (INPX,INPY,INPZ)

!---------------------------------------------------------------------
    use size_m
    use geom
    real :: outfld (lx2,ly2,lz2,1)
    real :: inpx   (lx1,ly1,lz1,1)
    real :: inpy   (lx1,ly1,lz1,1)
    real :: inpz   (lx1,ly1,lz1,1)
    common /ctmp0/ work (lx2,ly2,lz2,lelv)

    iflg = 1

    ntot2 = nx2*ny2*nz2*nelv
    call multd (work,inpx,rxm2,sxm2,txm2,1,iflg)
    call copy  (outfld,work,ntot2)
    call multd (work,inpy,rym2,sym2,tym2,2,iflg)
    call add2  (outfld,work,ntot2)
    if (ndim == 3) then
        call multd (work,inpz,rzm2,szm2,tzm2,3,iflg)
        call add2  (outfld,work,ntot2)
    endif

    return
    end subroutine opdiv

!-----------------------------------------------------------------------
    subroutine opgradt(outx,outy,outz,inpfld)
!------------------------------------------------------------------------

!     Compute DTx, DTy, DTz of an input field INPFLD

!-----------------------------------------------------------------------
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
    real :: outx   (lx1,ly1,lz1,1)
    real :: outy   (lx1,ly1,lz1,1)
    real :: outz   (lx1,ly1,lz1,1)
    real :: inpfld (lx2,ly2,lz2,1)

    call cdtp (outx,inpfld,rxm2,sxm2,txm2,1)
    call cdtp (outy,inpfld,rym2,sym2,tym2,2)
    if (ndim == 3) &
    call cdtp (outz,inpfld,rzm2,szm2,tzm2,3)

    return
    end subroutine opgradt
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine set_PNDoi(Pt,P,LkNt,N,D)

    use size_m   ! for write stmt


!     Set up operators for overintegration and interpolation

    integer :: N,D
    real ::    Pt(N,D),P(D,N),LkNt(N,0:N-1)

    parameter(lx=80)
    real :: zN(lx),zD(lx),wN(lx),wD(lx)

!     Compute Lagrangian interpolant points

    call zwgll(zN,wN,N)
    call zwgll(zD,wD,D)

    if (nid == 0) write(6,*) 'dealias, pndoi:',N,D
    call IGLLM (P,Pt,ZN,ZD,N,D,N,D)

    do j=1,D
        do i=1,N
            Pt(i,j) = wD(j)*Pt(i,j)/wN(i)
        enddo
    enddo
    return
    end subroutine set_PNDoi
!-----------------------------------------------------------------------
    subroutine wgradm1(ux,uy,uz,u,nel) ! weak form of grad

!     Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)

    use size_m
    use dxyz
    use geom
    use input
    use tstep
    use wz_m

    parameter (lxyz=lx1*ly1*lz1)
    real :: ux(lxyz,1),uy(lxyz,1),uz(lxyz,1),u(lxyz,1)

    common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)

    integer :: e

    N = nx1-1
    do e=1,nel
        if (if3d) then
            call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
            do i=1,lxyz
                ux(i,e) = w3m1(i,1,1)*(ur(i)*rxm1(i,1,1,e) &
                + us(i)*sxm1(i,1,1,e) &
                + ut(i)*txm1(i,1,1,e) )
                uy(i,e) = w3m1(i,1,1)*(ur(i)*rym1(i,1,1,e) &
                + us(i)*sym1(i,1,1,e) &
                + ut(i)*tym1(i,1,1,e) )
                uz(i,e) = w3m1(i,1,1)*(ur(i)*rzm1(i,1,1,e) &
                + us(i)*szm1(i,1,1,e) &
                + ut(i)*tzm1(i,1,1,e) )
            enddo
        else
#if 0
            if (ifaxis) then
                call setaxdy (ifrzer(e))  ! reset dytm1
                call setaxw1 (ifrzer(e))  ! reset w3m1
            endif

            call local_grad2(ur,us,u,N,e,dxm1,dytm1)

            do i=1,lxyz
                ux(i,e) =w3m1(i,1,1)*(ur(i)*rxm1(i,1,1,e) &
                + us(i)*sxm1(i,1,1,e) )
                uy(i,e) =w3m1(i,1,1)*(ur(i)*rym1(i,1,1,e) &
                + us(i)*sym1(i,1,1,e) )
            enddo
#endif
        endif

    enddo

    return
    end subroutine wgradm1
!-----------------------------------------------------------------------
    subroutine wlaplacian(out,a,diff,ifld)

!     compute weak form of the laplacian operator including the boundary
!     contribution

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

    real :: out(1),a(1),diff(1)
    real :: wrk(lx1,ly1,lz1,lelt)
    real :: h2(lx1,ly1,lz1,lelt)

    ntot = nx1*ny1*nz1*nelfld(ifld)
    if ( .NOT. iftmsh(ifld)) imesh = 1
    if (     iftmsh(ifld)) imesh = 2

    call rzero(h2,ntot)

    ifield_ = ifield
    ifield = ifld

    call bcneusc(out,1)
    call axhelm(wrk,a,diff,h2,imesh,1)
    call sub2 (out,wrk,ntot)
     
    ifield = ifield_

    return
    end subroutine wlaplacian
!-----------------------------------------------------------------------
