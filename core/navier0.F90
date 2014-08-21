    SUBROUTINE ESOLVER (RES,H1,H2,H2INV,INTYPE)
!---------------------------------------------------------------------

!     Choose E-solver

!--------------------------------------------------------------------
    use ctimer
    use size_m
    use esolv
    use input

    REAL :: RES   (LX2,LY2,LZ2,LELV)
    REAL :: H1    (LX1,LY1,LZ1,LELV)
    REAL :: H2    (LX1,LY1,LZ1,LELV)
    REAL :: H2INV (LX1,LY1,LZ1,LELV)
    common /scruz/ wk1(lx2*ly2*lz2*lelv) &
    , wk2(lx2*ly2*lz2*lelv) &
    , wk3(lx2*ly2*lz2*lelv)

    real :: kwave2

    if (icalld == 0) teslv=0.0
    icalld=icalld+1
    neslv=icalld
    etime1=dnekclock()

!     write(6,*) solver_type,' solver type',iesolv
    if (iesolv == 1) then
        if (solver_type == 'fdm') then
            write(*,*) "Oops, gfdm"
#if 0
            ntot2 = nx2*ny2*nz2*nelv
            kwave2 = 0.
            call gfdm_pres_solv  (wk1,res,wk2,wk3,kwave2)
            call copy            (res,wk1,ntot2)
#endif
        else
            if (param(42) == 1 .OR. solver_type == 'pdm') then
                CALL UZAWA (RES,H1,H2,H2INV,INTYPE,ICG)
            else
                call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
            endif
        endif
    else
        WRITE(6,*) 'ERROR: E-solver does not exist',iesolv
        WRITE(6,*) 'Stop in ESOLVER'
        CALL EXITT
    ENDIF

    teslv=teslv+(dnekclock()-etime1)

    RETURN
    END SUBROUTINE ESOLVER
    SUBROUTINE ESTRAT
!---------------------------------------------------------------------------

!     Decide strategy for E-solver

!---------------------------------------------------------------------------
    use size_m
    INCLUDE 'TOTAL'

    IESOLV = 1
    if (ifsplit) iesolv=0

    solver_type='itr'
    if (param(116) /= 0) solver_type='fdm'
!     if (param(90).ne.0)  solver_type='itn'

!     The following change recognizes that geometry is logically
!     tensor-product, but deformed:  pdm = Preconditioner is fdm

    if (param(59) /= 0 .AND. solver_type == 'fdm') solver_type='pdm'

    if (istep < 2 .AND. nid == 0) write(6,10) iesolv,solver_type
    10 format(2X,'E-solver strategy: ',I2,1X,A)




    RETURN
    END SUBROUTINE ESTRAT
    SUBROUTINE EINIT
!-----------------------------------------------------------------------------

!     Initialize E-solver

!-----------------------------------------------------------------------------
    use size_m
    use esolv
    INCLUDE 'SOLN'
    INCLUDE 'TSTEP'
    COMMON /SCRHI/  H2INV (LX1,LY1,LZ1,LELV)
    LOGICAL :: IFNEWP

    CALL ESTRAT
    RETURN
    END SUBROUTINE EINIT
!-----------------------------------------------------------------------
    subroutine dmp_map(imap)

!     Dump map file and element center point

    use size_m
    include 'TOTAL'

    common /ivrtx/ vertex ((2**ldim)*lelt)
    common /scruz/ xbar(ldim,lelt),ibar(lelt)
    integer :: vertex
    integer :: imap(nelgt)

    integer :: e,eg

    nxb = (nx1+1)/2
    nyb = (ny1+1)/2
    nzb = (nz1+1)/2
          
    do e=1,nelt
        xbar(ndim,e) = zm1(nxb,nyb,nzb,e)
        xbar(1   ,e) = xm1(nxb,nyb,nzb,e)
        xbar(2   ,e) = ym1(nxb,nyb,nzb,e)
        eg           = lglel(e)
        ibar(e)      = imap(eg)
    enddo
    call p_outvec_ir(ibar,xbar,ndim,'mpxyz.dat')

    return
    end subroutine dmp_map
!-----------------------------------------------------------------------
    subroutine p_outvec_ir(ia,a,lda,name9)
    use size_m
    include 'TOTAL'

    integer :: ia(1)
    real ::    a(lda,1)
    character(9) :: name9


    parameter (lbuf=50)
    common /scbuf/ buf(lbuf)
    integer :: ibuf(10),e,eg
    equivalence (buf,ibuf)

    if (nid == 0) then
        open(unit=49,file=name9)
        write(6,*) 'Opening ',name9,' in p_outveci. lda=',lda
    endif

    len = wdsize*(lda+1)
    dum = 0.

    do eg=1,nelgt

        mid   = gllnid(eg)
        e     = gllel (eg)
        mtype = 2000+eg

        if (nid == 0) then
            if (mid == 0) then
                call icopy(buf(1),ia(e),1)
                call  copy(buf(2),a(1,e),lda)
            else
                call csend (mtype,dum,wdsize,mid,nullpid)
                call crecv (mtype,buf,len)
            endif
            write(49,49) mid,ibuf(1),(buf(k+1),k=1,lda)
            49 format(2i12,1p3e16.7)
        elseif (nid == mid) then
            call icopy(buf(1),ia(e),1)
            call  copy(buf(2),a(1,e),lda)
            call crecv (mtype,dum,wdsize)
            call csend (mtype,buf,len,node0,nullpid)
        endif
    enddo

    if (nid == 0) then
        close(49)
        write(6,*) 'Done writing to ',name9,' p_outveci.'
    endif

    call nekgsync()

    return
    end subroutine p_outvec_ir
!-----------------------------------------------------------------------
