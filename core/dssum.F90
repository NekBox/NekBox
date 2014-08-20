    subroutine setupds(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
    use size_m
    include 'INPUT'
    include 'PARALLEL'
    include 'NONCON'

    integer :: gs_handle
    integer :: vertex(1)
    integer*8 :: glo_num(1),ngv

    common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

    t0 = dnekclock()

!     Global-to-local mapping for gs
    call set_vert(glo_num,ngv,nx,nel,vertex, .FALSE. )

!     Initialize gather-scatter code
    ntot      = nx*ny*nz*nel
    call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

!     call gs_chkr(glo_num)

    t1 = dnekclock() - t0
    if (nid == 0) then
        write(6,1) t1,gs_handle,nx,ngv,melg
        1 format('   setupds time',1pe11.4,' seconds ',2i3,2i12)
    endif

    return
    end subroutine setupds
!-----------------------------------------------------------------------
    subroutine dssum(u,nx,ny,nz)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)

    ifldt = ifield
!     if (ifldt.eq.0)       ifldt = 1
    if (ifldt == ifldmhd) ifldt = 1
!     write(6,*) ifldt,ifield,gsh_fld(ifldt),imesh,' ifldt'

    if (ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) then
        tdsmx=0.
        tdsmn=0.
    endif
    icalld=icalld+1
    etime1=dnekclock()
#endif


!                 T         ~  ~T  T
!     Implement QQ   :=   J Q  Q  J


!                  T
!     This is the J  part,  translating child data

!      call apply_Jt(u,nx,ny,nz,nel)



!                 ~ ~T
!     This is the Q Q  part

    call gs_op(gsh_fld(ifldt),u,1,1,0)  ! 1 ==> +



!     This is the J  part,  interpolating parent solution onto child

!      call apply_J(u,nx,ny,nz,nel)


#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tdsum=tdsum+timee
    ndsum=icalld
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine dssum
!-----------------------------------------------------------------------
    subroutine dsop(u,op,nx,ny,nz)
    include  'CTIMER'
    use size_m
    include 'PARALLEL'
    include 'INPUT'
    include 'TSTEP'

    real :: u(1)
    character(3) :: op
    character(10) :: s1,s2

!     o gs recognized operations:

!             o "+" ==> addition.
!             o "*" ==> multiplication.
!             o "M" ==> maximum.
!             o "m" ==> minimum.
!             o "A" ==> (fabs(x)>fabs(y)) ? (x) : (y), ident=0.0.
!             o "a" ==> (fabs(x)<fabs(y)) ? (x) : (y), ident=MAX_DBL
!             o "e" ==> ((x)==0.0) ? (y) : (x),        ident=0.0.

!             o note: a binary function pointer flavor exists.


!     o gs level:

!             o level=0 ==> pure tree
!             o level>=num_nodes-1 ==> pure pairwise
!             o level = 1,...num_nodes-2 ==> mix tree/pairwise.


    ifldt = ifield
!     if (ifldt.eq.0)       ifldt = 1
    if (ifldt == ifldmhd) ifldt = 1

!     if (nid.eq.0)
!    $   write(6,*) istep,' dsop: ',op,ifield,ifldt,gsh_fld(ifldt)

    if(ifsync) call nekgsync()

    if (op == '+  ') call gs_op(gsh_fld(ifldt),u,1,1,0)
    if (op == 'sum') call gs_op(gsh_fld(ifldt),u,1,1,0)
    if (op == 'SUM') call gs_op(gsh_fld(ifldt),u,1,1,0)

    if (op == '*  ') call gs_op(gsh_fld(ifldt),u,1,2,0)
    if (op == 'mul') call gs_op(gsh_fld(ifldt),u,1,2,0)
    if (op == 'MUL') call gs_op(gsh_fld(ifldt),u,1,2,0)

    if (op == 'm  ') call gs_op(gsh_fld(ifldt),u,1,3,0)
    if (op == 'min') call gs_op(gsh_fld(ifldt),u,1,3,0)
    if (op == 'mna') call gs_op(gsh_fld(ifldt),u,1,3,0)
    if (op == 'MIN') call gs_op(gsh_fld(ifldt),u,1,3,0)
    if (op == 'MNA') call gs_op(gsh_fld(ifldt),u,1,3,0)

    if (op == 'M  ') call gs_op(gsh_fld(ifldt),u,1,4,0)
    if (op == 'max') call gs_op(gsh_fld(ifldt),u,1,4,0)
    if (op == 'mxa') call gs_op(gsh_fld(ifldt),u,1,4,0)
    if (op == 'MAX') call gs_op(gsh_fld(ifldt),u,1,4,0)
    if (op == 'MXA') call gs_op(gsh_fld(ifldt),u,1,4,0)

    return
    end subroutine dsop
!-----------------------------------------------------------------------
    subroutine vec_dssum(u,v,w,nx,ny,nz)

!     Direct stiffness summation of the face data, for field U.

!     Boundary condition data corresponds to component IFIELD of
!     the CBC array.

    use ctimer
    INCLUDE 'SIZE'
    INCLUDE 'TOPOL'
    INCLUDE 'INPUT'
    INCLUDE 'PARALLEL'
    INCLUDE 'TSTEP'

    REAL :: U(1),V(1),W(1)

    if(ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) tvdss=0.0d0
    if (icalld == 0) tgsum=0.0d0
    icalld=icalld+1
    nvdss=icalld
    etime1=dnekclock()
#endif


!============================================================================
!     execution phase
!============================================================================

    ifldt = ifield
!     if (ifldt.eq.0)       ifldt = 1
    if (ifldt == ifldmhd) ifldt = 1

    call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)

#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tvdss=tvdss+timee
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine vec_dssum

!-----------------------------------------------------------------------
    subroutine vec_dsop(u,v,w,nx,ny,nz,op)

!     Direct stiffness summation of the face data, for field U.

!     Boundary condition data corresponds to component IFIELD of
!     the CBC array.

    use ctimer
    INCLUDE 'SIZE'
    INCLUDE 'TOPOL'
    INCLUDE 'INPUT'
    INCLUDE 'PARALLEL'
    INCLUDE 'TSTEP'

    real :: u(1),v(1),w(1)
    character(3) :: op

!============================================================================
!     execution phase
!============================================================================

    ifldt = ifield
!     if (ifldt.eq.0)       ifldt = 1
    if (ifldt == ifldmhd) ifldt = 1

!     write(6,*) 'opdsop: ',op,ifldt,ifield
    if(ifsync) call nekgsync()

    if (op == '+  ' .OR. op == 'sum' .OR. op == 'SUM') &
    call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)


    if (op == '*  ' .OR. op == 'mul' .OR. op == 'MUL') &
    call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,2,0)


    if (op == 'm  ' .OR. op == 'min' .OR. op == 'mna' &
     .OR. op == 'MIN' .OR. op == 'MNA') &
    call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,3,0)


    if (op == 'M  ' .OR. op == 'max' .OR. op == 'mxa' &
     .OR. op == 'MAX' .OR. op == 'MXA') &
    call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,4,0)


    return
    end subroutine vec_dsop
!-----------------------------------------------------------------------
    subroutine nvec_dssum(u,stride,n,gs_handle)

!     Direct stiffness summation of the array u for n fields

    use ctimer
    use size_m

    real :: u(1)
    integer :: n,stride,gs_handle

    if(ifsync) call nekgsync()

#ifndef NOTIMER
    icalld=icalld+1
    nvdss=icalld
    etime1=dnekclock()
#endif
    call gs_op_fields(gs_handle,u,stride,n,1,1,0)

#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tvdss=tvdss+timee
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine nvec_dssum

!----------------------------------------------------------------------
    subroutine matvec3(uout,Jmat,uin,iftrsp,n1,n2)

    use size_m

    real :: Jmat (n1,n1,2)
    real :: uin   (1)
    real :: uout  (1)
    logical :: iftrsp

    common /matvtmp/ utmp(lx1,ly1)

    if (ndim == 2) then
        call mxm (Jmat(1,1,1),n1,uin,n1,uout,n2)
    else
        if (iftrsp) then
            call transpose(uout,n2,uin,n1)
        else
            call copy     (uout,uin,n1*n2)
        endif
        call mxm (Jmat(1,1,1),n1,uout,n1,utmp,n2)
        call mxm (utmp,n2,Jmat(1,1,2),n1,uout,n1)
    endif

    return
    end subroutine matvec3
!-----------------------------------------------------------------------
    subroutine matvec3t(uout,Jmat,uin,iftrsp,n1,n2)

    use size_m

    real :: Jmat (n1,n1,2)
    real :: uin   (n1,n2)
    real :: uout  (n1,n2)
    logical :: iftrsp

    common /matvtmp/ utmp(lx1*ly1)

    call transpose(utmp,n2,uin,n1)
    call mxm (Jmat(1,1,2),n1,utmp,n1,uout,n2)
    call mxm (uout,n2,Jmat(1,1,1),n1,utmp,n1)
    if (iftrsp) then
        call copy     (uout,utmp,n1*n2)
    else
        call transpose(uout,n2,utmp,n1)
    endif

    return
    end subroutine matvec3t
!-----------------------------------------------------------------------
    subroutine matvect (out,d,vec,n1,n2)
    dimension d(n1,n2),out(1),vec(1)

!   handle non-square matrix in mat-vec mult -- TRANSPOSE
!    N1 is still the number of rows
!    N2 is still the number of cols


    call mxm(vec,1,d,n1,out,n2)

    return
    end subroutine matvect
!-----------------------------------------------------------------------
!      subroutine opq_in_place(a,b,c)
!      use size_m
!      real a(1),b(1),c(1)

!      call q_in_place(a)
!      call q_in_place(b)
!      if (ndim .eq.3) call q_in_place(c)

!      return
!      end
!-----------------------------------------------------------------------
    subroutine vectof_add(b,a,ie,iface,nx,ny,nz)

!     Copy vector A to the face (IFACE) of B
!     IFACE is the input in the pre-processor ordering scheme.

    DIMENSION A(NX,NY)
    DIMENSION B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                B(IX,IY,IZ,IE) = B(IX,IY,IZ,IE) + A(k,1)
    100 END DO
    return
    end subroutine vectof_add
!-----------------------------------------------------------------------
    subroutine zero_f(b,ie,iface,nx,ny,nz)

!     ZERO the face (IFACE) of B
!     IFACE is the input in the pre-processor ordering scheme.

    DIMENSION B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)

    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                B(IX,IY,IZ,IE) = 0.
    100 END DO
    return
    end subroutine zero_f
!-----------------------------------------------------------------------
    subroutine ftovec_0(a,b,ie,iface,nx,ny,nz)

!     Copy the face (IFACE) of B to vector A.
!     IFACE is the input in the pre-processor ordering scheme.

    DIMENSION A(NX,NY)
    DIMENSION B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                A(k,1)=B(IX,IY,IZ,IE)
                B(IX,IY,IZ,IE)=0.0
    100 END DO
    return
    end subroutine ftovec_0
!-----------------------------------------------------------------------
    subroutine ftovec(a,b,ie,iface,nx,ny,nz)

!     Copy the face (IFACE) of B to vector A.
!     IFACE is the input in the pre-processor ordering scheme.

    real :: A(NX,NY)
    real :: B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                A(k,1)=B(IX,IY,IZ,IE)
    100 END DO
    return
    end subroutine ftovec
!-----------------------------------------------------------------------
    subroutine vectof(b,a,ie,iface,nx,ny,nz)

!     Copy vector A to the face (IFACE) of B
!     IFACE is the input in the pre-processor ordering scheme.

    real :: A(NX,NY)
    real :: B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                B(IX,IY,IZ,IE) = A(k,1)
    100 END DO
    return
    end subroutine vectof
!-----------------------------------------------------------------------
    subroutine ftoveci(a,b,ie,iface,nx,ny,nz)

!     Copy the face (IFACE) of B to vector A.
!     IFACE is the input in the pre-processor ordering scheme.

    integer :: A(NX,NY)
    integer :: B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                A(k,1)=B(IX,IY,IZ,IE)
    100 END DO
    return
    end subroutine ftoveci
!-----------------------------------------------------------------------
    subroutine vectofi(b,a,ie,iface,nx,ny,nz)

!     Copy vector A to the face (IFACE) of B
!     IFACE is the input in the pre-processor ordering scheme.

    integer :: A(NX,NY)
    integer :: B(NX,NY,NZ,1)
    CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
    k = 0
    DO 100 IZ=KZ1,KZ2
        DO 100 IY=KY1,KY2
            DO 100 IX=KX1,KX2
                k = k + 1
                B(IX,IY,IZ,IE) = A(k,1)
    100 END DO
    return
    end subroutine vectofi
!-----------------------------------------------------------------------
    subroutine apply_Jt(u,nx,ny,nz,nel)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)


!                  T
!     This is the J  part,  translating child data

    do ie = 1 , nel
    !        Note, we zero out u() on this face after extracting, for
    !        consistency reasons discovered during Jerry's thesis.
    !        Thus,  "ftovec_0" rather than ftovec().   (iface -- Ed notation)
        do iface = 1 , 2*ndim
            im = mortar(iface,ie)
            if (im /= 0) then
                call ftovec_0(uin(1,iface),u,ie,iface,nx,ny,nz)
            endif
        enddo
        do iface=1,2*ndim
            im = mortar(iface,ie)
            if (im /= 0) then
                if (if3d) then
                    call matvec3t &
                    (uout,Jmat(1,1,1,im),uin(1,iface),ifJt(im),nx,nx)
                else
                    call matvect (uout,Jmat(1,1,1,im),uin(1,iface),nx,nx)
                endif
                call vectof_add(u,uout,ie,iface,nx,ny,nz)
            endif
        enddo
    enddo

    return
    end subroutine apply_Jt
!-----------------------------------------------------------------------
    subroutine apply_J(u,nx,ny,nz,nel)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)

!     This is the J  part,  interpolating parent solution onto child


    do ie = 1 , nel
        do iface = 1 , 2*ndim
            im = mortar(iface,ie)
            if (im /= 0) then
                call ftovec(uin(1,iface),u,ie,iface,nx,ny,nz)
            endif
        enddo
        do iface=1,2*ndim
            im = mortar(iface,ie)
            if (im /= 0) then
                call matvec3 &
                (uout,Jmat(1,1,1,im),uin(1,iface),ifJt(im),nx,nz)
                call vectof (u,uout,ie,iface,nx,ny,nz)
            endif
        enddo
    enddo

    return
    end subroutine apply_J
!-----------------------------------------------------------------------
    subroutine h1_proj(u,nx,ny,nz)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)

    if(ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) then
        tdsmx=0.
        tdsmn=0.
    endif
    icalld=icalld+1
    etime1=dnekclock()
#endif

    ifldt = ifield
!     if (ifldt.eq.0) ifldt = 1
    nel = nelv
    if (ifield >= 2) nel=nelt
    ntot = nx*ny*nz*nel




!                        ~  ~T
!     Implement   :=   J Q  Q  Mu


!                  T

    call col2  (u,umult,ntot)

!                 ~ ~T
!     This is the Q Q  part

    call gs_op(gsh_fld(ifldt),u,1,1,0)


!     This is the J  part,  interpolating parent solution onto child

    call apply_J(u,nx,ny,nz,nel)


#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tdsum=tdsum+timee
    ndsum=icalld
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine h1_proj
!-----------------------------------------------------------------------
    subroutine dssum_msk(u,mask,nx,ny,nz)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1),mask(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)

    if(ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) then
        tdsmx=0.
        tdsmn=0.
    endif
    icalld=icalld+1
    etime1=dnekclock()
#endif

    ifldt = ifield
!     if (ifldt.eq.0) ifldt = 1
    nel = nelv
    if (ifield >= 2) nel=nelt
    ntot = nx*ny*nz*nel



!                    T           ~  ~T  T
!     Implement Q M Q   :=   J M Q  Q  J


!                  T
!     This is the J  part,  translating child data

    call apply_Jt(u,nx,ny,nz,nel)



!                 ~ ~T
!     This is the Q Q  part

    call gs_op(gsh_fld(ifldt),u,1,1,0)
    call col2  (u,mask,ntot)


!     This is the J  part,  interpolating parent solution onto child

    call apply_J(u,nx,ny,nz,nel)


#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tdsum=tdsum+timee
    ndsum=icalld
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine dssum_msk
!-----------------------------------------------------------------------
    subroutine dssum_msk2(u,mask,binv,nx,ny,nz)
    use ctimer
    use size_m
    include 'INPUT'
    include 'NONCON'
    include 'PARALLEL'
    include 'TSTEP'
    real :: u(1),mask(1),binv(1)

    parameter (lface=lx1*ly1)
    common /nonctmp/ uin(lface,2*ldim),uout(lface)

    if(ifsync) call nekgsync()

#ifndef NOTIMER
    if (icalld == 0) then
        tdsmx=0.
        tdsmn=0.
    endif
    icalld=icalld+1
    etime1=dnekclock()
#endif


    ifldt = ifield
!     if (ifldt.eq.0) ifldt = 1
    nel = nelv
    if (ifield >= 2) nel=nelt
    ntot = nx*ny*nz*nel




!                    T           ~  ~T  T
!     Implement Q M Q   :=   J M Q  Q  J


!                  T
!     This is the J  part,  translating child data

    call apply_Jt(u,nx,ny,nz,nel)



!                 ~ ~T
!     This is the Q Q  part

    call gs_op(gsh_fld(ifldt),u,1,1,0)
    call col3  (u,mask,binv,ntot)


!     This is the J  part,  interpolating parent solution onto child

    call apply_J(u,nx,ny,nz,nel)


#ifndef NOTIMER
    timee=(dnekclock()-etime1)
    tdsum=tdsum+timee
    ndsum=icalld
    tdsmx=max(timee,tdsmx)
    tdsmn=min(timee,tdsmn)
#endif

    return
    end subroutine dssum_msk2
!-----------------------------------------------------------------------
