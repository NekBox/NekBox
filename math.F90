!-----------------------------------------------------------------------
subroutine tensor_product_transform(u, nu, v, nv, A, At, work1, work2)
  use kinds, only : DP
  implicit none
  integer, intent(in)   :: nu, nv
  real(DP), intent(in)  :: u(*)
  real(DP), intent(out) :: v(*)
  real(DP), intent(in)  :: A(*), At(*)
  real(DP), intent(out) :: work1(0:nu*nu*nv-1), work2(0:nu*nv*nv-1) ! scratch

  integer :: i
 
  call mxm(A,nv,u,nu,work1,nu*nu)
  do i=0,nu-1
      call mxm(work1(nv*nu*i),nv,At,nu,work2(nv*nv*i),nv)
  enddo
  call mxm(work2,nv*nv,At,nu,v,nv)
  return

end subroutine tensor_product_transform

!-----------------------------------------------------------------------
!> \brief  Tensor product application of v = (C x B x A) u .
!!  NOTE -- the transpose of B & C must be input, rather than B & C.
!!  -  scratch arrays: work1(nu*nu*nv), work2(nu*nv*nv)
subroutine tensor_product_multiply(u, nu, v, nv, A, Bt, Ct, work1, work2)
  use kinds, only : DP
  implicit none
  integer, intent(in)   :: nu, nv
  real(DP), intent(in)  :: u(*)
  real(DP), intent(out) :: v(*)
  real(DP), intent(in)  :: A(*), Bt(*), Ct(*)
  real(DP), intent(out) :: work1(0:nu*nu*nv-1), work2(0:nu*nv*nv-1) ! scratch

  integer :: i
 
  call mxm(A,nv,u,nu,work1,nu*nu)
  do i=0,nu-1
      call mxm(work1(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
  enddo
  call mxm(work2,nv*nv,Ct,nu,v,nv)
  return

end subroutine tensor_product_multiply



!-----------------------------------------------------------------------
!> \brief     Output: ur,us,ut         Input:u,N,e,D,Dt
subroutine local_grad3(ur,us,ut,u,N,e,D,Dt)
  use kinds, only : DP
  implicit none
  integer :: N,e
  real(DP) :: ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
  real(DP) :: u (0:N,0:N,0:N,*)
  real(DP) :: D (0:N,0:N),Dt(0:N,0:N)

  integer :: m1, m2, k
  m1 = N+1
  m2 = m1*m1

  call mxm(D ,m1,u(0,0,0,e),m1,ur,m2)
  do k=0,N
      call mxm(u(0,0,k,e),m1,Dt,m1,us(0,0,k),m1)
  enddo
  call mxm(u(0,0,0,e),m2,Dt,m1,ut,m1)

  return
end subroutine local_grad3

!-----------------------------------------------------------------------
subroutine helmholtz(h1, h2, nx, ny, nz, &
                     u, au, gx, gy, gz, b, &
                     work1, work2, work3)
  use kinds, only : DP
  use dxyz, only : wddx, wddyt, wddzt
  implicit none

  real(DP), intent(in) :: h1, h2
  integer, intent(in) :: nx, ny, nz
  real(DP), intent(in), dimension(nx, ny, nz) :: u, gx, gy, gz, b
  real(DP), intent(out), dimension(nx, ny, nz) :: au, work1, work2, work3

  integer :: iz

  call mxm   (wddx,nx,u(1,1,1),nx,work1,ny*nz)
  do iz=1,nz
      call mxm   (u(1,1,iz),nx,wddyt,ny,work2(1,1,iz),ny)
  END DO
  call mxm   (u(1,1,1),nx*ny,wddzt,nz,work3,nz)

  if (h2 /= 0._dp) then
    au(:,:,:) = h1* ( work1*gx + work2*gy + work3*gz ) + h2*b*u
  else
    au(:,:,:) = h1* ( work1*gx + work2*gy + work3*gz ) 
  endif

  return

end subroutine helmholtz

!----------------------------------------------------------------------
!> \brief clobbers r
subroutine hsmg_do_fast(e,r,s,d,nl)
  use kinds, only : DP
  use size_m, only : ndim, nelv
  use size_m, only : lx1, ly1, lz1
  use ctimer, only : schw_flop, schw_mop
  implicit none

  integer :: nl
  real(DP) :: e(nl**ndim,nelv)
  real(DP) :: r(nl**ndim,nelv)
  real(DP) :: s(nl*nl,2,ndim,nelv)
  real(DP) :: d(nl**ndim,nelv)
        
  integer :: ie,nn,i

  integer, parameter :: lwk=(lx1+2)*(ly1+2)*(lz1+2)
  real(DP) :: work1(nl*nl*nl),work2(nl*nl*nl)

  nn=nl**ndim
  schw_flop = schw_flop + nelv*nn
  schw_mop  = schw_mop  + nelv*nn ! r and e should be in cache
  schw_flop = schw_flop + 3*nn*(2*nl-1)*nelv
  schw_mop  = schw_mop + (nn + 3*nl*nl)*nelv
  schw_flop = schw_flop + 3*nn*(2*nl-1)*nelv
  schw_mop  = schw_mop + (nn + 3*nl*nl)*nelv

  do ie=1,nelv

    call tensor_product_multiply(r(1,ie), nl, r(1,ie), nl, s(1,2,1,ie), s(1,1,2,ie), s(1,1,3, ie), work1, work2)

    do i=1,nn
        r(i,ie)=d(i,ie)*r(i,ie)
    enddo

    call tensor_product_multiply(r(1,ie), nl, e(1,ie), nl, s(1,1,1,ie), s(1,2,2,ie), s(1,2,3, ie), work1, work2)

  enddo

  return
end subroutine hsmg_do_fast

!-----------------------------------------------------------------------
!> \brief Compute curl of U.
!!
!! \f$ (w_1, w_2, w_3) = \nabla \times (u_1, u_2, u_3) \f$
subroutine op_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lelv, nx1, ny1, nz1, nelv
  use geom, only : rxm1, rym1, rzm1, sxm1, sym1, szm1, txm1, tym1, tzm1
  use dxyz, only : dztm1, dytm1, dxm1
  use geom, only : jacm1, bm1, binvm1, jacmi
  use input, only : ifaxis, ifcyclic
  use tstep, only : ifield
  use mesh, only : if_ortho
  implicit none

  real(DP), intent(out) :: w1(lx1,ly1,lz1,lelv) !>!< 1st component of curl U
  real(DP), intent(out) :: w2(lx1,ly1,lz1,lelv) !>!< 2nd component of curl U
  real(DP), intent(out) :: w3(lx1,ly1,lz1,lelv) !>!< 3rd component of curl U
  real(DP), intent(in)  :: u1(lx1,ly1,lz1,lelv) !>!< 1st component of U
  real(DP), intent(in)  :: u2(lx1,ly1,lz1,lelv) !>!< 2nd component of U
  real(DP), intent(in)  :: u3(lx1,ly1,lz1,lelv) !>!< 3rd component of U
  real(DP), intent(out) :: work1(lx1,ly1,lz1,lelv) !>!< work array
  real(DP), intent(out) :: work2(lx1,ly1,lz1,lelv) !>!< work array
  logical, intent(in)   :: ifavg !>!< Average at boundary? 

  integer :: ntot, nxyz, ifielt, nxy1, nyz1, iel, iz
  real(DP), allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)

  allocate(tmp1(nx1,ny1,nz1), tmp2(nx1,ny1,nz1), tmp3(nx1, ny1, nz1))

  ntot  = nx1*ny1*nz1*nelv
  nxyz  = nx1*ny1*nz1
  NXY1  = NX1*NY1
  NYZ1  = NY1*NZ1

  if (if_ortho) then
    do iel = 1, nelv

      if (ifavg .and. .not. ifcyclic) then
        tmp3 = jacmi(:,:,:,iel) * bm1(:,:,:,iel)
      else
        tmp3 = jacmi(:,:,:,iel)
      endif

      ! work1=dw/dy ; work2=dv/dz
      do iz = 1, nz1
        CALL MXM  (U3(1,1,iz,iel),NX1,DYTM1,NY1,tmp1(1,1,iz),NY1)
      enddo
      CALL MXM  (U2(1,1,1,iel),NXY1,DZTM1,NZ1,tmp2,NZ1) 
      w1(:,:,:,iel) = (tmp1*sym1(:,:,:,iel) - tmp2*tzm1(:,:,:,iel)) * tmp3

      ! work1=du/dz ; work2=dw/dx
      CALL MXM  (U1(1,1,1,iel),NXY1,DZTM1,NZ1,tmp1,NZ1) 
      CALL MXM  (DXM1,NX1,U3(1,1,1,iel),NX1,tmp2,NYZ1)
      w2(:,:,:,iel) = (tmp1*tzm1(:,:,:,iel) - tmp2*rxm1(:,:,:,iel)) * tmp3 

      ! work1=dv/dx ; work2=du/dy
      CALL MXM   (DXM1,NX1,U2(1,1,1,iel),NX1,tmp1,NYZ1)
      do iz = 1, nz1
        CALL MXM  (U1(1,1,iz,iel),NX1,DYTM1,NY1,tmp2(1,1,iz),NY1)
      enddo
      w3(:,:,:,iel) = (tmp1*rxm1(:,:,:,iel) - tmp2*sym1(:,:,:,iel)) * tmp3 
    enddo
  else
    ! work1=dw/dy ; work2=dv/dz
    call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2)
    call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
    w1 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv) 

    ! work1=du/dz ; work2=dw/dx
    call dudxyz(work1,u1,rzm1,szm1,tzm1,jacm1,1,3)
    call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
    w2 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)

    ! work1=dv/dx ; work2=du/dy
    call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1)
    call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
    w3 = work1(:,:,:,1:nelv) - work2(:,:,:,1:nelv)
  
    if (ifavg .AND. .NOT. ifcyclic) then
      do iel = 1, nelv
        w1(:,:,:,iel) = w1(:,:,:,iel) * bm1(:,:,:,iel)
        w2(:,:,:,iel) = w2(:,:,:,iel) * bm1(:,:,:,iel)
        w3(:,:,:,iel) = w3(:,:,:,iel) * bm1(:,:,:,iel)
      enddo
    endif
  endif

  !  Avg at bndry
  if (ifavg .AND. .NOT. ifcyclic) then
    ifielt = ifield
    ifield = 1           
    call opdssum (w1,w2,w3)
    ifield = ifielt

    do iel = 1, nelv
      w1(:,:,:,iel) = w1(:,:,:,iel) * binvm1(:,:,:,iel)
      w2(:,:,:,iel) = w2(:,:,:,iel) * binvm1(:,:,:,iel)
      w3(:,:,:,iel) = w3(:,:,:,iel) * binvm1(:,:,:,iel)
    enddo
  endif

  return
end subroutine op_curl

!-----------------------------------------------------------------------
subroutine div_diag(alpha, beta, nx, ny, nz, prefactor, &
                    u, rx, v, sy, w, tz, res, work1, work2)
  use kinds, only : DP
  use dxyz, only : dxtm12, dym12, dzm12
  implicit none
  real(DP), intent(in) :: alpha, beta
  integer, intent(in) :: nx, ny, nz
  real(DP), intent(in), dimension(nx, ny, nz) :: prefactor, u, v, w, rx, sy, tz
  real(DP), intent(inout), dimension(nx, ny, nz) :: res
  real(DP) , intent(out), dimension(nx, ny, nz) :: work1, work2

  integer :: iz 

  ! X 
  work1 = u * rx * prefactor
  call mxm  (dxtm12,nx,work1,nx,work2,ny*nz)
  res = alpha * work2 + res * beta
  ! Y 
  work1 = v * sy * prefactor
  do iz=1,nz
      call mxm  (work1(:,:,iz),nx,dym12,ny,work2(:,:,iz),ny)
  enddo
  res = alpha * work2 + res
  ! Z
  work1 = w * tz * prefactor
  call mxm  (work1,nx*ny,dzm12,nz,work2,nz)
  res = alpha * work2 + res

  return
end subroutine div_diag

!-------------------------------------------------------------
!> \brief Compute DT*X (entire field)
!-------------------------------------------------------------
subroutine cdtp (dtx,x,rm2,sm2,tm2,isd)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, ly2, lz2, lelv
  use size_m, only : nx1, ny1, nz1, nx2, ny2, nz2, nelv, ndim
  use ctimer, only : icalld, tcdtp, ncdtp, etime1, dnekclock
  use dxyz, only : dym12, dam12, dcm12, dxtm12, dzm12
  use geom, only : ifrzer, jacm2, ym2, jacm1
  use input, only : ifaxis, ifsplit
  use ixyz, only : iym12, iam12, icm12
  use geom, only : bm1, bm2
  use wz_m, only : w3m2, w2am2, w2cm2
  implicit none

  integer :: isd
  real(DP) :: dtx  (lx1,ly1,lz1,lelv)
  real(DP) :: x    (lx2,ly2,lz2,lelv)
  real(DP) :: rm2  (lx2,ly2,lz2,lelv)
  real(DP) :: sm2  (lx2,ly2,lz2,lelv)
  real(DP) :: tm2  (lx2,ly2,lz2,lelv)

  real(DP) ::  wx  (lx1,ly1,lz1) &
  ,             ta1 (lx1,ly1,lz1) &
  ,             ta2 (lx1,ly1,lz1)

  integer :: e
  integer :: nxyz1, nxyz2, nxy1, nyz2, n1, n2, ny12, i1, i2, iz

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
  !      Collocate with weights
    wx = bm1(:,:,:,e) * x(:,:,:,e) / jacm1(:,:,:,e)

    ta1 = wx * rm2(:,:,:,e)
    call mxm  (dxtm12,nx1,ta1,nx2,dtx(:,:,:,e),nyz2)
    ta1 = wx * sm2(:,:,:,e)
    i1 = 1
    i2 = 1
    do iz=1,nz2
        call mxm  (ta1(:,:,iz),nx1,dym12,ny2,ta2(:,:,iz),ny1)
        i1 = i1 + n1
        i2 = i2 + n2
    enddo
    dtx(:,:,:,e) = dtx(:,:,:,e) + ta2
    ta1 = wx * tm2(:,:,:,e)
    call mxm  (ta1,nxy1,dzm12,nz2,ta2,nz1)
    dtx(:,:,:,e) = dtx(:,:,:,e) + ta2
  enddo

end subroutine cdtp 


!> \brief local inner product, with weight
real(DP) FUNCTION VLSC3(X,Y,B,N)
  use kinds, only : DP
  use opctr, only : isclld, nrout, myrout, rname, dct, ncall, dcount
  implicit none

  integer, intent(in) :: n
  real(DP), intent(in) :: X(n),Y(n),B(n)

  REAL(DP) :: DT, T
  integer :: isbcnt, i

#ifndef NOTIMER
  if (isclld == 0) then
      isclld=1
      nrout=nrout+1
      myrout=nrout
      rname(myrout) = 'VLSC3 '
  endif
  isbcnt = 3*n
  dct(myrout) = dct(myrout) + float(isbcnt)
  ncall(myrout) = ncall(myrout) + 1
  dcount      =      dcount + float(isbcnt)
#endif

  DT = 0.0
  DO 10 I=1,N
      T = X(I)*Y(I)*B(I)
      DT = DT+T
  10 END DO
  T=DT
  VLSC3 = T
  RETURN
END FUNCTION VLSC3

!-----------------------------------------------------------------------
!> \brief blank a string
SUBROUTINE BLANK(A,N)
  implicit none
  CHARACTER(1) :: A(*)
  integer :: n
  CHARACTER(1) :: BLNK = ' '
  integer :: i

  DO 10 I=1,N
      A(I)=BLNK
  10 END DO
  RETURN
END SUBROUTINE BLANK

!-----------------------------------------------------------------------
subroutine copy(a,b,n)
  use kinds, only : DP
  implicit none
  integer :: n
  real(DP) :: a(n),b(n)
  a = b

  return
end subroutine copy

!-----------------------------------------------------------------------
subroutine chcopy(a,b,n)
  implicit none
  integer :: n, i
  CHARACTER(1) :: A(n), B(n)

  DO 100 I = 1, N
      A(I) = B(I)
  100 END DO
  return
end subroutine chcopy

!-----------------------------------------------------------------------
!> \brief vector local max(abs( ))
real(DP) function vlamax(vec,n)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  REAL(DP), intent(in) :: VEC(n)
  integer :: i

  real(DP) :: TAMAX = 0.0

  DO I=1,N
      TAMAX = MAX(TAMAX,ABS(VEC(I)))
  END DO

  VLAMAX = TAMAX
  return
end function vlamax

!-----------------------------------------------------------------------
!> \brief Compute a Cartesian vector cross product.
subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)
  use kinds, only : DP
  implicit none

  integer,  intent(in)  :: n
  real(DP), intent(out) :: U1(n), U2(n), U3(n)
  real(DP), intent(in)  :: V1(n), V2(n), V3(n)
  real(DP), intent(in)  :: W1(n), W2(n), W3(n)
  integer :: i

  DO I=1,N
      U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
      U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
      U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  END DO

  return
end subroutine vcross

!-----------------------------------------------------------------------
!> \brief Yields MOD(I,N) with the exception that if I=K*N, result is N.
integer function mod1(i,n)
  implicit none
  integer, intent(in) :: i, n
  integer :: ii
  MOD1=0
  IF (I == 0) THEN
      return
  ENDIF
  IF (N == 0) THEN
      WRITE(6,*) &
      'WARNING:  Attempt to take MOD(I,0) in function mod1.'
      return
  ENDIF
  II = I+N-1
  MOD1 = MOD(II,N)+1
  return
end function mod1

!-----------------------------------------------------------------------
integer function log2(k)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: k
  real(DP) :: rk, rlog, rlog2

  RK=(K)
  RLOG=LOG10(RK)
  RLOG2=LOG10(2.0)
  RLOG=RLOG/RLOG2+0.5
  LOG2=INT(RLOG)
  return
end function log2

!-----------------------------------------------------------------------
!> \brief SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ) into item(i)
!! where JJ = ind(i)
subroutine iswap(b,ind,n,temp)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: ind(n)
  integer, intent(inout) :: b(n)
  integer, intent(out) :: temp(n) ! scratch
  integer :: i, jj

  DO I=1,N
      JJ=IND(I)
      TEMP(I)=B(JJ)
  END DO
  DO I=1,N
      B(I)=TEMP(I)
  END DO
  return
end subroutine iswap

!-----------------------------------------------------------------------
!     Vector reduction routines which require communication
!     on a parallel machine. These routines must be substituted with
!     appropriate routines which take into account the specific architecture.

!----------------------------------------------------------------------------
!> \brief Perform inner-product in double precision
real(DP) function glsc3(a,b,mult,n)
  use kinds, only : DP
  use opctr
  implicit none
  integer, intent(in) :: n
  REAL(DP), intent(in) :: A(n),B(n),MULT(n)
  REAL(DP) :: TMP,WORK(1)
  integer :: i, isbcnt

#ifndef NOTIMER
  if (isclld == 0) then
      isclld=1
      nrout=nrout+1
      myrout=nrout
      rname(myrout) = 'glsc3 '
  endif
  isbcnt = 3*n
  dct(myrout) = dct(myrout) + (isbcnt)
  ncall(myrout) = ncall(myrout) + 1
  dcount      =      dcount + (isbcnt)
#endif

  TMP = 0.0
  DO I=1,N
      TMP = TMP + A(I)*B(I)*MULT(I)
  END DO
  CALL GOP(TMP,WORK,'+  ',1)
  GLSC3 = TMP
  return
end function glsc3

!-----------------------------------------------------------------------
!> \brief Perform inner-product in double precision
real(DP) function glsc2(x,y,n)
  use kinds, only : DP
  use opctr
  integer, intent(in) :: n 
  real(DP), intent(in) :: x(n), y(n)
  real(DP) :: tmp,work(1)
  integer :: i, isbcnt

#ifndef NOTIMER
    if (isclld == 0) then
        isclld=1
        nrout=nrout+1
        myrout=nrout
        rname(myrout) = 'glsc2 '
    endif
    isbcnt = 2*n
    dct(myrout) = dct(myrout) + (isbcnt)
    ncall(myrout) = ncall(myrout) + 1
    dcount      =      dcount + (isbcnt)
#endif

  tmp=0.0
  do i=1,n
      tmp = tmp+ x(i)*y(i)
  END DO
  CALL GOP(TMP,WORK,'+  ',1)
  GLSC2 = TMP
  return
end function glsc2

!-----------------------------------------------------------------------
!> \brief Perform inner-product  x*x*y*z
real(DP) function glsc23(x,y,z,n)
  use kinds, only : DP
  implicit none
  integer :: n
  real(DP), intent(in) :: x(n), y(n),z(n)
  real(DP) :: tmp,work(1), ds
  integer :: i

  ds = 0.0
  do i=1,n
      ds=ds+x(i)*x(i)*y(i)*z(i)
  END DO
  tmp=ds
  call gop(tmp,work,'+  ',1)
  glsc23 = tmp
  return
end function glsc23

!-----------------------------------------------------------------------
real(DP) function glsum (x,n)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  real(DP), intent(in) :: X(n)
  real(DP) :: TMP(1),WORK(1), tsum
  integer :: i
  TSUM = 0._dp
  DO I=1,N
      TSUM = TSUM+X(I)
  END DO
  TMP(1)=TSUM
  CALL GOP(TMP,WORK,'+  ',1)
  GLSUM = TMP(1)
  return
end function glsum

!-----------------------------------------------------------------------
real(DP) function glamax(a,n)
  use kinds, only : DP
  implicit none
  integer :: n
  REAL(DP) :: A(n)
  real(DP) :: TMP(1),WORK(1), tmax
  integer :: i
  TMAX = 0.0
  DO I=1,N
      TMAX = MAX(TMAX,ABS(A(I)))
  END DO
  TMP(1)=TMAX
  CALL GOP(TMP,WORK,'M  ',1)
  GLAMAX=ABS(TMP(1))
  return
end function glamax

!-----------------------------------------------------------------------
integer function iglmin(a,n)
  implicit none
  integer, intent(in) :: n, a(n)
  integer :: tmp(1),work(1), tmin, i
  tmin=  999999999
  do i=1,n
      tmin=min(tmin,a(i))
  enddo
  tmp(1)=tmin
  call igop(tmp,work,'m  ',1)
  iglmin=tmp(1)
  return
end function iglmin

!-----------------------------------------------------------------------
integer function iglmax(a,n)
  implicit none
  integer, intent(in) :: n, a(n)
  integer :: tmp(1),work(1), tmax, i
  tmax= -999999999
  do i=1,n
      tmax=max(tmax,a(i))
  enddo
  tmp(1)=tmax
  call igop(tmp,work,'M  ',1)
  iglmax=tmp(1)
  return
end function iglmax

!-----------------------------------------------------------------------
integer function iglsum(a,n)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: a(n)
  integer :: tmp(1),work(1),tsum, i
  tsum= 0
  do i=1,n
      tsum=tsum+a(i)
  enddo
  tmp(1)=tsum
  call igop(tmp,work,'+  ',1)
  iglsum=tmp(1)
  return
end function iglsum

!-----------------------------------------------------------------------
!> \brief global sum (long integer)
integer(i8) function i8glsum(a,n)
  use kinds, only : i8
  integer,     intent(in) :: n
  integer(i8), intent(in) :: a(n)

  integer(i8) :: tsum, tmp(1),work(1)
  integer :: i

  tsum= 0
  do i=1,n
      tsum=tsum+a(i)
  enddo
  tmp(1)=tsum
  call i8gop(tmp,work,'+  ',1)
  i8glsum=tmp(1)
  return
END function

!-----------------------------------------------------------------------
real(DP) function glmax(a,n)
  use kinds, only : DP
  implicit none

  integer, intent(in) :: n
  REAL(DP), intent(in) :: A(n)
  real(DP) :: TMP(1),WORK(1), tmax
  integer :: i
 
  TMAX=-99.0e20
  DO I=1,N
      TMAX=MAX(TMAX,A(I))
  END DO
  TMP(1)=TMAX
  CALL GOP(TMP,WORK,'M  ',1)
  GLMAX=TMP(1)
  return
end function glmax

!-----------------------------------------------------------------------
real(DP) function glmin(a,n)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  REAL(DP), intent(in) :: A(n)
  real(DP) :: TMP(1),WORK(1), tmin
  integer :: i
  TMIN=99.0e20
  DO I=1,N
      TMIN=MIN(TMIN,A(I))
  END DO
  TMP(1)=TMIN
  CALL GOP(TMP,WORK,'m  ',1)
  GLMIN = TMP(1)
  return
end function glmin

!-----------------------------------------------------------------------
!> \brief If ANY LA=LB, then ALL LA=LB.
subroutine gllog(la,lb)
  use kinds, only : DP
  implicit none
  LOGICAL :: LA,LB
  real(DP) :: TMP(1),WORK(1)

  TMP(1)=1._dp
  IF (LB) THEN
      IF (LA) TMP(1)=0._dp
  ELSE
      IF ( .NOT. LA) TMP(1)=0._dp
  ENDIF
  CALL GOP(TMP,WORK,'*  ',1)
  IF (TMP(1) == 0._dp) LA=LB
  return
end subroutine gllog

!-----------------------------------------------------------------------
!> \brief   Use Heap Sort (p 231 Num. Rec., 1st Ed.)
subroutine isort(a,ind,n)
  implicit none
  integer :: n
  integer :: a(n),ind(n)
  integer :: aa, j, i, ii, ir, l

  dO 10 j=1,n
      ind(j)=j
  10 END DO

  if (n <= 1) return
  L=n/2+1
  ir=n
  100 continue
  if (l > 1) then
      l=l-1
      aa  = a  (l)
      ii  = ind(l)
  else
      aa =   a(ir)
      ii = ind(ir)
      a(ir) =   a( 1)
      ind(ir) = ind( 1)
      ir=ir-1
      if (ir == 1) then
          a(1) = aa
          ind(1) = ii
          return
      endif
  endif
  i=l
  j=l+l
  200 continue
  if (j <= ir) then
      if (j < ir) then
          if ( a(j) < a(j+1) ) j=j+1
      endif
      if (aa < a(j)) then
          a(i) = a(j)
          ind(i) = ind(j)
          i=j
          j=j+j
      else
          j=ir+1
      endif
      GOTO 200
  endif
  a(i) = aa
  ind(i) = ii
  GOTO 100

end subroutine isort

!-------------------------------------------------------
!> \brief Use Heap Sort (p 231 Num. Rec., 1st Ed.)
subroutine sort(a,ind,n)
  use kinds, only : DP
  implicit none
  integer :: n
  real(DP) :: a(n),aa
  integer :: ind(n)
  integer :: j, i, ii, ir, l

  dO 10 j=1,n
      ind(j)=j
  10 END DO

  if (n <= 1) return
  L=n/2+1
  ir=n
  100 continue
  if (l > 1) then
      l=l-1
      aa  = a  (l)
      ii  = ind(l)
  else
      aa =   a(ir)
      ii = ind(ir)
      a(ir) =   a( 1)
      ind(ir) = ind( 1)
      ir=ir-1
      if (ir == 1) then
          a(1) = aa
          ind(1) = ii
          return
      endif
  endif
  i=l
  j=l+l
  200 continue
  if (j <= ir) then
      if (j < ir) then
          if ( a(j) < a(j+1) ) j=j+1
      endif
      if (aa < a(j)) then
          a(i) = a(j)
          ind(i) = ind(j)
          i=j
          j=j+j
      else
          j=ir+1
      endif
      GOTO 200
  endif
  a(i) = aa
  ind(i) = ii
  GOTO 100
  end subroutine sort

!-----------------------------------------------------------------------
!> \brief Global maximum of long integer array
integer(i8) function i8glmax(a,n)
  use kinds, only : i8
  integer, intent(in) :: n
  integer(i8), intent(inout) :: a(n)
  integer(i8) :: tmp(1),work(1),tmax
  integer :: i

  tmax= -999999
  do i=1,n
      tmax=max(tmax,a(i))
  enddo

  tmp(1)=tmax
  call i8gop(tmp,work,'M  ',1)

  i8glmax=tmp(1)
  if (i8glmax == -999999) i8glmax=0

  return
END function

!-----------------------------------------------------------------------
!> \brief Construct A = I_n (identity matrix)
subroutine ident(a,n)
  use kinds, only : DP
  implicit none
  integer, intent(in) :: n
  real(DP), intent(out) ::  a(n,n)
  integer :: i

  a = 0._dp
  do i=1,n
      a(i,i) = 1._dp
  enddo
  return
end subroutine ident

!-----------------------------------------------------------------------
subroutine iswapt_ip(x,p,n)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: x(n)
  integer, intent(inout) :: p(n)

  integer :: j, k, loop_start, next, nextp, t1, t2

!   In-place permutation: x'(p) = x

  do k=1,n
      if (p(k) > 0) then   ! not swapped
          loop_start = k
          next       = p(loop_start)
          t1         = x(loop_start)
          do j=1,n
              if (next < 0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt
              elseif (next == loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
              else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
              endif
          enddo
          10 continue
      endif
  enddo

  do k=1,n
      p(k) = -p(k)
  enddo
  return
end subroutine iswapt_ip

