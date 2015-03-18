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
