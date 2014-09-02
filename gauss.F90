!> \brief the first subroutine to compute the matrix inverse
SUBROUTINE LU(A,N,NDIM,IR,IC)
  use kinds, only : DP
  implicit none
  real(DP) :: A(NDIM,1)
  integer :: n, ndim,IR(1),IC(1)

  integer :: i, j, k, l, m, irl, icm, k1
  real(DP) :: xmax, y, b, c

  DO 10 I=1,N
      IR(I)=I
      IC(I)=I
  10 END DO
  K=1
  L=K
  M=K
  XMAX=ABS(A(K,K))
  DO 100 I=K,N
      DO 100 J=K,N
          Y=ABS(A(I,J))
          IF(XMAX >= Y) GOTO 100
          XMAX=Y
          L=I
          M=J
  100 END DO
  DO 1000 K=1,N
      IRL=IR(L)
      IR(L)=IR(K)
      IR(K)=IRL
      ICM=IC(M)
      IC(M)=IC(K)
      IC(K)=ICM
      IF(L == K) GOTO 300
      DO 200 J=1,N
          B=A(K,J)
          A(K,J)=A(L,J)
          A(L,J)=B
      200 END DO
      300 IF(M == K) GOTO 500
      DO 400 I=1,N
          B=A(I,K)
          A(I,K)=A(I,M)
          A(I,M)=B
      400 END DO
      500 C=1./A(K,K)
      A(K,K)=C
      IF(K == N) GOTO 1000
      K1=K+1
      XMAX=ABS(A(K1,K1))
      L=K1
      M=K1
      DO 600 I=K1,N
          A(I,K)=C*A(I,K)
      600 END DO
      DO 800 I=K1,N
          B=A(I,K)
          DO 800 J=K1,N
              A(I,J)=A(I,J)-B*A(K,J)
              Y=ABS(A(I,J))
              IF(XMAX >= Y) GOTO 800
              XMAX=Y
              L=I
              M=J
      800 END DO
  1000 END DO
  RETURN
  END SUBROUTINE LU

!> \brief second part of the matrix inverse
SUBROUTINE SOLVE(F,A,K,N,NDIM,IR,IC)
  use kinds, only : DP
  implicit none

  integer :: n, k, ndim
  real(DP) :: A(NDIM,1),F(NDIM,1)
  integer :: IR(1),IC(1) 

  real(DP) :: G(2000)

  integer :: n1, kk, i, iri, i1, j, it, ici
  real(DP) :: B

  IF (N > 2000) THEN
      write(6,*) 'Abort IN Subrtouine SOLVE: N>2000, N=',N
      call exitt
  ENDIF

  N1=N+1
  DO 1000 KK=1,K
      DO 100 I=1,N
          IRI=IR(I)
          G(I)=F(IRI,KK)
      100 END DO
      DO 400 I=2,N
          I1=I-1
          B=G(I)
          DO 300 J=1,I1
              B=B-A(I,J)*G(J)
          300 END DO
          G(I)=B
      400 END DO
      DO 700 IT=1,N
          I=N1-IT
          I1=I+1
          B=G(I)
          IF(I == N) GOTO 701
          DO 600 J=I1,N
              B=B-A(I,J)*G(J)
          600 END DO
          701 G(I)=B*A(I,I)
      700 END DO
      DO 900 I=1,N
          ICI=IC(I)
          F(ICI,KK)=G(I)
      900 END DO
  1000 END DO
  RETURN
END SUBROUTINE SOLVE
