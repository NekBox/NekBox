!> \brief the first subroutine to compute the matrix inverse
SUBROUTINE LU(A,N,NDIM,IR,IC)
  use kinds, only : DP
  implicit none
  integer :: ndim, n, ir(1), ic(1)
  real(DP) :: A(NDIM,1)

  integer :: i, j, k, l, m, irl, icm, k1
  real(DP) :: xmax, y, b, c

  DO I=1,N
      IR(I)=I
      IC(I)=I
  END DO
  K=1
  L=K
  M=K
  XMAX=ABS(A(K,K))
  DO I=K,N
      DO J=K,N
          Y=ABS(A(I,J))
          IF(XMAX >= Y) cycle
          XMAX=Y
          L=I
          M=J
      END DO
  enddo 
  DO K=1,N
      IRL=IR(L)
      IR(L)=IR(K)
      IR(K)=IRL
      ICM=IC(M)
      IC(M)=IC(K)
      IC(K)=ICM
      IF(L /= K) then 
         DO J=1,N
             B=A(K,J)
             A(K,J)=A(L,J)
             A(L,J)=B
         END DO
      endif
      IF(M /= K) then 
         DO I=1,N
             B=A(I,K)
             A(I,K)=A(I,M)
             A(I,M)=B
         END DO
      endif
      C=1./A(K,K)
      A(K,K)=C
      IF(K == N) cycle
      K1=K+1
      XMAX=ABS(A(K1,K1))
      L=K1
      M=K1
      DO I=K1,N
          A(I,K)=C*A(I,K)
      END DO
      DO I=K1,N
          B=A(I,K)
          DO J=K1,N
              A(I,J)=A(I,J)-B*A(K,J)
              Y=ABS(A(I,J))
              IF(XMAX >= Y) cycle
              XMAX=Y
              L=I
              M=J
          END DO
       END DO 
  END DO 
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
