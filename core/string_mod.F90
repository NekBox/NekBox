module string
  implicit none

!  integer, parameter :: * = 132

contains

!-----------------------------------------------------------------------
integer function i1_from_char(s1)
  implicit none
  character(*) :: s1

  character(10) :: n10 = '0123456789'

  i1_from_char = indx2(n10,10,s1,1)-1

  return
end function i1_from_char

!-----------------------------------------------------------------------
integer function indx2(s1,l1,s2,l2)
  implicit none
  character(*) :: s1,s2
  integer :: l1, l2
  integer :: n1, i, i2

  n1=l1-l2+1
  indx2=0
  if (n1 < 1) return

  do i=1,n1
      i2=i+l2-1
      if (s1(i:i2) == s2(1:l2)) then
          indx2=i
          return
      endif
  enddo

  return
end function indx2

!-----------------------------------------------------------------------
INTEGER FUNCTION INDX1(S1,S2,L2)
  implicit none
  CHARACTER(*) :: S1,S2
  integer :: l2
  integer :: n1, i, i2

  N1=len_trim(s1)-L2+1
  INDX1=0
  IF (N1 < 1) return

  DO 100 I=1,N1
      I2=I+L2-1
      IF (S1(I:I2) == S2(1:L2)) THEN
          INDX1=I
          return
      ENDIF
  100 END DO

  return
END FUNCTION INDX1

!-----------------------------------------------------------------------
!> \brief INDX_CUT is returned with the location of S2 in S1 (0 if not found)
!!   S1 is returned with 1st occurance of S2 removed.
INTEGER FUNCTION INDX_CUT(S1,S2,L2)
  implicit none
  CHARACTER(*) :: S1,S2
  integer :: l2

  integer :: i1, n1, i, i2, n2
  I1=INDX1(S1,S2,L2)

  IF (I1 /= 0) THEN
  
      N1=len_trim(s1)-L2
      DO 100 I=I1,N1
          I2=I+L2
      !           remove the 1st occurance of S2 from S1.
          S1(I:I)=S1(I2:I2)
      100 END DO
      N2=N1+1
      DO 200 I=N2,len_trim(s1)
          S1(I:I)=' '
      200 END DO
  ENDIF

  INDX_CUT=I1
  return
END FUNCTION INDX_CUT

!-----------------------------------------------------------------------
!> \brief split string S1 into two parts, delimited by S2.
subroutine csplit(s0,s1,s2,l0)
  implicit none
  CHARACTER(*) :: S0,S1,S2
  integer :: l0

  integer :: i2, i1

  I2=INDX_CUT(S1,S2,L0)
  IF (I2 == 0) return

  I1=I2-1
  CALL BLANK(S0,len(s0))
  S0(1:I1)=S1(1:I1)
  CALL LSHFT(S1,I2)

  return
end subroutine csplit
!-----------------------------------------------------------------------
!> \brief shift string from IPT to the left
!!  INPUT : "abcde......    test    "
!!  OUTPUT: "e......    test        "     if ipt.eq.5
subroutine lshft(string,ipt)
  implicit none
  CHARACTER(1) :: STRING(*)
  integer :: ipt
  integer :: j, ij

  DO 20 J=1,133-IPT
      IJ=IPT+J-1
      STRING(J)=STRING(IJ)
  20 END DO
  DO 30 J=134-IPT,len(string)
      STRING(J)=' '
  30 END DO
  return
end subroutine lshft

!-----------------------------------------------------------------------
!> \brief left justify string
subroutine ljust(string)
  implicit none
  CHARACTER(1) :: STRING(*)
  integer :: i, j, ij

  IF (STRING(1) /= ' ') return

  DO 100 I=2,132
  
      IF (STRING(I) /= ' ') THEN
          DO 20 J=1,133-I
              IJ=I+J-1
              STRING(J)=STRING(IJ)
          20 END DO
          DO 30 J=134-I,132
              STRING(J)=' '
          30 END DO
          return
      ENDIF
  
  100 END DO

  return
end subroutine ljust

!-----------------------------------------------------------------------
!> \brief Capitalizes string of length n
subroutine capit(lettrs,n)
  implicit none
  integer :: n
  CHARACTER LETTRS(N)
  integer :: i, int

  DO 5 I=1,N
      INT=ICHAR(LETTRS(I))
      IF(INT >= 97 .AND. INT <= 122) THEN
          INT=INT-32
          LETTRS(I)=CHAR(INT)
      ENDIF
  5 END DO
  return
end subroutine capit

!-----------------------------------------------------------------------
!> \brief Read VALUE from LINE and set IFGTRL to .TRUE. if successful,
!!                                IFGTRL to .FALSE. otherwise.
!!   This complicated function is necessary thanks to the Ardent,
!!   which won't allow free formatted reads (*) from internal strings!
LOGICAL FUNCTION IFGTRL(VALUE,LINE)
  use kinds, only : DP
  implicit none
  real(DP) :: value
  CHARACTER(*) :: LINE

  CHARACTER(len(line)) :: WORK
  CHARACTER(8) ::  FMAT
  integer :: ifldw
  real(DP) :: TVAL

!   Note that the format Fn.0 is appropriate for fields of type:
!        34   34.0  34.0e+00
!   The only difficulty would be with '34' but since we identify
!   the field width exactly, there is no problem.

  IFGTRL= .FALSE. 
  VALUE=0.0

  WORK=LINE
  CALL LJUST(WORK)
  IFLDW=INDX1(WORK,' ',1)-1

  IF (IFLDW > 0) THEN
      WRITE(FMAT,10) IFLDW
      10 FORMAT('(F',I3.3,'.0)')
      READ(WORK,FMAT,ERR=100,END=100) TVAL
      VALUE=TVAL
      IFGTRL= .TRUE. 
      return
  ENDIF

  100 CONTINUE
  return
END FUNCTION IFGTRL

!-----------------------------------------------------------------------
!> \brief Read IVALUE from LINE and set IFGTIL to .TRUE. if successful,
!!                                IFGTIL to .FALSE. otherwise.
!!  This complicated function is necessary thanks to the Ardent,
!!  which won't allow free formatted reads (*) from internal strings!
LOGICAL FUNCTION IFGTIL(IVALUE,LINE)
  use kinds, only : DP
  implicit none
  integer :: ivalue
  CHARACTER(*) :: LINE
  CHARACTER(len(line)) :: WORK
  CHARACTER(8) ::  FMAT

  integer :: ifldw
  real(DP) :: tval

  IFGTIL= .FALSE. 
  IVALUE=0

  WORK=LINE
  CALL LJUST(WORK)
  IFLDW=INDX1(WORK,' ',1)-1
  write(*,*) "line:", line, " work:", work, " ifldw:", ifldw

  IF (IFLDW > 0) THEN
      WRITE(FMAT,10) IFLDW
      10 FORMAT('(F',I3.3,'.0)')
      READ(WORK,FMAT,ERR=100,END=100) TVAL
      IVALUE=INT(TVAL)
      IFGTIL= .TRUE. 
      return
  ENDIF

  100 CONTINUE
  return
END FUNCTION IFGTIL

integer function ltrunc(string,l)
  implicit none
  integer :: l
  CHARACTER(1) :: STRING(L)
  integer :: i, l1
  CHARACTER(1) :: BLNK=' '

  DO 100 I=L,1,-1
      L1=I
      IF (STRING(I) /= BLNK) GOTO 200
  100 END DO
  L1=0
  200 CONTINUE
  LTRUNC=L1
  return
end function ltrunc

!----------------------------------------------------------------------- 
subroutine cscan(sout,key,nk) 
  implicit none 
  character(*) :: sout,key 
  character(len(sout)) :: tmp_string 
  integer :: i, nk 
 
  do i=1,100000000 
      call blank(tmp_string,len(sout)) 
      read (nk,80,end=100,err=100) tmp_string 
      call chcopy(sout, tmp_string,len(sout))
  !        write (6,*) tmp_string
      if (indx1(tmp_string,key,nk) /= 0) return
  enddo
  100 continue

  80 format(a132)
  return

end subroutine cscan
!-----------------------------------------------------------------------


end module string
