!-------------------------------------------------------------------
!> \brief Set up deformed element logical switches
!-------------------------------------------------------------------
subroutine setdef()
  use kinds, only : DP
  use size_m, only : nid, nelt, lelt
  use input, only : param, ifmvbd, if3d, ccurve, xc, yc, zc
  use mesh, only : ifdfrm
  implicit none

  real(DP) :: XCC(8),YCC(8),ZCC(8)
  integer :: INDX(8)
  REAL :: VEC(3,12)
  LOGICAL :: IFVCHK

  integer :: ie, iedg, i, i1, j, i2
  real(DP) :: veclen

! Corner notation:

!                4+-----+3    ^ Y
!                /     /|     |
!               /     / |     |
!             8+-----+7 +2    +----> X
!              |     | /     /
!              |     |/     /
!             5+-----+6    Z


  DO IE=1,NELT
      IFDFRM(IE)= .FALSE. 
  END DO

  IF (IFMVBD) return

!   Force IFDFRM=.true. for all elements (for timing purposes only)

  IF (param(59) /= 0 .AND. nid == 0) &
  write(6,*) 'NOTE: All elements deformed , param(59) ^=0'
  IF (param(59) /= 0) return

!   Check against cases which won't allow for savings in HMHOLTZ

  INDX(1)=1
  INDX(2)=2
  INDX(3)=4
  INDX(4)=3
  INDX(5)=5
  INDX(6)=6
  INDX(7)=8
  INDX(8)=7

!   Check for deformation (rotation is acceptable).

  DO 500 IE=1,NELT
  
      call rzero(vec,36)
      IF (IF3D) THEN
          DO 100 IEDG=1,8
              IF(CCURVE(IEDG,IE) /= ' ') THEN
                  IFDFRM(IE)= .TRUE. 
                  GOTO 500
              ENDIF
          100 END DO
      
          DO 105 I=1,8
              XCC(I)=XC(INDX(I),IE)
              YCC(I)=YC(INDX(I),IE)
              ZCC(I)=ZC(INDX(I),IE)
          105 END DO
      
          DO 110 I=1,4
              VEC(1,I)=XCC(2*I)-XCC(2*I-1)
              VEC(2,I)=YCC(2*I)-YCC(2*I-1)
              VEC(3,I)=ZCC(2*I)-ZCC(2*I-1)
          110 END DO
      
          I1=4
          DO I=0,1
              DO J=0,1
                  I1=I1+1
                  I2=4*I+J+3
                  VEC(1,I1)=XCC(I2)-XCC(I2-2)
                  VEC(2,I1)=YCC(I2)-YCC(I2-2)
                  VEC(3,I1)=ZCC(I2)-ZCC(I2-2)
              enddo
          enddo
      
          I1=8
          DO 130 I=5,8
              I1=I1+1
              VEC(1,I1)=XCC(I)-XCC(I-4)
              VEC(2,I1)=YCC(I)-YCC(I-4)
              VEC(3,I1)=ZCC(I)-ZCC(I-4)
          130 END DO
      
          DO 140 I=1,12
              VECLEN = VEC(1,I)**2 + VEC(2,I)**2 + VEC(3,I)**2
              VECLEN = SQRT(VECLEN)
              VEC(1,I)=VEC(1,I)/VECLEN
              VEC(2,I)=VEC(2,I)/VECLEN
              VEC(3,I)=VEC(3,I)/VECLEN
          140 END DO
      
      !        Check the dot product of the adjacent edges to see that it is zero.
      
          IFDFRM(IE)= .FALSE. 
          IF (  IFVCHK(VEC,1,5, 9)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,1,6,10)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,2,5,11)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,2,6,12)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,3,7, 9)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,3,8,10)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,4,7,11)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,4,8,12)  ) IFDFRM(IE)= .TRUE. 
      
      !      Check the 2D case....
      
      ELSE
      
          DO 200 IEDG=1,4
              IF(CCURVE(IEDG,IE) /= ' ') THEN
                  IFDFRM(IE)= .TRUE. 
                  GOTO 500
              ENDIF
          200 END DO
      
          DO 205 I=1,4
              XCC(I)=XC(INDX(I),IE)
              YCC(I)=YC(INDX(I),IE)
          205 END DO
      
          VEC(1,1)=XCC(2)-XCC(1)
          VEC(1,2)=XCC(4)-XCC(3)
          VEC(1,3)=XCC(3)-XCC(1)
          VEC(1,4)=XCC(4)-XCC(2)
          VEC(1,5)=0.0
          VEC(2,1)=YCC(2)-YCC(1)
          VEC(2,2)=YCC(4)-YCC(3)
          VEC(2,3)=YCC(3)-YCC(1)
          VEC(2,4)=YCC(4)-YCC(2)
          VEC(2,5)=0.0
      
          DO 220 I=1,4
              VECLEN = VEC(1,I)**2 + VEC(2,I)**2
              VECLEN = SQRT(VECLEN)
              VEC(1,I)=VEC(1,I)/VECLEN
              VEC(2,I)=VEC(2,I)/VECLEN
          220 END DO
      
      !        Check the dot product of the adjacent edges to see that it is zero.
      
          IFDFRM(IE)= .FALSE. 
          IF (  IFVCHK(VEC,1,3,5)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,1,4,5)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,2,3,5)  ) IFDFRM(IE)= .TRUE. 
          IF (  IFVCHK(VEC,2,4,5)  ) IFDFRM(IE)= .TRUE. 
      ENDIF
  500 END DO
  return
end subroutine setdef

!> \brief Take the dot product of the three components of VEC to see if it's zero.
LOGICAL FUNCTION IFVCHK(VEC,I1,I2,I3)
  use kinds, only : DP
  implicit none

  real(DP) :: vec(3,12)
  integer :: i1, i2, i3
  LOGICAL :: IFTMP

  real(DP) :: epsm, dot1, dot2, dot3, dot

  IFTMP= .FALSE. 
  EPSM=1.0E-06

  DOT1=VEC(1,I1)*VEC(1,I2)+VEC(2,I1)*VEC(2,I2)+VEC(3,I1)*VEC(3,I2)
  DOT2=VEC(1,I2)*VEC(1,I3)+VEC(2,I2)*VEC(2,I3)+VEC(3,I2)*VEC(3,I3)
  DOT3=VEC(1,I1)*VEC(1,I3)+VEC(2,I1)*VEC(2,I3)+VEC(3,I1)*VEC(3,I3)

  DOT1=ABS(DOT1)
  DOT2=ABS(DOT2)
  DOT3=ABS(DOT3)
  DOT=DOT1+DOT2+DOT3
  IF (DOT > EPSM) IFTMP= .TRUE. 

  IFVCHK=IFTMP
  return
END FUNCTION IFVCHK

!-----------------------------------------------------------------------
!> \brief Generate xyz coordinates  for all elements.
!! Velocity formulation : mesh 3 is used
!! Stress   formulation : mesh 1 is used
!-----------------------------------------------------------------------
subroutine gencoor ()!xm3,ym3,zm3)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, lx3, ly3, lz3
  use geom, only : ifgmsh3, xm1, ym1, zm1
  implicit none

!max  real(DP) :: XM3(LX3,LY3,LZ3,1),YM3(LX3,LY3,LZ3,1),ZM3(LX3,LY3,LZ3,1)

!   Select appropriate mesh
  IF ( IFGMSH3 ) THEN
    write(*,*) "Oops: IFGMSH3"
!max        CALL GENXYZ (XM3,YM3,ZM3,NX3,NY3,NZ3)
  ELSE
      CALL GENXYZ (XM1,YM1,ZM1,NX1,NY1,NZ1)
  ENDIF

  return
end subroutine gencoor

!-----------------------------------------------------------------------
subroutine genxyz (xml,yml,zml,nxl,nyl,nzl)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, nelt, ndim
  use input, only : ccurve, ifaxis
  implicit none

  integer :: nxl, nyl, nzl
  real(DP) :: xml(nxl,nyl,nzl,1),yml(nxl,nyl,nzl,1),zml(nxl,nyl,nzl,1)

  real(DP) :: h(lx1,3,2), zgml(lx1,3)

  integer, parameter :: ldw=2*lx1*ly1*lz1

  character(1) :: ccv
  integer :: ie, nfaces, iface, isid


#ifdef MOAB
! already read/initialized vertex positions
  if (ifmoab) return
#endif

!   Initialize geometry arrays with bi- triquadratic deformations
  call linquad(xml,yml,zml,nxl,nyl,nzl)

  do ie=1,nelt

      call setzgml (zgml,ie,nxl,nyl,nzl,ifaxis)
      call sethmat (h,zgml,nxl,nyl,nzl)

  !        Deform surfaces - general 3D deformations
  !                        - extruded geometry deformations
      nfaces = 2*ndim
      do iface=1,nfaces
          ccv = ccurve(iface,ie)
          if (ccv == 's') then
            write(*,*) "Oops: ccv = 's'"
!max            call sphsrf(xml,yml,zml,iface,ie,nxl,nyl,nzl,work)
          endif
          if (ccv == 'e') then
            write(*,*) "Oops: ccv = 'e'" 
!max            call gensrf(xml,yml,zml,iface,ie,nxl,nyl,nzl,zgml)
          endif
      enddo

      do isid=1,8
          ccv = ccurve(isid,ie)
          if (ccv == 'C') then
            write(*,*) "Oops: ccv = 'C'"
!max          call arcsrf(xml,yml,zml,nxl,nyl,nzl,ie,isid)
          endif
      enddo

  enddo

  return
end subroutine genxyz

!-----------------------------------------------------------------------
subroutine sethmat(h,zgml,nxl,nyl,nzl)
  use kinds, only : DP
  use size_m, only : lx1
  use input, only : if3d
  implicit none

  real(DP) :: h(lx1,3,2),zgml(lx1,3)
  integer :: nxl, nyl, nzl

  integer :: ix, iy, iz

  do 10 ix=1,nxl
      h(ix,1,1)=(1.0-zgml(ix,1))*0.5
      h(ix,1,2)=(1.0+zgml(ix,1))*0.5
  10 END DO
  do 20 iy=1,nyl
      h(iy,2,1)=(1.0-zgml(iy,2))*0.5
      h(iy,2,2)=(1.0+zgml(iy,2))*0.5
  20 END DO
  if (if3d) then
      do 30 iz=1,nzl
          h(iz,3,1)=(1.0-zgml(iz,3))*0.5
          h(iz,3,2)=(1.0+zgml(iz,3))*0.5
      30 END DO
  else
      call rone(h(1,3,1),nzl)
      call rone(h(1,3,2),nzl)
  endif

  return
end subroutine sethmat

!-----------------------------------------------------------------------
subroutine setzgml (zgml,e,nxl,nyl,nzl,ifaxl)
  use kinds, only : DP
  use size_m, only : lx1, nx1, nx3, ny1, nz1
  use geom, only : ifgmsh3, ifrzer
  use wz_m, only : zgm1
  implicit none

  real(DP) :: zgml(lx1,3)
  integer :: e, nxl, nyl, nzl
  logical :: ifaxl

  integer :: k
  real(DP) :: zam1

  call rzero (zgml,3*nx1)


  if (nxl == 3 .AND. .NOT. ifaxl) then
      do k=1,3
          zgml(1,k) = -1
          zgml(2,k) =  0
          zgml(3,k) =  1
      enddo
  elseif (ifgmsh3 .AND. nxl == nx3) then
    write(*,*) "Oops: IFGMSH3"
#if 0
      call copy(zgml(1,1),zgm3(1,1),nx3)
      call copy(zgml(1,2),zgm3(1,2),ny3)
      call copy(zgml(1,3),zgm3(1,3),nz3)
      if (ifaxl .AND. ifrzer(e)) call copy(zgml(1,2),zam3,ny3)
#endif
  elseif (nxl == nx1) then
      call copy(zgml(1,1),zgm1(1,1),nx1)
      call copy(zgml(1,2),zgm1(1,2),ny1)
      call copy(zgml(1,3),zgm1(1,3),nz1)
      if (ifaxl .AND. ifrzer(e)) call copy(zgml(1,2),zam1,ny1)
  else
      call exitti('ABORT setzgml! $',nxl)
  endif

  return
  end subroutine setzgml

!-----------------------------------------------------------------------
!> \brief Compute Cartesian vector dot product.
REAL(DP) FUNCTION DOT(V1,V2,N)
  use kinds, only : DP
  implicit none

  integer :: n
  real(DP) :: V1(N),V2(N)
 
  real(DP) :: sum
  integer :: i

  SUM = 0
  DO 100 I=1,N
      SUM = SUM + V1(I)*V2(I)
  100 END DO
  DOT = SUM
  return
END FUNCTION DOT

!-----------------------------------------------------------------------
subroutine linquad(xl,yl,zl,nxl,nyl,nzl)
  use kinds, only : DP
  use size_m, only : ndim, nelt, lx1
  use input, only : ccurve, ifaxis
  implicit none

  integer :: nxl, nyl, nzl
  real(DP) :: xl(nxl*nyl*nzl,1),yl(nxl*nyl*nzl,1),zl(nxl*nyl*nzl,1)

  integer :: e, nedge, k
  logical :: ifmid

  nedge = 4 + 8*(ndim-2)

  do e=1,nelt ! Loop over all elements

      ifmid = .FALSE. 
      do k=1,nedge
          if (ccurve(k,e) == 'm') ifmid = .TRUE. 
      enddo

      if (lx1 == 2) ifmid = .FALSE. 
      if (ifmid) then
        write(*,*) "Oops: ifmid"
!max            call xyzquad(xl(1,e),yl(1,e),zl(1,e),nxl,nyl,nzl,e)
      else
          call xyzlin (xl(1,e),yl(1,e),zl(1,e),nxl,nyl,nzl,e,ifaxis)
      endif
  enddo

  return
end subroutine linquad

!-----------------------------------------------------------------------
!> \brief Generate bi- or trilinear mesh
subroutine xyzlin(xl,yl,zl,nxl,nyl,nzl,e,ifaxl)
  use kinds, only : DP
  use size_m
  use input
  implicit none

  integer :: nxl, nyl, nzl, e
  real(DP) :: xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
  logical :: ifaxl ! local ifaxis specification

! Preprocessor Corner notation:      Symmetric Corner notation:

!         4+-----+3    ^ s                    3+-----+4    ^ s
!         /     /|     |                      /     /|     |
!        /     / |     |                     /     / |     |
!      8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
!       |     | /     /                     |     | /     /
!       |     |/     /                      |     |/     /
!      5+-----+6    t                      5+-----+6    t

  integer, save :: indx(8) = (/ 1,2,4,3,5,6,8,7 /)

  integer, parameter :: ldw=4*lx1*ly1*lz1
  real(DP) :: xcb(2,2,2),ycb(2,2,2),zcb(2,2,2),w(ldw)

!  real(DP) :: zgml, jx,jy,jz,jxt,jyt,jzt, zlin
  real(DP) :: zgml(lx1,3),jx (lx1*2)
  real(DP) :: jxt(lx1*2),jyt(lx1*2),jzt(lx1*2),zlin(2)

  integer :: i, k, ndim2

  call setzgml (zgml,e,nxl,nyl,nzl,ifaxl)

  zlin(1) = -1
  zlin(2) =  1

  k = 1
  do i=1,nxl
      call fd_weights_full(zgml(i,1),zlin,1,0,jxt(k))
      call fd_weights_full(zgml(i,2),zlin,1,0,jyt(k))
      call fd_weights_full(zgml(i,3),zlin,1,0,jzt(k))
      k=k+2
  enddo
  call transpose(jx,nxl,jxt,2)

  ndim2 = 2**ndim
  ! Convert prex notation to lexicographical
  xcb = reshape(xc(indx(1:ndim2),e), (/2,2,2/))
  ycb = reshape(yc(indx(1:ndim2),e), (/2,2,2/))
  zcb = reshape(zc(indx(1:ndim2),e), (/2,2,2/))

!   Map R-S-T space into physical X-Y-Z space.

! NOTE:  Assumes nxl=nyl=nzl !

  call tensr3(xl,nxl,xcb,2,jx,jyt,jzt,w)
  call tensr3(yl,nxl,ycb,2,jx,jyt,jzt,w)
  call tensr3(zl,nxl,zcb,2,jx,jyt,jzt,w)

  return
end subroutine xyzlin
!-----------------------------------------------------------------------
