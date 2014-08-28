!-----------------------------------------------------------------------
!> \brief filter vx,vy,vz, and p by simple interpolation
subroutine q_filter(wght)
  use kinds, only : DP
  use size_m, only : lx1, ly1, lz1, lx2, lelt, lelv
  use size_m, only : nx1, ny1, nz1, nx2, nid, npert, ndim
  use input, only : ifflow, ifbase, if3d, ifsplit, ifldmhd, ifpert, ifheat
  use input, only : ifcvode, param, ifmhd, iflomach, npscal
  use soln, only : vx, vy, vz, pr, bx, by, bz, vxp, vyp, vzp, tp, t
  use tstep, only : ifield, istep
  use wz_m, only : zgm1, zgm2
  implicit none

  real(DP) :: wght

!     These are the dimensions that we interpolate onto for v and p:
  integer, parameter :: lxv=lx1-1
  integer, parameter :: lxp=lx2-1

  real(DP), save :: intdv(lx1,lx1)
  real(DP), save :: intuv(lx1,lx1)
  real(DP), save :: intdp(lx1,lx1)
  real(DP), save :: intup(lx1,lx1)
  real(DP), save :: intv(lx1,lx1)
  real(DP), save :: intp(lx1,lx1)

  real(DP) :: intw(lx1,lx1)
  real(DP) :: intt(lx1,lx1)
  real(DP) :: wk1  (lx1,lx1,lx1,lelt)
  real(DP) :: wk2  (lx1,lx1,lx1)
  real(DP) :: zgmv(lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)
  real(DP) :: tmax, omax
  common /ctmp0/ intw,intt &
  , wk1,wk2 &
  , zgmv,wgtv,zgmp,wgtp,tmax(100),omax(103)

!   outpost arrays
  integer, parameter :: lt=lx1*ly1*lz1*lelv

  character(18) :: sfmt

  integer, save :: icalld = 0

  integer :: imax, jmax, ncut, ifldt, j, mmax, nfldt, ifld, k
  integer, external :: iglmax
  real(DP) :: w0, umax, vmax, wmax, pmax
  real(DP), external :: glmax
  logical :: if_fltv

  imax = nid
  imax = iglmax(imax,1)
  jmax = iglmax(imax,1)
  if (icalld == 0) then
      icalld = 1
      ncut = param(101)+1
      call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
  elseif (icalld < 0) then   ! old (std.) filter
      icalld = 1
      call zwgll(zgmv,wgtv,lxv)
      call igllm(intuv,intw,zgmv,zgm1,lxv,nx1,lxv,nx1)
      call igllm(intdv,intw,zgm1,zgmv,nx1,lxv,nx1,lxv)
  
      call zwgl (zgmp,wgtp,lxp)
      call iglm (intup,intw,zgmp,zgm2,lxp,nx2,lxp,nx2)
      call iglm (intdp,intw,zgm2,zgmp,nx2,lxp,nx2,lxp)
  
  !        Multiply up and down interpolation into single operator
  
      call mxm(intup,nx2,intdp,lxp,intp,nx2)
      call mxm(intuv,nx1,intdv,lxv,intv,nx1)
  
  !        Weight the filter to make it a smooth (as opposed to truncated)
  !        decay in wave space

      w0 = 1.-wght
      call ident(intup,nx2)
      call add2sxy(intp,wght,intup,w0,nx2*nx2)

      call ident   (intuv,nx1)
      call add2sxy (intv ,wght,intuv,w0,nx1*nx1)

  endif

  ifldt  = ifield
!   ifield = 1

  if_fltv = .FALSE. 
  if ( ifflow .AND. .NOT. ifmhd ) if_fltv = .TRUE. 
  if ( ifield == 1  .AND. ifmhd ) if_fltv = .TRUE. 

!   Adam Peplinski; to take into account freezing of base flow
  if ( .NOT. ifbase             ) if_fltv = .FALSE. ! base-flow frozen

  if ( if_fltv ) then
      call filterq(vx,intv,nx1,nz1,wk1,wk2,intt,if3d,umax)
      call filterq(vy,intv,nx1,nz1,wk1,wk2,intt,if3d,vmax)
      if (if3d) &
      call filterq(vz,intv,nx1,nz1,wk1,wk2,intt,if3d,wmax)
      if (ifsplit .AND. .NOT. iflomach) &
      call filterq(pr,intv,nx1,nz1,wk1,wk2,intt,if3d,pmax)
  endif

  if (ifmhd .AND. ifield == ifldmhd) then
      call filterq(bx,intv,nx1,nz1,wk1,wk2,intt,if3d,umax)
      call filterq(by,intv,nx1,nz1,wk1,wk2,intt,if3d,vmax)
      if (if3d) &
      call filterq(bz,intv,nx1,nz1,wk1,wk2,intt,if3d,wmax)
  endif

  if (ifpert) then
      do j=1,npert

          ifield = 1
          call filterq(vxp(1,j),intv,nx1,nz1,wk1,wk2,intt,if3d,umax)
          call filterq(vyp(1,j),intv,nx1,nz1,wk1,wk2,intt,if3d,vmax)
          if (if3d) &
          call filterq(vzp(1,j),intv,nx1,nz1,wk1,wk2,intt,if3d,wmax)

          ifield = 2
          if (ifheat .AND. .NOT. ifcvode) &
          call filterq(tp(1,j,1),intv,nx1,nz1,wk1,wk2,intt,if3d,wmax)

      enddo
  endif

  mmax = 0
  if (ifflow) then
  !        pmax    = glmax(pmax,1)
      omax(1) = glmax(umax,1)
      omax(2) = glmax(vmax,1)
      omax(3) = glmax(wmax,1)
      mmax = ndim
  endif
           

  nfldt = 1+npscal
  if (ifheat .AND. .NOT. ifcvode) then
      do ifld=1,nfldt
          ifield = ifld + 1
          call filterq(t(1,1,1,1,ifld),intv &
          ,nx1,nz1,wk1,wk2,intt,if3d,tmax(ifld))
          mmax = mmax+1
          omax(mmax) = glmax(tmax(ifld),1)
      enddo
  endif

  if (nid == 0) then
      if (npscal == 0) then
      !           write(6,101) mmax
      !           write(sfmt,101) mmax
      ! 101       format('''(i8,1p',i1,'e12.4,a6)''')
      !           write(6,sfmt) istep,(omax(k),k=1,mmax),' qfilt'
      !         write(6,'(i8,1p4e12.4,a6)') istep,(omax(k),k=1,mmax),' qfilt'
      else
          if (if3d) then
              write(6,1) istep,ifield,umax,vmax,wmax
          else
              write(6,1) istep,ifield,umax,vmax
          endif
          1 format(4x,i7,i3,' qfilt:',1p3e12.4)
          if(ifheat .AND. .NOT. ifcvode) &
          write(6,'(1p50e12.4)') (tmax(k),k=1,nfldt)
      endif
  endif

  ifield = ifldt   ! RESTORE ifield


  return
  end subroutine q_filter
!-----------------------------------------------------------------------
    subroutine filterq(v,f,nx,nz,w1,w2,ft,if3d,dmax)

    use size_m
    use tstep

    real :: v(nx*nx*nz,nelt),w1(1),w2(1)
    logical :: if3d

    real :: f(nx,nx),ft(nx,nx)

    integer :: e

    call transpose(ft,nx,f,nx)

    nxyz=nx*nx*nz
    dmax = 0.


    nel = nelfld(ifield)


    if (if3d) then
        do e=1,nel
        !           Filter
            call copy(w2,v(1,e),nxyz)
            call mxm(f,nx,w2,nx,w1,nx*nx)
            i=1
            j=1
            do k=1,nx
                call mxm(w1(i),nx,ft,nx,w2(j),nx)
                i = i+nx*nx
                j = j+nx*nx
            enddo
            call mxm (w2,nx*nx,ft,nx,w1,nx)
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
        enddo
    
    else
        do e=1,nel
        !           Filter
            call copy(w1,v(1,e),nxyz)
            call mxm(f ,nx,w1,nx,w2,nx)
            call mxm(w2,nx,ft,nx,w1,nx)
        
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
        enddo
    endif

    return
    end subroutine filterq
!-----------------------------------------------------------------------
    subroutine local_grad3(ur,us,ut,u,N,e,D,Dt)
!     Output: ur,us,ut         Input:u,N,e,D,Dt
    real :: ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
    real :: u (0:N,0:N,0:N,1)
    real :: D (0:N,0:N),Dt(0:N,0:N)
    integer :: e

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
    subroutine gaujordf(a,m,n,indr,indc,ipiv,ierr,rmult)

!     Gauss-Jordan matrix inversion with full pivoting

!     Num. Rec. p. 30, 2nd Ed., Fortran


!     a     is an m x n matrix
!     rmult is a  work array of dimension m


    real :: a(m,n),rmult(m)
    integer :: indr(m),indc(n),ipiv(n)

!     call outmat(a,m,n,'ab4',n)
!     do i=1,m
!        write(6,1) (a(i,j),j=1,n)
!     enddo
! 1   format('mat: ',1p6e12.4)

    ierr = 0
    eps = 1.e-9
    call izero(ipiv,m)

    do k=1,m
        amx=0.
        do i=1,m                    ! Pivot search
            if (ipiv(i) /= 1) then
                do j=1,m
                    if (ipiv(j) == 0) then
                        if (abs(a(i,j)) >= amx) then
                            amx = abs(a(i,j))
                            ir  = i
                            jc  = j
                        endif
                    elseif (ipiv(j) > 1) then
                        ierr = -ipiv(j)
                        return
                    endif
                enddo
            endif
        enddo
        ipiv(jc) = ipiv(jc) + 1
    
    !       Swap rows
        if (ir /= jc) then
            do j=1,n
                tmp     = a(ir,j)
                a(ir,j) = a(jc,j)
                a(jc,j) = tmp
            enddo
        endif
        indr(k)=ir
        indc(k)=jc
    !       write(6 ,*) k,' Piv:',jc,a(jc,jc)
    !       write(28,*) k,' Piv:',jc,a(jc,jc)
        if (abs(a(jc,jc)) < eps) then
            write(6 ,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
            write(28,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
            ierr = jc
            call exitt
            return
        endif
        piv = 1./a(jc,jc)
        a(jc,jc)=1.
        do j=1,n
            a(jc,j) = a(jc,j)*piv
        enddo
    
        do j=1,n
            work    = a(jc,j)
            a(jc,j) = a(1 ,j)
            a(1 ,j) = work
        enddo
        do i=2,m
            rmult(i) = a(i,jc)
            a(i,jc)  = 0.
        enddo
    
        do j=1,n
            do i=2,m
                a(i,j) = a(i,j) - rmult(i)*a(1,j)
            enddo
        enddo
    
        do j=1,n
            work    = a(jc,j)
            a(jc,j) = a(1 ,j)
            a(1 ,j) = work
        enddo
    
    !       do i=1,m
    !          if (i.ne.jc) then
    !             rmult   = a(i,jc)
    !             a(i,jc) = 0.
    !             do j=1,n
    !                a(i,j) = a(i,j) - rmult*a(jc,j)
    !             enddo
    !          endif
    !       enddo
    
    enddo

!     Unscramble matrix
    do j=m,1,-1
        if (indr(j) /= indc(j)) then
            do i=1,m
                tmp=a(i,indr(j))
                a(i,indr(j))=a(i,indc(j))
                a(i,indc(j))=tmp
            enddo
        endif
    enddo

    return
    end subroutine gaujordf
!-----------------------------------------------------------------------
    subroutine legendre_poly(L,x,N)

!     Evaluate Legendre polynomials of degrees 0-N at point x

    real :: L(0:N)

    L(0) = 1.
    L(1) = x

    do j=2,N
        L(j) = ( (2*j-1) * x * L(j-1) - (j-1) * L(j-2) ) / j
    enddo

    return
    end subroutine legendre_poly
!-----------------------------------------------------------------------
    subroutine build_new_filter(intv,zpts,nx,kut,wght,nid)

!     This routing builds a 1D filter with a transfer function that
!     looks like:


!        ^
!    d_k |
!        |                 |
!     1  |__________      _v_
!        |          -_
!        |            \  wght
!        |             \  ___
!        |             |   ^
!     0  |-------------|---|>

!        0         c   N   k-->

!        Where c := N-kut is the point below which d_k = 1.



!      Here, nx = number of points

    real :: intv(nx,nx),zpts(nx)

    parameter (lm=40)
    parameter (lm2=lm*lm)
    real ::      phi(lm2),pht(lm2),diag(lm2),rmult(lm),Lj(lm)
    integer ::   indr(lm),indc(lm),ipiv(lm)

    if (nx > lm) then
        write(6,*) 'ABORT in build_new_filter:',nx,lm
        call exitt
    endif

    kj = 0
    n  = nx-1
    do j=1,nx
        z = zpts(j)
        call legendre_poly(Lj,z,n)
        kj = kj+1
        pht(kj) = Lj(1)
        kj = kj+1
        pht(kj) = Lj(2)
        do k=3,nx
            kj = kj+1
            pht(kj) = Lj(k)-Lj(k-2)
        enddo
    enddo
    call transpose (phi,nx,pht,nx)
    call copy      (pht,phi,nx*nx)
    call gaujordf  (pht,nx,nx,indr,indc,ipiv,ierr,rmult)

!     Set up transfer function

    call ident   (diag,nx)

    k0 = nx-kut
    do k=k0+1,nx
        kk = k+nx*(k-1)
        amp = wght*(k-k0)*(k-k0)/(kut*kut)   ! quadratic growth
        diag(kk) = 1.-amp
    enddo

    call mxm  (diag,nx,pht,nx,intv,nx)      !          -1
    call mxm  (phi ,nx,intv,nx,pht,nx)      !     V D V
    call copy (intv,pht,nx*nx)

    do k=1,nx*nx
        pht(k) = 1.-diag(k)
    enddo
    np1 = nx+1
    if (nid == 0) then
        write(6,6) 'filt amp',(pht (k),k=1,nx*nx,np1)
        write(6,6) 'filt trn',(diag(k),k=1,nx*nx,np1)
        6 format(a8,16f7.4,6(/,8x,16f7.4))
    endif

    return
    end subroutine build_new_filter
!-----------------------------------------------------------------------
    function ran1(idum)

    integer :: idum,ia,im,iq,ir,ntab,ndiv
    real ::    ran1,am,eps,rnmx

    parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836)
    parameter (ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)

!     Numerical Rec. in Fortran, 2nd eD.  P. 271

    integer :: j,k
    integer :: iv(ntab),iy
    save    iv,iy
    data    iv,iy /ntab*0,0/

    if (idum <= 0 .OR. iy == 0) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
            k    = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if(idum < 0) idum = idum+im
            if (j <= ntab) iv(j) = idum
        enddo
        iy = iv(1)
    endif
    k    = idum/iq
    idum = ia*(idum-k*iq)-ir*k
    if(idum < 0) idum = idum+im
    j     = 1+iy/ndiv
    iy    = iv(j)
    iv(j) = idum
    ran1  = min(am*iy,rnmx)
!     ran1  = cos(ran1*1.e8)

    return
    end function ran1
!-----------------------------------------------------------------------
    subroutine rand_fld_h1(x)

    use size_m
    real :: x(1)

    n=nx1*ny1*nz1*nelt
    id = n
    do i=1,n
        x(i) = ran1(id)
    enddo
    call dsavg(x)

    return
    end subroutine rand_fld_h1
!-----------------------------------------------------------------------
