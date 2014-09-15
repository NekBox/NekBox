!     Dimension file to be included

!     HCUBE array dimensions
module size_m

    integer :: ldim, lx1, ly1, lz1
    parameter (ldim=3)
    parameter (lx1=8,ly1=lx1,lz1=lx1)
    integer, parameter :: lxd=12,lyd=lxd,lzd=lxd
    integer, parameter :: lelx=1,lely=1,lelz=1

    integer :: lzl
    parameter (lzl=3 + 2*(ldim-3))

    integer :: lx2, ly2, lz2, lx3, ly3, lz3
    parameter (lx2=lx1)
    parameter (ly2=ly1)
    parameter (lz2=lz1)
    parameter (lx3=lx2)
    parameter (ly3=ly2)
    parameter (lz3=lz2)

    integer :: lelg, lp, lelt, lelv
    parameter (lelg = 2048)
    parameter (lp = 4)
    parameter (lelt= (lelg + 4 - 1)/lp, lelv=lelt)


!     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
!     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
!     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)

    integer :: lpelv, lpelt, lpert, lpx1, lpy1, lpz1, lpx2, lpy2, lpz2
    parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
    parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
    parameter (lpx2=1,lpy2=1,lpz2=1)

!     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
!     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
!     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)

    integer :: lbelv, lbelt, lbx1, lby1, lbz1, lbx2, lby2, lbz2
    parameter (lbelv=1,lbelt=1)                ! MHD
    parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
    parameter (lbx2=1,lby2=1,lbz2=1)

!     LX1M=LX1 when there are moving meshes; =1 otherwise
    integer :: lx1m, ly1m, lz1m, ldimt, ldimt1, ldimt3
    parameter (lx1m=1,ly1m=1,lz1m=1)
!    parameter (lx1m=lx1,ly1m=ly1,lz1m=lz1)
    parameter (ldimt= 4)                       ! 3 passive scalars + T
    parameter (ldimt1=ldimt+1)
    parameter (ldimt3=ldimt+3)

!     Note:  In the new code, LELGEC should be about sqrt(LELG)

    integer :: lelgec, lxyz2, lxz21
    PARAMETER (LELGEC = 1)
    PARAMETER (LXYZ2  = 1)
    PARAMETER (LXZ21  = 1)

    integer :: lmaxv, lmaxt, lmaxp, lxz, lorder, maxobj, maxmbr, lhis
    PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
    PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
    PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
    PARAMETER (LXZ=LX1*LZ1)
    PARAMETER (LORDER=3)
    PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
    PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
! Note: lhis*np > npoints in hpts.in

!     Common Block Dimensions

    integer :: lctmp0, lctmp1
    PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
    PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)

!     The parameter LVEC controls whether an additional 42 field arrays
!     are required for Steady State Solutions.  If you are not using
!     Steady State, it is recommended that LVEC=1.

    integer :: lvec
    PARAMETER (LVEC=1)

!     Uzawa projection array dimensions

    integer :: mxprev, lgmres
    parameter (mxprev = 40)
    parameter (lgmres = 10)

!     Split projection array dimensions

    integer :: lmvec, lsvec, lstore
    parameter(lmvec = 1)
    parameter(lsvec = 1)
    parameter(lstore=lmvec*lsvec)

!     NONCONFORMING STUFF

    integer :: maxmor
    parameter (maxmor = lelt)

!     Array dimensions

    integer :: NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2 &
    ,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID &
    ,NXD,NYD,NZD

! automatically added by makenek
    integer, parameter :: lxo   = lx1 ! max output grid size (lxo>=lx1)

! automatically added by makenek
    integer, parameter :: lpart = 1   ! max number of particles

! automatically added by makenek
    integer :: ax1,ay1,az1,ax2,ay2,az2
    parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! running averages

! automatically added by makenek
    integer, parameter :: lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1 !New Pressure Preconditioner

! automatically added by makenek
    integer, parameter :: lfdm=0  ! == 1 for fast diagonalization method

end module size_m
