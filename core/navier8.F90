!-----------------------------------------------------------------------
!> \brief Given global array, vertex, pointing to hex vertices, set up
!! a new array of global pointers for an nx^ndim set of elements.
subroutine set_vert(glo_num,ngv,nx,nel,vertex,ifcenter)
  use kinds, only : i8
  use size_m, only : nid
  use input, only : if3d
  implicit none

  integer(i8) :: glo_num(1),ngv
  integer :: vertex(*),nx, nel
  logical :: ifcenter

  integer :: nz

  if (if3d) then
      call setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
  else
!max        call setvert2d(glo_num,ngv,nx,nel,vertex,ifcenter)
  endif

!   Check for single-element periodicity 'p' bc
  nz = 1
  if (if3d) nz = nx
  call check_p_bc(glo_num,nx,nx,nz,nel)

  if(nid == 0) write(6,*) 'call usrsetvert'
  call usrsetvert(glo_num,nel,nx,nx,nx)
  if(nid == 0) write(6,'(A,/)') ' done :: usrsetvert'

  return
end subroutine set_vert

!> ?
subroutine set_up_h1_crs
  use kinds, only : DP, i8
  use size_m, only : nx1, ny1, nz1, nelv, ndim, nid
  use size_m, only : lx1, ly1, lz1, lelv
  use domain, only : nx_crs, nxyz_c, se_to_gcrs, lcr
  use geom, only : ifvcor, ifbcor
  use input, only : if3d, ifldmhd, cbc
  use mesh, only : vertex
  use parallel, only : xxth, mp=>np, nekcomm
  use tstep, only : ifield
  implicit none

  integer :: gs_handle
  integer :: null_space

  character(3) :: cb
  real(DP) :: cmlt(lcr,lelv),mask(lcr,lelv)
  integer, allocatable :: ia(:,:,:), ja(:,:,:)
  real :: z

  real(DP), allocatable :: a(:)
 
  real(DP), allocatable :: h1(:), h2(:), w1(:), w2(:)

  integer(i8) :: ngv

  real(DP) :: t0
  integer :: nxc, ncr, ntot, nzc, nfaces, ie, iface, nz, n
  real(DP), external :: dnekclock

  allocate(ia(lcr,lcr,lelv), ja(lcr,lcr,lelv))
  t0 = dnekclock()

!     nxc is order of coarse grid space + 1, nxc=2, linear, 3=quad,etc.
!     nxc=param(82)
!     if (nxc.gt.lxc) then
!        nxc=lxc
!        write(6,*) 'WARNING :: coarse grid space too large',nxc,lxc
!     endif
!     if (nxc.lt.2) nxc=2

  nxc     = 2
  nx_crs  = nxc

  if(nid == 0) write(6,*) 'setup h1 coarse grid, nx_crs=', nx_crs

  ncr     = nxc**ndim
  nxyz_c  = ncr

!   Set SEM_to_GLOB

  call get_vertex
  call set_vert(se_to_gcrs,ngv,nxc,nelv,vertex, .TRUE. )

!   Set mask
  z=0
  ntot=nelv*nxyz_c
  nzc=1
  if (if3d) nzc=nxc
  call rone(mask,ntot)
  call rone(cmlt,ntot)
  nfaces=2*ndim
!   ifield=1			!c? avo: set in set_overlap through 'TSTEP'?

  if (ifield == 1) then
      do ie=1,nelv
          do iface=1,nfaces
              cb=cbc(iface,ie,ifield)
              if (cb == 'O  '  .OR.  cb == 'ON '  .OR.  cb == 'MM '  .OR. &
              cb == 'mm '  .OR.  cb == 'ms '  .OR.  cb == 'MS ') &
              call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
          enddo
      enddo
  elseif (ifield == ifldmhd) then   ! no ifmhd ?avo?
      do ie=1,nelv
          do iface=1,nfaces
              cb=cbc(iface,ie,ifield)
              if (cb == 'ndd'  .OR.  cb == 'dnd'  .OR.  cb == 'ddn') &
              call facev(mask,ie,iface,z,nxc,nxc,nzc)
          enddo
      enddo
  endif

!   Set global index of dirichlet nodes to zero; xxt will ignore them

  call gs_setup(gs_handle,se_to_gcrs,ntot,nekcomm,mp)
  call gs_op   (gs_handle,mask,1,2,0)  !  "*"
  call gs_op   (gs_handle,cmlt,1,1,0)  !  "+"
  call gs_free (gs_handle)
  call set_jl_crs_mask(ntot,mask,se_to_gcrs)

  call invcol1(cmlt,ntot)

!   Setup local SEM-based Neumann operators (for now, just full...)

!    if (param(51).eq.1) then     ! old coarse grid
!       nxyz1=nx1*ny1*nz1
!       lda = 27*nxyz1*lelt
!       ldw =  7*nxyz1*lelt
!       call get_local_crs(a,lda,nxc,h1,h2,w,ldw)
!    else
!      NOTE: a(),h1,...,w2() must all be large enough
  n = nx1*ny1*nz1*nelv
  allocate(h1(nx1*ny1*nz1*nelv),h2(nx1*ny1*nz1*nelv))
  allocate(a(27*lx1*ly1*lz1*lelv))
  allocate(w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelv))
  h1=1._dp; h2 = 0._dp
  call get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)
  deallocate(h1,h2,w1,w2)
!    endif

  call set_mat_ij(ia,ja,ncr,nelv)
  null_space=0
  if (ifield == 1) then
      if (ifvcor)  null_space=1
  elseif (ifield == ifldmhd) then
      if (ifbcor)  null_space=1
  endif

  nz=ncr*ncr*nelv
  call crs_setup(xxth(ifield),nekcomm,mp, ntot,se_to_gcrs, &
  nz,ia,ja,a, null_space)
  deallocate(a)
!   call crs_stats(xxth(ifield))

  t0 = dnekclock()-t0
  if (nid == 0) then
      write(6,*) 'done :: setup h1 coarse grid ',t0, ' sec'
      write(6,*) ' '
  endif

  return
end subroutine set_up_h1_crs

!-----------------------------------------------------------------------
subroutine set_jl_crs_mask(n, mask, se_to_gcrs)
  use kinds, only : DP, i8
  implicit none
  integer :: n
  real(DP) :: mask(*)
  integer(i8) :: se_to_gcrs(*)

  integer :: i
  do i=1,n
      if(mask(i) < 0.1) se_to_gcrs(i)=0
  enddo
  return
end subroutine set_jl_crs_mask

!-----------------------------------------------------------------------
subroutine set_mat_ij(ia,ja,n,ne)
  implicit none
  integer :: n,ne
  integer :: ia(n,n,ne), ja(n,n,ne)

  integer :: i,j,ie
  do ie=1,ne
      do j=1,n
          do i=1,n
              ia(i,j,ie)=(ie-1)*n+i-1
              ja(i,j,ie)=(ie-1)*n+j-1
          enddo
      enddo
  enddo
  return
end subroutine set_mat_ij

!-----------------------------------------------------------------------
!> \brief Use Heap Sort (p 231 Num. Rec., 1st Ed.)
subroutine ituple_sort(a,lda,n,key,nkey,ind,aa)
  implicit none

  integer :: lda, n, nkey
  integer :: a(lda,n),aa(lda)
  integer :: ind(*),key(nkey)
  logical :: iftuple_ialtb

  integer :: j, L, ir, ii, i

  dO 10 j=1,n
      ind(j)=j
  10 END DO

  if (n <= 1) return
  L=n/2+1
  ir=n
  100 continue
  if (l > 1) then
      l=l-1
  !           aa  = a  (l)
      call icopy(aa,a(1,l),lda)
      ii  = ind(l)
  else
  !           aa =   a(ir)
      call icopy(aa,a(1,ir),lda)
      ii = ind(ir)
  !           a(ir) =   a( 1)
      call icopy(a(1,ir),a(1,1),lda)
      ind(ir) = ind( 1)
      ir=ir-1
      if (ir == 1) then
      !              a(1) = aa
          call icopy(a(1,1),aa,lda)
          ind(1) = ii
          return
      endif
  endif
  i=l
  j=l+l
  200 continue
  if (j <= ir) then
      if (j < ir) then
          if (iftuple_ialtb(a(1,j),a(1,j+1),key,nkey)) j=j+1
      endif
      if (iftuple_ialtb(aa,a(1,j),key,nkey)) then
      !              a(i) = a(j)
          call icopy(a(1,i),a(1,j),lda)
          ind(i) = ind(j)
          i=j
          j=j+j
      else
          j=ir+1
      endif
      GOTO 200
  endif
!       a(i) = aa
  call icopy(a(1,i),aa,lda)
  ind(i) = ii
  GOTO 100
end subroutine ituple_sort

!-----------------------------------------------------------------------
logical function iftuple_ialtb(a,b,key,nkey)
  implicit none
  integer :: nkey
  integer :: a(*),b(*)
  integer :: key(nkey)

  integer :: i, k
 
  do i=1,nkey
      k=key(i)
      if (a(k) < b(k)) then
          iftuple_ialtb = .TRUE. 
          return
      elseif (a(k) > b(k)) then
          iftuple_ialtb = .FALSE. 
          return
      endif
  enddo
  iftuple_ialtb = .FALSE. 
  return
end function iftuple_ialtb

!-----------------------------------------------------------------------
logical function iftuple_ianeb(a,b,key,nkey)
  implicit none
  integer :: nkey
  integer :: a(*),b(*)
  integer :: key(nkey)

  integer :: i, k
  do i=1,nkey
      k=key(i)
      if (a(k) /= b(k)) then
          iftuple_ianeb = .TRUE. 
          return
      endif
  enddo
  iftuple_ianeb = .FALSE. 
  return
end function iftuple_ianeb

!-----------------------------------------------------------------------
!> \brief  Spectral interpolation from A to B via tensor products
!!     -  scratch arrays: w(na*na*nb + nb*nb*na)
!!     5/3/00  -- this routine replaces specmp in navier1.f, which
!!                has a potential memory problem
subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
  use kinds, only : DP
  implicit none

  logical :: if3d
  integer :: nb, na, ldw
  real(DP) :: b(nb,nb,nb),a(na,na,na)
  real(DP) :: w(ldw)

  integer :: ltest, nab, nbb, k, l, iz
  real(DP) :: ba, ab

  ltest = na*nb
  if (if3d) ltest = na*na*nb + nb*na*na
  if (ldw < ltest) then
      write(6,*) 'ERROR specmp:',ldw,ltest,if3d
      call exitt
  endif

  if (if3d) then
      nab = na*nb
      nbb = nb*nb
      call mxm(ba,nb,a,na,w,na*na)
      k=1
      l=na*na*nb + 1
      do iz=1,na
          call mxm(w(k),nb,ab,na,w(l),nb)
          k=k+nab
          l=l+nbb
      enddo
      l=na*na*nb + 1
      call mxm(w(l),nbb,ab,na,b,nb)
  else
      call mxm(ba,nb,a,na,w,na)
      call mxm(w,nb,ab,na,b,nb)
  endif
  return
end subroutine specmpn

!-----------------------------------------------------------------------
!> \brief Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
subroutine irank(A,IND,N)
  implicit none
  integer :: A(*),IND(*), n

  integer :: j, L, ir, i, q, indx

  if (n <= 1) return
  DO 10 J=1,N
      IND(j)=j
  10 END DO

  if (n == 1) return
  L=n/2+1
  ir=n
  100 continue
  IF (l > 1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
  ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
          ind(1)=indx
          return
      endif
  ENDIF
  i=l
  j=l+l
  200 continue
  IF (J <= IR) THEN
      IF (J < IR) THEN
          IF ( A(IND(j)) < A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q < A(IND(j))) THEN
          IND(I)=IND(J)
          I=J
          J=J+J
      ELSE
          J=IR+1
      ENDIF
      GOTO 200
  ENDIF
  IND(I)=INDX
  GOTO 100
end subroutine irank

!-----------------------------------------------------------------------
!> \brief This routine generates Nelv submatrices of order ncl using
!! Galerkin projection
subroutine get_local_crs_galerkin(a,ncl,nxc,h1,h2,w1,w2)
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1, nelv, lx1, ldim
  implicit none

  integer :: ncl, nxc
  real ::    a(ncl,ncl,*),h1(*),h2(*)
  real ::    w1(nx1*ny1*nz1,nelv),w2(nx1*ny1*nz1,nelv)

  integer, parameter :: lcrd=lx1**ldim
  real(DP) :: b(lcrd,8)

  integer :: e, i, j, isd, imsh, nxyz
  real(DP), external :: vlsc2

  do j=1,ncl
      call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
  enddo

  isd  = 1
  imsh = 1

  nxyz = nx1*ny1*nz1
  do j = 1,ncl
      do e = 1,nelv
          call copy(w1(1,e),b(1,j),nxyz)
      enddo

      call axhelm (w2,w1,h1,h2,imsh,isd)        ! A^e * bj

      do e = 1,nelv
          do i = 1,ncl
              a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * A^e * bj
          enddo
      enddo

  enddo

  return
end subroutine get_local_crs_galerkin

!-----------------------------------------------------------------------
subroutine gen_crs_basis(b,j) ! bi- tri-linear
  use kinds, only : DP
  use size_m
  implicit none

  real(DP) :: b(nx1,ny1,nz1)
  integer :: j

  real(DP) :: z0(lx1),z1(lx1)
  real(DP) :: zr(lx1),zs(lx1),zt(lx1)

  integer :: p,q,r
  integer :: i

  call zwgll(zr,zs,nx1)

  do i=1,nx1
      z0(i) = .5*(1-zr(i))  ! 1-->0
      z1(i) = .5*(1+zr(i))  ! 0-->1
  enddo

  call copy(zr,z0,nx1)
  call copy(zs,z0,nx1)
  call copy(zt,z0,nx1)

  if (mod(j,2) == 0)                        call copy(zr,z1,nx1)
  if (j == 3 .OR. j == 4 .OR. j == 7 .OR. j == 8) call copy(zs,z1,nx1)
  if (j > 4)                               call copy(zt,z1,nx1)

  if (ndim == 3) then
      do r=1,nx1
          do q=1,nx1
              do p=1,nx1
                  b(p,q,r) = zr(p)*zs(q)*zt(r)
              enddo
          enddo
      enddo
  else
      do q=1,nx1
          do p=1,nx1
              b(p,q,1) = zr(p)*zs(q)
          enddo
      enddo
  endif

  return
end subroutine gen_crs_basis

!-----------------------------------------------------------------------
subroutine get_vertex
  use zper, only : ifgtp
  implicit none

  integer, save :: icalld = 0

  if (icalld > 0) return
  icalld = 1

  if (ifgtp) then
      write(*,*) "Oops: ifgtp"
!max        call gen_gtp_vertex    (vertex, ncrnr)
  else
      call get_vert
  endif

  return
end subroutine get_vertex

!-----------------------------------------------------------------------
subroutine assign_gllnid(gllnid,iunsort,nelgt,nelgv,np)
  implicit none
  integer :: gllnid(*),iunsort(*),nelgt,nelgv, np
  integer :: eg

  integer :: nelgs, nnpstr, npstar, log2p, np2
  integer, external :: ivlmax, log2

  log2p = log2(np)
  np2   = 2**log2p
  if (np2 == np .AND. nelgv == nelgt) then   ! std power of 2 case

      npstar = ivlmax(gllnid,nelgt)+1
      nnpstr = npstar/np
      do eg=1,nelgt
          gllnid(eg) = gllnid(eg)/nnpstr
      enddo

      return

  elseif (np2 == np) then   ! std power of 2 case, conjugate heat xfer

  !        Assign fluid elements
      npstar = max(np,ivlmax(gllnid,nelgv)+1)
      nnpstr = npstar/np
      do eg=1,nelgv
          gllnid(eg) = gllnid(eg)/nnpstr
      enddo

  !        Assign solid elements
      nelgs  = nelgt-nelgv  ! number of solid elements
      npstar = max(np,ivlmax(gllnid(nelgv+1),nelgs)+1)
      nnpstr = npstar/np
      do eg=nelgv+1,nelgt
          gllnid(eg) = gllnid(eg)/nnpstr
      enddo

      return

  elseif (nelgv /= nelgt) then
      call exitti &
      ('Conjugate heat transfer requires P=power of 2.$',np)
  endif


!  Below is the code for P a non-power of two:

!  Split the sorted gllnid array (read from .map file)
!  into np contiguous partitions.

!  To load balance the partitions in case of mod(nelgt,np)>0
!  add 1 contiguous entry out of the sorted list to NODE_i
!  where i = np-mod(nelgt,np) ... np

#if 0
  nel   = nelgt/np       ! number of elements per processor
  nmod  = mod(nelgt,np)  ! bounded between 1 ... np-1
  npp   = np - nmod      ! how many paritions of size nel
     
! sort gllnid
  call isort(gllnid,iunsort,nelgt)

! setup partitions of size nel
  k   = 0
  do ip = 0,npp-1
      do e = 1,nel
          k = k + 1
          gllnid(k) = ip
      enddo
  enddo
! setup partitions of size nel+1
  if(nmod > 0) then
      do ip = npp,np-1
          do e = 1,nel+1
              k = k + 1
              gllnid(k) = ip
          enddo
      enddo
  endif

! unddo sorting to restore initial ordering by
! global element number
  call iswapt_ip(gllnid,iunsort,nelgt)
#endif

  return
end subroutine assign_gllnid

!-----------------------------------------------------------------------
subroutine get_vert
  use size_m, only : ndim
  use input, only : ifmoab
  use mesh, only : vertex
  use parallel, only : nelgt
  use zper, only : ifgfdm
  implicit none

  integer :: ncrnr
  integer, save :: icalld = 0

  if (icalld > 0) return
  icalld = 1

  ncrnr = 2**ndim

  if (ifmoab) then
#ifdef MOAB
      call nekMOAB_loadConn (vertex, nelgt, ncrnr)
#endif
  else
      call get_vert_map(vertex, ncrnr, nelgt, '.map', ifgfdm)
  endif

  return
end subroutine get_vert

!-----------------------------------------------------------------------
subroutine get_vert_map(vertex, nlv, nel, suffix, ifgfdm)
  use size_m, only : lx1, ly1, lz1, lelv, ldim, nid, nelt, nelv, ndim
  use input, only : reafle
  use parallel, only : np, gllnid, isize, gllel, nelgt, nelgv, cr_h
  use string, only : ltrunc
  implicit none

  logical :: ifgfdm

  integer :: nlv, nel
  integer :: vertex(nlv,*)
  character(4) :: suffix

  integer, parameter :: mdw=2+2**ldim
  integer, parameter :: ndw=7*lx1*ly1*lz1*lelv/mdw

  integer, allocatable :: wk(:,:)   ! room for long ints, if desired

  integer :: e,eg,eg0,eg1, iok, lfname, neli, nnzi, npass, len, msg_id
  integer :: ipass, m, k, ntuple, lng, i, key, nkey, iflag, nv, mid
  integer, external :: iglmax, irecv

  character(132) :: mapfle
  character(1) ::   mapfle1(132)
  equivalence  (mapfle,mapfle1)

  allocate(wk(mdw,ndw))

  iok = 0
  if (nid == 0) then
      lfname = ltrunc(reafle,132) - 4
      call blank (mapfle,132)
      call chcopy(mapfle,reafle,lfname)
      call chcopy(mapfle1(lfname+1),suffix,4)
      open(unit=80,file=mapfle,status='old',err=99)
      read(80,*,err=99) neli,nnzi
      iok = 1
  endif
  99 continue
  iok = iglmax(iok,1)
  if (iok == 0) goto 999     ! Mapfile not found

  if (nid == 0) then
      neli = iglmax(neli,1)   ! communicate to all procs
  else
      neli = 0
      neli = iglmax(neli,1)   ! communicate neli to all procs
  endif

  npass = 1 + (neli/ndw)
  if (npass > np) then
      if (nid == 0) write(6,*) npass,np,neli,ndw,'Error get_vert_map'
      call exitt
  endif

  len = 4*mdw*ndw
  if (nid > 0 .AND. nid < npass) msg_id=irecv(nid,wk,len)
  call nekgsync

  if (nid == 0) then
      eg0 = 0
      do ipass=1,npass
          eg1 = min(eg0+ndw,neli)
          m   = 0
          do eg=eg0+1,eg1
              m = m+1
              read(80,*,end=998) (wk(k,m),k=2,mdw)
              if( .NOT. ifgfdm)  gllnid(eg) = wk(2,m)  !proc map,  must still be divided
              wk(1,m)    = eg
          enddo
          if (ipass < npass) call csend(ipass,wk,len,ipass,0) !send to ipass
          eg0 = eg1
      enddo
      close(80)
      ntuple = m
  elseif (nid < npass) then
    write(*,*) "Oops: nid < npass"
#if 0
      call msgwait(msg_id)
      ntuple = ndw
#endif
  else
      ntuple = 0
  endif

!   Distribute and assign partitions
  if ( .NOT. ifgfdm) then             ! gllnid is already assigned for gfdm
      lng = isize*neli
      call bcast(gllnid,lng)
      call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch

  !       if(nid.eq.0) then
  !         write(99,*) (gllnid(i),i=1,nelgt)
  !       endif
  !       call exitt
  endif

  nelt=0 !     Count number of elements on this processor
  nelv=0
  do eg=1,neli
      if (gllnid(eg) == nid) then
          if (eg <= nelgv) nelv=nelv+1
          if (eg <= nelgt) nelt=nelt+1
      endif
  enddo
  if (np <= 64) write(6,*) nid,nelv,nelt,nelgv,nelgt,' NELV'

!   NOW: crystal route vertex by processor id

  do i=1,ntuple
      eg=wk(1,i)
      wk(2,i)=gllnid(eg)        ! processor id for element eg
  enddo

  key = 2  ! processor id is in wk(2,:)
  call crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)

  if ( .NOT. ifgfdm) then            ! no sorting for gfdm?
      key = 1  ! Sort tuple list by eg := wk(1,:)
      nkey = 1
      call crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
  endif

  iflag = 0
  if (ntuple /= nelt) then
      write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
      write(6,*) 'Check that .map file and .rea file agree'
      iflag=1
  else
      nv = 2**ndim
      do e=1,nelt
          call icopy(vertex(1,e),wk(3,e),nv)
      enddo
  endif

  iflag = iglmax(iflag,1)
  if (iflag > 0) then
      do mid=0,np-1
          call nekgsync
          if (mid == nid) &
          write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FB'
          call nekgsync
      enddo
      call nekgsync
      call exitt
  endif

  return

  999 continue
  if (nid == 0) write(6,*) 'ABORT: Could not find map file ',mapfle
  call exitt

  998 continue
  if (nid == 0) write(6,*)ipass,npass,eg0,eg1,mdw,m,eg,'get v fail'
  call exitt0  ! Emergency exit

  return
end subroutine get_vert_map
!-----------------------------------------------------------------------
!> \brief Compute rank of each unique entry a(1,i)
!!   Output:   ind(i)  i=1,...,n    (global) rank of entry a(*,i)
!!             nn  = max(rank)
!!             a(j,i) is permuted
!!   Input:    a(j,i) j=1,...,m;  i=1,...,n
!!             m      :   leading dim. of v  (ldv must be .ge. m)
!!             key    :   sort key
!!             nkey   :
!!   Although not mandatory, this ranking procedure is probably
!!   most effectively employed when the keys are pre-sorted. Thus,
!!   the option is provided to sort vi() prior to the ranking.
subroutine irank_vecn(ind,nn,a,m,n,key,nkey,aa)
  implicit none
  integer :: nn, m, n, nkey
  integer :: ind(n),a(m,n)
  integer :: key(nkey),aa(m)

  logical :: iftuple_ianeb,a_ne_b
  integer :: nk, i

  nk = min(nkey,m)
  call ituple_sort(a,m,n,key,nk,ind,aa)

!   Find unique a's
  call icopy(aa,a,m)
  nn     = 1
  ind(1) = nn

  do i=2,n
      a_ne_b = iftuple_ianeb(aa,a(1,i),key,nk)
      if (a_ne_b) then
          call icopy(aa,a(1,i),m)
          nn = nn+1
      endif
      ind(i) = nn ! set ind() to rank
  enddo

  return
end subroutine irank_vecn

!-----------------------------------------------------------------------
!> \brief Return a unique rank for each matched tuple set. Global.  Balanced.
!! tuple is destroyed.
!! By "balanced" we mean that none of the tuple entries is likely to
!! be much more uniquely populated than any other, so that any of
!! the tuples can serve as an initial (parallel) sort key
!! First two slots in tuple(:,i) assumed empty
subroutine gbtuple_rank(tuple,m,n,nmax,cr_h,nid,np,ind)
  implicit none
  integer :: m, n, nmax, nid, np
  integer :: ind(nmax),tuple(m,nmax),cr_h

  integer, parameter :: mmax=40
  integer :: key(mmax),wtuple(mmax)
  integer :: nk, i, k, nkey, ni, ky, nimx, nu, nu_tot, nu_prior
  integer, external :: iglmax, igl_running_sum

  if (m > mmax) then
      write(6,*) nid,m,mmax,' gbtuple_rank fail'
      call exitt
  endif

  do i=1,n
      tuple(1,i) = mod(tuple(3,i),np) ! destination processor
      tuple(2,i) = i                  ! return location
  enddo

  ni= n
  ky=1  ! Assumes crystal_new already called
  call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

  nimx = iglmax(ni,1)
  if (ni > nmax)   write(6,*) ni,nmax,n,'cr_xfer problem, A'
  if (nimx > nmax) call exitt

  nkey = m-2
  do k=1,nkey
      key(k) = k+2
  enddo

  call irank_vecn(ind,nu,tuple,m,ni,key,nkey,wtuple)! tuple re-ordered,
! but contents same

  nu_tot   = igl_running_sum(nu) ! running sum over P processors
  nu_prior = nu_tot - nu

  do i=1,ni
      tuple(3,i) = ind(i) + nu_prior  ! global ranking
  enddo

  call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

  nk = 1  ! restore to original order, local rank: 2; global: 3
  key(1) = 2
  call ituple_sort(tuple,m,n,key,nk,ind,wtuple)

  return
end subroutine gbtuple_rank

!-----------------------------------------------------------------------
!> \brief setup unique ids for dssum.
!!  note:
!!  total number of unique vertices, edges and faces has to be smaller
!!  than 2**31 (integer-4 limit).
!!  if nelgt < 2**31/12 we're ok for sure (independent of N)!
subroutine setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
  use kinds, only : i8
  use size_m, only : lelt, ndim, nid, nelt
  use parallel, only : cr_h, np, nelgt, lglel
  use topol, only : icface, skpdat
  implicit none

  integer(i8) :: glo_num(*),ngv
  integer :: vertex(0:1,0:1,0:1,*),nx
  logical :: ifcenter

  integer ::  edge(0:1,0:1,0:1,3,lelt),enum(12,lelt),fnum(6,lelt)

  integer, parameter :: nsafe=8   ! OFTEN, nsafe=2 suffices

  integer :: vertex_flat(8)
  integer :: etuple(4,12*lelt*nsafe)
  integer :: ftuple(5,6,lelt*nsafe)
  integer :: ind(12*lelt*nsafe)
  equivalence  (etuple,ftuple)

  integer :: gvf(4),facet(4),key(3),e
  logical :: ifij
        
  integer(i8) :: igv,ig0
  integer(i8) :: ngvv,ngve,ngvs,ngvi,ngvm
  integer(i8) :: n_on_edge,n_on_face,n_in_interior
  integer(i8) :: i8glmax

  integer :: ny, nz, nxyz, nlv, nel, k, j, i, il, ile, kswap, m, n, nmax
  integer :: n_unique_edges, iedg_loc, i0, i0e, nfaces, ncrnr, ifac, icrn
  integer :: n_unique_faces, iface, i1, is, j0, j1, js, idir, jdir, it, jt
  integer :: nxx, l
  integer, external :: iglmax

  ny   = nx
  nz   = nx
  nxyz = nx*ny*nz

  key(1)=1
  key(2)=2
  key(3)=3

  ngvs = -1

!   Assign hypercube ordering of vertices
!   -------------------------------------

!   Count number of unique vertices
  nlv  = 2**ndim
  ngvv = iglmax(vertex,nlv*nel)

  do e=1,nel
      do k=0,1
          do j=0,1
              do i=0,1
              !           Local to global node number (vertex)
                  il  = 1 + (nx-1)*i + nx*(nx-1)*j + nx*nx*(nx-1)*k
                  ile = il + nx*ny*nz*(e-1)
                  glo_num(ile)   = vertex(i,j,k,e)
              enddo
          enddo
      enddo
  enddo
  ngv  = ngvv

  if (nx == 2) return

!   Assign global vertex numbers to SEM nodes on each edge
!   ------------------------------------------------------

!   Assign edge labels by bounding vertices.
  do e=1,nel
      do k=0,1
          do j=0,1
              do i=0,1
                  edge(i,j,k,1,e) = vertex(i,j,k,e)  ! r-edge
                  edge(j,i,k,2,e) = vertex(i,j,k,e)  ! s-edge
                  edge(k,i,j,3,e) = vertex(i,j,k,e)  ! t-edge
              enddo
          enddo
      enddo
  enddo

!   Sort edges by bounding vertices.
  do e=1,nel
    do i = 1,3
      do k=0,1
        do j=0,1
          if (edge(0,j,k,i,e) > edge(1,j,k,i,e)) then
              kswap = edge(0,j,k,i,e)
              edge(0,j,k,i,e) = edge(1,j,k,i,e)
              edge(1,j,k,i,e) = kswap
          endif
          etuple(3,1+j+k*2+(i-1)*4+12*(e-1)) = edge(0,j,k,i,e)
          etuple(4,1+j+k*2+(i-1)*4+12*(e-1)) = edge(1,j,k,i,e)
        enddo
      enddo
    enddo
  enddo

!   Assign a number (rank) to each unique edge
  m    = 4
  n    = 12*nel
  nmax = 12*lelt*nsafe  ! nsafe for crystal router factor of safety
  call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
  do i=1,nel
    do j = 1,12
      enum(j,i) = etuple(3,j+(i-1)*12)
    enddo
  enddo
  n_unique_edges = iglmax(enum,12*nel)

  n_on_edge = nx-2
  ngve      = n_unique_edges*n_on_edge
  do e=1,nel
      iedg_loc = 0
  
  !        Edges 1-4
      do k=0,1
          do j=0,1
              igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
              i0  = nx*(nx-1)*j + nx*nx*(nx-1)*k
              i0e = i0 + nxyz*(e-1)
              if (glo_num(i0e+1) < glo_num(i0e+nx)) then
                  do i=2,nx-1                                   ! std forward case
                      glo_num(i0e+i) = igv + i-1
                  enddo
              else
                  do i=2,nx-1                                   ! backward case
                      glo_num(i0e+i) = igv + 1 + n_on_edge-(i-1)
                  enddo
              endif
              iedg_loc = iedg_loc + 1
          enddo
      enddo
  
  !        Edges 5-8
      do k=0,1
          do i=0,1
              igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
              i0  = 1+(nx-1)*i + nx*nx*(nx-1)*k
              i0e = i0 + nxyz*(e-1)
              if (glo_num(i0e) < glo_num(i0e+nx*(nx-1))) then
                  do j=2,nx-1                                   ! std forward case
                      glo_num(i0e+(j-1)*nx) = igv + j-1
                  enddo
              else
                  do j=2,nx-1                                   ! backward case
                      glo_num(i0e+(j-1)*nx) = igv + 1 + n_on_edge-(j-1)
                  enddo
              endif
              iedg_loc = iedg_loc + 1
          enddo
      enddo
  
  !        Edges 9-12
      do j=0,1
          do i=0,1
              igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
              i0  = 1 + (nx-1)*i + nx*(nx-1)*j
              i0e = i0 + nxyz*(e-1)
              if (glo_num(i0e) < glo_num(i0e+nx*nx*(nx-1))) then
                  do k=2,nx-1                                   ! std forward case
                      glo_num(i0e+(k-1)*nx*nx) = igv + k-1
                  enddo
              else
                  do k=2,nx-1                                   ! backward case
                      glo_num(i0e+(k-1)*nx*nx) = igv + 1 + n_on_edge-(k-1)
                  enddo
              endif
              iedg_loc = iedg_loc + 1
          enddo
      enddo
  enddo
  ngv   = ngv + ngve

!   Asign global node numbers on the interior of each face
!   ------------------------------------------------------

!   Assign faces by 3-tuples

!   (The following variables all take the symmetric
!   notation of IFACE as arguments:)

!   ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
!                       as depicted below, e.g. ICFACE(i,2)=2,4,6,8.

!                      3+-----+4    ^ Y
!                      /  2  /|     |
!   Edge 1 extends    /     / |     |
!     from vertex   7+-----+8 +2    +----> X
!     1 to 2.        |  4  | /     /
!                    |     |/     /
!                   5+-----+6    Z
!                       3

  nfaces=ndim*2
  ncrnr =2**(ndim-1)
  do e=1,nel
      vertex_flat = reshape(vertex(:,:,:,e), (/8/))
      do ifac=1,nfaces
          do icrn=1,ncrnr
              i                  = icface(icrn,ifac)
              facet(icrn)        = vertex_flat(i)
          enddo
          call isort(facet,ind,ncrnr)
          call icopy(ftuple(3,ifac,e),facet,ncrnr-1)
      enddo
  enddo

!   Assign a number (rank) to each unique face
  m    = 5
  n    = 6*nel
  nmax = 6*lelt*nsafe  ! nsafe for crystal router factor of safety
  call gbtuple_rank(ftuple,m,n,nmax,cr_h,nid,np,ind)
  do e = 1,nel
    do i=1,6
      fnum(i,e) = ftuple(3,i,e)
    enddo
  enddo
  n_unique_faces = iglmax(fnum,6*nel)

  call dsset (nx,ny,nz)
  do e=1,nel
      do iface=1,nfaces
          i0 = skpdat(1,iface)
          i1 = skpdat(2,iface)
          is = skpdat(3,iface)
          j0 = skpdat(4,iface)
          j1 = skpdat(5,iface)
          js = skpdat(6,iface)
      
      !        On each face, count from minimum global vertex number,
      !        towards smallest adjacent vertex number.  e.g., suppose
      !        the face is defined by the following global vertex numbers:
      
      
      !                    11+--------+81
      !                      |c      d|
      !                      |        |
      !                      |        |
      !                      |a      b|
      !                    15+--------+62
      
      !        We would count from c-->a, then towards d.
      
          gvf(1) = int(glo_num(i0+nx*(j0-1)+nxyz*(e-1)))
          gvf(2) = int(glo_num(i1+nx*(j0-1)+nxyz*(e-1)))
          gvf(3) = int(glo_num(i0+nx*(j1-1)+nxyz*(e-1)))
          gvf(4) = int(glo_num(i1+nx*(j1-1)+nxyz*(e-1)))
      
          call irank(gvf,ind,4)
      
      !        ind(1) tells which element of gvf() is smallest.
      
          ifij = .FALSE. 
          if (ind(1) == 1) then
              idir =  1
              jdir =  1
              if (gvf(2) < gvf(3)) ifij = .TRUE. 
          elseif (ind(1) == 2) then
              idir = -1
              jdir =  1
              if (gvf(1) < gvf(4)) ifij = .TRUE. 
          elseif (ind(1) == 3) then
              idir =  1
              jdir = -1
              if (gvf(4) < gvf(1)) ifij = .TRUE. 
          elseif (ind(1) == 4) then
              idir = -1
              jdir = -1
              if (gvf(3) < gvf(2)) ifij = .TRUE. 
          else
            idir = 0
            jdir = 0
            write(*,*) "Wasn't supposed to get here!"
          endif
      
          if (idir < 0) then
              it=i0
              i0=i1
              i1=it
              is=-is
          endif
      
          if (jdir < 0) then
              jt=j0
              j0=j1
              j1=jt
              js=-js
          endif
      
          nxx = nx*nx
          n_on_face = (nx-2)*(ny-2)
          ngvs  = n_unique_faces*n_on_face
          ig0 = ngv + n_on_face*(fnum(iface,e)-1)
          if (ifij) then
              k=0
              l=0
              do j=j0,j1,js
                  do i=i0,i1,is
                      k=k+1
                  !              this is a serious kludge to stay on the face interior
                      if (k > nx .AND. k < nxx-nx .AND. &
                      mod(k,nx) /= 1 .AND. mod(k,nx) /= 0) then
                      !                 interior
                          l = l+1
                          glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
                      endif
                  enddo
              enddo
          else
              k=0
              l=0
              do i=i0,i1,is
                  do j=j0,j1,js
                      k=k+1
                  !              this is a serious kludge to stay on the face interior
                      if (k > nx .AND. k < nxx-nx .AND. &
                      mod(k,nx) /= 1 .AND. mod(k,nx) /= 0) then
                      !                 interior
                          l = l+1
                          glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
                      endif
                  enddo
              enddo
          endif
      enddo
  enddo
  ngv   = ngv + ngvs

!   Finally,  number interiors (only ifcenter=.true.)
!   -------------------------------------------------

  n_in_interior = (nx-2)*(ny-2)*(nz-2)
  ngvi = n_in_interior*nelgt
  if (ifcenter) then
      do e=1,nel
          ig0 = ngv + n_in_interior*(lglel(e)-1)
          l = 0
          do k=2,nz-1
              do j=2,ny-1
                  do i=2,nx-1
                      l = l+1
                      glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = ig0+l
                  enddo
              enddo
          enddo
      enddo
      ngv = ngv + ngvi
  else
      do e=1,nel
          l = 0
          do k=2,nz-1
              do j=2,ny-1
                  do i=2,nx-1
                      l = l+1
                      glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = 0
                  enddo
              enddo
          enddo
      enddo
  endif

!   Quick check on maximum #dofs:
  m    = nxyz*nelt
  ngvm = i8glmax(glo_num,m)
  ngvv = ngvv + ngve + ngvs  ! number of unique ids w/o interior
  ngvi = ngvi + ngvv         ! total number of unique ids
  if (nid == 0) write(6,1) nx,ngvv,ngvi,ngv,ngvm
  1 format('   setvert3d:',i4,4i12)

  return
end subroutine setvert3d

!-----------------------------------------------------------------------
subroutine check_p_bc(glo_num,nx,ny,nz,nel)
  use kinds, only : i8
  use size_m, only : ndim, nelt
  use input, only : ifflow, cbc
  implicit none

  integer :: nx, ny, nz, nel
  integer(i8) :: glo_num(nx,ny,nz,nel)

  integer(i8) :: gmn
  integer :: e,f,fo,ef,efo, ifld, nface, k, j, i
  integer, save :: eface0(6) = (/ 4,2,1,3,5,6 /)

  ifld = 2
  if (ifflow) ifld = 1

  nface=2*ndim
  do e=1,nelt
      do f=1,nface,2
          fo  = f+1
          ef  = eface0(f)
          efo = eface0(fo)
          if (cbc(ef,e,ifld) == 'p  ' .AND. cbc(efo,e,ifld) == 'p  ') then
              if (f == 1) then  ! r=-/+1
                  do k=1,nz
                      do j=1,ny
                          gmn = min(glo_num(1,j,k,e),glo_num(nx,j,k,e))
                          glo_num(1 ,j,k,e) = gmn
                          glo_num(nx,j,k,e) = gmn
                      enddo
                  enddo
              elseif (f == 3) then  ! s=-/+1
                  do k=1,nz
                      do i=1,nx
                          gmn = min(glo_num(i,1,k,e),glo_num(i,ny,k,e))
                          glo_num(i,1 ,k,e) = gmn
                          glo_num(i,ny,k,e) = gmn
                      enddo
                  enddo
              else
                  do j=1,ny
                      do i=1,nx
                          gmn = min(glo_num(i,j,1,e),glo_num(i,j,nz,e))
                          glo_num(i,j,1 ,e) = gmn
                          glo_num(i,j,nz,e) = gmn
                      enddo
                  enddo
              endif
          endif
      enddo
  enddo

  return
end subroutine check_p_bc

!-----------------------------------------------------------------------
