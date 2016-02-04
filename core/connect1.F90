!-----------------------------------------------------------------------
!> \brief Parallel compatible routine to find connectivity of element structure.
!! On Processor 0:
!! .Verify right-handedness of elements.
!! .Verify element-to-element reciprocity of BC's
!! .Verify correlation between E-E BC's and physical coincidence
!! .Set rotations
!! .Determine multiplicity
!! .Set up direct stiffness summation arrays.
!! All Processors:
!! .Disperse/Receive BC and MULT temporary data read from preprocessor.
subroutine setup_topo()
  use kinds, only : DP, i8
  use size_m, only : nid, ndim, nx1, ny1, nz1, nelv, nelt, nfield
  use size_m, only : lx1, ly1, lz1, lelv
  use ds, only : setupds, dssum
  use input, only : ifflow, ifmvbd, ifheat, param
  use mesh, only : vertex
  use mvgeom, only : wmult
  use parallel, only : nelgv, nelgt, gsh_fld, nelg
  use soln, only : vmult, tmult
  use tstep, only : ifield
  use zper, only : ifgtp
  implicit none

  integer(i8), allocatable :: glo_num(:)

  integer :: nxl, nyl, nzl, mfield, ncrnr, ntotv, ntott
  real(DP) :: vmltmax, ivmltmax
  real(DP), external :: glmax

  allocate(glo_num(lx1*ly1*lz1*lelv))

  if(nid == 0) write(6,*) 'setup mesh topology'

!   Initialize key arrays for Direct Stiffness SUM.

  NXL=3
  NYL=3
  NZL=1+2*(NDIM-2)

  call initds
  call dsset (nx1,ny1,nz1)
  call setedge

!=================================================
!     Establish (global) domain topology
!=================================================
!     .Generate topologically correct mesh data.
!     .Set up element centers, face centers, etc.
!     .Check  right handedness of elements.
!     .Check  element boundary conditions.
!     .Establish Element-Element rotations
!     .Construct the element to processor map and

  call genxyzl
  call setside
  if (param(75) < 1) call verify

  CALL SETCDOF
!max  IF (IFAXIS            ) CALL SETRZER
!max    IF (IFMVBD            ) CALL CBCMESH
#if 0
  IF (IFMODEL .AND. IFKEPS) CALL CBCTURB
#endif
  if (param(75) < 1) CALL CHKAXCB

!========================================================================
!     Set up element-processor mapping and establish global numbering
!========================================================================

  mfield=2
  if (ifflow) mfield=1
  if (ifmvbd) mfield=0

  ncrnr = 2**ndim

  if (nelgv == nelgt) then
      if (ifgtp) then
        write(*,*) "Oops: ifgtp"
!max          call gen_gtp_vertex    (vertex, ncrnr)
      else
          call get_vert
      endif
      call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
      gsh_fld(2)=gsh_fld(1)

  !        call gs_counter(glo_num,gsh_fld(1))

  else

  !        For conjugate heat transfer, it is assumed that fluid
  !        elements are listed both globally and locally with lower
  !        element numbers than the solid elements.
  !        We currently assume that there is at least one fluid elem.
  !        per processor.
  

      call get_vert
  !        call outmati(vertex,4,nelv,'vrtx V')
      call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)

  !        call get_vert  (vertex, ncrnr, nelgt, '.mp2')  !  LATER !
  !        call outmati(vertex,4,nelt,'vrtx T')
      call setupds(gsh_fld(2),nx1,ny1,nz1,nelt,nelgt,vertex,glo_num)

  
  !        Feb 20, 2012:  It appears that we do not need this restriction: (pff)
  !        check if there is a least one fluid element on each processor
  !        do iel = 1,nelt
  !           ieg = lglel(iel)
  !           if (ieg.le.nelgv) goto 101
  !        enddo
  !        if(nid.eq.0) write(6,*)
  !    &     'ERROR: each domain must contain at least one fluid element!'
  !        call exitt
  ! 101   continue

  endif

!max     if (ifmvbd) call setup_mesh_dssum ! Set up dssum for mesh

!========================================================================
!     Set up multiplicity and direct stiffness arrays for each IFIELD
!========================================================================

  ntotv = nx1*ny1*nz1*nelv
  ntott = nx1*ny1*nz1*nelt

  if (ifheat) then
      ifield = 2
      allocate(tmult(nx1,ny1,nz1,nelt,1))
      tmult(:,:,:,:,1) = 1._dp
      call dssum   (tmult(:,1,1,1,1))
      tmult(:,:,:,:,1) = 1._dp / tmult(:,:,:,:,1)
  endif

  if (ifflow) then
    if (ifheat) then
      vmult => tmult(:,:,:,:,1)
    else
      ifield = 1
      allocate(vmult(nx1,ny1,nz1,nelv))
      vmult = 1._dp
      call dssum   (vmult(:,1,1,1))
      vmltmax=glmax(vmult,ntotv)
      ivmltmax=vmltmax
      if (nid == 0) write(6,*) ivmltmax,' max multiplicity'
      vmult = 1._dp / vmult
    endif
  endif
  if ( .NOT. ifflow) call copy(vmult,tmult,ntott)
  if (ifmvbd)  call copy (wmult,vmult,ntott)
  do ifield=3,nfield                  ! Additional pass. scalrs.
      if (nelg(ifield) == nelgv) then
          gsh_fld(ifield) = gsh_fld(1)
          call copy (tmult(1,1,1,1,ifield-1),vmult,ntotv)
      else
          gsh_fld(ifield) = gsh_fld(2)
          call copy (tmult(1,1,1,1,ifield-1),tmult,ntott)
      endif
  enddo

  if(nid == 0) then
      write(6,*) 'done :: setup mesh topology'
      write(6,*) ' '
  endif

  return
end subroutine setup_topo

!-----------------------------------------------------------------------
!> \brief     -- Direct Stiffness Initialization Routine --
!!    Set up required data for packing data on faces of spectral cubes.
subroutine initds
  use size_m
  use topol, only : group, eface1, eface, nomlis
  implicit none

  integer :: j, idim, iface

!   Nominal ordering for direct stiffness summation of faces
  J=0
  DO IDIM=1,NDIM
      DO IFACE=1,2
          J=J+1
          NOMLIS(IFACE,IDIM)=J
      enddo
  END DO

!   Assign Ed's numbering scheme to PF's scheme.

  EFACE(1)=4
  EFACE(2)=2
  EFACE(3)=1
  EFACE(4)=3
  EFACE(5)=5
  EFACE(6)=6

!   Assign inverse of Ed's numbering scheme to PF's scheme.

  EFACE1(1)=3
  EFACE1(2)=2
  EFACE1(3)=4
  EFACE1(4)=1
  EFACE1(5)=5
  EFACE1(6)=6

!   Assign group designation to each face to determine ordering of indices.

  GROUP(1)=0
  GROUP(2)=1
  GROUP(3)=1
  GROUP(4)=0
  GROUP(5)=0
  GROUP(6)=1

  RETURN
end subroutine initds

!-----------------------------------------------------------------------
!> \brief Initialize EDGE arrays for face and edge specific tasks.
!!     .NOTE: Sevaral arrays in common are initialized via
!!            BLOCKDATA EDGEC
!!     Computed arrays:
!!     IEDGE  -  Minimal list of wire frame nodes.
!!               Used to search for all physical
!!               coincidences.
!!     IEDGEF - .Ordered list of wire frame nodes
!!               associated with faces 1 through 6.
!!              .Each of 4 sides of square frame stored
!!               individually so that rotations are
!!               readily handled.
!!              .Two types of node orderings stored -
!!               (0) is clockwise marching
!!               (1) is counter-clockwise marching
!!                   for image face.
!!     IFACE         - indicates the face number.  Two notations
!!                     are currently in use:
!!                     i) Preprocessor notation:
!!                                       +--------+     ^ S
!!                                      /        /|     |
!!                                     /    3   / |     |
!!                               4--> /        /  |     |
!!                                   +--------+ 2 +     +----> R
!!                                   |        |  /     /
!!                                   |    6   | /     /
!!                                   |        |/     /
!!                                   +--------+     T
!!                                       1
!!                    ii) Symmetric notation:
!!                                       +--------+     ^ S
!!                                      /        /|     |
!!                                     /    4   / |     |
!!                               1--> /        /  |     |
!!                                   +--------+ 2 +     +----> R
!!                                   |        |  /     /
!!                                   |    6   | /     /
!!                                   |        |/     /
!!                                   +--------+     T
!!                                       3
!!     EFACE(IFACE)  -   Given face number IFACE in symmetric notation,
!!                       returns preprocessor notation face number.
!!     EFACE1(IFACE) -   Given face number IFACE in preprocessor notation,
!!                       returns symmetric notation face number.
!!  The following variables all take the symmetric
!!  notation of IFACE as arguments:
!!     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
!!                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.
!!                        3+-----+4    ^ Y
!!                        /  2  /|     |
!!     Edge 1 extends    /     / |     |
!!       from vertex   7+-----+8 +2    +----> X
!!       1 to 2.        |  4  | /     /
!!                      |     |/     /
!!                     5+-----+6    Z
!!                         3
!!     IEDGFC(i,IFACE) -   Gives the 4 edges which border the face IFACE
!!                         Edge numbering is as follows:
!!                            Edge = 1,2,3,4     run in +r direction
!!                            Edge = 5,6,7,8     run in +s direction
!!                            Edge = 9,10,11,12  run in +t direction
!!                         Ordering of each edge is such that a monotonically
!!                         increasing sequence of vertices is associated with
!!                         the start point of a corresponding set of
!!                         monotonically increasing edge numbers, e.g.,
!!     ICEDG(i,IEDGE)  -   Gives 3 variables for determining the stride along
!!                         a given edge, IEDGE;  i=1 gives the starting vertex
!!                                               i=2 gives the stopping vertex
!!                                               i=3 gives the stride size.
subroutine setedge()
  use size_m, only : ndim
  use topol, only : iedgef, iedge, invedg, skpdat, group
  implicit none

  integer :: ITMP(3,3,3)
  INTEGER :: ORDER

  integer :: nxl, nyl, nzl, nxy, nxyz, nfaces
  integer :: i3d, i, i3, iz, i2, iy, i1, ix
  integer :: j, j1, j2, jskip2
  integer :: jf2, js2, jskip1, jf1, js1, iface, image, ii
  integer :: iedgef_flat(48)
  NXL=3
  NYL=3
  NZL=1+2*(NDIM-2)
  NXY   =NXL*NYL
  NXYZ  =NXL*NYL*NZL
  NFACES=2*NDIM

!----------------------------------------------------------------------
!     Set up edge arrays (temporary - required only for defining DS)
!----------------------------------------------------------------------
!     Fill corners - 1 through 8.

  I3D=1
  IF (NDIM == 2) I3D=0

  I=0
  DO I3=0,I3D
      IZ=1+(NZL-1)*I3
      DO I2=0,1
          IY=1+(NYL-1)*I2
          DO I1=0,1
              IX=1+(NXL-1)*I1
              I=I+1
              IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
          enddo
      enddo
  END DO

!   Fill X-direction edges.

  DO I3=0,I3D
      IZ=1+(NZL-1)*I3
      DO I2=0,1
          IY=1+(NYL-1)*I2
          DO IX=2,NXL-1
              I=I+1
              IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
          enddo
      enddo
  enddo

!   Fill Y-direction edges.

  DO I3=0,I3D
      IZ=1+(NZL-1)*I3
      DO I1=0,1
          IX=1+(NXL-1)*I1
          DO IY=2,NYL-1
              I=I+1
              IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
          enddo
      enddo
  enddo

!   Fill Z-direction edges.

  IF (NDIM == 3) THEN
      DO I2=0,1
          IY=1+(NYL-1)*I2
          DO I1=0,1
              IX=1+(NXL-1)*I1
              DO IZ=2,NZL-1
                  I=I+1
                  IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
              enddo
          enddo
      enddo
  ENDIF

  invedg = 0
  DO II=1,20
      IX=IEDGE(II)
      INVEDG(IX)=II
  END DO


!   GENERAL FACE, GENERAL ROTATION EDGE NUMBERS.
#if 0

  IF (NDIM == 3) THEN
  
  !        Pack 3-D edge numbering:
  
  !        Fill temporary array with local index numbers:
      itmp = reshape( (/(i,i=1,nxyz)/), (/nxl, nyl, nzl/)) 
  
  !        Two sets are required, the base cube and the image cube
  !        which is being summed with it.
  
      DO 1000 IMAGE=0,1
      
      !        Pack edges for each face, no rotation.
      
          DO 500 IFACE=1,NFACES
              JS1    = SKPDAT(1,IFACE)
              JF1    = SKPDAT(2,IFACE)
              JSKIP1 = SKPDAT(3,IFACE)
              JS2    = SKPDAT(4,IFACE)
              JF2    = SKPDAT(5,IFACE)
              JSKIP2 = SKPDAT(6,IFACE)
          
          !           Choose proper indexing order according to face type and image.
          
              ORDER = (-1)**(GROUP(IFACE)+IMAGE)
              IF (ORDER == 1) THEN
              
              !              Forward ordering:
              
              !            +-------------+    ^ v1
              !            | --------->| |    |
              !            | ^    2    | |    +-->
              !            | |         | |      v2
              !            | |1       3| |
              !            | |    4    V |
              !            | |<--------- |
              !            F-------------I     F is fiducial node.
              
              !                                I is location of fiducial node for
              !                                     image face.
              
              !           Load edge 1:
              
                  J=0
                  J2=JS2
                  DO 100 J1=JS1,JF1-JSKIP1,JSKIP1
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  100 END DO
              
              !           Load edge 2:
              
                  J=0
                  J1=JF1
                  DO 200 J2=JS2,JF2-JSKIP2,JSKIP2
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  200 END DO
              
              
              !           Load edge 3:
              
                  J=0
                  J2=JF2
                  DO 300 J1=JF1,JS1+JSKIP1,-JSKIP1
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  300 END DO
              
              !           Load edge 4:
              
                  J=0
                  J1=JS1
                  DO 400 J2=JF2,JS2+JSKIP2,-JSKIP2
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  400 END DO
              ! insert
              iedgef(:,:,:,image) = reshape(iedgef_flat, (/2,4,6/)) 

              ELSE
              
              !           Reverse ordering:
              
              !            +-------------+
              !            | |<--------- |       ^ v2
              !            | |    2    ^ |       |
              !            | |         | |    <--+
              !            | |3       1| |     v1
              !            | V    4    | |
              !            | --------->| |
              !            I-------------F     F is fiducial node.
              
              !                                I is location of fiducial node for
              !                                     image face.
              
              !           Load edge 1:
              
                  J=0
                  J1=JS1
                  write(*,*) IFACE, JS2, JF2, JSKIP2, IMAGE
                  DO 105 J2=JS2,JF2-JSKIP2,JSKIP2
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  105 END DO
              
              !           Load edge 2:
              
                  J=0
                  J2=JF2
                  DO 205 J1=JS1,JF1-JSKIP1,JSKIP1
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  205 END DO
              
              !           Load edge 3:
              
                  J=0
                  J1=JF1
                  DO 305 J2=JF2,JS2+JSKIP2,-JSKIP2
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  305 END DO
              
              !           Load edge 4:
              
                  J=0
                  J2=JS2
                  DO 405 J1=JF1,JS1+JSKIP1,-JSKIP1
                      J=J+1
                      IEDGEF_flat(J+(IFACE-1)*8)=j1 + 3*(j2-1) !ITMP(J1,J2,1)
                  405 END DO
              ENDIF
              ! insert
              iedgef(:,:,:,image) = reshape(iedgef_flat, (/2,4,6/)) 
          
          500 END DO
      1000 END DO
  ELSE
  
  !        Load edge information for 2-D case
  
      IEDGEF(1,1,1,0) = NXY - NXL + 1
      IEDGEF(1,2,1,0) = 1
      IEDGEF(1,1,2,0) = NXL
      IEDGEF(1,2,2,0) = NXY
      IEDGEF(1,1,3,0) = 1
      IEDGEF(1,2,3,0) = NXL
      IEDGEF(1,1,4,0) = NXY
      IEDGEF(1,2,4,0) = NXY - NXL + 1
  
      IEDGEF(1,1,1,1) = 1
      IEDGEF(1,2,1,1) = NXY - NXL + 1
      IEDGEF(1,1,2,1) = NXY
      IEDGEF(1,2,2,1) = NXL
      IEDGEF(1,1,3,1) = NXL
      IEDGEF(1,2,3,1) = 1
      IEDGEF(1,1,4,1) = NXY - NXL + 1
      IEDGEF(1,2,4,1) = NXY
  ENDIF
#endif

  RETURN
end subroutine setedge

!-----------------------------------------------------------------------
!> \brief Set up arrays IXCN,ESKIP,SKPDAT,NEDG,NOFFST for new NX,NY,NZ
subroutine dsset(nx,ny,nz)
  use topol, only : ixcn, skpdat, eskip, nedg
  implicit none

  integer :: nx, ny, nz

  INTEGER, save :: NXO = 0, NYO = 0, NZO = 0
  integer :: ic, icz, icy, icx, nxy, ied, iedm

!   Check if element surface counters are already set from last call...
  IF (NXO == NX .AND. NYO == NY .AND. NZO == NZ) RETURN

!   else, proceed....

  NXO = NX
  NYO = NY
  NZO = NZ

!   Establish corner to elemental node number mappings

  IC=0
  DO ICZ=0,1
      DO ICY=0,1
          DO ICX=0,1
          !       Supress vectorization to
          !        IF(ICX.EQ.0)DUMMY=0
          !        IF(ICX.EQ.1)DUMMY=1
          !        DUMMY2=DUMMY2+DUMMY
              IC=IC+1
              IXCN(IC)= 1 + (NX-1)*ICX + NX*(NY-1)*ICY + NX*NY*(NZ-1)*ICZ
          enddo
      enddo
  enddo

!   Assign indices for direct stiffness summation of arbitrary faces.


!   Y-Z Planes (Faces 1 and 2)

  SKPDAT(1,1)=1
  SKPDAT(2,1)=NX*(NY-1)+1
  SKPDAT(3,1)=NX
  SKPDAT(4,1)=1
  SKPDAT(5,1)=NY*(NZ-1)+1
  SKPDAT(6,1)=NY

  SKPDAT(1,2)=1             + (NX-1)
  SKPDAT(2,2)=NX*(NY-1)+1   + (NX-1)
  SKPDAT(3,2)=NX
  SKPDAT(4,2)=1
  SKPDAT(5,2)=NY*(NZ-1)+1
  SKPDAT(6,2)=NY

!   X-Z Planes (Faces 3 and 4)

  SKPDAT(1,3)=1
  SKPDAT(2,3)=NX
  SKPDAT(3,3)=1
  SKPDAT(4,3)=1
  SKPDAT(5,3)=NY*(NZ-1)+1
  SKPDAT(6,3)=NY

  SKPDAT(1,4)=1           + NX*(NY-1)
  SKPDAT(2,4)=NX          + NX*(NY-1)
  SKPDAT(3,4)=1
  SKPDAT(4,4)=1
  SKPDAT(5,4)=NY*(NZ-1)+1
  SKPDAT(6,4)=NY

!   X-Y Planes (Faces 5 and 6)

  SKPDAT(1,5)=1
  SKPDAT(2,5)=NX
  SKPDAT(3,5)=1
  SKPDAT(4,5)=1
  SKPDAT(5,5)=NY
  SKPDAT(6,5)=1

  SKPDAT(1,6)=1           + NX*NY*(NZ-1)
  SKPDAT(2,6)=NX          + NX*NY*(NZ-1)
  SKPDAT(3,6)=1
  SKPDAT(4,6)=1
  SKPDAT(5,6)=NY
  SKPDAT(6,6)=1

!   Set up skip indices for each of the 12 edges

!       Note that NXY = NX*NY even for 2-D since
!       this branch does not apply to the 2D case anyway.

!   ESKIP(*,1) = start location
!   ESKIP(*,2) = end
!   ESKIP(*,3) = stride

  NXY=NX*NY
  ESKIP( 1,1) = IXCN(1) + 1
  ESKIP( 1,2) = IXCN(2) - 1
  ESKIP( 1,3) = 1
  ESKIP( 2,1) = IXCN(3) + 1
  ESKIP( 2,2) = IXCN(4) - 1
  ESKIP( 2,3) = 1
  ESKIP( 3,1) = IXCN(5) + 1
  ESKIP( 3,2) = IXCN(6) - 1
  ESKIP( 3,3) = 1
  ESKIP( 4,1) = IXCN(7) + 1
  ESKIP( 4,2) = IXCN(8) - 1
  ESKIP( 4,3) = 1
  ESKIP( 5,1) = IXCN(1) + NX
  ESKIP( 5,2) = IXCN(3) - NX
  ESKIP( 5,3) = NX
  ESKIP( 6,1) = IXCN(2) + NX
  ESKIP( 6,2) = IXCN(4) - NX
  ESKIP( 6,3) = NX
  ESKIP( 7,1) = IXCN(5) + NX
  ESKIP( 7,2) = IXCN(7) - NX
  ESKIP( 7,3) = NX
  ESKIP( 8,1) = IXCN(6) + NX
  ESKIP( 8,2) = IXCN(8) - NX
  ESKIP( 8,3) = NX
  ESKIP( 9,1) = IXCN(1) + NXY
  ESKIP( 9,2) = IXCN(5) - NXY
  ESKIP( 9,3) = NXY
  ESKIP(10,1) = IXCN(2) + NXY
  ESKIP(10,2) = IXCN(6) - NXY
  ESKIP(10,3) = NXY
  ESKIP(11,1) = IXCN(3) + NXY
  ESKIP(11,2) = IXCN(7) - NXY
  ESKIP(11,3) = NXY
  ESKIP(12,1) = IXCN(4) + NXY
  ESKIP(12,2) = IXCN(8) - NXY
  ESKIP(12,3) = NXY

!   Load reverse direction edge arrays for reverse mappings...

  DO 20 IED=1,12
      IEDM=-IED
      ESKIP(IEDM,1) =  ESKIP(IED,2)
      ESKIP(IEDM,2) =  ESKIP(IED,1)
      ESKIP(IEDM,3) = -ESKIP(IED,3)
  20 END DO

!   Compute offset for global edge vector given current element
!   dimensions.....

!   NGSPED(ITE,ICMP) = number of global (ie, distinct) special edges
!                      of type ITE (1,2, or 3)  for field ICMP.

!                      ITE = 1 implies an "X" edge
!                      ITE = 2 implies an "Y" edge
!                      ITE = 3 implies an "Z" edge

!   Set up number of nodes along each of the 3 types of edges
!   (endpoints excluded).

  NEDG(1)=NX-2
  NEDG(2)=NY-2
  NEDG(3)=NZ-2

  RETURN
end subroutine dsset

!-----------------------------------------------------------------------
!> \brief Generate xyz coordinates
subroutine genxyzl()
  use kinds, only : DP
  use size_m, only : ndim, nelt
  use input, only : xc, yc, zc
  use scratch, only : xml, yml, zml
  implicit none

  real(DP) :: XCB(2,2,2),YCB(2,2,2),ZCB(2,2,2), H(3,3,2)
  integer :: indx(8)

  integer :: nxl, nyl, nzl, ntot3
  integer :: ix, iy, iz, ixt, iyt, izt, ie, ndim2
  NXL=3
  NYL=3
  NZL=1+2*(NDIM-2)
  NTOT3=NXL*NYL*NZL*NELT

! Preprocessor Corner notation:      Symmetric Corner notation:

!         4+-----+3    ^ s                    3+-----+4    ^ s
!         /     /|     |                      /     /|     |
!        /     / |     |                     /     / |     |
!      8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
!       |     | /     /                     |     | /     /
!       |     |/     /                      |     |/     /
!      5+-----+6    t                      5+-----+6    t

  DO 10 IX=1,NXL
      H(IX,1,1)=0.5*FLOAT(3-IX)
      H(IX,1,2)=0.5*FLOAT(IX-1)
  10 END DO
  DO 20 IY=1,NYL
      H(IY,2,1)=0.5*FLOAT(3-IY)
      H(IY,2,2)=0.5*FLOAT(IY-1)
  20 END DO
  DO 30 IZ=1,NZL
      H(IZ,3,1)=0.5*FLOAT(3-IZ)
      H(IZ,3,2)=0.5*FLOAT(IZ-1)
  30 END DO

  INDX(1)=1
  INDX(2)=2
  INDX(3)=4
  INDX(4)=3
  INDX(5)=5
  INDX(6)=6
  INDX(7)=8
  INDX(8)=7

  xml = 0._dp; yml = 0._dp; zml = 0._dp
  xcb = 0._dp; ycb = 0._dp; zcb = 0._dp

  DO 5000 IE=1,NELT
  
      NDIM2 = 2**NDIM
      xcb = reshape(xc(indx(1:ndim2), ie), (/2,2,2/))
      ycb = reshape(yc(indx(1:ndim2), ie), (/2,2,2/))
      zcb = reshape(zc(indx(1:ndim2), ie), (/2,2,2/))
 
  !        Map R-S-T space into physical X-Y-Z space.
  
      DO IZT=1,ndim-1
          DO IYT=1,2
              DO IXT=1,2
              
                  DO IZ=1,NZL
                      DO IY=1,NYL
                          DO IX=1,NXL
                              XML(IX,IY,IZ,IE)=XML(IX,IY,IZ,IE)+ &
                              H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*XCB(IXT,IYT,IZT)
                              YML(IX,IY,IZ,IE)=YML(IX,IY,IZ,IE)+ &
                              H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*YCB(IXT,IYT,IZT)
                              ZML(IX,IY,IZ,IE)=ZML(IX,IY,IZ,IE)+ &
                              H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*ZCB(IXT,IYT,IZT)
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
  
  5000 END DO
  RETURN
end subroutine genxyzl

!-----------------------------------------------------------------------
!> \brief Verify right-handedness of elements.
!! .Verify element-to-element reciprocity of BC's
!! .Verify correlation between E-E BC's and physical coincidence
subroutine verify
  implicit none 
  call verrhe
  return
end subroutine verify

!-----------------------------------------------------------------------
subroutine setside
  use kinds, only : DP
  use size_m, only : nelt, ndim
  use input, only : if3d, xc, yc, zc
  use scratch, only : xyz, side
  use topol, only : indx, icface
  implicit none

  integer :: ie, j, ivtx, nfaces, ncrnr, icrn, ifac, icr1, ivt1, idim
  integer, external :: mod1
  real(DP) :: avwght

!   SIDE(i,IFACE,IE) -  Physical (xyz) location of element side midpoint.
!                       i=1,2,3 gives x,y,z value, respectively.
!                       i=4  gives average dimension of face for setting
!                            tolerances.
  INDX(1)=1
  INDX(2)=2
  INDX(3)=4
  INDX(4)=3
  INDX(5)=5
  INDX(6)=6
  INDX(7)=8
  INDX(8)=7

!   Flip vertex array structure

!   write(6,*) nelv,nelt,if3d
  xyz = 0._dp
  if (if3d) then
      do ie=1,nelt
          do j=1,8
              ivtx = indx(j)
              xyz(1,ivtx,ie) = xc(j,ie)
              xyz(2,ivtx,ie) = yc(j,ie)
              xyz(3,ivtx,ie) = zc(j,ie)
          !           write(6,1) ie,j,ivtx,xc(j,ie),yc(j,ie),zc(j,ie),' xcz'
          !           write(6,1) ie,j,ivtx,(xyz(k,ivtx,ie),k=1,3),' vtx'
          !   1       format(3i5,1p3e12.4,a4)
          enddo
      enddo
  else
      do ie=1,nelt
          do j=1,4
              ivtx = indx(j)
              xyz(1,ivtx,ie) = xc(j,ie)
              xyz(2,ivtx,ie) = yc(j,ie)
              xyz(3,ivtx,ie) = 0.0
          enddo
      enddo
  endif

!   Compute location of center and "diameter" of each element side.

  NFACES=NDIM*2
  NCRNR =2**(NDIM-1)
  side = 0._dp
  DO ICRN=1,NCRNR
      DO IFAC=1,NFACES
          IVTX = ICFACE(ICRN,IFAC)
          ICR1 = NCRNR+(ICRN-1)
          ICR1 = MOD1(ICR1,NCRNR)
          IVT1 = ICFACE(ICR1,IFAC)
          DO 400 IE=1,NELT
              DO 300 IDIM=1,NDIM
                  SIDE(IDIM,IFAC,IE)=SIDE(IDIM,IFAC,IE)+XYZ(IDIM,IVTX,IE)
                  SIDE(   4,IFAC,IE)=SIDE(   4,IFAC,IE)+ &
                  ( XYZ(IDIM,IVTX,IE)-XYZ(IDIM,IVT1,IE) )**2
              300 END DO
              SIDE(4,IFAC,IE)=SQRT( SIDE(4,IFAC,IE) )
          400 END DO
      enddo
  END DO
  AVWGHT=1.0/FLOAT(NCRNR)
  side = side*avwght

!   call exitt
  RETURN
end subroutine setside

!-----------------------------------------------------------------------
!> \brief Verify right-handedness of given elements.
!!     8 Mar 1989 21:58:26   PFF
subroutine verrhe()
  use kinds, only : DP
  use size_m, only : nid, nelt
  use input, only : if3d
  use parallel, only : lglel, nelgt
  use scratch, only : xyz
  implicit none

  integer :: ie, ieg
  real(DP) :: v1, v2, v3, v4, v5, v6, v7, v8
  LOGICAL :: IFCSTT
  real(DP), external :: volum0

  IFCSTT= .TRUE. 
  IF ( .NOT. IF3D) THEN
    write(*,*) "Oops: not if3d"
#if 0
      DO 1000 IE=1,NELT
      
      !        CRSS2D(A,B,O) = (A-O) X (B-O)
      
          C1=CRSS2D(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,1,IE))
          C2=CRSS2D(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,2,IE))
          C3=CRSS2D(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,3,IE))
          C4=CRSS2D(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,4,IE))
      
          IF (C1 <= 0.0 .OR. C2 <= 0.0 .OR. &
          C3 <= 0.0 .OR. C4 <= 0.0 ) THEN
          
              ieg=lglel(ie)
              WRITE(6,800) IEG,C1,C2,C3,C4
              call exitt
              800 FORMAT(/,2X,'WARNINGa: Detected non-right-handed element.', &
              /,2X,'Number',I8,'  C1-4:',4E12.4)
              IFCSTT= .FALSE. 
          !           CALL QUERY(IFYES,'Proceed                                 ')
          !           IF (.NOT.IFYES) GOTO 9000
          ENDIF
      1000 END DO
#endif 
  !     Else 3-D:
  ELSE
      DO 2000 IE=1,NELT
      
      !        VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
      
          V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
          V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
          V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
          V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
          V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
          V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
          V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
          V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
      
          IF (V1 <= 0.0 .OR. V2 <= 0.0 .OR. &
          V3 <= 0.0 .OR. V4 <= 0.0 .OR. &
          V5 <= 0.0 .OR. V6 <= 0.0 .OR. &
          V7 <= 0.0 .OR. V8 <= 0.0    ) THEN
          
              ieg=lglel(ie)
              WRITE(6,1800) IEG,V1,V2,V3,V4,V5,V6,V7,V8
              call exitt
              1800 FORMAT(/,2X,'WARNINGb: Detected non-right-handed element.', &
              /,2X,'Number',I8,'  V1-8:',4E12.4 &
              /,2X,'      ',4X,'       ',4E12.4)
              IFCSTT= .FALSE. 
          ENDIF
      2000 END DO
  ENDIF

!   Print out results from right-handed check

  IF ( .NOT. IFCSTT) WRITE(6,2001)

!   Check consistency accross all processors.

  CALL GLLOG(IFCSTT, .FALSE. )

  IF ( .NOT. IFCSTT) THEN
      IF (NID == 0) WRITE(6,2003) NELGT
      call exitt
  ELSE
      IF (NID == 0) WRITE(6,2002) NELGT
  ENDIF

  2001 FORMAT(//,'  Elemental geometry not right-handed, ABORTING' &
  ,' in routine VERRHE.')
  2002 FORMAT('   Right-handed check complete for',I8,' elements. OK.')
  2003 FORMAT('   Right-handed check failed for',I8,' elements.' &
  ,'   Exiting in routine VERRHE.')
  RETURN
end subroutine verrhe

!-----------------------------------------------------------------------
!> \brief Given four points in R , (P1,P2,P3,P0), VOLUM0 returns
!! the volume enclosed by the parallelagram defined by the
!! vectors { (P1-P0),(P2-P0),(P3-P0) }.  This routine has
!! the nice feature that if the 3 vectors so defined are
!! not right-handed then the volume returned is negative.
real(DP) FUNCTION VOLUM0(P1,P2,P3,P0)
  use kinds, only : DP
  implicit none
  REAL(DP) :: P1(3),P2(3),P3(3),P0(3)
  real(DP) :: u1, u2, u3, v1, v2, v3, w1, w2, w3, cross1, cross2, cross3
  U1=P1(1)-P0(1)
  U2=P1(2)-P0(2)
  U3=P1(3)-P0(3)

  V1=P2(1)-P0(1)
  V2=P2(2)-P0(2)
  V3=P2(3)-P0(3)

  W1=P3(1)-P0(1)
  W2=P3(2)-P0(2)
  W3=P3(3)-P0(3)

  CROSS1 = U2*V3-U3*V2
  CROSS2 = U3*V1-U1*V3
  CROSS3 = U1*V2-U2*V1

  VOLUM0  = W1*CROSS1 + W2*CROSS2 + W3*CROSS3
           
  RETURN
END FUNCTION VOLUM0

!-----------------------------------------------------------------------
!> \brief ifcase in preprocessor notation
subroutine facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
  implicit none
  integer :: kx1, ky1, kz1, kx2, ky2, kz2, nx, ny, nz, iface

  KX1=1
  KY1=1
  KZ1=1
  KX2=NX
  KY2=NY
  KZ2=NZ
  IF (IFACE == 1) KY2=1
  IF (IFACE == 2) KX1=NX
  IF (IFACE == 3) KY1=NY
  IF (IFACE == 4) KX2=1
  IF (IFACE == 5) KZ2=1
  IF (IFACE == 6) KZ1=NZ
  return
end subroutine facind

!-----------------------------------------------------------------------
!> \brief Assign the value VAL to face(IFACE,IE) of array A.
!! IFACE is the input in the pre-processor ordering scheme.
subroutine facev(a,ie,iface,val,nx,ny,nz)
  use kinds, only : DP
  use size_m, only : lelt
  implicit none
  integer :: ie, iface, nx, ny, nz
  real(DP) :: a(NX,NY,NZ,LELT), val

  integer :: kx1, ky1, kz1, kx2, ky2, kz2
  integer :: ix, iy, iz

  CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
  DO IZ=KZ1,KZ2
      DO IY=KY1,KY2
          DO IX=KX1,KX2
              A(IX,IY,IZ,IE)=VAL
          enddo
      enddo
  enddo
  RETURN
end subroutine facev
!-----------------------------------------------------------------------
