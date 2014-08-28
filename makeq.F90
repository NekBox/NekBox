!-----------------------------------------------------------------------
!> \brief Generate forcing function for the solution of a passive scalar.
!! !! NOTE: Do not change the content of the array BQ until the current
subroutine makeq()
  use kinds, only : DP
  use size_m, only : nx1, ny1, nz1
  use size_m, only : lx1, ly1, lz1, lelt
  use input, only : ifmhd, ifaxis, ifadvc, ifchar, iftran, ifcvode, iftmsh
  use input, only : ifmvbd
  use mass, only : bintm1, bm1
  use soln, only : t, vdiff, bq
  use tstep, only : ifield, nelfld
  implicit none

  logical ::  if_conv_std
  integer :: ntot

  if_conv_std = .TRUE. 
  if (ifmhd .AND. ifaxis) if_conv_std = .FALSE. ! conv. treated in induct.f

  call makeq_aux ! nekuq, etc.

  if (ifadvc(ifield) .AND. .NOT. ifchar .AND. if_conv_std)  call convab

  if (iftran) then
      if(ifcvode) then
        write(*,*) "Oops: ifcvode"
#if 0
          ntot = nx1*ny1*nz1*nelfld(ifield)
          call wlaplacian(w1,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield), &
          ifield)
          call add2(bq(1,1,1,1,ifield-1),w1,ntot)
          if(iftmsh(ifield)) then
              call dssum(bq,nx1,ny1,nz1)
              call col2(bq,bintm1,ntot)
              call col2(bq,bm1,ntot)
          endif
#endif
      else
          if (ifmvbd) then       ! ifchar is false
          write(*,*) "Oops: ifmvbd"
#if 0
              call admesht
              call makeabq
              call makebdq
#endif
          elseif (ifchar .AND. ifadvc(ifield)) then
          write(*,*) "Oops: ifchar and ifadvc"
#if 0
              call makeabq
              call convch
#endif
          else
              call makeabq
              call makebdq
          endif
      endif
  endif

  return
end subroutine makeq
