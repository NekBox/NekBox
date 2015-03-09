!> \file makeq.F90 \copybrief makeq

!-----------------------------------------------------------------------
!> \brief Generate forcing function for the solution of a passive scalar.
!! !! NOTE: Do not change the content of the array BQ until the current
subroutine makeq()
  use kinds, only : DP
  use input, only : ifmhd, ifaxis, ifadvc, ifchar, iftran, ifcvode
  use input, only : ifmvbd
  use tstep, only : ifield
  use ctimer, only : tmakeq, nmakeq, dnekclock
  implicit none

  logical ::  if_conv_std
  real(DP) :: etime

  nmakeq = nmakeq + 1
  etime = dnekclock()

  if_conv_std = .TRUE. 
  if (ifmhd .AND. ifaxis) if_conv_std = .FALSE. ! conv. treated in induct.f

  call makeq_aux ! nekuq, etc.

  etime = etime - dnekclock()
  if (ifadvc(ifield) .AND. .NOT. ifchar .AND. if_conv_std)  call convab
  etime = etime + dnekclock()

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

  tmakeq = tmakeq + (dnekclock() - etime)

  return
end subroutine makeq
