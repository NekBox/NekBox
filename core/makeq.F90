!-----------------------------------------------------------------------
    subroutine makeq

!     Generate forcing function for the solution of a passive scalar.
!     !! NOTE: Do not change the content of the array BQ until the current

    use size_m
    use dealias
  use dxyz
  use eigen
  use esolv
  use geom
  use input
  use ixyz
  use mass
  use mvgeom
  use parallel
  use soln
  use steady
  use topol
  use tstep
  use turbo
  use wz_m
  use wzf

    logical ::  if_conv_std
    common /SCRUZ/ w1(lx1,ly1,lz1,lelt)

    if_conv_std = .TRUE. 
    if (ifmhd .AND. ifaxis) if_conv_std = .FALSE. ! conv. treated in induct.f

    call makeq_aux ! nekuq, etc.

    if (ifadvc(ifield) .AND. .NOT. ifchar .AND. if_conv_std)  call convab

    if (iftran) then
        if(ifcvode) then
            ntot = nx1*ny1*nz1*nelfld(ifield)
            call wlaplacian(w1,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield), &
            ifield)
            call add2(bq(1,1,1,1,ifield-1),w1,ntot)
            if(iftmsh(ifield)) then
                call dssum(bq,nx1,ny1,nz1)
                call col2(bq,bintm1,ntot)
                call col2(bq,bm1,ntot)
            endif
        else
            if (ifmvbd) then       ! ifchar is false
#if 0
                call admesht
                call makeabq
                call makebdq
#endif
            elseif (ifchar .AND. ifadvc(ifield)) then
                call makeabq
                call convch
            else
                call makeabq
                call makebdq
            endif
        endif
    endif

    return
    end subroutine makeq
