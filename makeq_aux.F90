!> \file makeq_aux.F90 \copybrief makeq_aux

subroutine makeq_aux
  use input, only : ifmhd, ifaxis, ifcvode
  use tstep, only : ifield
  implicit none

  logical ::  ifturb,if_conv_std

  if_conv_std = .TRUE. 
  if (ifmhd .AND. ifaxis) if_conv_std = .FALSE. ! conv. treated in induct.f

  if(ifcvode .AND. ifield == 2) call setprop

  ifturb = .false.
#if 0
  call whatfld (ifturb)
  if (ifturb) call maketq ! zero bq
#endif

  if ( .NOT. ifturb .AND. if_conv_std)  call makeuq !zero bq
        
  return
end subroutine makeq_aux
