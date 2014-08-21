    subroutine makeq_aux

    use size_m
    use input
    use tstep

    logical ::  ifturb,if_conv_std

    if_conv_std = .TRUE. 
    if (ifmhd .AND. ifaxis) if_conv_std = .FALSE. ! conv. treated in induct.f

    if(ifcvode .AND. ifield == 2) call setprop

#if 0
    call whatfld (ifturb)
    if (ifturb) call maketq ! zero bq
#endif

    if ( .NOT. ifturb .AND. if_conv_std)  call makeuq !zero bq
        
    return
    end subroutine makeq_aux
