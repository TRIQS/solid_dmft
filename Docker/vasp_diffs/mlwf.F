476a477
>    REAL(q) :: FERMILEVEL
479a481
>    FERMILEVEL = W%EFERMI(1)
491a494,502
> 
>    ! The Vasp-W90 interface does not want to be executed multiple times
>    ! We therefore need to reset a few things manually in ICHARG=5 mode
>    IF (MOD(INFO%ICHARG,10)==5) THEN
>        CALL MLWF_FREE(P_MLWF)
>        P_MLWF%PROJ_MODE = PROJ_UNKNOWN
>        IF (NUM_WANN == 0) NUM_WANN = P_MLWF%NUM_WANN
>    ENDIF
> 
563c574
<       CALL LPRJ_WRITE(MLWF%LPRJ_functions,MLWF%LPRJ_COVL,IO%IU6,IO%IU0,W)
---
>       CALL LPRJ_WRITE(MLWF%LPRJ_functions,MLWF%LPRJ_COVL,IO%IU6,IO%IU0,W, FERMILEVEL)
576,579d586
<    ! or the one-shot method
<    ! or unitary matrix contrained optimization (umco)
<    ELSE
<       CALL MLWF_TRAFO_RUN(P_MLWF,W,P,CQIJ,LATT_CUR,T_INFO,INFO,IO)
