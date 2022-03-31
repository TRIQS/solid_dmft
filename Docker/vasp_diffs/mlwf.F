476a477
>    REAL(q) :: FERMILEVEL
479c480,481
< 
---
>     
>    FERMILEVEL = W%EFERMI(1)
563c565
<       CALL LPRJ_WRITE(MLWF%LPRJ_functions,MLWF%LPRJ_COVL,IO%IU6,IO%IU0,W)
---
>       CALL LPRJ_WRITE(MLWF%LPRJ_functions,MLWF%LPRJ_COVL,IO%IU6,IO%IU0,W, FERMILEVEL)
578,579c580,581
<    ELSE
<       CALL MLWF_TRAFO_RUN(P_MLWF,W,P,CQIJ,LATT_CUR,T_INFO,INFO,IO)
---
>    !ELSE
>    !   CALL MLWF_TRAFO_RUN(P_MLWF,W,P,CQIJ,LATT_CUR,T_INFO,INFO,IO)
1125c1127
<       MLWF%LPRJ_functions(IS)%spin_qaxis = proj_s_qaxisx(3,IS)
---
>       MLWF%LPRJ_functions(IS)%spin_qaxis = proj_s_qaxisx(:,IS)
