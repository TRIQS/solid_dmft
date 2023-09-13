69a70
>       USE mlwf
656a658,659
>       flush(17)
>       print *, ' '
706,708c709,717
<          CALL LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IO%IU6, IO%IU0, W)
<          CALL LPRJ_LDApU(IO, W)   ! write a LDA+U GAMMA file
<          CALL LPRJ_DEALLOC
---
>          CALL LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IO%IU6, IO%IU0, W, EFERMI)
>          ! CALL LPRJ_LDApU(IO, W)   ! write a LDA+U GAMMA file only for debug
>          CALL LPRJ_DEALLOC
> 
>          IF (WANNIER90()) THEN
>               CALL KPAR_SYNC_ALL(WDES,W)
>               CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
>               CALL WRITE_EIGENVAL(WDES, W, IO%IU6)
>          ENDIF
