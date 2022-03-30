69a70
>       USE mlwf
159a161,162
> ! existence of file vasp.suppress_projs, only needed in ICHARG=5 mode
>       LOGICAL :: LSUPPRESS_PROJS_EXISTS
656a660,661
>       flush(17)
>       print *, ' '
697,704c702,720
<          IF (IO%LORBIT==14) CALL SPHPRO_FAST( &
<               GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
<               INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,INNER_STEP, IO%LORBIT,PAR, &
<               EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
<          CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,LPRJ_functions,LPRJ_COVL,T_INFO,INFO,IO)
<          CALL LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IO%IU6, IO%IU0, W)
<          CALL LPRJ_LDApU(IO, W)   ! write a LDA+U GAMMA file
<          CALL LPRJ_DEALLOC
---
>          INQUIRE(FILE='vasp.suppress_projs',&
>                  EXIST=LSUPPRESS_PROJS_EXISTS)
> 
>          IF (.NOT. LSUPPRESS_PROJS_EXISTS) THEN
>             IF (IO%LORBIT==14) CALL SPHPRO_FAST( &
>                  GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
>                  INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,INNER_STEP, IO%LORBIT,PAR, &
>                  EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
>             CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,LPRJ_functions,LPRJ_COVL,T_INFO,INFO,IO)
>             CALL LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IO%IU6, IO%IU0, W, EFERMI)
>             ! CALL LPRJ_LDApU(IO, W)   ! write a LDA+U GAMMA file
>             CALL LPRJ_DEALLOC
> 
>             IF (WANNIER90()) THEN
>                 CALL KPAR_SYNC_ALL(WDES,W)
>                 CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
>                 CALL WRITE_EIGENVAL(WDES, W, IO%IU6)
>             ENDIF
>          ENDIF
