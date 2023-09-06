2800c2800
<       SUBROUTINE LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IU6,IU0,W)
---
>       SUBROUTINE LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IU6,IU0,W,EFERMI)
2807a2808
>       REAL(q) :: EFERMI
2939c2940
<          WRITE(99,'(4I6,"  # of spin, # of k-points, # of bands, # of proj" )') NS,NK,NB,NF
---
>          WRITE(99,'(4I6,F12.7,"  # of spin, # of k-points, # of bands, # of proj, Efermi" )') W%WDES%NCDIJ,NK,NB,NF,EFERMI
