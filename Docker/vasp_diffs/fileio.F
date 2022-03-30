1726c1726,1728
<         WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
---
>         IF (IO%IU0>=0) THEN
>           WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
>         ENDIF
