1755c1755,1757
<         WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
---
>         IF (IO%IU0>=0) THEN
>           WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
>         ENDIF
