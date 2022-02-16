#!/bin/sh
awk '/varray name="epsilon_diag"/ { on=1 } on == 1 && /<v>/ { print $2, $3+0 } /<\/varray>/ { on = 0 } ' vasprun.xml | sort -n -k 1  > screen_G.dat


