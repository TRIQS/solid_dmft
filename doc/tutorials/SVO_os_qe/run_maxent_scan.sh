#!/bin/bash
export WORKDIR=`pwd`

########################################################
##
##CHANGE PATH TO MAXENT TO THE RIGHT FOLDER!
##
########################################################
export PATH_TO_MAXENT="${WORKDIR}/../solid_dmft/python/solid_dmft/postprocessing/"



echo "Bash version ${BASH_VERSION}..."
for J in 0.0; do
  for U in 3.5 4 4.5 5 5.5 6 6.5; do
    echo "Moving to directory for U=$U J=$J eV"
    cd  "./J$J/U$U/"
    echo "Running the maxent script"
  
    mpirun -n 30 python3 "${PATH_TO_MAXENT}maxent_gf_imp.py" "./out/svo.h5"
  
    cd $WORKDIR
  done
done
