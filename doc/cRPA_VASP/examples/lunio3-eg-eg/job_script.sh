#!/bin/bash
#
#SBATCH --job-name=lu-pbnm-rel-eg-eg-high
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --account=mr6
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --exclusive
#======START=====

module load slurm

echo "Current directory is $PWD"
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_NNODES nodes ($SLURM_NODELIST), $SLURM_NTASKS_PER_NODE tasks/node"

export OMP_NUM_THREADS=1
VASP="srun /apps/monch/p504/daint/vasp/5.4.3/mc/bin/vasp_std"

#cp INCAR.DFT INCAR
#$VASP
#cat INCAR OUTCAR > OUTCAR.DFT

#cp INCAR.EXACT INCAR
#$VASP
#cat INCAR OUTCAR > OUTCAR.EXACT

cp INCAR.CRPA INCAR
$VASP
cat INCAR OUTCAR > OUTCAR.CRPA_target.static

