#!/bin/bash
#
#SBATCH --job-name=srvo3-crpa
#SBATCH --partition=debug
#SBATCH --constraint=mc
#SBATCH --account=mr6
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
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

cp KPOINTS.dos KPOINTS
cp INCAR.DFT INCAR
$VASP
cat INCAR OUTCAR > OUTCAR.DFT

cp KPOINTS.dos KPOINTS
cp INCAR.EXACT INCAR
$VASP
cat INCAR OUTCAR > OUTCAR.EXACT

cp INCAR.CRPA INCAR
$VASP
cat INCAR OUTCAR > OUTCAR.CRPA_target.static

