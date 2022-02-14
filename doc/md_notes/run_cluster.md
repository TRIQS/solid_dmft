# Running solid_dmft on a cluster

## Running on CSCS daint

in some directories one can also find example job files to run everything on
daint. Note, that one has to first load the desired docker images with sarus
on daint: https://user.cscs.ch/tools/containers/sarus/, see the README.md in the `/Docker` folder.

## one shot job on daint

one shot is quite straight forward. Just get the newest version of these
scripts, go to a working directory and then create job file that looks like
this:
```
#!/bin/bash
#SBATCH --job-name="svo-test"
#SBATCH --time=1:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --account=eth3
#SBATCH --ntasks-per-core=1
#SBATCH --constraint=mc
#SBATCH --partition=normal
#SBATCH --output=out.%j
#SBATCH --error=err.%j

#======START=====

srun sarus run --mpi --mount=type=bind,source=$SCRATCH,destination=$SCRATCH --mount=type=bind,source=/apps,destination=/apps load/library/triqs-2.1-vasp bash -c "cd $PWD ; python /apps/ethz/eth3/dmatl-theory-git/solid_dmft/solid_dmft.py"
```
thats it. This line automatically runs the docker image and executes the
`solid_dmft.py` script. Unfortunately the new sarus container enginge does not mounts automatically user directories. Therefore, one needs to specify with `--mount` to mount the scratch and apps folder manually. Then, one executes in the container bash to first go into the current dir and then executes python and the dmft script.

## CSC calculations on daint

CSC calculations need the parameter `csc = True` and the mandatory parameters from the group `dft`.
Then, solid_dmft automatically starts VASP on as many cores as specified.
Note that VASP runs on cores that are already used by solid_dmft.
This minimizes the time that cores are idle while not harming the performance because these two processes are never active at the same time.

For the latest version in the Dockerfile_MPICH, we need the sarus version >= 1.3.2, which can be loaded from the daint modules as `sarus/1.3.2` but isn't the default version.
The reason for this is that only from this sarus version on, having more than one version of libgfortran in the docker image is supported, which comes from Vasp requiring the use of gfortran7 and everything else using gfortran9.

A slurm job script should look like this:
```
#!/bin/bash
#SBATCH --job-name="svo-csc-test"
#SBATCH --time=4:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=36
#SBATCH --account=eth3
#SBATCH --ntasks-per-core=1
#SBATCH --constraint=mc
#SBATCH --partition=normal
#SBATCH --output=out.%j
#SBATCH --error=err.%j

# path to solid_dmft.py script
SCRIPTDIR=/apps/ethz/eth3/dmatl-theory-git/solid_dmft/solid_dmft.py
# Sarus image that is utilized
IMAGE=load/library/triqs_mpich

srun --cpu-bind=none --mpi=pmi2 sarus run --mount=type=bind,source=/apps,destination=/apps --mount=type=bind,source=$SCRATCH,destination=$SCRATCH --workdir=$PWD $IMAGE python3 $SCRIPTDIR
```
Note that here the mpi option is given to the `srun` command and not the sarus command, as for one-shot calculations.
This is important for the python to be able to start VASP.

In general I found 1 node for Vasp is in most cases enough, which means that we set `n_cores` in the dmft\_config.ini to 36 here.
Using more than one node results in a lot of MPI communication, which in turn slows down the calculation significantly.
For a 80 atom unit cell 2 nodes are useful, but for a 20 atom unit cell not at all!
