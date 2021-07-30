![logo_soliDMFT](doc/logos/logo_solid_dmft.png)

This program allows to perform DFT+DMFT ''one-shot'' and CSC
calculations from h5 archives or VASP input files for multiband systems using
the TRIQS package, in combination with the CThyb solver and SumkDFT from
DFT-tools. Runs with triqs 3.x.x

Copyright (C) 2021: A. Hampel, M. Merkel, and S. Beck
(see LICENSE.txt for details)


## Source code files and their use

- __bin/solid_dmft:__ main file that runs the DMFT calculation and starts a CSC flow 
- __main.py:__ main function that invokes `csc_flow_control` or a one shot 
  calculation directly by invoking `dmft_cycle` on a given h5 archives
- __read_config.py:__ contains the functions to read the dmft config file. Take a
  look in `read_config_doc.md` for a detailed list of parameters
- __dmft_cycle.py:__ contains the `dmft_cycle` function that run a predefined
  number of DMFT iterations
- __csc_flow.py:__ contains the `csc_flow_control` function to steer the CSC
  calculation and call then ones per DFT+DMFT cycle the `dmft_cycle` function
- __/dmft_tools__
    - __convergence.py:__ contains functions to calculate convergence properties and
      defining definitions for Gf differences. Results will be stored in a dictionary:
      `convergence_obs`
    - __formatter.py:__
    - __interaction_hamiltonian.py:__
    - __legendre_filter.py:__
    - __manipulate_chemical_potential.py:__
    - __observables.py:__ contains all functions necessary to calculate and write the
      observables during the run, which are stored in a general dictionary: `observables`
    - __solver.py:__
- __/dft_managers.py__
    - __vasp_manager.py:__ contains all functions to communicate with vasp,
    which in a CSC calculation in running the whole time (in the background,
    if not needed)
    - __qe_manager.py:__ contains all functions to start quantum espresso
- __/util:__ non-essential scripts that are independent of the rest of the code
- __/postprocessing:__ different scripts to perform postprocessing steps like analytical 
  continuations with triqs's MaxEnt


## Getting started

To start take a look in the `example` directory. There one finds several
examples to run. Best start with the svo-one-shot example. The
`dmft_config.ini` file contains the configuration for the DMFT run, which is
explained in the read\_config method in the main script. The `svo.h5` is the DMFT
input data, which is obtained from projection on localized Wannier functions
(see folder `svo-one-shot/converter`).

Furthermore, there is a `read_config_doc.md` file containing the docstrings from
the main script in a readable format. If one wishes to do CSC calculations the
docker container must contain also a installed VASP version >5.4.4 that
understands the ICHARG=5 flag.

To test triqs, you can use the official image from https://hub.docker.com/r/flatironinstitute/triqs/.
For a full installation including Vasp run our own Docker image (see folder `/docker`).

### CSC calculations locally

Here one needs a special docker image with vasp included. This can be done by
building the Dockerfile in `/Docker`.
Then start this docker image as done above and go to the directory with all
necessary input files (start with `svo-csc` example). You need a pre-converged
CHGCAR and preferably a WAVECAR, a set of INCAR, POSCAR, KPOINTS and POTCAR
files, the PLO cfg file `plo.cfg` and the usual DMFT input file
`dmft_config.ini`, which specifies the number of ranks for the DFT code and the DFT code executable in the `[dft]` section.

The whole machinery is started by calling `solid_dmft.py` as for one-shot calculations. Importantly the flag `csc = True` has to be set in the general section in the config file. Then:
```
mpirun -n 12 /work/solid_dmft.py
```
The programm will then run the `csc_flow_control` routine, which starts VASP accordingly by spawning a new child process. After VASP is finished it will run the converter, run the dmft_cycle, and then VASP again until the given
limit of DMFT iterations is reached. This should also work on most HPC systems (tested on slurm with OpenMPI), as the the child mpirun call is performed without the slurm environment variables. This tricks slrum into starting more ranks than it has available. Note, that maybe a slight adaption of the environment variables is needed to make sure VASP is found on the second node. The variables are stored `args` in the function `start_vasp_from_master_node` of the module `csc_flow.py`

One remark regarding the number of iterations per DFT cycle. Since VASP uses a
block Davidson scheme for minimizing the energy functional not all eigenvalues
of the Hamiltonian are updated simultaneously therefore one has to make several
iterations before the changes from DMFT in the charge density are completely
considered. The default value are __6__ DFT iterations, which is very
conservative, and can be changed by changing the config parameter `n_iter` in the `[dft]` section. In general one should use `IALGO=90` in VASP, which performs an exact diagonalization rather than a partial diagonalization scheme, but this is very slow for larger systems.


## Running on CSCS daint

in some directories one can also find example job files to run everything on
daint. Note, that one has to first load the desired docker images with sarus
on daint: https://user.cscs.ch/tools/containers/sarus/, see the README.md in the `/Docker` folder.

### one shot job on daint

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

### CSC calculations on daint

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




## orbital order in the W90 converter

Some interaction Hamiltonians are sensitive to the order of orbitals (i.e. density-density or Slater Hamiltonian), others are invariant under rotations in orbital space (i.e. the Kanamori Hamiltonian).
For the former class and W90-based DMFT calculations, we need to be careful because the order of W90 (z^2, xz, yz, x^2-y^2, xy) is different from the order expected by TRIQS (xy, yz, z^2, xz, x^2-y^2).
Therefore, we need to specify the order of orbitals in the projections block (example for Pbnm or P21/n cell, full d shell):
```
begin projections
# site 0
f=0.5,0.0,0.0:dxy
f=0.5,0.0,0.0:dyz
f=0.5,0.0,0.0:dz2
f=0.5,0.0,0.0:dxz
f=0.5,0.0,0.0:dx2-y2
# site 1
f=0.5,0.0,0.5:dxy
f=0.5,0.0,0.5:dyz
f=0.5,0.0,0.5:dz2
f=0.5,0.0,0.5:dxz
f=0.5,0.0,0.5:dx2-y2
# site 2
f=0.0,0.5,0.0:dxy
f=0.0,0.5,0.0:dyz
f=0.0,0.5,0.0:dz2
f=0.0,0.5,0.0:dxz
f=0.0,0.5,0.0:dx2-y2
# site 3
f=0.0,0.5,0.5:dxy
f=0.0,0.5,0.5:dyz
f=0.0,0.5,0.5:dz2
f=0.0,0.5,0.5:dxz
f=0.0,0.5,0.5:dx2-y2
end projections
```
Warning: simply using `Fe:dxy,dyz,dz2,dxz,dx2-y2` does not work, VASP/W90 brings the d orbitals back to W90 standard order.

The 45-degree rotation for the sqrt2 x sqrt2 x 2 cell can be ignored because the interaction Hamiltonian is invariant under swapping x^2-y^2 and xy.



