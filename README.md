![logo_soliDMFT](doc/logo_solid_dmft.png)

This program allows to perform DFT+DMFT ''one-shot'' and CSC
calculations from h5 archives or VASP input files for multiband systems using
the TRIQS package, in combination with the CThyb solver and SumkDFT from
DFT-tools. Runs with triqs 3.x.x

Copyright (C) 2021: A. Hampel, M. Merkel, and S. Beck
(see LICENSE.txt for details)


## Source code files and their use

- __solid_dmft:__ main file that runs the DMFT calculation and starts a CSC flow 
- - __main.py:__ main function that invokes `csc_flow_control` or a one shot 
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

## LOCPROJ bug for individual projections:

Example use of LOCPROJ for t2g manifold of SrVO3 (the order of the orbitals seems to be mixed up... this example leads to x^2 -y^2, z^2, yz... )
In the current version there is some mix up in the mapping between selected orbitals in the INCAR and actual selected in the LOCPROJ. This is 
what the software does (left side is INCAR, right side is resulting in the LOCPROJ)

* xy -> x2-y2
* yz -> z2
* xz -> yz
* x2-y2 -> xz
* z2 -> xy

```
LOCPROJ = 2 : dxz : Pr 1
LOCPROJ = 2 : dx2-y2 : Pr 1
LOCPROJ = 2 : dz2 : Pr 1
```
However, if the complete d manifold is chosen, the usual VASP order (xy, yz, z2, xz, x2-y2) is obtained in the LOCPROJ. This is done as shown below
```
LOCPROJ = 2 : d : Pr 1
```

## convergence of projectors with Vasp

for a good convergence of the projectors it is important to convergence the wavefunctions to high accuracy. Otherwise this often leads to off-diagonal elements in the the local Green's function. To check convergence pay attention to the rms and rms(c) values in the Vasp output. The former specifies the convergence of the KS wavefunction and the latter is difference of the input and out charge density. Note, this does not necessarily coincide with good convergence of the total energy in DFT! Here an example of two calculations for the same system, both converged down to `EDIFF= 1E-10` and Vasp stopped. First run:

```
       N       E                     dE             d eps       ncg     rms          rms(c)
...
DAV:  25    -0.394708006287E+02   -0.65893E-09   -0.11730E-10 134994   0.197E-06  0.992E-05
...
```
second run with different smearing:
```
...
DAV:  31    -0.394760088659E+02    0.39472E-09    0.35516E-13 132366   0.110E-10  0.245E-10
...
```
The total energy is lower as well. But more importantly the second calculation produces well converged projectors preserving symmetries way better, with less off-diagonal elements in Gloc, making it way easier for the solver. Always pay attention to rms.


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


## remarks on the Vasp version

### General remarks

One can use the official Vasp 5.4.4 patch 1 version with a few modifications:

- there is a bug in `fileio.F` around line 1710 where the code tries print out something like "reading the density matrix from Gamma", but this should be done only by the master node. Adding a `IF (IO%IU0>=0) THEN ... ENDIF` around it fixes this
- in the current version of the dft_tools interface the file `LOCPROJ` should contain the fermi energy in the header. Therefore  one should replace the following line in `locproj.F`:
```
WRITE(99,'(4I6,"  # of spin, # of k-points, # of bands, # of proj" )') NS,NK,NB,NF
```
by
```
WRITE(99,'(4I6,F12.7,"  # of spin, # of k-points, # of bands, # of proj, Efermi" )') W%WDES%NCDIJ,NK,NB,NF,EFERMI
```
and add the variable `EFERMI` accordingly in the function call.
- Vasp gets sometimes stuck and does not write the `OSZICAR` file correctly due to a stuck buffer. Adding a flush to the buffer to have a correctly written `OSZICAR` to extract the DFT energy helps, by adding in `electron.F` around line 580 after
```
CALL STOP_TIMING("G",IO%IU6,"DOS")
```
two lines:
```
flush(17)
print *, ' '
```
- this one is __essential__ for the current version of the DMFT code. Vasp spends a very long time in the function `LPRJ_LDApU` and this function is not needed! It is used for some basic checks and a manual LDA+U implementation. Removing the call to this function in `electron.F` in line 644 speeds up the calculation by up to 30%! If this is not done, Vasp will create a GAMMA file each iteration which needs to be removed manually to not overwrite the DMFT GAMMA file!
- make sure that mixing in VASP stays turned on. Don't set IMIX or the DFT steps won't converge!

### Enabling CSC calculations with Wannier90 projectors

You basically only need to add two things to have W90 run in Vasp's CSC mode, all in `electron.F`:

- the line `USE mlwf` at the top of the `SUBROUTINE ELMIN` together with all the other `USE ...` statements.
- right below where you removed the call to `LPRJ_LDApU` (see above, around line 650), there is the line `CALL LPRJ_DEALLOC_COVL`. Just add the following block right below, inside the same "if" as the `CALL LPRJ_DEALLOC_COVL`:
```
IF (WANNIER90()) THEN
   CALL KPAR_SYNC_ALL(WDES,W)
   CALL MLWF_WANNIER90(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
ENDIF
```
Then, the only problem you'll have is the order of compilation in the `.objects` file. It has to change because now electron.F references mlwf. For that move the entries `twoelectron4o.o` and `mlwf.o` (in this order) up right behind `linear_optics.o`. Then, move the lines from `electron.o` to `stm.o` behind the new position of `mlwf.o`.

Remarks:

- W90-CSC requires Wannier90 v3, in v2 the tag write_u_matrices does not work correctly. Until now, linking W90 v3 to Vasp with the `DVASP2WANNIER90v2` has worked without any problems even though it is not officially supported
- symmetries in Vasp should remain turned on, otherwise the determination of rotation matrices in dft_tools' wannier converter will most likely fail

### Speeding up the DFT part by not writing projectors at every step
This is very important for CSC calculations with W90 but also speeds up the PLO-based ones.

Writing the Wannier projectors is a time consuming step (and to a lesser extent, the PLO projectors) and basically needs only to be done in the DFT iteration right before a DMFT iteration. Therefore, solid_dmft writes the file `vasp.suppress_projs` that tells Vasp when __not__ to compute/write the projectors. This requires two small changes in `electron.F` in the Vasp source code:

- adding the definition of a logical variable where all other variables are defined for `SUBROUTINE ELMIN`, e.g. around line 150, by inserting `LOGICAL :: LSUPPRESS_PROJS_EXISTS`
- go to the place where you removed the call to `LPRJ_LDApU` (see above, around line 650). This is inside a `IF (MOD(INFO%ICHARG,10)==5) THEN ... ENDIF` block. This whole block has to be disabled when the file `vasp.suppress_projs` exists. So, right under this block's "IF", add the lines
```
INQUIRE(FILE='vasp.suppress_projs',&
        EXIST=LSUPPRESS_PROJS_EXISTS)

IF (.NOT. LSUPPRESS_PROJS_EXISTS) THEN
```
and right before this block's "ENDIF", add another `ENDIF`.
