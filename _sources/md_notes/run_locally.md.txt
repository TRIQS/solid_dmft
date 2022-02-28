# Run on your machine

## CSC calculations locally

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
