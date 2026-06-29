# Interface to VASP


## General remarks

solid_dmft officially supports VASP 6.5.0 and newer out of the box, without any source-code modifications. The previously required manual patches to the VASP 5.4.4 source are no longer needed and have been removed from these notes. See the VASP [6.5.0 changelog](https://vasp.at/wiki/Changelog#6.5.0) for details.

To use the interface you need a VASP build with HDF5 support enabled (`-DVASP_HDF5` / the `vaspout.h5` and `vaspgamma.h5` output). VASP then communicates with solid_dmft through these HDF5 files instead of the legacy `LOCPROJ`/`GAMMA` text files. Relevant VASP documentation:

- [GAMMA](https://vasp.at/wiki/GAMMA)
- [vaspgamma.h5](https://vasp.at/wiki/Vaspgamma.h5)
- [DFT+DMFT calculations tutorial](https://vasp.at/wiki/DFT%2BDMFT_calculations)

A few things to keep in mind:

- make sure to set `LSYNCH5 = .TRUE.` in the `INCAR` so that VASP syncs `vaspout.h5` to disk and solid_dmft can read it while VASP holds the file open during CSC.
- keep an eye on the charge mixing. The DFT charge density has to converge within each CSC cycle, so watch the `rms(c)` (charge convergence) value in the VASP std out. For difficult cases, e.g. systems with charge order, it usually helps to slow down the mixing, as done in the `PrNiO3_csc_vasp_plo_cthyb` tutorial:
```
! slow mixing for charge order
IMIX = 1
AMIX = 0.04
BMIX = 0.3
```
Don't turn mixing off entirely, or the DFT steps won't converge.

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

## Enabling CSC calculations with Wannier90 projectors

The current VASP release 6.6.0 needs some changes to be able run CSC DFT+DMFT calculations.

### VASP source changes

**1. `src/electron.F`, `SUBROUTINE ELMIN`** — run the Wannier90 interface inside the
charge-update loop. Add `USE mlwf` to the `USE` statements, and in the
`IF (MOD(INFO%ICHARG,10)==5)` block that writes the projectors:

```
! gate to the non-delay iterations (same condition as the GAMMA/lock block below)
IF (MOD(INFO%ICHARG,10)==5 .AND. N >= ABS(INFO%NELMDL)) THEN
   IF (IO%LORBIT==14) CALL SPHPRO_FAST( ... )
   ! only write LOCPROJ-tag projectors if a LOCPROJ tag is actually set
   IF (ALLOCATED(LPRJ_functions)) THEN
      CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,LPRJ_functions,LPRJ_COVL,T_INFO,INFO,IO)
      CALL LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IO%IU6,IO%IU0,W)
   ENDIF
   ...
   CALL LPRJ_DEALLOC
   IF (WANNIER90()) THEN
      CALL KPAR_SYNC_ALL(WDES,W)
      CALL MLWF_FREE(MLWF_GLOBAL)
      CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
   ENDIF
ENDIF
```

**2. `src/mlwf.F`** — two fixes are needed because `MLWF_MAIN` is now called once per
charge-update cycle on the saved `MLWF_GLOBAL`:

- In `MLWF_MAIN`, reset the projection mode so it is re-determined on every call
  (add right before the `LSCDM`/`LPRJ_functions`/`LOCPROJ_AUTO` mode checks):
```
P_MLWF%PROJ_MODE = PROJ_UNKNOWN
```
- In `MLWF_WANNIER90_SETUP`, drop the stray `num_wann=0`; it resets the module-level
  `NUM_WANN` (which holds the user's tag). The count returned by `wannier_setup` is
  stored in `READ_NUM_WANN`:
```
! MLWF%NNTOT=0; num_bands=0; num_wann=0    ->
MLWF%NNTOT=0; num_bands=0
```

With the CMake build nothing else is required. For the classic makefile build,
`electron.F` now references `mlwf`, so the `.objects` order must change: move the
entries `twoelectron4o.o` and `mlwf.o` (in this order) up right behind
`linear_optics.o`, then move the lines from `electron.o` to `stm.o` behind the new
position of `mlwf.o`.

### INCAR

```
ICHARG  = 5
NELM    = 1000
NELMIN  = 1000          ! keep VASP iterating; solid_dmft drives termination
LSYNCH5 = .TRUE.        ! sync vaspout.h5 so solid_dmft can read it while VASP holds it open
LWANNIER90     = .TRUE.
LWRITE_MMN_AMN = .TRUE. ! write mmn amn files so that we can call manually wannier90.x
NUM_WANN = 3
WANNIER90_WIN = "
begin projections
V:dxy,dyz,dxz
end projections

exclude_bands = 1-16, 22-32
dis_win_min = 4.1
dis_win_max = 6.9
write_hr = true
write_u_matrices = true
"
```

Do **not** set `LWANNIER90_RUN = .TRUE.`: solid_dmft runs Wannier90 itself; VASP only
writes the `.amn`/`.mmn`/`.eig` and `.win` files each cycle.

### solid_dmft

```
[general]
seedname = wannier90
[dft]
dft_code = vasp
projector_type = w90
```

Remarks:

- Symmetries in Vasp should remain turned on; the wannier converter reduces the
  density correction to the irreducible k-points and needs the symmetry information to
  determine the rotation matrices.
