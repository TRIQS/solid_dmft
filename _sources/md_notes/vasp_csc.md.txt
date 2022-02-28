# Interface to VASP


## General remarks

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

## Speeding up by not writing projectors at every step


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
