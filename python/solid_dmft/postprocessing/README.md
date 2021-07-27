This folder contains scripts for post-processing with MaxEnt as well as a few other utility scripts. All of these scripts are completely independent of the main code of solid-dmft (except of maxent_gf_latt.py, which is optionally used to set the chemical potential to the middle of the band gap).

# MaxEnt-Scripts

This is a folder containing stand-alone scripts for analytic continuations.
They all use maxent (https://triqs.github.io/maxent/master/index.html) and are generally designed to read DMFT results from the h5 archive, perform the analytic continuation and then write the information back to the h5 archive.
Available now are functions to get the

- impurity spectral function A\_imp(omega): `maxent_gf_imp.py`
- lattice spectral function A\_latt(omega): `maxent_gf_latt.py`
- impurity self-energy Sigma(omega): `maxent_sigma.py`

With the impurity self-energy Sigma(omega), you can use `make_spaghetti.py` to generate

- the lattice spectral function A\_latt(omega)
- the k-resolved lattice spectral function A\_latt(k, omega) (i.e. DMFT bands, spaghetti).
For names of high-symmetry points, you can use the k coordinates in `dft_bands_input/kpoints` in the h5 archive. For the DMFT bands, you need the Kohn-Sham Hamiltonian along a k path in the same basis as the DMFT calculations, i.e., `dft_input/hopping`.

## WARNING

Analytic continuations are by their nature rather inexact.
For the continuation of Sigma(i omega), we need to construct an auxiliary Green's function.
This makes the procedure even less exact than it already is so use the results with care!

## Usage
These files can be run either from the terminal with
```
python <path to script>/<script_name> <h5_name> (<iteration>)
```
or from within another python program by
```
import sys
sys.path.append('<path to script>')
import <script_name>
<script_name>.main(h5_name (, iteration))
```
If you run it from within another python program, the `main()` method not only writes the results back but also returns them directly.
If you run it from the terminal, you can use wildcards (like `*`) in the h5 name, where this path with wildcards has to be surrounded by a single inverted comma `'`.
Then, the program will loop over all possible paths, with a parallelization of up to 8 maxent runs at the same time.
In both ways, if you leave out the iteration, the script will read and write from the group 'DMFT\_results/last\_iter'.

## Notes

- if you get the error `ImportError: No module named functools_lru_cache` when running the scripts from the terminal, you need to either update your docker container with the newest Dockerfile or import the script to a jupyter notebook, as detailed above.
- it is normal that in `maxent_sigma.py`, the step of extracting Sigma(omega) takes much longer than the maxent run
- `maxent_gf_latt.py` can only calculate the trace of the lattice GF.
- in `make_spaghetti.py`, the calculation of the k-integrated spectral function is not necessary.
Removing it will speed up this script by a factor of ~2 but you lose a way to compare your results to normal spectral functions.

# Other scripts

To be written...