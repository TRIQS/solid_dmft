# Examples

Here one finds multiple examples for running the DFT+DMFT code. For the CSC
examples a working, adapted, VASP version 5.4.4 and suitable POTCAR files are
needed. There are two example systems: SrVO3 and LuNiO3. For the LuNiO3 CSC examples
please use the following PAW-PBE POTCARs: Lu_3, Ni_pv, O. For the SrVO3 CSC
example please use the PAW-PBE POTCARs: Sr_sv, V, O.

In the example for svo-one-shot is also a folder 'converter' which contains the
work flow for creating the h5 archive.

To run any of the one-shot examples execute:
```
mpirun solid_dmft
```
within one of the directories. For the CSC calculations follow the instructions
in the documentation of this repository.
