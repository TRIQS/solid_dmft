# Examples

This folder contains examples of solid_dmft calculations that are updated on a semi-regular
basis. The current version works with solid_dmft 3.1.0 and, for CSC calculations,
with VASP 6.3.0.

For a better documented introduction into solid_dmft, please look in the folder ../tutorials.

There are two example systems: SrVO3 and LuNiO3.

## POTCARs

For the CSC examples, the Vasp POTCAR files are needed. For the LuNiO3 CSC examples
please use the following PAW-PBE POTCARs: Lu_3, Ni_pv, O. For the SrVO3 CSC
example please use the PAW-PBE POTCARs: Sr_sv, V, O.

## How to run

In the example for svo-one-shot is also a folder 'converter' which contains the
work flow for creating the h5 archive.

To run any of the one-shot examples execute:
```
mpirun solid_dmft
```
within one of the directories. For the CSC calculations follow the instructions
in the documentation of this repository.

