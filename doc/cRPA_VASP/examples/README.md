# examples for VASP cRPA calculations

Three examples for a typical VASP cRPA calculation. In each folder is a
`job_script.sh` that shows the workflow. Note, that the wannier90.win file is
already complete and wannier90 needs not to be run in between manually. However,
if you change the number of bands or kpoints during any of the VASP runs the
wannier90 call will fail or produce spurious output. Therefore, check always the
wannier90.wout files for weird results.

For all calculations one can compare the results of cutting bands vs. using the
wannier functions as target states. For example for SrVO3 this would correspond
to these two different options. To select target states via bands:
```
NCRPAHIGH= 23
NCRPALOW= 21
```
or wannier target states:
```
NTARGET_STATES = 3*1
```
