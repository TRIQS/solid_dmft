# Wannier90 interface

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
