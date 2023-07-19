
# run the scf calculation
mpirun -n 36 pw.x -i svo.scf.in | tee vo.scf.out


# OPTIONAL: run the DFT bandstructure calculation
mpirun -n 36 pw.x -i svo.bnd.in | tee svo.bnd.out
mpirun -n 36 bands.x -i svo.bands.in | tee svo.bands.out


# run the nscf calculation
mpirun -n 36 pw.x -i svo.nscf.in | tee svo.nscf.out

# Wannierize

## pre-process wannier90
wannier90.x -pp svo
## run interface wannier90-quantum espresso
mpirun -n 36 pw2wannier90.x -i svo.pw2wan.in | tee svo.pw2wan.out
## run wannier90
wannier90.x svo

# Interface with TRIQS, this will create the svo.h5 archive to start the calculation
python3 ./convert_wannier.py
