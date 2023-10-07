import sys

import numpy as np
from h5 import HDFArchive
from scipy.constants import physical_constants

###
# This script read bdft output and dump g0w0 eigenvalues to si.eig for wannier90 interpolation
###

# first arg is the h5 to use
if len(sys.argv) < 2:
    print ("Usage: python qp_evs_to_eig.py <h5>")
    quit()
print('h5 archive:', str(sys.argv[1]))
bdft_output = str(sys.argv[1])

Hartree_eV = physical_constants['Hartree energy in eV'][0]

###### params ######
# bdft output is defined by "ouptut" in "evgw0" block.
# number of bands used in wannier90.
# It should be consistent with  "num_bands" in si.win
nbnd_pp = 3
# number of bands to include in the beginning
excl = 20
# iteration of evGW0
it = None
# seed for w90
seed = 'svo'

############
############

# Read evGW0 eigenvalues
with HDFArchive(bdft_output, 'r') as ar:
    if not it:
        it = ar['scf']['final_iter']
    eig = ar[f'scf/iter{it}/E_ska'].real * Hartree_eV


# Write eigenvalues in the format of Quantum Espresso .eig file
ns, nkpts, nbnd_gw = eig.shape
with open(f'{seed}.eig', 'w') as f:
    for k in range(1, nkpts+1):
        for i in range(excl+1, excl+nbnd_pp+1):
            f.write("   {}   {}   {}\n".format(i-excl, k, eig[0, k-1, i-1]))
