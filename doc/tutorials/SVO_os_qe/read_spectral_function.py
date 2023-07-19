#%%

import numpy as np
import matplotlib.pyplot as plt
from h5 import HDFArchive


def read_h5(external_path, h5_internal_path):
    """Reads the A spectral function from h5 archive."""

    with HDFArchive(external_path, 'r') as archive:
        read_arr = archive[h5_internal_path]
    return read_arr
    
# %%
U = 4
J = 0.75
filepath = f'./J{J}/U{U}'

filename = f'{filepath}/out/svo.h5'
internal_prefix =  'DMFT_results/last_iter/Aimp_maxent_0'

# %%
A_imp = np.array(read_h5(filename, f'{internal_prefix}/Aimp_w_line_fit/up_0'))
omega = np.array(read_h5(filename, f'{internal_prefix}/mesh'))
# %%
plt.plot(omega, A_imp[0,0][:].real)

plt.xlabel('Energy [eV]')
plt.xlim([-7,7])
plt.ylabel('Spectral functions [eV^-1]')
plt.title(f'J = {J} eV U = {U} eV')
plt.savefig(f'A_func_{J=}_{U=}.jpg')
# %%
