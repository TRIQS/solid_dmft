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
U = 4.25
J = 0.0

filename = f'./J{J}/U{U}/out/svo.h5'
internal_prefix =  'DMFT_results/last_iter/Aimp_w_0/'

# %%
a = np.array(read_h5(filename, f'{internal_prefix}Aimp_w_line_fit/up_0'))
#b = np.array(file['DMFT_results']['last_iter']['Aimp_w_0']['Aimp_w_chi2_curvature']['down_0'])
en = np.array(read_h5(filename, f'{internal_prefix}mesh/up_0'))
print(a.shape)
# %%
plt.plot(en, a[0,0][:])
# plt.plot(en, b[0,0][:,0])

plt.xlabel('Energy [eV]')
plt.xlim([-7,7])
plt.ylabel('Spectral functions [eV^-1]')
plt.title(f'Correlated metal Phase (J = {J} eV U = {U} eV)')
plt.savefig(f'A_func_{J=}_{U=}.jpg')
# %%
