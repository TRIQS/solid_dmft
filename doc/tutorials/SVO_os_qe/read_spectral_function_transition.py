#%%

import numpy as np
#avoid connecting to the X window manager, needs to be here before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from h5 import HDFArchive


colorcycle = [
 '#004c6d',
 '#006d92',
 '#0090b7',
 '#00a3ca',
 '#00b5dc',
 '#00c9ee',
 '#00dcff']

def read_h5(external_path, h5_internal_path):
    """Reads the A spectral function from h5 archive."""

    with HDFArchive(external_path, 'r') as archive:
        read_arr = archive[h5_internal_path]
    return read_arr
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_vec = np.arange(4, 7.5, 0.5)
J=0.0
print(U_vec)

for u in range(len(U_vec)):
    U = np.format_float_positional(float(U_vec[u]), trim='-')
    filepath = f'./J{J}/U{U}/out'
    print(f'I am currently reading the folder: {filepath}')
        
    filename = f'./J{J}/U{U}/out/svo.h5'
    internal_prefix = 'DMFT_results/last_iter/Aimp_w_0/'
    a = np.array(read_h5(filename, f'{internal_prefix}Aimp_w_line_fit/up_0'))
    en = np.array(read_h5(filename, f'{internal_prefix}mesh/up_0'))
    plt.plot(en, a[0,0][:]-u*1, label=f'U={U}eV ', color=colorcycle[u])


plt.legend()
plt.grid(axis='y')
plt.xlabel('Energy [eV]')
plt.xlim([-7,7])
plt.ylabel('Spectral functions [eV^-1]')
plt.title(f'Mott-Transition at J = 0 eV')
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])
plt.savefig(f'A_func_transition.jpg', dpi=600 )
# %%
