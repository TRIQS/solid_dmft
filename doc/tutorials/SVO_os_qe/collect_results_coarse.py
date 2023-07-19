import numpy as np
#avoid connecting to the X window manager, needs to be here before importing pyplot
import matplotlib
matplotlib.use('Agg')
from heatmap_func import heatmap, annotate_heatmap
import matplotlib.pyplot as plt
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_vec = np.array([0.0, 0.25, 0.5, 0.75, 1])
U_vec = np.array(range(2, 11, 2))
G_array = np.zeros((len(J_vec), len(U_vec)))
for j in range(len(J_vec)) :
    J = J_vec[j] 
    for u in range(len(U_vec)):
        U = U_vec[u]
        filepath = f'./J{J}/U{U}/out'
        print(f'I am currently reading the folder: {filepath}')
        
        with open(f'{filepath}/observables_imp0.dat') as f:
            for line in f:
                pass
            last_line = line
            G_bet2 = line.split("|")[2]
            G_bet2 = list( map(float, G_bet2.split()) )
            G_array[j,u] = np.abs(G_bet2[0])
            print(G_bet2)

#plotting routine
fig, ax = plt.subplots(figsize=(8, 8))
im, cbar = heatmap(G_array.transpose(),U_vec,J_vec, ax=ax,
                   cmap="viridis", cbarlabel="|G(beta/2)| [eV^-1]", origin='lower', aspect='auto')
texts = annotate_heatmap(im, valfmt="{x:.2f}")
plt.ylabel("Hubbard U [eV]", size =16)
plt.xlabel("Hund J [eV]", size=16)

plt.savefig("MIT_coarse.jpg", dpi=300, bbox_inches='tight')
plt.show()

