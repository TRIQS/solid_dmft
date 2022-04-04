import numpy as np

#avoid connecting to the X window manager, needs to be here before importing pyplot
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_vec = np.array([0.0, 0.25, 0.5, 0.65, 0.75, 1])
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
fig, ax = plt.subplots()
im = ax.imshow(X=G_array.transpose(), cmap="Purples", extent=[0, 1, 2, 10], aspect='auto', origin='lower')
fig.colorbar(im, ax=ax)


ax.set_xticklabels(J_vec)
#ax.set_xticks( range(len(J_vec)) )
ax.set_yticklabels(U_vec)
ax.set_yticks(U_vec)
#ax.set_xticklabels(piv.columns, rotation=90)
#ax.set_yticklabels(piv.index)
ax.set_xlabel("J [eV]")
ax.set_ylabel("U [eV]")
ax.set_title("G(beta/2) [eV^-1]")

plt.tight_layout()
plt.gcf()
plt.savefig('MIT_transition.jpg', dpi=300)
plt.show()

