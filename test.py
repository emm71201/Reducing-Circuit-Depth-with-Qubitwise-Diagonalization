# Author: Edison Murairi
# Date: Feb. 3rd, 2023
# I use the algorithm in E. Van de Berg to generate commuting Paulis randomly

#%% importing
import helpers
import diagonalize as dg
from tableau_operations import *
import brute_force_optimization as bfo
# %%
pstringsLiH = helpers.read_hamiltonian("misc/LiH2.txt")
pstrings = helpers.read_hamiltonian("misc/miller.txt")

cluster = []
i = 0
j  = 0
while i < len(pstrings):

    curr = pstrings[i]
    flag = True

    for j in range(len(cluster)):

        if not curr.commute(cluster[j]):
            flag = False
    
    if flag:
        cluster.append(curr)
    
    i += 1

x,z = np.load("random_paulis_sets/tableau_n_5_j_0.npy")
#%%
# bfo.heuristics(None, x,z)

# %%
print("without swaps: ")
X,Z,S,U=dg.main_diagonalizer(pstrings, reduced_h=False, connectivity="full", optimize=False)
print(U)



# %%
