# Author Edison
# Date: Feb. 3rd, 2023
# Benchmarking the diagonalization algorithm

#%%
import sys
sys.path.append("../Simultaneous-diagonalization-main")
import cl
import helpers as hlp
from diagonalize import *
import galois
GF = galois.GF(2)
import numpy as np
#%%[markdown]
# GET the number of CNOT and the Circuit Depth for Diagonalizing the Set
def get_count(X,Z):

    pstrings = hlp.tableau_to_pstrings(X,Z)
    X,Z,S,U = main_diagonalizer(pstrings)
    depth = U.depth()

    if "cx" in U.count_ops():
        return U.count_ops()["cx"], depth
    
    return 0, depth
# %%[markdown]
# Generate and save sets of commuting Paulis randomly using Ewout's codes
def get_paulis(n):

    res = cl.random_commuting_matrix(n, n)[:,:-2]

    X = res[:,:n]
    Z = res[:, n:]

    return X,Z
    
nlist = [5,10,15,20,25]

# for n in nlist:

#     print("n = ", n)

#     for j in range(20):

#         x,z = get_paulis(n)

#         np.save("./random_paulis_sets/tableau_n_{0}_j_{1}.npy".format(n,j), (x,z))
# %%
nlist = [5,10,15,20]
for n in nlist:
    print("n = ", n)

    with open(f"benchmark_result/benchmark{n}.csv", "w") as f:

        for j in range(20):

            print(f"j = {j+1}")
            x,z = np.load("./random_paulis_sets/tableau_n_{0}_j_{1}.npy".format(n,j))
            count,depth = get_count(x,z)

            f.write(f"{n},{count},{depth}\n")
# %%
