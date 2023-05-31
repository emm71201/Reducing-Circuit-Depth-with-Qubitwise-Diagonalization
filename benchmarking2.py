#%%[markdown]
# Author Edison <br>
# Date: Feb. 3rd, 2023 <br>
# Benchmarking the diagonalization algorithm
#%%
import sys
sys.path.append("../Simultaneous-diagonalization-main")
import cl
import helpers as hlp
from diagonalize import *
import diagonalize as dg
import galois
GF = galois.GF(2)
import numpy as np

#%%

def get_random_paulis(n, k):

    T = cl.random_commuting_matrix(n,k,flagSign=False, flagComplex=False, init=None)

    X = T[:, 0:n]
    Z = T[:, n+1:2*n+1]

    return X, Z

n = 25
for j in range(20):

    print(f"j = {j+1}")
    fname = f"./random_paulis_sets/tableau_n_{n}_j_{j}.npy"
    np.save(fname, get_random_paulis(n,n))
    
#%%[markdown]
# GET the number of CNOT and the Circuit Depth for Diagonalizing the Set
def get_count(X,Z):

    pstrings = hlp.tableau_to_pstrings(X,Z)
    X,Z,S,U = dg.main_diagonalizer(pstrings, reduced_h=False)
    if X.any():
        print("Did not work")

    depth = U.depth()

    if "cx" in U.count_ops():
        return U.count_ops()["cx"], depth
    
    return 0, depth
# %%[markdown]
# Diagonalize each set and record the CNOT count and the depth

# nlist = [20]
for n in [20, 25]:
    print("n = ", n)
    with open(f"benchmark_result/result_heuristics_{n}.csv", "w") as f:

        for j in range(20):

            print(f"j = {j+1}")
            x,z = np.load("./random_paulis_sets/tableau_n_{0}_j_{1}.npy".format(n,j))
            count,depth = get_count(x,z)

            f.write(f"{n},{count},{depth}\n")

# %%[markdown]
# Fix n and vary r. 
# Randomly choose which Pauli to remove
import random as rd
j=0
n = 10
x,z = np.load("./random_paulis_sets/tableau_n_{0}_j_{1}.npy".format(n,j))
def remove_rows(x,z,rows):

    return np.delete(x,rows, axis=0), np.delete(z,rows, axis=0)

def compile_r(x,z,r):

    """ use only r Paulis in the diagonalization """
    nremoved = x.shape[0] - r
    
    try:
        rows = rd.sample(range(x.shape[0]), nremoved)
    except ValueError:
        print("Something is wrong")
        return

    xfiltered, zfiltered = remove_rows(x,z,rows)

    count,depth = get_count(xfiltered,zfiltered)

    return count, depth


def experiment(n, j):

    x,z = np.load("./random_paulis_sets/tableau_n_{0}_j_{1}.npy".format(n,j))

    res = []

    for r in range(1,n+1):
        count, depth = compile_r(x,z,r)

        res.append([n,r,count,depth])
    
    return res

def run_experiments(n):

    res = []
    for j in range(20):

        res.append(experiment(n,j))
    
    return res
# %%
for n in [15]:
    np.save(f"benchmark_fixed_n/fixed_n_{n}", run_experiments(n))
# %%
