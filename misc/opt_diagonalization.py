# Author: Edison Murairi
# Date: March 14th, 2023
# Given a stabilizer S, find a Pauli opertor in N(S) with miniamal weight d

#%%
import helpers
import itertools as it
from diag_optimization import weight
import numpy as np
import galois
GF = galois.GF(2)

from tableau_operations import makeTableauMatrix
from tableau_operations import getIndependentPauliStrings

# %%
def pstring_to_vector(mypstr):
    n = len(mypstr)
    res = np.zeros(2*n, dtype=int)
    for i in range(n):
        if mypstr[i] == "1":
            res[i] = 1
        if mypstr[i] == "2":
            res[i] = 1
            res[i + n] = 1
        if mypstr[i] == "3":
            res[i+n] = 1
    return GF(res)

def commute(p1, p2):
 
    assert len(p1) == len(p2)
    n = len(p1)//2

    return int(p1[:n].dot(p2[n:]) + p1[n:].dot(p2[:n])) == 0


def commute_with_others(p, others):
    for other in others:
        if not commute(p, other):
            return False
    
    return True

def insertpmatrices(n, qbits, pmatrices):

    pstr = ['0']*n

    for q in range(len(qbits)):
        pstr[qbits[q]] = str(pmatrices[q])
    
    return "".join(pstr)


def candidate_weight_w(x,z,wght):

    n = x.shape[1]
    T = makeTableauMatrix(z,x)
    qbitsList = it.combinations(range(n), r=wght)

    for qbits in qbitsList:

        pmatricesList = it.product(range(1,4), repeat=wght)

        for pmatrices in pmatricesList:

            pstr = insertpmatrices(n, qbits, pmatrices)
            p_op = pstring_to_vector(pstr)
            if commute_with_others(p_op, T) and weight(p_op) != 0:
                return p_op
    
    return None

# %%
def search_for_vector(x,z):

    r, n = x.shape

    for wght in range(r//2 + 2):

        cand = candidate_weight_w(x,z,wght)

        if not cand is None:

            return cand
    
    print("we should never get here")

    return 
# %%
n = 15
x,z = np.load(f"sets/tableau_n_{n}_j_3.npy")

#%% 
helpers.tab