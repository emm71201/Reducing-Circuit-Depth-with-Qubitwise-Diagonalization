#%% Setup
import sys
sys.path.append("../")
import tableau
import pstring 
import numpy as np
# %%
weight2 = lambda pstr: sum(pstr != "0" for c in pstr)
commute2 = lambda p1,p2: bool((sum ( p1[c] != p2[c] and (p1[c] != "0" and p2[c] != "0") for c in range(len(p1))) + 1) % 2)

#%%
def compute_swap_cost_linear(cand):

    """ compute the number of swap gates assuming linear connectivity """
    n = len(cand)//2

    pivots = [j for j in range(n) if cand[j] == 1 or cand[j+n] == 1 ]

    return sum( abs(pivots[j] - pivots[j+1]) for j in range(len(pivots) - 1)) - 1


def find_normalizer(nullspace)