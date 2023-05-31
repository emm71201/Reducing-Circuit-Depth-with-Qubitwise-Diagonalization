#%% importing
from tableau_operations import makeTableauMatrix
import itertools as it
from diag_optimization import weight
#%% My second implementation of Brute force. This time, I search through the whole normalizer
def commute(p1, p2):
 
    assert len(p1) == len(p2)
    n = len(p1)//2
    
    return int(p1[:n].dot(p2[n:]) + p1[n:].dot(p2[:n])) == 0

def commute_with_others(p, others):
    for other in others:
        if not commute(p, other):
            return False
    
    return True
# %%
def normalizer(T):

    n = T.shape[1]//2
    candidates = it.product([0,1], repeat = 2*n)

    for candidate in candidates:

        candidate = GF(candidate)

        if commute_with_others(candidate, T):
            yield candidate
        
# %%
def search_normalizer(x,z):
    T = makeTableauMatrix(x,z)
    normal = normalizer(T)
    #normal = list(normal)

    currp = T[0]
    currw = weight(currp)


    for nextp in normal:

        nextw = weight(nextp)
        if nextw < currw and nextw !=0:
            currp = nextp
            currw = nextw

    
    return currp

# condidering linear connectivity

# weight2 = lambda pstr: sum(pstr != 0 for c in pstr)
# commute2 = lambda p1,p2: bool((sum ( p1[c] != p2[c] and (p1[c] != "0" and p2[c] != "0") for c in range(len(p1))) + 1) % 2)

#%%
def compute_swap_cost_linear(cand):

    """ compute the number of swap gates assuming linear connectivity """
    n = len(cand)//2

    pivots = [j for j in range(n) if cand[j] == 1 or cand[j+n] == 1 ]

    return sum( abs(pivots[j] - pivots[j+1]) for j in range(len(pivots) - 1)) - 1

def search_normalizer2(x,z):
    T = makeTableauMatrix(x,z)
    normal = normalizer(T)
    #normal = list(normal)

    currp = T[0]
    currswap = compute_swap_cost_linear(T[0])

    if currswap == 0:
        return currp


    for nextp in normal:

        next_swap_cost = compute_swap_cost_linear(nextp)

        if next_swap_cost == 0:
            return nextp

        if next_swap_cost < currswap:

            currswap = next_swap_cost
            currp = next_swap_cost 

    return currp