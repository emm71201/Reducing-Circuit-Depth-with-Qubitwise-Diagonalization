
# Author: Edison Murairi
# Date: February 21st, 2023
# Given a stabilizer S, find a Pauli opertor in N(S) with miniamal weight d

#%%
import diag_helpers as helpers
import itertools as it
#from diag_optimization import weight
import numpy as np
import galois
GF = galois.GF(2)

weight = lambda pstr,n: sum(int(pstr[c])==1 or int(pstr[c+n]) == 1 for c in range(n))

from tableau_operations import makeTableauMatrix

def generators_strings(x,z):

    return helpers.tableau_to_pstrings(x,z)

def template(k,n,r):

    return it.combinations(range(r), r=k)

def apply_template(nullspace, template):

    res = GF([0]*nullspace.shape[1])

    for t in template:

        res += nullspace[t]

    return res

def linear_swap_cost(vector):

    n = len(vector)//2

    pivots = [j for j in range(n) if (vector[j] == 1 or vector[j+n] == 1)]

    cost = 0
    for j in range(len(pivots)-1, 0, -1):


        cost += (pivots[j] - pivots[j-1]) - 1
    
    return cost

def linear_swap_cost_optimize(vectors):

    curr_vector = vectors[0]
    curr_cost = linear_swap_cost(curr_vector)

    for vector in vectors:

        tmp_cost = linear_swap_cost(vector)

        if tmp_cost < curr_cost:

            curr_vector = vector
            curr_cost = tmp_cost
    
    return curr_vector


def brute_force_optimize(nullspace, connectivity="full"):

    r = nullspace.shape[0]
    n = nullspace.shape[1]//2
    k = 1

    # curr_vector = nullspace[0]
    # curr_weight = weight(curr_vector, n)

    curr_vectors = [nullspace[0]]
    curr_weight = weight(curr_vectors[0], n)

    while k <= r//2 + 1:

        templts = template(k, n, r)

        for templ in templts:

            tmp_vector = apply_template(nullspace, templ)
            tmp_weight = weight(tmp_vector, n)

            if tmp_weight == 1:
                return tmp_vector

            # if tmp_weight < curr_weight:

            #     curr_vector = tmp_vector
            #     curr_weight = tmp_weight
            
            if tmp_weight < curr_weight:

                curr_vectors = [tmp_vector]
            
            if tmp_weight == curr_weight:

                curr_vectors.append(tmp_vector)

        k += 1
    

    if connectivity == "linear":

        return linear_swap_cost_optimize(curr_vectors)

    return curr_vectors[0]

#%%
def heuristics(nullspace):

    r = nullspace.shape[0]
    n = nullspace.shape[1]//2

    curr_vector = nullspace[0]
    curr_weight = weight(curr_vector, n)

    for j in range(r):

        tmp_vector = nullspace[j]
        tmp_weight = weight(tmp_vector, n)

        if tmp_weight < curr_weight:

            curr_vector = tmp_vector
            curr_weight = tmp_weight
    
    return curr_vector


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
    T = makeTableauMatrix(z,x)
    normal = normalizer(T)
    #normal = list(normal)

    n = T.shape[1]//2

    currp = T[0]
    currw = weight(currp, n)


    for nextp in normal:

        nextw = weight(nextp, n)
        if nextw < currw and nextw !=0:
            currp = nextp
            currw = nextw

    return currp

