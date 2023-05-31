
#%% importing
from tableau_operations import makeTableauMatrix
import itertools as it
from diag_optimization import weight
import galois
GF = galois.GF(2)
# %%
def weight(vector):

    n = len(vector)//2
    res = 0

    for j in range(n):
        if vector[j] == 1 or vector[j+n] == 1:
            return res
    
    return res
    
def search_normalizer(null_sapce):

    r = null_sapce.shape[0]
    n = null_sapce.shape[1]//2

    curr_vector = null_sapce[0]
    curr_weight = weight(curr_vector)

    for curr_r in range(1, r//2 + 1 + 1):

        templates = it.combinations(null_sapce, r=curr_r)

        for template in templates:

            vector = GF([0]*(2*n))
            for t in template:
                vector += t

            tmp_weight = weight(vector)

            if tmp_weight < curr_weight:

                curr_weight = tmp_weight
                curr_vector = vector
    
    return curr_vector

