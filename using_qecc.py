#%%
import qecc as q
import numpy as np
import helpers
#%%

def to_pstring(pvector):
    n = len(pvector)//2
    res = ""
    for i in range(n):
        if pvector[i] == 0 and pvector[i+n] == 0:
            res += "I"
        if pvector[i] == 0 and pvector[i+n] == 1:
            res += "Z"
        if pvector[i] == 1 and pvector[i+n] == 0:
            res += "X"
        if pvector[i] == 1 and pvector[i+n] == 1:
            res += "Y"
    
    return res
        
def get_logical(null_space):
    """Return one logical operator as a string of Pauli matrices 
    null space basis is actually the generator of the normalizer group.
     We just return the first pauli in this basis """
    
    pvector = null_space[0]

    return to_pstring(pvector)

def find_normalizer(x,z, null_space):

    pstrings = [p.pauli_matrix_form() for p in helpers.tableau_to_pstrings(x,z)]

    logicalx = get_logical(null_space)

    stab = q.StabilizerCode(pstrings, [logicalx], [])

    return stab.normalizer_group()

def minimize_weight(x,z,null_space):

    r,n = x.shape

    normalizer = find_normalizer(z,x,null_space)
    # for n in normalizer:
    #     print(n)

    curr_p = next(normalizer)
    if curr_p.wt == 0:
        curr_p = next(normalizer)
    curr_wt = curr_p.wt

    for tmp_p in normalizer:

        tmp_wt = tmp_p.wt
        if tmp_wt < curr_wt and tmp_wt != 0:
            curr_p = tmp_p
            curr_wt = tmp_wt

        # if curr_wt < r//4 + 1:

        #     xpart = curr_p.as_bsv().x
        #     zpart = curr_p.as_bsv().z

        #     print(f"weight = {curr_wt}")

        #     return np.concatenate((xpart, zpart))
    
    xpart = curr_p.as_bsv().x
    zpart = curr_p.as_bsv().z

    # print(f"weight = {curr_wt}")

    return np.concatenate((xpart, zpart))

# %%
