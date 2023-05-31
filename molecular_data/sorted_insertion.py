# implementing sorted insertion in Crawford et. al.
# By Edison Murairi

#%% importing
import sys
sys.path.append("../pdiag/")
import tableau
import helpers
import diagonalize
# %%
def get_count(pstrings):

    X,Z,S,U = diagonalize.main_diagonalizer(pstrings)
    depth = U.depth()

    if "cx" in U.count_ops():
        return U.count_ops()["cx"], depth
    
    return 0, depth

def commute_with_all(p, others):

    """check if pauli string p commutes with a list of pauli strings other"""

    for other in others:

        if not p.commute(other):

            return False
    
    return True

def grouping(pstrings):

    res = {}
    maxkey = 0
    res[maxkey] = [pstrings[0]]

    indx = 1
    while indx < len(pstrings):

        p = pstrings[indx]
        inserted = False

        for grp in range(maxkey + 1):

            if commute_with_all(p, res[grp]) and not inserted:
                
                res[grp].append(p)
                inserted = True
                break
        
        if not inserted:
            maxkey += 1
            res[maxkey] = [p]
        
        indx += 1
    
    return res
            

def sorted_insertion(pstrings):

    psorted = sorted(pstrings, key=lambda p: abs(p.coef), reverse=True)

    return grouping(psorted)

#%% running the experiments
files = ["LiH", "HF", "H2O", "BH3", "CH4"]
for f in files:

    print("Molecule: ", f)

    with open(f"{f}_results.csv", "w") as infile:
        ham = helpers.read_hamiltonian(f"../molecular_data/{f}.csv")
        n = len(ham[0].string)

        grouped = sorted_insertion(ham)

        infile.write("key, #Paulis, #qbits, cxcount, depth\n")

        for key, item in grouped.items():

            m = len(item)

            cxcount, depth = get_count(item)

            infile.write(f"{key}, {m}, {n}, {cxcount}, {depth}\n")
        
# %%
