# Author: Eidson Murairi
# Date: January 4th, 2023

# Implement the Hamiltonian diagonalization algorithm developed by the authors
# The Hamiltonian is represented as a linear combination of Pauli operators (Pauli strings)
# We are given commuting Pauli operstors to diagonalize
# We represent the Paulis in Tableau
# We select the independent Pauli strings
# A basis that diagonalizes the independent Pauli strings also diagonalizes the whole set of commuting strings
# The columns of the tableau for independent commuting Pauli strings are linearly dependent.
# Therefore, the tableau for commuting Pauli strings has a non-trivial null space.
# We use the basis of the null space as instructions for diagonalization.
# We also implement a modification to reduce the number of single qubit gates
# See the paper for more details.


from tableau import *
from tableau_operations import *
import brute_force_optimization as bfo
import copy as cp

def mask_column(x,z,col, qmap):

    """This function mask a column, meaning that we will not consider such a column
    during the diagonalization.
    We will use it to mask columns that are already diagonal """
    n = x.shape[1]

    #update qmap
    for c in range(col, n-1):

        qmap[c] = qmap[c+1]
    del qmap[n-1]

    x = numpy.delete(x, col, 1)
    z = numpy.delete(z, col, 1)
    ind_pstrings = getIndependentPauliStrings(x,z)
    x,z,_,_ = tableau(ind_pstrings)

    return x,z,qmap

def mask_diagonal_columns(x,z,qmap):

    """This function finds all the columns that are already diagonal and mask them """

    n = x.shape[1]
    col = 0
    while col < n:

        if not x[:,col].any():

            x,z,qmap = mask_column(x,z,col,qmap)
            col -= 1
            n -= 1

        col += 1

    return x,z,qmap

#%%
getkey = lambda mydict, val: ([k for k in mydict if mydict[k]==val])

def apply_linear_connectivity(pivots, qmap, x,z,s,u, X,Z,S,U):

    # print("pivots = ", pivots)

    # start from last pivot
    for j in range(len(pivots)-1, 0, -1):

        # print(f"j = {j}")

        p1 = pivots[j]
        p2 = pivots[j-1]
        P2 = qmap[p2]
        

        P1tmp = qmap[p1]
        while abs(P1tmp - P2) > 1:

            # print(P1tmp, P2)
            X,Z,S,U = swapgate(X,Z,S,U, P1tmp, P1tmp - 1)

            tmpKey1, tmpKey2 = None, None
            
            try:
                tmpKey1 = getkey(qmap, P1tmp)[0]
            except:
                pass
            try:
                tmpKey2 = getkey(qmap, P1tmp - 1)[0]
            except:
                pass

            
            if not(tmpKey1 is None) and not(tmpKey2 is None):

                x,z,s,u = swapgate(x,z,s,u,tmpKey1, tmpKey2)


            if not(tmpKey1 is None) and tmpKey2 is None:

                qmap[tmpKey1] = P1tmp - 1
 
            
            if tmpKey1 is None and not(tmpKey2 is None):


                qmap[tmpKey2] = P1tmp

            P1tmp -= 1
        

        X,Z,S,U = cxgate(X,Z,S,U, P1tmp, P2)
        x,z,s,u = cxgate(x,z,s,u, p2+1, p2)
        
    
    return X,Z,S,U, x,z,s,u, pivots[0]


def step1(basis, qmap, x,z, X,Z,S,U):

    ntmp = len(basis)//2
    u = QuantumCircuit(ntmp)
    s = GF(numpy.zeros(x.shape[0], dtype=int))

    pivots = []
    for j in range(ntmp):

        a = basis[j]
        b = basis[j + ntmp]
        if a == 1 or b == 1:

            pivots.append(j)

        if a == 0 and b == 1:

            # apply h gate on the main circuit
            X,Z,S,U = hgate(X,Z,S,U, qmap[j])

            # apply h gate on the dummy circuit
            x,z,s,u = hgate(x,z,s,u, j)

        if a == 1 and b == 1:

            # apply s and h on the main circuit
            X,Z,S,U = sgate(X,Z,S,U, qmap[j])
            X,Z,S,U = hgate(X,Z,S,U, qmap[j])

            # apply s and h on the dummy circuit
            x,z,s,u = sgate(x,z,s,u, j)
            x,z,s,u = hgate(x,z,s,u, j)

    return X,Z,S,U, x,z,s,u, pivots, qmap

def step2(qmap, pivots, x,z,s,u, X,Z,S,U):

    # p1 = pivots[0]
    # for p2 in pivots[1:]:

    #     # cnot on the main circuit
    #     X,Z,S,U = cxgate(X,Z,S,U, qmap[p2], qmap[p1])
    #     # cnot on the dummy citcuit
    #     x,z,s,u = cxgate(x,z,s,u, p2,p1)

    
    while len(pivots) != 1:
        newpivots = []
        for j in range(0, len(pivots), 2):
            p1 = pivots[j]
            if j + 1 < len(pivots):
                p2 = pivots[j+1]

                #cnot on the main circuit
                X,Z,S,U = cxgate(X,Z,S,U, qmap[p2], qmap[p1])
                # cnot on the dummy circuit
                x,z,s,u = cxgate(x,z,s,u, p2,p1)

            newpivots.append(p1)
        
        pivots = newpivots

    return X,Z,S,U, x,z,s,u, pivots




def step12(basis, qmap, x,z, X,Z,S,U, connectivity="full"):

    """An implementation that combines both step 1 and step 2. The advantage is that it reduces 
    the number of single qubit gates """

    # process all the cases basis[i] = basis[i + n] = 1
    # process all the cases basis[i] = 0 and basis[i + n] = 1
    # perform the cnot remaining (the original step 2)

    ntmp = len(basis)//2
    u = QuantumCircuit(ntmp)
    s = GF(numpy.zeros(x.shape[0], dtype=int))

    zpivots = []
    xpivots = []
    for j in range(ntmp):

        a = basis[j]
        b = basis[j + ntmp]

        if a == 1 and b == 1:

            # apply s and h on the main circuit
            X,Z,S,U = sgate(X,Z,S,U, qmap[j])
            X,Z,S,U = hgate(X,Z,S,U, qmap[j])

            # apply s and h on the dummy circuit
            x,z,s,u = sgate(x,z,s,u, j)
            x,z,s,u = hgate(x,z,s,u, j)

            xpivots.append(j)
        
        if a == 1 and b == 0:
            xpivots.append(j)
        
        if a == 0 and b == 1:
            zpivots.append(j)


    #print(xpivots)
    if len(zpivots) != 0:
        # zp1 = zpivots[0]
        # for zp in zpivots[1:]:

        #     X,Z,S,U = cxgate(X,Z,S,U, qmap[zp1], qmap[zp])
        #     x,z,s,u = cxgate(x,z,s,u, zp1, zp)

        while len(zpivots) != 1:

            newzpivots = []
            for j in range(0, len(zpivots), 2):

                zp1 = zpivots[j]
                if j + 1 < len(zpivots):
                    zp2 = zpivots[j+1]
                    X,Z,S,U = cxgate(X,Z,S,U, qmap[zp1], qmap[zp2])
                    x,z,s,u = cxgate(x,z,s,u, zp1, zp2)
                newzpivots.append(zp1)
            zpivots = newzpivots
        


        zp1 = zpivots[0]
        X,Z,S,U = hgate(X,Z,S,U, qmap[zp1])
        x,z,s,u = hgate(x,z,s,u,zp1)
        xpivots.append(zp1)

    if len(xpivots) != 0:
        # xp1 = xpivots[0]
        # for xp in xpivots[1:]:

        #     X,Z,S,U = cxgate(X,Z,S,U, qmap[xp], qmap[xp1])
        #     x,z,s,u = cxgate(x,z,s,u, xp, xp1)
        
        while len(xpivots) != 1:

            newxpivots = []
            for j in range(0, len(xpivots), 2):

                xp1 = xpivots[j]
                if j + 1 < len(xpivots):
                    xp2 = xpivots[j+1]

                    X,Z,S,U = cxgate(X,Z,S,U, qmap[xp2], qmap[xp1])
                    x,z,s,u = cxgate(x,z,s,u, xp2, xp1)
                newxpivots.append(xp1)
            xpivots = newxpivots
    
    return X,Z,S,U, x,z,s,u, xpivots

   
def reduce_column(basis, qmap, x,z, X,Z,S,U, reduced_h = False, connectivity="full"):

    """We use the basis vector to diagonalize one column """

    ntmp = len(basis)//2
    u = QuantumCircuit(ntmp)
    s = GF(numpy.zeros(x.shape[0], dtype=int))


    if connectivity == "linear":

        X,Z,S,U, x,z,s,u, pivots, qmap = step1(basis, qmap, x,z, X,Z,S,U)
        
        X,Z,S,U, x,z,s,u, pp = apply_linear_connectivity(pivots, qmap, x,z,s,u, X,Z,S,U)

        return X,Z,S,U,x,z, pp
    
    if reduced_h:

        X,Z,S,U, x,z,s,u, pivots = step12(basis, qmap, x,z, X,Z,S,U)

        return X,Z,S,U, x,z, pivots[0]


    X,Z,S,U, x,z,s,u, pivots, qmap = step1(basis, qmap, x,z, X,Z,S,U)

    X,Z,S,U, x,z,s,u, pivots = step2(qmap, pivots, x,z,s,u, X,Z,S,U)

    return X,Z,S,U, x,z, pivots[0]
        

    # p1 = pivots[0]
    # for p2 in pivots[1:]:

    #     # cnot on the main circuit
    #     X,Z,S,U = cxgate(X,Z,S,U, qmap[p2], qmap[p1])
    #     # cnot on the dummy citcuit
    #     x,z,s,u = cxgate(x,z,s,u, p2,p1)
    

    # return X,Z,S,U, x,z, pivots[0]

def main_diagonalizer(pstrings, reduced_h = False, connectivity="full", optimize=True):

    """This is the main diagonalization loop """

    X,Z,S,Coefs = tableau(pstrings)
    # check that the dimensions are right
    assert X.shape == Z.shape
    assert S.shape == Coefs.shape
    assert X.shape[0] == S.shape[0]

    # initialize
    n = X.shape[1]
    U = QuantumCircuit(n)
    qmap = {j:j for j in range(n)}

    if not X.any():
        return X,Z,S,U

    x = cp.deepcopy(X)
    z = cp.deepcopy(Z)
    x,z, qmap = mask_diagonal_columns(x,z,qmap)
    #x,z = getIndependentPauliStrings(x,z)

    #get independent
    ind_pstrings = getIndependentPauliStrings(x,z)
    x,z,_,_ = tableau(ind_pstrings)

    # start the loop
    for bb in range(n):

        # if there is only one column (qubit), I should diagonalize immediately
        if x.shape[1] == 1:

            if x[:,0].any():

                if not z[:,0].any():

                    # x is 1 and z is 0 ---> The pauli matrix is X
                    X,Z,S,U = hgate(X,Z,S,U, qmap[0])

                else:
                    # x is 1 and z is 1 ---> The pauli matrix is Y
                    X,Z,S,U = sgate(X,Z,S,U, qmap[0])
                    X,Z,S,U = hgate(X,Z,S,U, qmap[0])

            return X,Z,S,U

        nullspace = makeTableauMatrix(x,z).null_space()
        if nullspace.shape[0] == 0:
            return X,Z,S,U
        if not optimize:
            #print("Use the Heuristics guaranteeing a vector with weight at most r + 1")
            #basis = bfo.heuristics(x,z)
            print("Not implemented")
        else:
            basis = bfo.brute_force_optimize(nullspace, connectivity=connectivity)
        X,Z,S,U, x,z, pivot = reduce_column(basis, qmap, x,z, X,Z,S,U, reduced_h, connectivity)
        x,z,qmap = mask_column(x,z,pivot,qmap)


    return X, Z, S, U

