from pstring import *

def pauli_to_numerical(pauli):

    """convert a Pauli string from I,X,Y,Z notation to integer notation
    For example 'IXYZ' is converted to '123' """

    res = ""
    for p in pauli:
        if p == "I":
            res += "0"
        if p == "X":
            res += "1"
        if p == "Y":
            res += "2"
        if p == "Z":
            res += "3"

    return res

def numerical_to_pauli(pauli):

    """inverse of the function above """

    res = ""
    for p in pauli:
        if int(p) == 0:
            res += "I"
        elif int(p) == 1:
            res += "X"
        elif int(p) == 2:
            res += "Y"
        elif int(p) == 3:
            res += "Z"
    
    return res

def tableau_to_pstrings(X,Z):

    """ convert tableau to Pauli strings """

    result = []
    assert X.shape == Z.shape
    m,n = X.shape

    for i in range(m):
        pauli = ""
        
        for j in range(n):
            
            if int(X[i,j]) == 0 and int(Z[i,j]) == 0:

                pauli += "I"
            
            if int(X[i,j]) == 0 and int(Z[i,j]) == 1:

                pauli += "Z"
            
            if int(X[i,j]) == 1 and int(Z[i,j]) == 0:

                pauli += "X"
            
            if int(X[i,j]) == 1 and int(Z[i,j]) == 1:

                pauli += "Y"
        
        result.append(pstring(pauli_to_numerical(pauli), 1))
    
    return result

def read_hamiltonian(filepath):
    """Reads the Hamiltonian
    filepath : The path to the file containing the Hamiltonian"""

    pstrings = []
    with open(filepath, "r") as f:
        for line in f:
            tmp = line.split(",")
            pauli = tmp[0]
            coef = tmp[1]
            try:
                coef = float(coef)
            except:
                print("The coefficient of the Pauli strings must be numerical.")
                print("Check the coefficients from the Hamiltonian file")
                return

            if pauli != "I"*len(pauli):
                pstrings.append(pstring(pauli_to_numerical(pauli),coef))

    return pstrings
