import numpy
from qiskit import *
from tableau import *
from pstring import *
import galois
GF = galois.GF(2)

def hgate(X,Z,S, U, a):

    xtmp = numpy.copy(X[:,a])
    ztmp = numpy.copy(Z[:,a])

    S = (S + X[:,a] * Z[:,a])

    X[:,a] = ztmp

    Z[:,a] = xtmp

    U.h(a)

    return X,Z,S, U

def sgate(X,Z,S, U, a):

    S = (S + X[:,a] * Z[:,a])

    tmp = (Z[:,a] + X[:,a])

    Z[:,a] = tmp

    U.s(a)

    return X, Z, S, U

def cxgate(X,Z,S, U, a, b):

    ones = GF(numpy.ones(len(X), dtype=int))

    stmp = S + X[:,a] * Z[:,b] * (X[:,b] + Z[:,a] + ones)

    S = stmp

    ztmp = (Z[:,a] + Z[:,b])

    Z[:,a] = ztmp

    xtmp = (X[:,b] + X[:,a])

    X[:,b] = xtmp

    U.cx(a,b)

    return X, Z, S, U

def czgate(X,Z,S,U,a,b):

    X,Z,S,U = hgate(X,Z,S,U,b)
    X,Z,S,U = cxgate(X,Z,S,U,a,b)
    X,Z,S,U = hgate(X,Z,S,U,b)

    return X, Z, S, U

def swapgate(X,Z,S,U,a,b):

    X[:,[a,b]] = X[:,[b,a]]
    Z[:,[a,b]] = Z[:,[b,a]]

    U.swap(a,b)

    return X, Z, S, U


def swaprows(X,Z,S,a,b):

    X[[a,b]] = X[[b,a]]
    Z[[a,b]] = Z[[b,a]]
    S[[a,b]] = S[[b,a]]

    return X,Z,S

def swapcolumns(X,Z,S,U,a,b):

    X[:,[a,b]] = X[:,[b,a]]
    Z[:,[a,b]] = Z[:,[b,a]]

    return X,Z,S,U


# In[144]:


#### other helper functions

def makeTableauMatrix(X,Z):
    """The tableau matrix is defined as [X,Z] """

    return GF(numpy.concatenate((X,Z), axis=1,dtype=int))

def rank(M):
    """ return the rank of the rank of the Matrix M """

    return numpy.linalg.matrix_rank(M)

def tableauMatrixToPauliStrings(T,S,Coefs):

    T = numpy.array(T)
    S = numpy.array(S)
    Coefs = numpy.array(Coefs)

    numrows = T.shape[0]
    numcols = T.shape[1]//2
    pstrs = []
    for r in range(numrows):

        row = T[r]
        coef = Coefs[r] * (-1)**S[r]

        strform = ""

        for bit in range(numcols):

            if row[bit] == 0:

                if row[bit + numcols] == 0:

                    strform += "0"

                if row[bit + numcols] == 1:
                    strform += "3"

            if row[bit] == 1:

                if row[bit + numcols] == 0:

                    strform += "1"

                if row[bit + numcols] == 1:
                    strform += "2"

        pstrs.append(pstring(strform, coef))

    return pstrs

def getIndependentPauliStrings(X,Z):

    """ return the independent pauli strings
    Their coefficients are 1 """

    rowspace = makeTableauMatrix(X,Z).row_space()
    numrows = rowspace.shape[0]

    row_reduced = []
    for row in rowspace:
        if row.any():
            row_reduced.append(row)
    
    rowspace = GF(numpy.array(row_reduced))
    numrows = rowspace.shape[0]

    S = numpy.zeros(numrows, dtype=int)
    Coefs = numpy.ones(numrows, dtype=int)


    return tableauMatrixToPauliStrings(rowspace,S,Coefs)