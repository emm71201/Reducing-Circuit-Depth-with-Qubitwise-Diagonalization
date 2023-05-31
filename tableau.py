import numpy
import galois
## all operations are in GF2 except the coefficients of the pauli strings
GF = galois.GF(2)

def tableau(pstrs):

    """Convert a list of Pauli strings to Tableau representation
    Input: a list of pauli strings objects
    Return: X,Z,S, and Coeficients"""

    n = len(pstrs[0].string)
    k = len(pstrs)

    X = numpy.zeros(shape = (k, n),dtype=int)
    Z = numpy.empty(shape = (k, n),dtype=int)
    S = numpy.empty(shape = (k), dtype=int)

    Coefs = numpy.empty(shape = (k))

    for i in range(k):

        p = pstrs[i]

        tmpx = numpy.empty(shape=(n),dtype=int)
        tmpy = numpy.empty(shape=(n), dtype = int)
        tmpz = numpy.empty(shape=(n), dtype = int)

        for j in range(n):

            if p.string[j]=="0":
                tmpx[j]=0
                tmpz[j]=0

            elif p.string[j] == "1":

                tmpx[j] = 1
                tmpz[j]= 0

            elif p.string[j] == "2":

                tmpx[j] = 1
                tmpz[j] = 1

            elif p.string[j] == "3":

                tmpx[j] = 0
                tmpz[j] = 1

        if p.coef >= 0:
            S[i] = 0

        else:
            S[i] = 1

        Coefs[i] = abs(p.coef)

        X[i] = GF(tmpx)
        Z[i] = GF(tmpz)

    return GF(X), GF(Z), GF(S), Coefs
