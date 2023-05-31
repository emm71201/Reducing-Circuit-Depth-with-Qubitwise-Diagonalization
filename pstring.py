import numpy as np

def commute(b1, b2):

    """check if two pauli matrices b1 and b2 commute
    The pauli matrices I, X,Y,Z are represented as integers from 0 to 3 respectively
    """

    if b1 == "0" or b2 == "0":
        return True

    if b1 == b2:
        return True

    return False

class pstring:

    """The Pauli string Class """

    def __init__(self, string, coef):

        self.string = string
        self.coef = coef

    def __le__(self, other):

        return self.string <= other.string

    def __lt__(self, other):

        return self.string <= other.string


    def weight(self):
        res = 0
        for c in self.string:
            if int(c) != 0:
                res += 1

        return res 

    def pauli_matrix_form(self):

        """ return the Pauli string in Pauli matrix form """

        result = ""
        for p in self.string:
            if p == "0":
                result += "I"
            elif p == "1":
                result += "X"
            elif p == "2":
                result += "Y"
            else:
                result += "Z"

        return result

    def commute(self, other):

        """Return True if the Pauli string self commutes with other
        Otherwise return False """

        str1 = self.string
        str2 = other.string
        assert len(str1) == len(str2)

        count = 0

        for j in range(len(str1)):

            if str1[j] == str2[j] or str1[j] == '0' or str2[j] == '0':
                pass
            else:
                count += 1

        if count %2 == 0:
            return True

        return False
    
    def to_vector(self):

        n = len(self.string)
        res = np.zeros(2*n, dtype=int)

        for j in range(n):
            if self.string[j] == "1":
                res[j] = 1
            if self.string[j] == "2":
                res[j] = 1
                res[j + n] = 1
            if self.string[j] == "3":
                res[j+n] = 1
        
        return res


    def __str__(self):

        return "{0} * {1}".format(self.coef, self.pauli_matrix_form())
