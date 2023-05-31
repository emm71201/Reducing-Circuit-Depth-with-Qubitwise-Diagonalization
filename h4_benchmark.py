# Author: Edison Murairi
# Date: May. 18th, 2023

#%% importing
import helpers
import diagonalize as dg
from tableau_operations import *
from pstring import *
#%%
def parser(fname, n):

    hamiltonian = {}

    with open(fname, "r") as infile:

        cluster = 0
        hamiltonian[cluster] = []

        for line in infile:

            if line.strip() != "end":

                line = line.strip()
                line += ";"

                j = 0
                pstr = ""
                while line[j] != ";":

                    if line[j] != " ":
                        pstr += line[j]
                    
                    if len(pstr) == n:

                        mypstr = helpers.pauli_to_numerical(pstr)

                        hamiltonian[cluster].append(pstring(mypstr, 1))

                        pstr = ""
                    
                    j += 1
            
            if line.strip() == "end":

                cluster += 1

                hamiltonian[cluster] = []
            
    
    return hamiltonian

h4 = parser("misc/h4.csv", 8)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[0], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[1], reduced_h=False, connectivity="linear")
print(U)
#%%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[2], reduced_h=False, connectivity="linear")
print(U)
#%%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[3], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[4], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[5], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[6], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[7], reduced_h=False, connectivity="linear")
print(U)
# %%
print("with linear connectivity: ")
X,Z,S,U=dg.main_diagonalizer(h4[8], reduced_h=False, connectivity="linear")
print(U)
# %%
