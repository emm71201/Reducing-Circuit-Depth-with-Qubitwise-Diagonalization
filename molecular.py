# Author: Edison Murairi
# Date: May. 18th, 2023

#%% importing
import helpers
import diagonalize as dg
from tableau_operations import *
from grouping import *
# %%
# names = [ "BH3", "CH4", "H2O","HF", "LiH"]
names = ["HeH+", "LiH", "BeH2", "NH3", "BH3"]
#%%
def process(molecule_name):

    hamiltonian = helpers.read_hamiltonian("molecular_data/" + molecule_name + ".csv")
    n = len(hamiltonian[0].string)

    print(f"Molecule processing: {molecule_name}")

    with open("molecular_data/" + molecule_name + "_results.csv", "a") as ofile:

        ofile.write("cluster,n,N,cx,depth\n")

        clusters = make_clusters(hamiltonian, strategy="independent_set")
        nclusters = max(clusters.keys())

        rs = []
        
        print("\t Diagonalizing clusters:")

        for cluster in clusters:

            #X,Z,S,U = dg.main_diagonalizer(clusters[cluster], optimize=True, connectivity="full")

            rind = dg.main_diagonalizer(clusters[cluster], optimize=True, connectivity="full")
            rs.append(rind)

            # if X.any():
            #     print("Did not work")
            #     return
            
            # try:
            #     cxcount = U.count_ops()["cx"]
            # except KeyError:
            #     cxcount = 0

            # ofile.write(f"{cluster},{n},{len(clusters[cluster])},{cxcount},{U.depth()}\n")

            # print(f"\t\t Finished diagonalizing cluster: {cluster} out of {nclusters}")
        
        for r in rs[:-1]:
            ofile.write(f"{r},")
        ofile.write(f"{rs[-1]}")


    return 

# for molecule in names:
#     process(molecule)

with open("molecular_data/" + "ngenerators.csv", "w") as ofile:
    for molecule in names:

        hamiltonian = helpers.read_hamiltonian("molecular_data/" + molecule + ".csv")
        n = len(hamiltonian[0].string)

        ostring = f"{molecule},"

        print(f"Molecule processing: {molecule}")

        clusters = make_clusters(hamiltonian, strategy="independent_set")
        nclusters = max(clusters.keys())

        for cluster in clusters:

            #X,Z,S,U = dg.main_diagonalizer(clusters[cluster], optimize=True, connectivity="full")
            print("\t Diagonalizing clusters:")
            rind = dg.main_diagonalizer(clusters[cluster], optimize=True, connectivity="full")

            ostring += f"{rind},"
        
        ostring = ostring[:len(ostring)-1] + "\n"
        ofile.write(ostring)
        

# %% Analysis
import pandas as pd
dfs = {}
for molecule in names:
    try:
        dfs[molecule] = pd.read_csv("molecular_data/" + molecule + "_results.csv")
    except:
        print(molecule)

# %%
generators = {}
with open("molecular_data/ngenerators.csv", "r") as infile:

    for line in infile:

        data = line.strip().split(",")
        generators[data[0]] = sum([int(item) for item in data[1:]])/len([int(item) for item in data[1:]])
generators
#%%
number_clusters = {}
for molecule in names:

    with open("molecular_data/" + molecule + "_results.csv", "r") as ofile:

        number_clusters[molecule] = len(ofile.readlines()) - 1

number_clusters

#%%
agg_data = {"Molecule":[], "n":[], "r":[], "N":[], "kappa":[], "CNOT":[], "CNOT_STD":[], "Depth":[], "Depth_STD":[]}
for molecule in names:

    data = dfs[molecule].mean(axis=0)
    agg_data["Molecule"].append(molecule)
    agg_data["n"].append(data["n"])
    agg_data["r"].append(generators[molecule])
    agg_data["N"].append(dfs[molecule]["N"].sum())
    agg_data["kappa"].append(number_clusters[molecule])
    agg_data["CNOT"].append(data["cx"])
    agg_data["CNOT_STD"].append(dfs[molecule]["cx"].std(ddof=0))
    agg_data["Depth"].append(data["depth"])
    agg_data["Depth_STD"].append(dfs[molecule]["depth"].std(ddof=0))
# %%
print(pd.DataFrame(agg_data).round(decimals=2).to_latex(index=False))
#%%
agg_data = pd.DataFrame(agg_data)
agg_data
# %%
agg_data["n"]*np.log2(agg_data["r"])
# %%
