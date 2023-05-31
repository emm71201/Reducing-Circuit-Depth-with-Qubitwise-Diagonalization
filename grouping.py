from pstring import *
#from numpy import *
import networkx as nx
import numpy as np

def common_entry(arr1, arr2):

    """checks if two arrays have at least one common entry"""

    for x in arr1:

        for y in arr2:

            if x == y:

                return True

    return False

def no_commute(pstrings):
    """make a graph of pstrings such that \n
    there is an edge between two pstring if they do not commute"""
    G = nx.Graph()

    for i in range(len(pstrings)):
        p1 = pstrings[i]
        G.add_node(p1)

        for j in range(len(pstrings)):

            p2 = pstrings[j]

            if not p1.commute(p2):
                #G.add_edge(p1.string, p2.string)
                G.add_edge(p1, p2)


    return G

def no_commute_v2(pstrings, filename="nocommute_graph.dimacs"):
    """make a graph of pstrings such that \n
    there is an edge between two pstring if they do not commute"""

    keyTable = {}
    for j in range(len(pstrings)):
        keyTable[pstrings[j].string] = j+1
    np.save("keyTable.npy", keyTable)
    vertexCount = len(pstrings)
    edgeCount = 0

    with open(filename, "w") as f:
        for i in range(len(pstrings)):
            p1 = pstrings[i]

            for j in range(i+1,len(pstrings)):

                p2 = pstrings[j]

                if not p1.commute(p2):
                    f.write("e {0} {1}\n".format(keyTable[p1.string], keyTable[p2.string]))
                    edgeCount += 1

            clear_output(wait=True)
            os.system("clear")
            print("i = {0}".format(i))
        f.write("p edge {0} {1}".format(vertexCount, edgeCount))


def commute_v2(pstrings, filename="commute_graph.dimacs"):
    """make a graph of pstrings such that \n
    there is an edge between two pstring if they commute"""
    #gr = {}
    keyTable = {}
    for j in range(len(pstrings)):
        keyTable[pstrings[j].string] = j+1
    save("keyTable.npy", keyTable)
    vertexCount = len(pstrings)
    edgeCount = 0

    with open(filename, "w") as f:
        for i in range(len(pstrings)):
            p1 = pstrings[i]

            for j in range(i+1,len(pstrings)):

                p2 = pstrings[j]

                if p1.commute(p2):
                    f.write("e {0} {1}\n".format(keyTable[p1.string], keyTable[p2.string]))
                    edgeCount += 1

            clear_output(wait=True)
            os.system("clear")
            print("i = {0}".format(i))
        f.write("p edge {0} {1}".format(vertexCount, edgeCount))

def make_clusters(pstrings, strategy="DSATUR"):
    """Provide clusters of commuting pauli strings by graph coloring
    The default strategy is 'DSATUR' as implemented in networkx"""

    strategies = ['largest_first','random_sequential', 'smallest_last',\
    'independent_set', 'connected_sequential_bfs', 'connected_sequential_dfs', \
    'connected_sequential', 'saturation_largest_first',\
    'DSATUR']

    G = no_commute(pstrings)

    if strategy not in strategies:
        print("The strategy is not available. We proceed with the default, DSATUR")
        strategy="DSATUR"

    colors = nx.coloring.greedy_color(G, strategy=strategy)
    maxcol = max(colors.values())
    result = {col:[] for col in range(maxcol + 1)}
    for pauli in colors.keys():
        col = colors[pauli]
        result[col].append(pauli)

    return result
