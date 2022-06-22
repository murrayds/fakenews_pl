import networkx as nx
import pandas as pd
from collections import Counter
import numpy as np
import json

txt_file = open(snakemake.input[0], "r")
file_content = txt_file.read().splitlines()

mapping = {}
index = 0
for key in list(Counter(file_content).keys()):
    mapping[key] = index
    index = index + 1


index = 0
for i in range(len(file_content)):
    file_content[i] = mapping[file_content[i]]

G = nx.read_gml(snakemake.input[1])

print("Number of nodes: ", len(G.nodes()))
print("Number of edges: ", G.number_of_edges())

attrib = {}
comm = {}
counter = 0
for n in G.nodes():
    c = file_content[counter]
    attrib[n] = {"community": c}
    if c not in comm:
        comm[c] = []
    comm[c].append(n)
    counter += 1

nx.set_node_attributes(G, attrib)

M = len(np.unique(file_content))

blocks = np.zeros(shape = (M, M))

for edge in G.edges():
    source = G.nodes[edge[0]]['community']
    target = G.nodes[edge[1]]['community']
    blocks[source][target] += 1
    blocks[target][source] += 1

CoarseG = nx.from_numpy_matrix(blocks)
nx.write_gml(CoarseG, snakemake.output[0])

# Now save the community mapping, which is needed later
with open(snakemake.output[1], "w") as f:
    comm_json = json.dumps(comm)
    # write json object to file
    f.write(comm_json)
