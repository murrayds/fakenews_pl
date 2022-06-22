import pandas as pd
import networkx as nx

# The minimum degree a node needs to have unless it gets removed
MINIMUM_DEGREE = 10
# Load the network

G = nx.read_gml(snakemake.input[0])

to_remove = []
for node in G.nodes:
   ndict = G.nodes[node]
   if int(ndict.get('numfake')) < 1:
        to_remove.append(node)

G.remove_nodes_from(to_remove)

nx.write_gml(G, snakemake.output[0])
