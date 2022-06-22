import pandas as pd
import networkx as nx

# The minimum degree a node needs to have unless it gets removed
MINIMUM_DEGREE = 20
# Load the network

G = nx.read_gml(snakemake.input[0])

# Remove nodes with to low of a degree
to_remove = [node for node,degree in dict(G.degree()).items() if degree <= MINIMUM_DEGREE]
G.remove_nodes_from(to_remove)

# Also remove nodes that do not have any tweets, i.e., they don't
# even have a '0' in the antivax data
to_remove = []
for node in G.nodes:
   ndict = G.nodes[node]
   if ndict.get('numfake') is None:
        to_remove.append(node)

G.remove_nodes_from(to_remove)

nx.write_gml(G, snakemake.output[0])
