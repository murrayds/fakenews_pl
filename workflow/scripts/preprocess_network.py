import networkx as nx
import json
import pandas as pd

G = nx.read_gml(snakemake.input[0])

print('Original data: ', nx.info(G))

components = list(G.subgraph(c) for c in nx.connected_components(G))
components.sort(key=lambda c: c.size(), reverse=True)

print('--------------')
print('Number of components: ', len(components))

# Get the giant component
G_new = components[0]

print('--------------')
print('Giant connected component: ', nx.info(G_new))
print('Is connected: ', nx.is_connected(G_new))

#convert the data into a json file format - as an input to NeuraLayout
G_json = nx.readwrite.json_graph.node_link_data(G_new)

with open(snakemake.output[0], "w") as write_file:
    json.dump(G_json, write_file)
