"""
    basic_spring_layout.py

    Performs a really quick & dirty "spring layout" on the network based on the
    communities. Nothing fancy, just meant for testing.
"""

import networkx as nx
import pandas as pd
from collections import Counter
import numpy as np
import json

seed = 1111
comm_scale = 7
total_scale = 15

G = nx.read_gml(snakemake.input[0])
CoarseG = nx.read_gml(snakemake.input[1])

# Now load the node mapping
with open(snakemake.input[2]) as f:
    nodemap = json.loads(f.read())

# Load the anti-vax dataframe
antivax = pd.read_csv(snakemake.input[3], sep = "\t")
antivax = antivax.drop('rank', axis = 1)
antivax = antivax.rename(columns = {'num_fake': "numfake"})

# First, we layout the coarse-grained network
pos = nx.spring_layout(CoarseG,
                       dim = 3,
                       scale = total_scale,
                       seed = seed)

# Next, layout each individual network
records = []
for i in range(len(nodemap)):
    # get sub-network
    subnodes = nodemap['{}'.format(i)]

    subnet = nx.induced_subgraph(G, subnodes)

    # get the total number of fake news tweets
    total_numfake = antivax[antivax['userid'].isin(subnodes)].numfake.sum()

    sublayout = nx.spring_layout(subnet, dim = 3, scale = comm_scale, seed = seed)
    records.append(sublayout)


# Calculate the center of mass of the communities
mass = []
for i in range(len(records)):
    masses = np.array(list(records[i].values()))

    center_of_mass = np.average(masses, axis=0)
    mass.append(center_of_mass)

print(pos)
# Now tesselate the points
tesselated = records.copy()
for i in range(len(records)):
    subnet_center = np.array(pos["{}".format(i)])
    nodes_center = np.array(mass[i])

    diff = subnet_center - nodes_center
    for node in records[i]:
        tesselated[i][node] = np.add(records[i][node], diff)

# Now iterate through each row of 'tesselated', creating a single matrix
coords = []
for i in range(len(tesselated)):
    nodes = tesselated[i]
    for node in nodes:
        record = {
            'userid': node,
            'x': nodes[node][0],
            'y': nodes[node][1],
            'z': nodes[node][2],
            'community': i
        }
        coords.append(record)

df = pd.DataFrame(coords)

df_merged = df.merge(antivax, on = "userid", how = "left")

# Now save the final file
df_merged.to_csv(snakemake.output[0])
