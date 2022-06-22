import networkx as nx
import json
import pandas as pd
import numpy as np


with open(snakemake.input[0]) as f:
    df = json.load(f)

# Populate the list of starting & ending coordinates for edges, based on information
# stored in the node dataframe
node_ids = list(df['links'].keys())
records = []
for i in range(len(node_ids)):
    node_id = node_ids[i]

    xstart = df['links'][node_id]['points'][0][0]
    ystart = df['links'][node_id]['points'][0][1]
    zstart = df['links'][node_id]['points'][0][2]

    xend = df['links'][node_id]['points'][1][0]
    yend = df['links'][node_id]['points'][1][1]
    zend = df['links'][node_id]['points'][1][2]

    record = {
        'xstart': xstart,
        'ystart': ystart,
        'zstart': zstart,
        'xend': xend,
        'yend': yend,
        'zend': zend
    }

    records.append(record)


# Now, calculate the direction & distance vectors for each link
links = []
for r in records:
    start = np.array((r['xstart'], r['ystart'], r['zstart']))
    end = np.array((r['xend'], r['yend'], r['zend']))
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]

    dist = np.sqrt(np.square(dx) + np.square(dy) + np.square(dz))
    dirvec = np.array([dx / dist, dy / dist, dz / dist])
    link = {
        'x': start[0],
        'y': start[1],
        'z': start[2],
        'dist': dist,
        'dirvec_x': dirvec[0],
        'dirvec_y': dirvec[1],
        'dirvec_z': dirvec[2],
        'xend': end[0],
        'yend': end[1],
        'zend': end[2],
    }
    links.append(link)

# Convert to a pandas dataframe and save the output
links_df = pd.DataFrame(links)
links_df = links_df.sort_values(by = ['x'])
links_df.to_csv(snakemake.output[0], index = False)
