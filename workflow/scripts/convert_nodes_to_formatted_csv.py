import networkx as nx
import json
import pandas as pd
import numpy as np


with open(snakemake.input[0]) as f:
    df = json.load(f)


# Iterate over the nodes in the network, extracting the key information, and
# output to a csv
node_records = []
for i in range(len(df['nodes']['positions'])):

    record = {
        'x': df['nodes']['positions'][i][0],
        'y': df['nodes']['positions'][i][1],
        'z': df['nodes']['positions'][i][2],
        'numfake': df['nodes']['numfake'][i],
        'verified': df['nodes']['verified'][df['nodes']['labels'][i]],
        'userid': df['nodes']['labels'][i]
    }
    node_records.append(record)


# Convert to a dataframe and save to a csv
node_df = pd.DataFrame(node_records)
node_df = node_df.sort_values(by = ['x'])

node_df.to_csv(snakemake.output[0], index=False)
