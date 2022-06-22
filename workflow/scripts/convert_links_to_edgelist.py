import pandas as pd
import networkx as nx
import json

with open(snakemake.input[0]) as f:
    df = json.load(f)


all_links = list(df['links'].keys())
records = []
for link in all_links:
    spl = link.split(";")
    record = {
        "Source": spl[0],
        "Target": spl[1],
    }
    records.append(record)


links_df = pd.DataFrame(records)
links_df.to_csv(snakemake.output[0])
