import pandas as pd
import networkx as nx

print("Loading edgelist")
# Get the original network edge list
edgelist = pd.read_csv(snakemake.input[0],
                       sep='\s+',
                       skiprows = 2,
                       header=None,
                       names=['Source', 'Target'])

print("Loading hashed node names")
# Get the hashed nodes
with open(snakemake.input[1], 'r') as f:
    hashed = f.read().splitlines()

print("Relabeling nodes with hashed labels")
def match_hash(x):
    return [hashed[x.Source - 1], hashed[x.Target - 1]]

new_edgelist = edgelist.apply(match_hash, axis = 1)

# Convert to the correct format of dataframe
new_edgelist_df = pd.DataFrame(list(new_edgelist), columns = ["Source", "Target"])

print("Converting to network object")
# Now convert to a networkx format
G = nx.from_pandas_edgelist(new_edgelist_df, source = "Source", target = "Target")


print(len(G.nodes()), " nodes in network")
print(G.number_of_edges(), " edges in network")
print("Assigning metadata")

#
# Now, we add the metadata we have to the network
#

# First, add how much antivax info they share
antivax = pd.read_csv(snakemake.input[2], sep = "\t")

attrib = {}
for index, row in antivax.iterrows():
    attrib[row['userid']] = {"numfake": row['num_fake']}

nx.set_node_attributes(G, attrib)


# Add the verified status of the user
verified = pd.read_csv(snakemake.input[3], sep = "\t")
attrib = {}
for index, row in verified.iterrows():
    attrib[row['userid']] = {"verified": row['verified']}

nx.set_node_attributes(G, attrib)

# Finally, add the remaining demographic characteristics
demographics = pd.read_csv(snakemake.input[4], sep = "\t")

attrib = {}
for index, row in demographics.iterrows():
    attrib[row['userid']] = {"sex": row['sex'],
                             "party": row['party_reg'],
                             "race": row['race'],
                             "age": row['age']}

nx.set_node_attributes(G, attrib)

# Save the file
print("Saving network to file")
nx.write_gml(G, snakemake.output[0])
