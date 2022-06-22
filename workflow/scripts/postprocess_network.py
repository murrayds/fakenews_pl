import networkx as nx
import json
import pandas as pd
from networkx.readwrite import json_graph
from numpy import *

def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return json_graph.node_link_graph(js_graph)


#outputs - positions
dp = pd.read_csv(snakemake.input[0]).values
G = read_json_file(snakemake.input[1])

label = [i for i in G.nodes]
label_dict = {  label[i]:i for i in range(len(label))}

verified = {k:v['verified']  for k,v in G.nodes(data = True)}
numfake = [v['numfake'] for k,v in G.nodes(data = True)]
colors = [ log(v['numfake'])+2  for k,v in G.nodes(data = True)]
sizes = [ log(v['numfake'])+2  for k,v in G.nodes(data = True)]

pts = [[dp[i][0], dp[i][1], dp[i][2]] for i in range(len(label))]

label_pts = { i: [dp[i][0], dp[i][1], dp[i][2]] for i in range(len(label))}

end_points = [ [ label_dict[i[0]], label_dict[i[1]] ] for i in G.edges]

link_labels = [ str(i[0]) +';' + str(i[1])  for i in G.edges]

end_position = []
for i in range(len(end_points)):
    end_position += [[label_pts[end_points[i][0]], label_pts[end_points[i][1]]] ]


newlink= dict()
for i in range(len(link_labels)):
    newlink[link_labels[i]] = {'end_points': end_points[i], 'points': end_position[i], 'radius': 0.03, 'weight':1, 'color': '5e6366' }

data = {
    'scale': 1,
    'nodes': {
        'positions': pts,
        'labels': label,
        'numfake': numfake,
        'colors': colors,
        'verified': verified,
        'sizes':sizes,
    },
    'info': {'A': 100, 'k': 10, 'pow': 2, 'n_radius': 1,'links': {'T0': 0, 'ce': 0, 'weighted': False, 'thickness': 0.03, 'lables': link_labels}, 'nodes': {'radius': 1, 'weight': False, 'labels': label, 'colors': colors,'groups': verified, 'sizes':sizes,}},
    'links': newlink,
     }


with open(snakemake.output[0], "w") as write_file:
    json.dump(data, write_file)
