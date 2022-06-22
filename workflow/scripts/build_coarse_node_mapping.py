"""
    basic_spring_layout.py

    Performs a really quick & dirty "spring layout" on the network based on the
    communities. Nothing fancy, just meant for testing.
"""

import networkx as nx
import pandas as pd
from collections import Counter
import numpy as np


seed = 1111
comm_scale = 7
total_scale = 15


G = nx.read_gml("fakenews_community_network.gml")

pos = nx.spring_layout(G,
                       dim = 3,
                       scale = total_scale,
                       seed = seed
                      )
