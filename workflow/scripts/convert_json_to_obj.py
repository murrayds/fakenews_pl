import pyobj_v2 as pj
import os
from glob import glob
from numpy import *

params = {
    'tube_cross' : 3, # cross section of links (tubes)
    'tu_rad' :0.01,#0.5, # link thickness relative to real value
    'rad_n' : 0.5, #  Node radius
    'segs' : 2, # link segments
}

file_dict = {
    'node': os.path.dirname(snakemake.output[0]) + '/',
    'link': os.path.dirname(snakemake.output[1]) + '/',
}

out_name = os.path.basename(snakemake.output[0]).replace("-nodes.obj", "")

pj.save_obj(snakemake.input[0],
            out_name = out_name,
            pth = file_dict,
            **params)
