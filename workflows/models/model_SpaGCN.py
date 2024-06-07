import os

# cluster R path
os.environ['R_HOME'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/Rnonchev/lib/R'
os.environ['R_USER'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/nonchev/lib/python3.9/site-packages/rpy2'

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2,40))

import warnings
warnings.filterwarnings('ignore')


import socket
node_name = socket.gethostname()
warnings.warn(node_name)

import sys
from tqdm import tqdm
import pandas as pd
import scanpy as sc
import numpy as np
import random
import torch
import yaml
import json


import SpaGCN as spg
import cv2
import sys

# ### 3. Input parameters and data paths


sample = str(sys.argv[1])
param_name = str(sys.argv[2])
img_format = str(sys.argv[3])
cross_validation_combination = str(sys.argv[4])
model = str(sys.argv[5])
mode = str(sys.argv[6])
out_folder = str(sys.argv[7])


adata_in = f"{out_folder}/data/h5ad/{sample}.h5ad"
img_path = f"{out_folder}/data/image/{sample}.{img_format}"
json_path = f"{out_folder}/data/meta/{sample}.json"
parameter_path = f"{out_folder}/{cross_validation_combination}/{model}_fine_tune/parameters/{param_name}.yaml"
cluster_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/clusters/model-{sample}-{param_name}.csv"

# my parameters 
with open(parameter_path, "r") as stream:
    parameters = yaml.safe_load(stream)
parameters

s=parameters["s"]
histology=parameters["histology"]


del parameters["s"]
del parameters["histology"]



adata = sc.read(adata_in)
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

img=cv2.imread(img_path)


x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].astype(int).tolist()
y_pixel=adata.obs["y_pixel"].astype(int).tolist()

adj=spg.calculate_adj_matrix(x=x_array,y=y_array, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=49, alpha=s, histology=histology)

p=0.5 
#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

n_clusters=adata.obs.ground_truth.unique().size
#Set seed
r_seed=t_seed=n_seed=100
#Seaech for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)


clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, **parameters)
pred, prob=clf.predict()
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=0)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=pred, dis=adj_2d, shape="hexagon")


# ### 7. Save clusters


adata.obs["spagcn_cluster"] = np.array(refined_pred).astype(str)
cluster_out = adata.obs[["spagcn_cluster"]].reset_index().rename({'index': "spot_barcode"}, axis=1)

cluster_out.to_csv(cluster_out_path, index=False)


