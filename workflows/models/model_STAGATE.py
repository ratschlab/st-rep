import os


os.environ['R_HOME'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/Rnonchev/lib/R'
os.environ['R_USER'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/nonchev/lib/python3.9/site-packages/rpy2'

import warnings
warnings.filterwarnings('ignore')
import socket
node_name = socket.gethostname()
warnings.warn(node_name)


# In[ ]:

import sys
from tqdm import tqdm
import pandas as pd
import scanpy as sc
import numpy as np
import random
import yaml
import json


# In[ ]:


SEED = 2023


# In[ ]:


import STAGATE


# ### 3. Input parameters and data paths

# In[ ]:


sample = str(sys.argv[1])
param_name = str(sys.argv[2])
img_format = str(sys.argv[3])
cross_validation_combination = str(sys.argv[4])
model = str(sys.argv[5])
mode = str(sys.argv[6])
out_folder = str(sys.argv[7])


# In[ ]:


adata_in = f"{out_folder}/data/h5ad/{sample}.h5ad"
img_path = f"{out_folder}/data/image/{sample}.{img_format}"
json_path = f"{out_folder}/data/meta/{sample}.json"
parameter_path = f"{out_folder}/{cross_validation_combination}/{model}_fine_tune/parameters/{param_name}.yaml"
cluster_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/clusters/model-{sample}-{param_name}.csv"


# In[ ]:


# my parameters 
with open(parameter_path, "r") as stream:
    parameters = yaml.safe_load(stream)
parameters

k_cutoff = parameters["k_cutoff"]
del parameters["k_cutoff"]

adata = sc.read(adata_in)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)



# In[ ]:


STAGATE.Cal_Spatial_Net(adata, k_cutoff=k_cutoff, model="KNN")


STAGATE.Stats_Spatial_Net(adata)

adata = STAGATE.train_STAGATE(adata, alpha=0, random_seed=SEED, **parameters)

sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=adata.obs.ground_truth.unique().size)

adata.obs["stagate_cluster"] = adata.obs["mclust"].astype(str)


cluster_out = adata.obs[["stagate_cluster"]].reset_index().rename({'index': "spot_barcode"}, axis=1)

cluster_out.to_csv(cluster_out_path, index=False)




