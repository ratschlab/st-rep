import os


os.environ['R_HOME'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/Rnonchev/lib/R'
os.environ['R_USER'] = '/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/nonchev/lib/python3.9/site-packages/rpy2'


# In[ ]:


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
import torch
import yaml
import json


# In[ ]:


from GraphST import GraphST
import numpy as np
# clustering
from GraphST.utils import clustering
SEED = 2023

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


# In[ ]:


radius = int(parameters["radius"])
del parameters["radius"]


# In[ ]:


adata = sc.read(adata_in)


# ### 5. Load and train model

# In[ ]:


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# In[ ]:


model = GraphST.GraphST(adata, device=device, random_seed=SEED, **parameters)


# In[ ]:


adata = model.train()

clustering(adata, adata.obs.ground_truth.unique().size, radius=radius, method="mclust", refinement=radius > 0)

adata.obs["graphst_cluster"] = adata.obs.domain.values


# ### 7. Save clusters

# In[ ]:


cluster_out = adata.obs[["graphst_cluster"]].reset_index().rename({'index': "spot_barcode", 'graphst_cluster': param_name}, axis=1)
cluster_out.to_csv(cluster_out_path, index=False)



