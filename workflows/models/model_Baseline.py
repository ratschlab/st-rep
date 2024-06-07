import os
import sys
import warnings
warnings.filterwarnings('ignore')



import socket
node_name = socket.gethostname()

import logging
# Configure the logging module
logging.basicConfig(level=logging.INFO)  # Set the desired logging level
logging.getLogger("pyvips").setLevel(logging.CRITICAL)


import sys
from aestetik.utils.utils_transcriptomics import preprocess_adata
from aestetik.utils.utils_clustering import clustering
from tqdm import tqdm
import pandas as pd
import scanpy as sc
import numpy as np
import random
import torch
import yaml
import json


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
json_path = f"{out_folder}/data/meta/{sample}.json"
img_path = f"{out_folder}/data/image/{sample}.{img_format}"
#parameter_path = f"{out_folder}/{cross_validation_combination}/{model}_fine_tune/parameters/{param_name}.yaml"
cluster_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/clusters_default/model-{sample}-{param_name}.csv"


# In[ ]:


# my parameters 
#with open(parameter_path, "r") as stream:
#    parameters = yaml.safe_load(stream)


# In[ ]:


#spot_diameter_fullres = round(json.load(open(json_path))["spot_diameter_fullres"])


# In[ ]:


#dot_size = float(json.load(open(json_path))["dot_size"])


# In[ ]:


img_features_path = f"{out_folder}/data/image_features/{sample}_inception.npy"


# ### 4. Spot-expression data preprocessing

# In[ ]:


adata = sc.read(adata_in)
adata = preprocess_adata(adata)
n_components = adata.obsm["X_pca"].shape[1]

# In[ ]:


adata.obsm["X_pca_counts"] = adata.obsm["X_pca"]


# In[ ]:

img_features = np.load(img_features_path)
adata.obsm["X_pca_image"] = img_features[:,0:n_components]


# ### 5. Load and train model


# In[ ]:

adata.obsm["X_pca"] = np.concatenate((adata.obsm["X_pca_counts"], adata.obsm["X_pca_image"]), axis=1)

# In[ ]:

clustering(adata, adata.obs.ground_truth.unique().size, "X_pca", "bgm")

# ### 7. Save clusters

# In[ ]:


cluster_out = adata.obs[["X_pca_cluster"]].reset_index().rename({'index': "spot_barcode", 'baseline_cluster': param_name}, axis=1)
cluster_out


# In[ ]:


cluster_out.to_csv(cluster_out_path, index=False)


# In[ ]:




