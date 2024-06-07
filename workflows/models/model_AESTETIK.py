#!/usr/bin/env python
# coding: utf-8

# # Proof-of-concept model architecture

# ### 1. Import packages

# In[ ]:


import os
import sys

# In[ ]:

import warnings
warnings.filterwarnings('ignore')
import sys

# In[ ]:


import logging
# Configure the logging module
logging.basicConfig(level=logging.INFO)  # Set the desired logging level
logging.getLogger("pyvips").setLevel(logging.CRITICAL)


# In[ ]:


import pandas as pd
import scanpy as sc
import numpy as np
import yaml
import json


# In[ ]:


from aestetik.utils.utils_transcriptomics import preprocess_adata
from aestetik import AESTETIK

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
parameter_path = f"{out_folder}/{cross_validation_combination}/{model}_fine_tune/parameters/{param_name}.yaml"
cluster_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/clusters/model-{sample}-{param_name}.csv"
latent_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/latent/model-{sample}-{param_name}.csv"
loss_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/loss/model-{sample}-{param_name}.npy"


# In[ ]:


# my parameters 
with open(parameter_path, "r") as stream:
    parameters = yaml.safe_load(stream)


# In[ ]:


spot_diameter_fullres = round(json.load(open(json_path))["spot_diameter_fullres"])


# In[ ]:


dot_size = float(json.load(open(json_path))["dot_size"])


# In[ ]:


img_features_path = f"{out_folder}/data/image_features/{sample}_inception.npy"


# ### 4. Spot-expression data preprocessing

# In[ ]:


adata = sc.read(adata_in)
adata = preprocess_adata(adata)
n_components = adata.obsm["X_pca"].shape[1]

# In[ ]:


adata.obsm["X_pca_transcriptomics"] = adata.obsm["X_pca"]


# In[ ]:

img_features = np.load(img_features_path)
adata.obsm["X_pca_morphology"] = img_features[:,0:n_components]


# ### 5. Load and train model


# In[ ]:


model = AESTETIK(adata, 
                 nCluster=adata.obs.ground_truth.unique().size, 
                 img_path=img_path,
                 spot_diameter_fullres=spot_diameter_fullres,
                 **parameters)


# In[ ]:


model.prepare_input_for_model()


# In[ ]:


model.train()


# ### 6. Compute embeddings and clustering

# In[ ]:


model.compute_spot_representations()


# ### 7. Save clusters

# In[ ]:


cluster_out = adata.obs[[f"AESTETIK_cluster"]].reset_index().rename({'index': "spot_barcode", f'AESTETIK_cluster': param_name}, axis=1)
cluster_out


# In[ ]:


cluster_out.to_csv(cluster_out_path, index=False)
# save latent space
latent_space = pd.DataFrame(adata.obsm["AESTETIK"], index=adata.obs.index)
latent_space.to_csv(latent_out_path)
# save loss
np.save(loss_out_path, model.losses)