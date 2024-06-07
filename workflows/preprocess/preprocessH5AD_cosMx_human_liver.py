from aestetik.utils.utils_transcriptomics import preprocess_adata
from sklearn.metrics.cluster import adjusted_rand_score
from tqdm import tqdm
from PIL import Image
import pandas as pd
import scanpy as sc
import numpy as np
import pyvips
import json
import sys
import os 
import glob

import sys
sys.path.append('../')
from src.preprocess_utils.preprocess_image import get_low_res_image

sample = str(sys.argv[1])
downsample_factor = int(sys.argv[2])
out_folder = str(sys.argv[3])
print(sample)

adata_in = f"data/{sample}.h5ad"
img_path = f"data/{sample}_pca_rgb.jpg"
adata_out = f"{out_folder}/data/h5ad/{sample}.h5ad"

adata = sc.read(adata_in)


adata.obs["ground_truth"] = adata.obs.cellType.values
adata.obs["index"] = adata.obs.cell_id
adata.obs.set_index("index", inplace=True)


adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.var["gene_ids"] = adata.var.genename


adata.obs["x_array"]=adata.obs["x_array"]
adata.obs["y_array"]=adata.obs["y_array"]
adata.obs["x_pixel"]=adata.obs["x_pixel"].astype(np.float32)
adata.obs["y_pixel"]=adata.obs["y_pixel"].astype(np.float32)

img = get_low_res_image(img_path, downsample_factor=downsample_factor)


adata.obsm['spatial'] = adata.obs[["x_pixel", "y_pixel"]].values
# adjust coordinates to new image dimensions
adata.obsm['spatial'] = adata.obsm['spatial']/downsample_factor
# create 'spatial' entries
adata.uns['spatial'] = dict()
adata.uns['spatial']['library_id'] = dict()
adata.uns['spatial']['library_id']['images'] = dict()
adata.uns['spatial']['library_id']['images']['hires'] = img


selected_col = ["cell_id", "Run_Tissue_name", 
                "x_slide_mm", "y_slide_mm", 
                "cellType", "niche", "ground_truth",
                "x_array", "y_array", "x_pixel", "y_pixel"
               ]
adata.obs = adata.obs[selected_col]

adata.write_h5ad(adata_out)