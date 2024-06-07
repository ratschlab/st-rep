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

SAMPLE = str(sys.argv[1])
downsample_factor = int(sys.argv[2])
out_folder = str(sys.argv[3])
print(SAMPLE)

base_path = f"data/human_breast_cancers/"
adata_h5 = f"{base_path}/filtered_count_matrices/{SAMPLE}_filtered_count_matrix/"
adata_out = f"{out_folder}/data/h5ad/{SAMPLE}.h5ad"
img_path = f"{base_path}/spatial/{SAMPLE}_spatial/tissue_hires_image.png"
spatial_info = f"{base_path}/spatial/{SAMPLE}_spatial/tissue_positions_list.csv"
json_file = f"{base_path}/spatial/{SAMPLE}_spatial/scalefactors_json.json"
ground_truth_cluster_path = f"{base_path}/metadata/{SAMPLE}_metadata.csv"

adata = sc.read_10x_mtx(adata_h5)


## preprocess tutorial scanpy

sc.pp.filter_genes(adata, min_counts=1)


spatial=pd.read_csv(spatial_info,sep=",",header=None,na_filter=False,index_col=0) 
json_info = json.load(open(json_file))


adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]

adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=adata.obs["x4"] * json_info["tissue_hires_scalef"]
adata.obs["y_pixel"]=adata.obs["x5"] * json_info["tissue_hires_scalef"]



ground_truth = pd.read_csv(ground_truth_cluster_path, index_col=0)
ground_truth["Barcode"] = ground_truth.index
adata.obs.loc[ground_truth.Barcode, "ground_truth"] = ground_truth.Classification.values
adata = adata[~adata.obs.ground_truth.isna(), :]
adata = adata[adata.obs.ground_truth != "Artefact", :]
adata.obs.ground_truth = adata.obs.ground_truth.astype(str)


#Select captured samples
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")

img = get_low_res_image(img_path, downsample_factor=downsample_factor)

adata.obsm['spatial'] = adata.obs[["y_pixel", "x_pixel"]].values
# adjust coordinates to new image dimensions
adata.obsm['spatial'] = adata.obsm['spatial']/downsample_factor
# create 'spatial' entries
adata.uns['spatial'] = dict()
adata.uns['spatial']['library_id'] = dict()
adata.uns['spatial']['library_id']['images'] = dict()
adata.uns['spatial']['library_id']['images']['hires'] = img


adata.write_h5ad(adata_out)