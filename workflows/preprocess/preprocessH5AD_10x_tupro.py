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
replicate_new = SAMPLE
replicate = replicate_new if "-2" not in replicate_new else replicate_new.replace("-2", "_2")
replicate = replicate.replace("r", "rep") # we did this due to name length size
SAMPLE = replicate.split("-")[0]
print(SAMPLE)

base_path = f"data/projects2022-digpathways_datasets/10x_tupro"
adata_h5 = f"{base_path}/data_derived/{SAMPLE}/{replicate}/preprocessed.h5ad"
annotation_path = f"data/spot_annotations_kalin/{replicate}.csv"
adata_out = f"{out_folder}/data/h5ad/{replicate_new}.h5ad"
image_info_file = f"{base_path}/10x_image-paths_resolution.csv"

df = pd.read_csv(image_info_file)
replicate_info = df[df.replicate == replicate]
img_path = replicate_info.image_path.values[0]

resolution = replicate_info.resolution.values[0]
tile_size=160

if resolution == 0.3448:
    tile_size_sample = tile_size
elif resolution == 0.1736:
    tile_size_sample = int(tile_size * (0.3448 / 0.1736))
    downsample_factor = 15

adata = sc.read(adata_h5)

adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_counts=1)

adata.obs["index"] = adata.obs.cell_id
adata.obs.set_index("index", inplace=True)

adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.var["gene_ids"] = adata.var.ens

img = get_low_res_image(img_path, downsample_factor)
# make sure squidpy (or scanpy) can find the relevant var names
adata.var_names = adata.var['symbol'].to_list()
# adjust coordinates to new image dimensions
adata.obsm['spatial'] = adata.obsm['spatial']/downsample_factor
# flip X and Y axes
adata.obsm['spatial'] = np.flip(adata.obsm['spatial'], 1)
# create 'spatial' entries
adata.uns['spatial'] = dict()
adata.uns['spatial']['library_id'] = dict()
adata.uns['spatial']['library_id']['images'] = dict()
adata.uns['spatial']['library_id']['images']['hires'] = img

adata.obs["x_array"]=adata.obs["spot_X"]
adata.obs["y_array"]=adata.obs["spot_Y"]
adata.obs["x_pixel"]=adata.obs["X"]
adata.obs["y_pixel"]=adata.obs["Y"]

print("Data loaded...")

## preprocess tutorial scanpy

sc.pp.filter_genes(adata, min_counts=1)

## add labels

annotations = pd.read_csv(annotation_path)

adata.obs = adata.obs.merge(annotations).set_index("cell_id")
adata.obs["index"] = adata.obs.index
adata.obs.set_index("index", inplace=True)
adata.obs["ground_truth"] = adata.obs.label.values


# remove bad labels

adata = adata[adata.obs.label != "Whitespace",:]
adata = adata[adata.obs.label != "Mixed",:]

clusters = adata.obs.ground_truth.value_counts()[adata.obs.ground_truth.value_counts() / len(adata) >= 0.01].index.values
adata = adata[adata.obs.ground_truth.isin(clusters),:] # remove rare clusters 1% (noise)
print("ground_truth loaded...")

adata.write_h5ad(adata_out)