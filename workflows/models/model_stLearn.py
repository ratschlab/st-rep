from pathlib import Path
import stlearn as st
import scanpy as sc
import tempfile
import random
import json
import yaml
import sys

SEED = 2023
random.seed(SEED)

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
    
n_components = int(parameters["n_components"])
del parameters["n_components"]
cnn_base = parameters["cnn_base"]
del parameters["cnn_base"]


spot_diameter_fullres = round(json.load(open(json_path))["spot_diameter_fullres"])

org_adata = sc.read(adata_in)
adata = st.create_stlearn(count = org_adata.X, spatial = org_adata.obs[["y_pixel", "x_pixel"]].rename(columns={"x_pixel":"imagerow", "y_pixel": "imagecol"}), library_id="Sample_test", scale=1,
                 image_path = img_path)

adata.obs["ground_truth"] = org_adata.obs.ground_truth.values
adata.obs["array_row"] = org_adata.obs["x_array"].values
adata.obs["array_col"] = org_adata.obs["y_array"].values
adata.obs["index"] = org_adata.obs.index
adata.obs.set_index("index", inplace=True)

# pre-processing for gene count table
st.pp.filter_genes(adata,min_cells=1)
st.pp.normalize_total(adata)
st.pp.log1p(adata)

# run PCA for gene expression data
st.em.run_pca(adata, n_comps=n_components)


temp_dir = tempfile.gettempdir()
TILE_PATH = Path(f"{temp_dir}/tiles")
TILE_PATH.mkdir(parents=True, exist_ok=True)


# pre-processing for spot image
st.pp.tiling(adata, TILE_PATH, crop_size=spot_diameter_fullres)

# this step uses deep learning model to extract high-level features from tile images
# may need few minutes to be completed
st.pp.extract_feature(adata, cnn_base=cnn_base, n_components=n_components)


st.spatial.SME.SME_normalize(adata, use_data="raw", **parameters)

adata.X = adata.obsm['raw_SME_normalized']

st.pp.scale(adata)
st.em.run_pca(adata, n_comps=n_components)

st.tl.clustering.kmeans(adata, n_clusters=adata.obs.ground_truth.unique().size, use_data="X_pca", key_added="stlearn_cluster")


cluster_out = adata.obs[["stlearn_cluster"]].reset_index().rename({'index': "spot_barcode"}, axis=1)


cluster_out.to_csv(cluster_out_path, index=False)

