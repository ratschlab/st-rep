from aestetik.utils.utils_morphology import extract_morphology_embeddings
from torchvision.models import inception_v3, Inception_V3_Weights
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from torchvision import transforms
import tensorflow.compat.v1 as tf
import muse_sc as muse
from tqdm import tqdm
import scanpy as sc
import numpy as np
import torch
import scipy
import json
import yaml
import sys


config = tf.ConfigProto()
config.gpu_options.allow_growth = True

SEED = 2023
np.random.seed(SEED)

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
img_features_path = f"{out_folder}/data/image_features/{sample}_inception.npy"
parameter_path = f"{out_folder}/{cross_validation_combination}/{model}_fine_tune/parameters/{param_name}.yaml"
cluster_out_path = f"{out_folder}/{cross_validation_combination}/{model}_{mode}/clusters/model-{sample}-{param_name}.csv"


# my parameters 
with open(parameter_path, "r") as stream:
    parameters = yaml.safe_load(stream)
    
n_components = int(parameters["n_components"])
del parameters["n_components"]


spot_diameter_fullres = round(json.load(open(json_path))["spot_diameter_fullres"])
spot_diameter_fullres

adata = sc.read(adata_in)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=500)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
data_a = adata.X[ :, adata.var.highly_variable.values]
# check sparse
if type(data_a) == scipy.sparse._csr.csr_matrix:
    data_a = data_a.toarray()
    
    
weights = Inception_V3_Weights.DEFAULT
image_model_children = inception_v3(weights=weights)
image_model_children.fc = torch.nn.Identity()

image_model_children.eval()    
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
feature_dim = 2048

preprocess = transforms.Compose([
                transforms.ToTensor(),
                weights.transforms(antialias=True),
            ])


data_b = extract_morphology_embeddings(img_path, 
                                      image_model_children,
                                      x_pixel=adata.obs.y_pixel, 
                                      y_pixel=adata.obs.x_pixel, 
                                      spot_diameter=spot_diameter_fullres,
                                      device=device,
                                      feature_dim=feature_dim,
                                      preprocess=preprocess,
                                      apply_pca=False)

view_a_feature = PCA(n_components=n_components).fit_transform(data_a)
view_b_feature = PCA(n_components=n_components).fit_transform(data_b)

view_a_label = KMeans(n_clusters=adata.obs.ground_truth.unique().size, random_state=SEED, n_init="auto").fit(view_a_feature).predict(view_a_feature)
view_b_label = KMeans(n_clusters=adata.obs.ground_truth.unique().size, random_state=SEED, n_init="auto").fit(view_b_feature).predict(view_b_feature)


muse_feature, reconstruct_x, reconstruct_y, \
latent_x, latent_y = muse.muse_fit_predict(data_a,
                                           data_b,
                                           view_a_label,
                                           view_b_label,
                                           **parameters)

muse_cluster = KMeans(n_clusters=adata.obs.ground_truth.unique().size, random_state=SEED, n_init="auto").fit(muse_feature).predict(muse_feature)

adata.obs["muse_cluster"] = muse_cluster.astype(str)



cluster_out = adata.obs[["muse_cluster"]].reset_index().rename({'index': "spot_barcode"}, axis=1)


cluster_out.to_csv(cluster_out_path, index=False)

