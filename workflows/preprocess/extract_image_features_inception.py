from torchvision.models import mobilenet_v2, MobileNet_V2_Weights
from torchvision.models import inception_v3, Inception_V3_Weights
from aestetik.utils.utils_morphology import extract_morphology_embeddings
from collections import OrderedDict
from torchvision import transforms
import torch.nn.functional as F
import scanpy as sc
import numpy as np
import torch
import json
import glob
import sys
import os

sample = str(sys.argv[1])
out_folder = str(sys.argv[2])

adata_in = f"{out_folder}/data/h5ad/{sample}.h5ad"
json_path = f"{out_folder}/data/meta/{sample}.json"
img_path = glob.glob(f"{out_folder}/data/image/{sample}*")[0]
img_format = img_path.split(".")[-1]

features_inception_v3_out_path = f"{out_folder}/data/image_features/{sample}_inception.npy"

spot_diameter_fullres = round(json.load(open(json_path))["spot_diameter_fullres"])

adata = sc.read(adata_in)
n_components = 100


########################## Inception_V3_Weights FEATURE EXTRACTION ##########################

weights = Inception_V3_Weights.DEFAULT
morphology_model = inception_v3(weights=weights)
morphology_model.fc = torch.nn.Identity()

morphology_model.eval()    
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
feature_dim = 2048

preprocess = transforms.Compose([
                transforms.ToTensor(),
                weights.transforms(antialias=True),
            ])

features_inception_v3 = extract_morphology_embeddings(img_path, 
                                             morphology_model,
                                             x_pixel=adata.obs.y_pixel, 
                                             y_pixel=adata.obs.x_pixel, 
                                             spot_diameter=spot_diameter_fullres,
                                             device=device,
                                             n_components=n_components,
                                             feature_dim=feature_dim,
                                             preprocess=preprocess)



########################## SAVE FEATURES ##########################

np.save(features_inception_v3_out_path, features_inception_v3)