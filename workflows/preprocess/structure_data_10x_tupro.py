import json
import shutil
import sys
import pandas as pd
import glob

SAMPLE = str(sys.argv[1])
out_folder = str(sys.argv[2])
replicate_new = SAMPLE
replicate = replicate_new if "-2" not in replicate_new else replicate_new.replace("-2", "_2")
replicate = replicate.replace("r", "rep") # we did this due to name length size
SAMPLE = replicate.split("-")[0]

base_path = f"data/projects2022-digpathways_datasets/10x_tupro"
image_info_file = f"{base_path}/10x_image-paths_resolution.csv"

df = pd.read_csv(image_info_file)
replicate_info = df[df.replicate == replicate]
img_path = replicate_info.image_path.values[0]

img_name = replicate_info.image_name.values[0].replace(".tif", "")
annotations_file = f"{base_path}/region_annotations/HE_tissue_classification/{replicate}-{img_name}.annotations"
annotated_image_path = glob.glob(f"{base_path}/region_annotations/HE_tissue_classification/{replicate}-{img_name}*.tif")[0]

resolution = replicate_info.resolution.values[0]
resolution = replicate_info.resolution.values[0]
tile_size=160
dot_size = 3
if resolution == 0.3448:
    tile_size_sample = tile_size
elif resolution == 0.1736:
    tile_size_sample = int(tile_size * (0.3448 / 0.1736))
    dot_size = 12

img_format = img_path.split(".")[-1]
json_out_path = f"{out_folder}/data/meta/{replicate_new}.json"
img_out_path = f"{out_folder}/data/image/{replicate_new}.{img_format}"

# move image 
shutil.copy(img_path, img_out_path)


# general meta info about sample
json_info = {"SAMPLE": replicate_new, "spot_diameter_fullres": tile_size_sample, "dot_size": dot_size, "image_quality": "low"}
json_info

with open(json_out_path, 'w') as f:
    f.write(json.dumps(json_info))