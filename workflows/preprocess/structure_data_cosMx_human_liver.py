import json
import shutil
import sys
import glob 

sample = str(sys.argv[1])
out_folder = str(sys.argv[2])

img_path = f"data/{sample}_pca_rgb.jpg"
img_format = img_path.split(".")[-1]
img_out_path = f"{out_folder}/data/image/{sample}.{img_format}"
json_out_path = f"{out_folder}/data/meta/{sample}.json"

# move image 
shutil.copy(img_path, img_out_path)

spot_diameter_fullres = 100

# general meta info about sample
json_info = {"SAMPLE": sample, "spot_diameter_fullres": spot_diameter_fullres, "dot_size": 30}
json_info

with open(json_out_path, 'w') as f:
    f.write(json.dumps(json_info))