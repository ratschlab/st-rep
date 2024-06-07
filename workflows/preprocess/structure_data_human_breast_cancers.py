import json
import shutil
import sys

SAMPLE = str(sys.argv[1])
out_folder = str(sys.argv[2])

base_path = f"data/human_breast_cancers/"
img_path = f"{base_path}/spatial/{SAMPLE}_spatial/tissue_hires_image.png"
img_format = img_path.split(".")[-1]
img_out_path = f"{out_folder}/data/image/{SAMPLE}.{img_format}"
json_file = f"{base_path}/spatial/{SAMPLE}_spatial/scalefactors_json.json"
json_out_path = f"{out_folder}/data/meta/{SAMPLE}.json"

# move image 
shutil.copy(img_path, img_out_path)

json_info = json.load(open(json_file))
          
tissue_hires_scalef = float(json_info["tissue_hires_scalef"])
spot_diameter_fullres = float(json_info["spot_diameter_fullres"])
scaled_spot_diameter_fullres = max(spot_diameter_fullres, 75) # inception_v3 in stLearn requires >= 75


dot_size = 5 if SAMPLE in ["1142243F", "1160920F"] else 10
                      
# general meta info about sample
json_info = {"SAMPLE": SAMPLE, "spot_diameter_fullres": scaled_spot_diameter_fullres, "dot_size": dot_size}

with open(json_out_path, 'w') as f:
    f.write(json.dumps(json_info))