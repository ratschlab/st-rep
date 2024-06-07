import json
import shutil
import sys

SAMPLE = str(sys.argv[1])
out_folder = str(sys.argv[2])

img_path = f"data/HumanPilot10x/image_full/{SAMPLE}_full_image.tif"
json_path = f"data/spatial-datasets/data-raw/2020_maynard_prefrontal-cortex/{SAMPLE}/spatial/scalefactors_json.json"
json_out_path = f"{out_folder}/data/meta/{SAMPLE}.json"

img_format = img_path.split(".")[-1]
img_out_path = f"{out_folder}/data/image/{SAMPLE}.{img_format}"

spot_diameter_fullres = round(json.load(open(json_path))["fiducial_diameter_fullres"])
spot_diameter_fullres

# move image 
shutil.copy(img_path, img_out_path)


# general meta info about sample
json_info = {"SAMPLE": SAMPLE, "spot_diameter_fullres": spot_diameter_fullres, "dot_size": 5}
json_info

with open(json_out_path, 'w') as f:
    f.write(json.dumps(json_info))