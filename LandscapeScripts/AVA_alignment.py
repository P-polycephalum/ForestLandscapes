from arosics import COREG, COREG_LOCAL
import os
import shutil
import rasterio
import numpy as np
from shapely.geometry import box
from rasterio.mask import mask
import geopandas as gpd
from matplotlib import pyplot as plt


plot="BCI_50ha"
path_drone= r"\\stri-sm01\ForestLandscapes\LandscapeProducts\Drone"
path_reference= r"\\stri-sm01\ForestLandscapes\UAVSHARE\AVUELO_crownmap\BCI_25haplot\2024-07-16_orthoWhole_bci_resFull_clipped\2024-07-16_orthoWhole_bci_resFull_clipped.tif"
path_crownmap= r"\\stri-sm01\ForestLandscapes\LandscapeRaw\Crownmaps\Big_plots\lefo\BCI_ava_crownmap_2025.gpkg"

temp_path= r"\\stri-sm01\ForestLandscapes\UAVSHARE\BCI_ava_timeseries\orthomosaic"
os.makedirs(temp_path, exist_ok=True)

for landscape in os.listdir(path_drone):
    year_path = os.path.join(path_drone, landscape)
    if not os.path.isdir(year_path):
        continue
    for mission in os.listdir(year_path):
        if "BCI_ava" in mission:
            mission_path = os.path.join(year_path, mission, "orthophoto")
            if not os.path.isdir(mission_path):
                continue
            for item in os.listdir(mission_path):
                if "orthomosaic" in item.lower() and item.lower().endswith(".tif"):
                    orthomosaic_path = os.path.join(mission_path, item)
                    target_path = os.path.join(temp_path, item)
                    if not os.path.exists(target_path):
                        print(f"Found orthomosaic: {orthomosaic_path} -> copying to {target_path}")
                        shutil.copy2(orthomosaic_path, target_path)
                    else:
                        print(f"Target already exists, skipping: {target_path}")

orthomosaics= os.listdir(temp_path)
#lets get the extent of one of them

crownmap_ava_2025= gpd.read_file(path_crownmap)
bounds_crownmap= crownmap_ava_2025.total_bounds
box_crownmap= box(bounds_crownmap[0]-10, bounds_crownmap[1]-10, bounds_crownmap[2]+10, bounds_crownmap[3]+10)



path_mavic= r"\\stri-sm01\ForestLandscapes\UAVSHARE\BCI_50ha_timeseries\mavic"
#timeseries parameters
res=0.05
target_crs = "EPSG:32617"
left=625753.8462543194
bottom=1011723.3728508463
right=626809.573072755
top=1012295.6060808859
box_timeseries= box(left, bottom, right, top)

#now we need to pull the orthomosaics from the mavic only present in landscape 2024 2025 and 2026
for landscape in os.listdir(path_drone):
    if landscape in ["2024", "2025", "2026"]:
        year_path = os.path.join(path_drone, landscape)
        if not os.path.isdir(year_path):
            continue
        for mission in os.listdir(year_path):
            if "BCI_50ha" in mission and "M3E" in mission and any(date in mission for date in ["2024_08_07", "2024-07-22", "20241112"]):
                mission_path = os.path.join(year_path, mission, "orthophoto")
                if os.path.isdir(mission_path):
                    files = [f for f in os.listdir(mission_path) if f.lower().endswith('.tif') and 'orthomosaic' in f.lower()]
                    if not files:
                        continue
                
                    ortho_file = files[0]
                    out_path = os.path.join(path_mavic, ortho_file)
                    if not os.path.exists(out_path):
                        print(f"Orthomosaic already exists, skipping: {out_path}")                    
                        continue
                    ortho_path = os.path.join(mission_path, ortho_file)

                    # if "BCI_50ha_2025_01_03_orthomosaic.tif" in ortho_file or "BCI_50ha_2024_02_21_orthomosaic.tif" in ortho_file:
                    #     print(f"Skipping {ortho_file} due to known issues.")
                    #     continue

                    if os.path.exists(ortho_path):
                    
                        print(f"Found orthomosaic: {ortho_path} -> processing to {out_path}")
                        os.makedirs(os.path.dirname(out_path), exist_ok=True)
                        with rasterio.open(ortho_path) as src:
                            out_image, out_transform = mask(src, [box_timeseries], crop=True, all_touched=True)
                            out_image_selected = out_image[[0, 1, 2, 7], :, :]
                            out_meta = src.meta.copy()
                            out_meta.update({
                                "driver": "GTiff",
                                "height": out_image_selected.shape[1],
                                "width": out_image_selected.shape[2],
                                "transform": out_transform,
                                "count": 4
                            })
                            with rasterio.open(out_path, "w", **out_meta) as dest:
                                dest.write(out_image_selected)



path_cropped= temp_path.replace("orthomosaic", "cropped")
os.makedirs(path_cropped, exist_ok=True)
for landscape in os.listdir(temp_path):
    src_path = os.path.join(temp_path, landscape)
    print(f"Processing: {src_path}")

    if not os.path.exists(src_path):
        print(f"  Source not found, skipping: {src_path}")
        continue
    if not os.path.isfile(src_path):
        print(f"  Not a file, skipping: {src_path}")
        continue

    aligned_path = os.path.join(path_cropped, landscape)
    if os.path.exists(aligned_path):
        print(f"  Target already exists, skipping: {aligned_path}")
        continue

    try:
        with rasterio.open(src_path) as src:
            out_image, out_transform = mask(src, [box_crownmap], crop=True, all_touched=True)
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
                "count": out_image.shape[0]
            })

            os.makedirs(os.path.dirname(aligned_path), exist_ok=True)
            with rasterio.open(aligned_path, "w", **out_meta) as dest:
                dest.write(out_image)
        print(f"  Written: {aligned_path}")
    except Exception as e:
        print(f"  Failed processing {src_path}: {e}")
        continue














