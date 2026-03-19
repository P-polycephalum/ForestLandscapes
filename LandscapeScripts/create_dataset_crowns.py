import json
import os
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import box

plot="50ha"
spp="Cavanillesia platanifolia"

master_control_json=f"D:\\BCI_{plot}_timeseries\\master_control.json"
master_control_gdf= f"D:\\BCI_{plot}_timeseries\\master_control_polygons.gpkg"


crownmap=gpd.read_file(master_control_gdf)

with open(master_control_json, 'r') as f:
    master_control = json.load(f)

json_filtered = {k: v for k, v in master_control['crowns'].items() if v.get('latin') == spp}
keys_to_keep = set(json_filtered.keys())


path_orthomosaics=r"D:\BCI_50ha_timeseries\Product_local"
orthomosaics=os.listdir(path_orthomosaics)

path_crowns_out=r"D:\BCI_50ha_timeseries\crowns"
os.makedirs(path_crowns_out, exist_ok=True)
#only cavallinesian trees
subset=crownmap[crownmap['polygon_id'].isin(keys_to_keep)]
subset['global_id']= subset['polygon_id'].apply(lambda x: x.split("_")[0])
subset['date']= subset['polygon_id'].apply(lambda x: x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[3])

for ortho in orthomosaics:
    dat= "_".join(ortho.split("_")[2:5])
    print(f"Processing orthomosaic {ortho}...")

    #get all trees with that date
    subset2= subset[subset['date']== dat]

    with rasterio.open(os.path.join(path_orthomosaics, ortho)) as src:
        for tree in subset2.itertuples():
            global_id = tree.global_id
            dat= tree.date
            polygon_id = global_id + "_" + dat
            print(f"Processing tree {global_id}...")
            bounds = tree.geometry.bounds
            box_crown_5 = box(bounds[0] - 5, bounds[1] - 5, bounds[2] + 5, bounds[3] + 5)

            out_image, out_transform = mask(src, [box_crown_5], crop=True)
            out_meta = src.meta.copy()
            x_min, y_min = out_transform * (0, 0)
            xres, yres = out_transform[0], out_transform[4]
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })
            tree_folder= os.path.join(path_crowns_out, f"{global_id}")
            os.makedirs(tree_folder, exist_ok=True)
            with rasterio.open(os.path.join(tree_folder, f"{polygon_id}.tif"), "w", **out_meta) as dest:
                dest.write(out_image)
