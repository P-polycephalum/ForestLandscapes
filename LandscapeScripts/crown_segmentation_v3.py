import gc
import os
import torch
from LandscapeScripts.utils import tile_ortho, crown_segment, crown_avoid, crownmap_QC, crownmap_metrics
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box as box1
from sam2.build_sam import build_sam2
from sam2.sam2_image_predictor import SAM2ImagePredictor

#load the sam2 model
sam2_checkpoint = r"D:\BCI_ava_timeseries\sam2.1_hiera_large.pt"
model_cfg = r"D:\BCI_ava_timeseries\sam2.1_hiera_l.yaml"
device="cuda"


# Define paths, still monolothic
wd_path=r"D:\BCI_ava_timeseries"
tile_folder= os.path.join(wd_path,"tiles")
segmented_crowns_folder= os.path.join(wd_path,"segmented_crowns")
os.makedirs(segmented_crowns_folder, exist_ok=True)
os.makedirs(tile_folder, exist_ok=True)

#reference crownmap
ava_crownmap= os.path.join(wd_path,"crownmap","BCI_ava_crownmap_2025.gpkg")
ava_crownmap_gdf=gpd.read_file(ava_crownmap)

#reference orthomosaic
reference_image= os.path.join(wd_path,"AVA_plot_clipped.tif")

#orthomosaics
orthomosaics = [f for f in os.listdir(os.path.join(wd_path, "Product_local")) if f.endswith(".tif")]


def clear_tile_folder(folder_path: str) -> None:
    if not os.path.isdir(folder_path):
        return
    for entry in os.listdir(folder_path):
        entry_path = os.path.join(folder_path, entry)
        if os.path.isfile(entry_path):
            os.remove(entry_path)

#segment the reference
tile_ortho(reference_image,100,30,tile_folder)
segmented_crownmap= crown_segment(tile_folder,ava_crownmap_gdf,model_cfg,sam2_checkpoint,device=device)
segmented_crownmap= crownmap_metrics(ava_crownmap_gdf, segmented_crownmap)
segmented_crownmap= segmented_crownmap.sort_values(by="similarity", ascending=False).drop_duplicates('global_id')
crownavoid_gdf= crown_avoid(segmented_crownmap)
crownavoid_gdf=crownmap_metrics(ava_crownmap_gdf, crownavoid_gdf)
crown_qc= crownmap_QC(ava_crownmap_gdf, crownavoid_gdf, precision_threshold=0.6, recall_threshold=0.6, iou_threshold=0.6)
crown_qc.to_file(os.path.join(segmented_crowns_folder,"BCI_ava_reference_crownmap.gpkg"))


for orthomosaic in orthomosaics[128:]:
    print(f"Processing {orthomosaic}...")
    orthomosaic_path = os.path.join(wd_path, "Product_local", orthomosaic)
    output_name = f"BCI_ava_crownmap_{'_'.join(orthomosaic.split('_')[2:5])}.gpkg"
    output_path = os.path.join(segmented_crowns_folder, output_name)
    tile_ortho(orthomosaic_path, 100, 30, tile_folder)
    segmented_crownmap = crown_segment(tile_folder, ava_crownmap_gdf, checkpoint=sam2_checkpoint,model_cfg=model_cfg, device=device)
    segmented_crownmap = crownmap_metrics(ava_crownmap_gdf, segmented_crownmap)
    segmented_crownmap = segmented_crownmap.sort_values(by="similarity", ascending=False).drop_duplicates("global_id")
    crownavoid_gdf = crown_avoid(segmented_crownmap)
    crownavoid_gdf = crownmap_metrics(ava_crownmap_gdf, crownavoid_gdf)
    crown_qc = crownmap_QC(
        ava_crownmap_gdf,
        crownavoid_gdf,
        precision_threshold=0.6,
        recall_threshold=0.6,
        iou_threshold=0.6,
    )
    crown_qc.to_file(output_path)
    ava_crownmap_gdf = gpd.read_file(output_path)

    del segmented_crownmap, crownavoid_gdf, crown_qc
    clear_tile_folder(tile_folder)
    gc.collect()
    if device == "cuda":
        torch.cuda.empty_cache()



for orthomosaic in orthomosaics[126::-1]:
    print(f"Processing {orthomosaic}...")
    orthomosaic_path = os.path.join(wd_path, "Product_local", orthomosaic)
    output_name = f"BCI_ava_crownmap_{'_'.join(orthomosaic.split('_')[2:5])}.gpkg"
    output_path = os.path.join(segmented_crowns_folder, output_name)
    tile_ortho(orthomosaic_path, 100, 30, tile_folder)
    segmented_crownmap = crown_segment(tile_folder, ava_crownmap_gdf, checkpoint=sam2_checkpoint,model_cfg=model_cfg, device=device)
    segmented_crownmap = crownmap_metrics(ava_crownmap_gdf, segmented_crownmap)
    segmented_crownmap = segmented_crownmap.sort_values(by="similarity", ascending=False).drop_duplicates("global_id")
    crownavoid_gdf = crown_avoid(segmented_crownmap)
    crownavoid_gdf = crownmap_metrics(ava_crownmap_gdf, crownavoid_gdf)
    crown_qc = crownmap_QC(
        ava_crownmap_gdf,
        crownavoid_gdf,
        precision_threshold=0.6,
        recall_threshold=0.6,
        iou_threshold=0.6,
    )
    crown_qc.to_file(output_path)
    ava_crownmap_gdf = gpd.read_file(output_path)

    del segmented_crownmap, crownavoid_gdf, crown_qc
    clear_tile_folder(tile_folder)
    gc.collect()
    if device == "cuda":
        torch.cuda.empty_cache()






