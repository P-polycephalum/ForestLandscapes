"""Prepare master control JSON and pre-crop crowns from aligned tile mosaics."""

import json
import os

import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import box
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import rasterio
import zarr

ROOT_PATH = r"D:\BCI_50ha_timeseries"
MASTER_GDF_PATH = os.path.join(ROOT_PATH, "master_gdf.gpkg")
MASTER_JSON_PATH = os.path.join(ROOT_PATH, "master_control.json")
REFERENCE_CROWNMAP_PATH = os.path.join(ROOT_PATH, "crownmap", "BCI_50ha_2022_2023_crownmap_raw.shp")
TILES_ROOT = os.path.join(ROOT_PATH, "tiles", "aligned_local")
CROWNS_DIR = os.path.join(ROOT_PATH, "crowns")

TILE = "0_0"
CROP_BUFFER_METERS = 5.0
MAX_SHAPE_DELTA_PIXELS = 3


def coerce_array_shape(image: np.ndarray, expected_shape: tuple[int, int, int]) -> np.ndarray | None:
    """Pad/crop small edge differences so arrays can be stacked safely.

    Returns None when the mismatch is larger than MAX_SHAPE_DELTA_PIXELS.
    """
    count_exp, height_exp, width_exp = expected_shape
    count_in, height_in, width_in = image.shape

    if count_in != count_exp:
        return None

    if abs(height_in - height_exp) > MAX_SHAPE_DELTA_PIXELS:
        return None
    if abs(width_in - width_exp) > MAX_SHAPE_DELTA_PIXELS:
        return None

    out = image

    # Align top-left and reconcile trailing-edge pixel drift.
    if height_in > height_exp:
        out = out[:, :height_exp, :]
    elif height_in < height_exp:
        pad_h = height_exp - height_in
        out = np.pad(out, ((0, 0), (0, pad_h), (0, 0)), mode="constant", constant_values=0)

    if width_in > width_exp:
        out = out[:, :, :width_exp]
    elif width_in < width_exp:
        pad_w = width_exp - width_in
        out = np.pad(out, ((0, 0), (0, 0), (0, pad_w)), mode="constant", constant_values=0)

    return out


def build_master_control(root_path: str, crownmap_path: str) -> dict:
    """Create the JSON schema used by labeler workflows."""
    return {
        "paths": {
            "root": root_path,
            "crowns_dir": os.path.join(root_path, "crowns"),
            "crownmap_path": crownmap_path,
            "orthomosaics_dir": os.path.join(root_path, "tiles", "local_aligned"),
        },
        "variables": {
            "global_id": {"description": "Unique tree identifier"},
            "date": {
                "description": "Acquisition date for crown polygon. Exception for reference crowns where date is 'reference'"
            },
            "latin": {"description": "Latin species name"},
            "labels": {
                "description": "Variable labels for the tree crown including leaf coverage, flowering, fruiting, and new leaves"
            },
            "area": {"description": "Area of the tree crown polygon in square meters"},
            "score": {"description": "Score from SAM2"},
            "tag": {"description": "Tag from ForestGEO census data"},
            "IoU": {
                "description": "Intersection over Union similarity score between segmented crown and reference crown"
            },
            "Precision": {
                "description": "Precision score for the segmented crown compared to the reference crown"
            },
            "Recall": {
                "description": "Recall score for the segmented crown compared to the reference crown"
            },
            "F1": {
                "description": "F1 score for the segmented crown compared to the reference crown"
            },
            "similarity": {
                "description": "Overall similarity score between segmented crown and reference crown, calculated as the average of IoU, Precision, Recall, and F1 scores"
            },
        },
        "crowns": {},
    }

def save_json(data: dict, output_path: str) -> None:
    with open(output_path, "w") as file:
        json.dump(data, file, indent=4)


def populate_crowns(master_gdf: gpd.GeoDataFrame, reference_crownmap: gpd.GeoDataFrame, master_control: dict) -> gpd.GeoDataFrame:
    """Add derived fields and initialize crown-level metadata records."""
    latin_map = reference_crownmap.set_index("GlobalID")["latin"]
    master_gdf["latin"] = master_gdf["global_id"].map(latin_map)
    master_gdf["polygon_id"] = master_gdf["global_id"] + "_" + master_gdf["date"]
    master_gdf["area"] = master_gdf["geometry"].area

    for row in master_gdf.itertuples(index=False):
        master_control["crowns"][str(row.polygon_id)] = {
            "global_id": row.global_id,
            "date": row.date,
            "latin": row.latin,
            "labels": [
                {
                    "leafing": None,
                    "flowering": None,
                    "fruiting": None,
                    "new_leaves": None,
                }
            ],
            "area": row.area,
            "score": None,
            "tag": None,
            "IoU": None,
            "Precision": None,
            "Recall": None,
            "F1": None,
            "quality": None,
            "edited": False,
            "no_edits_needed": False,
        }

    return master_gdf


def crop_crowns_from_tile(master_gdf: gpd.GeoDataFrame, tile: str) -> None:
    """Crop crown-centered windows for each date and write per-tree tiles."""
    os.makedirs(CROWNS_DIR, exist_ok=True)

    tile_dir = os.path.join(TILES_ROOT, tile)
    tif_files = sorted(file for file in os.listdir(tile_dir) if file.lower().endswith(".tif"))

    for tif_name in tif_files:
        print(f"Processing date {tif_name}...")
        date_value = "_".join(tif_name.split("_")[2:5])
        subset_gdf = master_gdf[master_gdf["date"] == date_value]
        print(f"{len(subset_gdf)} crowns to process for this date")

        with rasterio.open(os.path.join(tile_dir, tif_name)) as src:
            for row in subset_gdf.itertuples(index=False):
                xmin, ymin, xmax, ymax = row.geometry.bounds
                crown_box = box(
                    xmin - CROP_BUFFER_METERS,
                    ymin - CROP_BUFFER_METERS,
                    xmax + CROP_BUFFER_METERS,
                    ymax + CROP_BUFFER_METERS,
                )

                out_image, out_transform = mask(src, [crown_box], crop=True)
                out_meta = src.meta.copy()
                out_meta.update(
                    {
                        "driver": "GTiff",
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "transform": out_transform,
                    }
                )

                crown_folder = os.path.join(CROWNS_DIR, row.global_id)
                os.makedirs(crown_folder, exist_ok=True)
                crown_path = os.path.join(crown_folder, f"{row.polygon_id}.tif")

                with rasterio.open(crown_path, "w", **out_meta) as dst:
                    dst.write(out_image)


def crop_crowns_to_zarr(master_gdf: gpd.GeoDataFrame, tile: str) -> None:
    """
    Crop crown-centered windows for each date and save one Zarr stack per tree.

    Assumes:
    - master_gdf has columns: global_id, date, geometry
    - tile folder contains one tif per date
    - dates in tif filename are parsed the same way as below

    Output:
    CROWNS_DIR/
        <global_id>/
            images.zarr
            dates.csv
    """
    os.makedirs(CROWNS_DIR, exist_ok=True)

    tile_dir = os.path.join(TILES_ROOT, tile)
    tif_files = sorted(
        file for file in os.listdir(tile_dir)
        if file.lower().endswith(".tif")
    )

    if len(tif_files) == 0:
        print(f"No tif files found in {tile_dir}")
        return

    # ------------------------------------------------------------------
    # 1. Build one fixed crop box per tree across all dates
    # ------------------------------------------------------------------
    tree_boxes = {}

    for global_id, group in master_gdf.groupby("global_id"):
        xmin, ymin, xmax, ymax = group.total_bounds
        tree_boxes[global_id] = box(
            xmin - CROP_BUFFER_METERS,
            ymin - CROP_BUFFER_METERS,
            xmax + CROP_BUFFER_METERS,
            ymax + CROP_BUFFER_METERS,
        )

    # ------------------------------------------------------------------
    # 2. Collect crops by tree across time
    # ------------------------------------------------------------------
    tree_images = defaultdict(list)
    tree_dates = defaultdict(list)
    tree_meta = {}

    for tif_name in tif_files:
        print(f"Processing date {tif_name}...")
        date_value = "_".join(tif_name.split("_")[2:5])

        subset_gdf = master_gdf[master_gdf["date"] == date_value]
        print(f"{len(subset_gdf)} crowns to process for this date")

        if subset_gdf.empty:
            continue

        tif_path = os.path.join(tile_dir, tif_name)

        with rasterio.open(tif_path) as src:
            for row in subset_gdf.itertuples(index=False):
                global_id = row.global_id
                crop_box = tree_boxes[global_id]

                out_image, out_transform = mask(src, [crop_box], crop=True)

                # Save metadata once per tree
                if global_id not in tree_meta:
                    tree_meta[global_id] = {
                        "crs": str(src.crs),
                        "transform": out_transform,
                        "count": out_image.shape[0],
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "dtype": str(out_image.dtype),
                    }

                # Enforce per-tree shape consistency for safe time stacking.
                expected_shape = (
                    tree_meta[global_id]["count"],
                    tree_meta[global_id]["height"],
                    tree_meta[global_id]["width"],
                )
                if out_image.shape != expected_shape:
                    adjusted = coerce_array_shape(out_image, expected_shape)
                    if adjusted is None:
                        print(
                            f"Skipping {global_id} on {date_value}: "
                            f"shape {out_image.shape} cannot be reconciled to {expected_shape}"
                        )
                        continue
                    out_image = adjusted

                if out_image.shape != expected_shape:
                    print(
                        f"Skipping {global_id} on {date_value}: "
                        f"shape {out_image.shape} != expected {expected_shape}"
                    )
                    continue

                tree_images[global_id].append(out_image)
                tree_dates[global_id].append(date_value)

    # ------------------------------------------------------------------
    # 3. Write one Zarr per tree
    # ------------------------------------------------------------------
    for global_id, image_list in tree_images.items():
        if len(image_list) == 0:
            continue

        print(f"Writing Zarr for tree {global_id} with {len(image_list)} dates")

        # stack shape: (time, band, y, x)
        stack = np.stack(image_list, axis=0)

        crown_folder = os.path.join(CROWNS_DIR, global_id)
        os.makedirs(crown_folder, exist_ok=True)

        zarr_path = os.path.join(crown_folder, "images.zarr")

        # Overwrite existing Zarr if present
        if os.path.exists(zarr_path):
            import shutil
            shutil.rmtree(zarr_path)

        z = zarr.open(
            zarr_path,
            mode="w",
            shape=stack.shape,
            chunks=(1, stack.shape[1], min(512, stack.shape[2]), min(512, stack.shape[3])),
            dtype=stack.dtype,
        )
        z[:] = stack

        # simple metadata
        z.attrs["global_id"] = global_id
        z.attrs["axis_order"] = "time,band,y,x"
        z.attrs["dates"] = tree_dates[global_id]
        z.attrs["crs"] = tree_meta[global_id]["crs"]
        z.attrs["transform"] = tuple(tree_meta[global_id]["transform"])  # affine tuple

        # save dates table too, easier to inspect later
        pd.DataFrame({"date": tree_dates[global_id]}).to_csv(
            os.path.join(crown_folder, "dates.csv"),
            index=False,
        )


def main() -> None:
    master_gdf = gpd.read_file(MASTER_GDF_PATH)
    reference_crownmap = gpd.read_file(REFERENCE_CROWNMAP_PATH)

    master_control = build_master_control(ROOT_PATH, REFERENCE_CROWNMAP_PATH)
    save_json(master_control, MASTER_JSON_PATH)

    master_gdf = populate_crowns(master_gdf, reference_crownmap, master_control)
    save_json(master_control, MASTER_JSON_PATH)

    crop_crowns_to_zarr(master_gdf, TILE)


if __name__ == "__main__":
    main()



path=r"C:\Users\vasquezvicente\repo\stem_summ.gpkg"
import geopandas as gpd

gp=gpd.read_file(path)
gp
gp['meanWD']