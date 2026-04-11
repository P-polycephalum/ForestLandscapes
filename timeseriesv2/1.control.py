import os
import subprocess
import geopandas as gpd
import pandas as pd
from shapely.geometry import box
from tqdm import tqdm
import geopandas as gpd
import json
from collections import defaultdict

def split_into_buckets(gdf, x_size=250, y_size=250, n_tiles=None, grid_shape=None):
    """
    Split polygons into buckets using centroid assignment.

    Modes:
    1) Explicit size mode (default): provide x_size and y_size.
    2) Fixed shape mode: provide grid_shape=(rows, cols), e.g. (2, 4).
    3) Fixed count mode: provide n_tiles; rows/cols are chosen from factor pairs
       to make bucket cells as close to square as possible.
    """
    if gdf.empty:
        return {}

    gdf = gdf.copy()
    minx, miny, maxx, maxy = gdf.total_bounds
    x_range = maxx - minx
    y_range = maxy - miny

    if x_range <= 0 or y_range <= 0:
        raise ValueError("Invalid bounds: width/height must be positive.")

    rows = None
    cols = None

    if grid_shape is not None:
        if len(grid_shape) != 2:
            raise ValueError("grid_shape must be (rows, cols).")
        rows, cols = int(grid_shape[0]), int(grid_shape[1])
        if rows <= 0 or cols <= 0:
            raise ValueError("grid_shape values must be > 0.")
        x_size = x_range / cols
        y_size = y_range / rows
    elif n_tiles is not None:
        n_tiles = int(n_tiles)
        if n_tiles <= 0:
            raise ValueError("n_tiles must be > 0.")

        # Choose factor pair that makes cells closest to square.
        best_rows, best_cols = 1, n_tiles
        best_score = float("inf")
        for r in range(1, int(np.sqrt(n_tiles)) + 1):
            if n_tiles % r != 0:
                continue
            c = n_tiles // r
            for rows_candidate, cols_candidate in ((r, c), (c, r)):
                cell_w = x_range / cols_candidate
                cell_h = y_range / rows_candidate
                score = abs(np.log(cell_w / cell_h))
                if score < best_score:
                    best_score = score
                    best_rows, best_cols = rows_candidate, cols_candidate

        rows, cols = best_rows, best_cols
        x_size = x_range / cols
        y_size = y_range / rows
    else:
        if x_size <= 0 or y_size <= 0:
            raise ValueError("x_size and y_size must be > 0.")
        cols = int(np.ceil(x_range / x_size))
        rows = int(np.ceil(y_range / y_size))

    centroids = gdf.geometry.centroid
    col_float = (centroids.x - minx) / x_size
    row_float = (centroids.y - miny) / y_size

    # Clamp to valid range so edge centroids land in the last bucket.
    gdf["col"] = np.clip(np.floor(col_float).astype(int), 0, cols - 1)
    gdf["row"] = np.clip(np.floor(row_float).astype(int), 0, rows - 1)

    gdf["bucket_id"] = gdf["row"].astype(str) + "_" + gdf["col"].astype(str)

    buckets = {
        bucket_id: group.drop(columns=["row", "col", "bucket_id"])
        for bucket_id, group in gdf.groupby("bucket_id")
    }
    print(
        f"Bucket grid: {rows} rows x {cols} cols ({rows * cols} tiles), "
        f"tile size ~ {x_size:.2f}m x {y_size:.2f}m"
    )
    return buckets

def _run_gdal(cmd):
    """Run GDAL command with captured output (avoid noisy raw stderr)."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"stdout: {p.stdout[-400:]}\n"
            f"stderr: {p.stderr[-400:]}"
        )

def is_valid_tiff(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    p = subprocess.run(["gdalinfo", path], capture_output=True, text=True)
    return p.returncode == 0

def build_combined_vrt(ortho_path, dsm_path, combined_vrt, bounds, res):
    xmin, ymin, xmax, ymax = bounds

    ortho_vrt = combined_vrt.replace(".vrt", "_ortho.vrt")
    dsm_vrt = combined_vrt.replace(".vrt", "_dsm.vrt")

    # ORTHO
    _run_gdal([
        "gdalbuildvrt",
        "-te", str(xmin), str(ymin), str(xmax), str(ymax),
        "-tr", str(res), str(res),
        "-tap",
        "-r", "bilinear",
        "-overwrite",
        ortho_vrt,
        ortho_path
    ])

    # DSM
    _run_gdal([
        "gdalbuildvrt",
        "-te", str(xmin), str(ymin), str(xmax), str(ymax),
        "-tr", str(res), str(res),
        "-tap",
        "-r", "nearest",
        "-overwrite",
        dsm_vrt,
        dsm_path
    ])

    # STACK
    _run_gdal([
        "gdalbuildvrt",
        "-separate",
        "-overwrite",
        combined_vrt,
        ortho_vrt,
        dsm_vrt
    ])

    return ortho_vrt, dsm_vrt

def extract_tile_from_vrt(combined_vrt, output_path, bounds):
    xmin, ymin, xmax, ymax = bounds

    _run_gdal([
        "gdal_translate",
        "-projwin", str(xmin), str(ymax), str(xmax), str(ymin),
        "-of", "GTiff",
        "-ot", "Byte",
        "-scale",
        "-co", "COMPRESS=LZW",
        combined_vrt,
        output_path
    ])
# timeseries parameters
res = 0.05
target_crs = "EPSG:32617"
left = 625753.8462543194
bottom = 1011723.3728508463
right = 626809.573072755
top = 1012295.6060808859
box_timeseries = box(left, bottom, right, top)

timeseries_orthomosaics = r"D:\BCI_50ha_timeseries\orthomosaic"
timeseries_dsms = r"D:\BCI_50ha_timeseries\dsm"
tiles_folder = r"D:\BCI_50ha_timeseries\tiles"
os.makedirs(tiles_folder, exist_ok=True)

############################################################################################
crownmap = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_2022_2023_crownmap_raw.shp"
crownmap_gdf = gpd.read_file(crownmap)

buck = split_into_buckets(crownmap_gdf, grid_shape=(2, 4))
bucket_attributes = {}
for bucket_id, bucket_gdf in buck.items():
    minx, miny, maxx, maxy = bucket_gdf.total_bounds
    buffered_box = box(minx - 5, miny - 5, maxx + 5, maxy + 5)
    bucket_attributes[bucket_id] = {
        "polygons": bucket_gdf,
        "box": buffered_box,
    }

# gdfs_to_concat = []
# for bucket_id, bucket_gdf in buck.items():
#     bucket_gdf_copy = bucket_gdf.copy()
#     bucket_gdf_copy["bucket_id"] = bucket_id
#     gdfs_to_concat.append(bucket_gdf_copy)

# master_gdf = gpd.GeoDataFrame(
#     pd.concat(gdfs_to_concat, ignore_index=True),
#     crs=crownmap_gdf.crs
# )

# master_gdf_path = r"D:\BCI_50ha_timeseries\master_gdf.geoparquet"
# master_gdf.to_parquet(master_gdf_path, index=False)
# print(f"Saved master GDF → {len(master_gdf)} rows")

# Get list of tif files to process
tif_files = [f for f in os.listdir(timeseries_orthomosaics) if f.endswith(".tif")]
print(f"Found {len(tif_files)} orthomosaics to process\n")

def collect_invalid_tiles(tiles_folder):
    tile_files = [f for f in os.listdir(tiles_folder) if f.endswith(".tif")]
    invalid_files = []
    for tf in tqdm(tile_files, desc="Validating tiles"):
        path = os.path.join(tiles_folder, tf)
        if not is_valid_tiff(path):
            invalid_files.append(tf)
    return invalid_files

invalid=collect_invalid_tiles(tiles_folder)

def delete_invalid_tiles(tiles_folder, invalid_files):
    for tf in tqdm(invalid_files, desc="Deleting invalid tiles"):
        path = os.path.join(tiles_folder, tf)
        if os.path.exists(path):
            os.remove(path)
            print(f"Deleted invalid file: {tf}")

max_passes = 5
prev_invalid_count = None
pass_num = 1

while pass_num <= max_passes:
    tif_files = [f for f in os.listdir(timeseries_orthomosaics) if f.endswith(".tif")]
    print(f"\n=== PASS {pass_num} ===")
    print(f"Found {len(tif_files)} orthomosaics to process")

    # Cleanup invalid tiles before processing this pass
    invalid_files = collect_invalid_tiles(tiles_folder)
    print(f"Invalid tiles before processing: {len(invalid_files)}")
    if invalid_files:
        delete_invalid_tiles(tiles_folder, invalid_files)

    # Process all available orthomosaics
    for combined_file in tqdm(tif_files, desc=f"Processing timeseries (pass {pass_num})"):
        ortho_path = os.path.join(timeseries_orthomosaics, combined_file)
        dsm_path = os.path.join(timeseries_dsms, combined_file.replace("orthomosaic", "dsm"))
        name = os.path.splitext(combined_file)[0]

        if not is_valid_tiff(ortho_path) or not is_valid_tiff(dsm_path):
            print(f"[WARNING] Skipping {combined_file} due to invalid orthomosaic or DSM.")
            continue

        combined_vrt = os.path.join(tiles_folder, f"{name}_combined.vrt")

        try:
            ortho_vrt, dsm_vrt = build_combined_vrt(
                ortho_path,
                dsm_path,
                combined_vrt,
                box_timeseries.bounds,
                res
            )
        except Exception as e:
            print(f"[ERROR] Failed building VRT for {name}: {e}")
            continue

        file_failed = False

        for bucket_id in buck:
            bounds = bucket_attributes[bucket_id]["box"].bounds
            output_path = os.path.join(tiles_folder, f"{name}_tile_{bucket_id}.tif")

            if os.path.exists(output_path):
                continue

            try:
                extract_tile_from_vrt(combined_vrt, output_path, bounds)
            except Exception as e:
                print(f"[ERROR] {name} bucket {bucket_id}: {e}")
                file_failed = True
                break

        # cleanup VRTs
        for f in [combined_vrt, ortho_vrt, dsm_vrt]:
            if os.path.exists(f):
                os.remove(f)

        # delete inputs only if all expected tiles exist
        if (not file_failed) and all(
            os.path.exists(os.path.join(tiles_folder, f"{name}_tile_{bid}.tif"))
            for bid in buck
        ):
            if os.path.exists(ortho_path):
                os.remove(ortho_path)
            if os.path.exists(dsm_path):
                os.remove(dsm_path)

    # Post-pass diagnostics
    invalid_files = collect_invalid_tiles(tiles_folder)
    invalid_count = len(invalid_files)
    print(f"Invalid tiles after pass {pass_num}: {invalid_count}")

    if invalid_count == 0:
        print("No invalid tiles remain. Stopping retry loop.")
        break

    delete_invalid_tiles(tiles_folder, invalid_files)

    # anti-infinite-loop protection
    if prev_invalid_count is not None and invalid_count >= prev_invalid_count:
        print("Invalid count did not improve. Stopping to avoid infinite loop.")
        break

    prev_invalid_count = invalid_count
    pass_num += 1
def find_missing_tiles_by_date(tiles_folder, buck):
    """
    Returns:
        dict[str, list[str]] where key = date base name (without _tile_*),
        value = missing bucket IDs for that date.
    """
    expected_bucket_ids = {str(bid) for bid in buck.keys()}  # e.g. {"0_0","0_1",...,"1_3"}
    tiles_by_date = {}

    tile_files = [f for f in os.listdir(tiles_folder) if f.endswith(".tif")]
    for tf in tile_files:
        if "_tile_" not in tf:
            continue
        stem = os.path.splitext(tf)[0]  # BCI_50ha_2024_08_07_orthomosaic_tile_0_0
        date_base, bid = stem.rsplit("_tile_", 1)  # (..._orthomosaic, 0_0)
        tiles_by_date.setdefault(date_base, set()).add(bid)

    missing_by_date = {}
    for date_base, present_ids in tiles_by_date.items():
        missing = sorted(expected_bucket_ids - present_ids)
        if missing:
            missing_by_date[date_base] = missing

    return missing_by_date

missing_by_date = find_missing_tiles_by_date(tiles_folder, buck)

if not missing_by_date:
    print("All dates have all expected tiles.")
else:
    print(f"Dates with missing tiles: {len(missing_by_date)}")
    for date_base, missing_ids in missing_by_date.items():
        print(f"[MISSING] {date_base}: {missing_ids}")

##################################################################################
dir_address = r"D:\BCI_50ha_timeseries"
json_path = os.path.join(dir_address, "bucket_attributes.json")
# List tile rasters once.
tile_files = sorted([f for f in os.listdir(tiles_folder) if f.endswith(".tif")])
print(f"Found {len(tile_files)} tiles to process\n")

# Add the per-bucket time-series file list to bucket attributes.
for bucket_id in buck.keys():
    tile_suffix = f"_tile_{bucket_id}.tif"
    bucket_files = [
        os.path.join(tiles_folder, tile_file)
        for tile_file in tile_files
        if tile_file.endswith(tile_suffix)
    ]
    bucket_attributes[bucket_id]["files"] = sorted(bucket_files)

bucket_attributes.keys()

#SAVE updated bucket attributes with file lists to json
with open(json_path, "w") as f:
    json.dump(bucket_attributes, f, indent=4)  