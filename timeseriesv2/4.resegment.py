import torch
import os
import numpy as np
import gc
import geopandas as gpd
import pandas as pd
import cv2
import zarr
from affine import Affine
from canopyrs.engine.config_parsers import SegmenterConfig
from canopyrs.engine.models.segmenter.sam3 import Sam3PredictorWrapper
from shapely.ops import transform as shp_transform
from shapely.affinity import translate
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, box
from PIL import Image

# ---------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------
def crowns_to_boxes(gdf_px, width, height):
    boxes = []
    for _, row in gdf_px.iterrows():
        xmin, ymin, xmax, ymax = row.geometry.bounds

        x0 = int(np.clip(xmin, 0, width))
        x1 = int(np.clip(xmax, 0, width))
        y0 = int(np.clip(ymin, 0, height))
        y1 = int(np.clip(ymax, 0, height))

        boxes.append([x0, y1, x1, y0])
    return boxes

def pixel_to_utm(geom, transform):
    return shp_transform(lambda x, y: transform * (x, y), geom)

def crownmap_metrics(original_crownmap, segmented_crownmap):
    segmented_ids = set(segmented_crownmap['GlobalID'].unique())
    original_ids = set(original_crownmap['GlobalID'].unique())
    missing_in_segmented = original_ids - segmented_ids
    extra_in_segmented = segmented_ids - original_ids
    print(f"Missing IDs in segmented: {missing_in_segmented if missing_in_segmented else 'None'}")
    print(f"Extra IDs in segmented: {extra_in_segmented if extra_in_segmented else 'None'}")

    original_polys = {row['GlobalID']: row['geometry'] for _, row in original_crownmap.iterrows()}
    segmented_polys = {row['GlobalID']: row['geometry'] for _, row in segmented_crownmap.iterrows()}

    common_ids = original_ids & segmented_ids
    for crown in common_ids:
        original_geom = original_polys.get(crown)
        segmented_geom = segmented_polys.get(crown)
        original_geom= original_geom.buffer(0)
        segmented_geom= segmented_geom.buffer(0)
        if original_geom is None or segmented_geom is None:
            print(f"Warning: Crown {crown} missing geometry in one of the inputs. Skipping.")
            continue

        intersection_area = original_geom.intersection(segmented_geom).area
        union_area = original_geom.union(segmented_geom).area
        precision = intersection_area / segmented_geom.area if segmented_geom.area > 0 else 0
        recall = intersection_area / original_geom.area if original_geom.area > 0 else 0
        iou = intersection_area / union_area if union_area > 0 else 0
        f1= 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        score = 0.6 * iou + (0.4 * f1)
        # print(f"GlobalID: {crown}, IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}, F1: {f1:.2f}, Score: {score:.2f}")
        segmented_crownmap.loc[segmented_crownmap['GlobalID'] == crown, 'IoU'] = iou
        segmented_crownmap.loc[segmented_crownmap['GlobalID'] == crown, 'Precision'] = precision
        segmented_crownmap.loc[segmented_crownmap['GlobalID'] == crown, 'Recall'] = recall
        segmented_crownmap.loc[segmented_crownmap['GlobalID'] == crown, 'F1'] = f1
        segmented_crownmap.loc[segmented_crownmap['GlobalID'] == crown, 'similarity'] = score
    
    return segmented_crownmap

def crown_avoid(crown_gdf):
    """Resolve overlapping crown polygons and keep valid polygon geometries.

    Takes a GeoDataFrame of crown polygons, converts multipolygons to their
    largest part, subtracts overlaps based on relative area, and drops any
    non-polygon geometries that remain.
    """
    
    crown_avoidance = crown_gdf.copy()
    crown_avoidance['geometry'] = crown_avoidance.geometry.buffer(0)

    for index, row in crown_avoidance.iterrows():
        if isinstance(row.geometry, MultiPolygon):
            multi_polygon = row.geometry
            polygons = [polygon for polygon in multi_polygon.geoms]
            largest_polygon = max(polygons, key=lambda polygon: polygon.area)
            crown_avoidance.at[index, 'geometry'] = largest_polygon
            #print(f"Converted MultiPolygon to largest Polygon for index {index}.")

    sindex = crown_avoidance.sindex
    modifications = {}  # Dictionary to collect modifications
    for idx, polygon in crown_avoidance.iterrows():
        possible_matches_index = list(sindex.intersection(polygon['geometry'].bounds))
        possible_matches = crown_avoidance.iloc[possible_matches_index]
        adjacents = possible_matches[possible_matches.geometry.intersects(polygon['geometry']) & (possible_matches.index != idx)]
        if adjacents.empty:
            continue
        else:
            for adj_idx, adj_polygon in adjacents.iterrows():
                if polygon['similarity'] > adj_polygon['similarity']:
                    # Adjacent loses overlap
                    modifications[adj_idx] = modifications.get(adj_idx, adj_polygon.geometry).difference(polygon.geometry)
                    #print(f"Polygon {idx} is more similar than {adj_idx}. Subtracting overlap from adjacent.")   
                elif polygon['similarity'] < adj_polygon['similarity']:
                    # Current polygon loses overlap
                    modifications[idx] = modifications.get(idx, polygon.geometry).difference(adj_polygon.geometry)
                    #print(f"Polygon {idx} is less similar than {adj_idx}. Subtracting overlap from polygon.")

    for idx, new_geom in modifications.items():
        crown_avoidance.at[idx, 'geometry'] = new_geom
    for index, row in crown_avoidance.iterrows():
        if isinstance(row.geometry, MultiPolygon):
            multi_polygon = row.geometry
            polygons = [polygon for polygon in multi_polygon.geoms]
            largest_polygon = max(polygons, key=lambda polygon: polygon.area)
            crown_avoidance.at[index, 'geometry'] = largest_polygon

    for index, row in crown_avoidance.iterrows():
        geom = row["geometry"]
        if isinstance(geom, GeometryCollection):
            polygons = [g for g in geom.geoms if isinstance(g, Polygon)]
            if polygons:
                crown_avoidance.at[index, "geometry"] = polygons[0]
            else:
                crown_avoidance.at[index, "geometry"] = pd.NA
        elif not isinstance(geom, Polygon):
            crown_avoidance.at[index, "geometry"] = pd.NA
    return crown_avoidance

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

    # Use centroid to assign each polygon to one grid cell.
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

RESEG_THRESHOLD = 0.7

def _pick_geometry_reseg(row, ref_geoms):
    """Keep new geometry only if it meets RESEG_THRESHOLD; otherwise fall back to seed."""
    sim = row["similarity"] if "similarity" in row.index and pd.notna(row["similarity"]) else 0.0
    if sim >= RESEG_THRESHOLD:
        return row["geometry"]
    gid = row["GlobalID"]
    if gid in ref_geoms.index:
        return ref_geoms.loc[gid]
    return row["geometry"]

def build_bucket_attributes(crownmap_path, tiles_folder, grid_shape=(2, 4), buffer_size=5):
    crownmap_gdf = gpd.read_file(crownmap_path)
    tile_files = sorted(
        file_name
        for file_name in os.listdir(tiles_folder)
        if file_name.lower().endswith(".tif")
    )
    buckets = split_into_buckets(crownmap_gdf, grid_shape=grid_shape)

    bucket_attributes = {}

    for bucket_id, bucket_gdf in buckets.items():
        minx, miny, maxx, maxy = bucket_gdf.total_bounds
        buffered_box = box(
            minx - buffer_size,
            miny - buffer_size,
            maxx + buffer_size,
            maxy + buffer_size,
        )

        tile_suffix = f"_tile_{bucket_id}.tif"
        bucket_files = sorted(
            os.path.join(tiles_folder, file_name)
            for file_name in tile_files
            if file_name.endswith(tile_suffix)
        )

        bucket_attributes[bucket_id] = {
            "polygons": bucket_gdf,
            "box": buffered_box,
            "files": bucket_files,
        }

    return crownmap_gdf, buckets, bucket_attributes

def _safe_time(time_str):
    """Turn a time string into a safe filename component."""
    return time_str.replace(":", "-").replace("/", "-")

def _part_path(bucket_id, time_str):
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    os.makedirs(bucket_dir, exist_ok=True)
    return os.path.join(bucket_dir, f"{_safe_time(time_str)}.parquet")

def _part_exists(bucket_id, time_str):
    return os.path.exists(_part_path(bucket_id, time_str))

def _read_part(bucket_id, time_str):
    return gpd.read_parquet(_part_path(bucket_id, time_str))

def _write_part(gdf, bucket_id, time_str):
    gdf.to_parquet(_part_path(bucket_id, time_str), index=False)

def _read_bucket(bucket_id):
    """Read all time slices for one bucket — cheap, only that bucket's files."""
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    if not os.path.isdir(bucket_dir):
        return gpd.GeoDataFrame()
    files = sorted(f for f in os.listdir(bucket_dir) if f.endswith(".parquet"))
    if not files:
        return gpd.GeoDataFrame()
    return gpd.GeoDataFrame(
        pd.concat([gpd.read_parquet(os.path.join(bucket_dir, f)) for f in files],
                  ignore_index=True),
        crs="EPSG:32617"
    )

def crowns_to_boxes_local(gdf):
    boxes = []
    for geom in gdf.geometry:
        minx, miny, maxx, maxy = geom.bounds
        boxes.append([minx, miny, maxx, maxy])
    return boxes

def _segment_crowns(seed_crowns, time_idx):
    """Run SAM on seed_crowns bounding boxes at time_idx. Returns list of row dicts (pixel space)."""
    crown_buckets = split_into_buckets(seed_crowns, grid_shape=(3, 3))
    time_rows = []
    for bucket_id, bucket_crowns in crown_buckets.items():
        print(f"  Sub-bucket {bucket_id} — crowns: {len(bucket_crowns)}")
        bucket_crowns_px = bucket_crowns.copy()
        bucket_crowns_px["geometry"] = bucket_crowns.geometry.apply(
            lambda geom: shp_transform(lambda x, y, z=None: att_inversed * (x, y), geom)
        )
        minx, miny, maxx, maxy = bucket_crowns_px.total_bounds

        buffer = 100
        minx -= buffer; miny -= buffer; maxx += buffer; maxy += buffer

        xmin_px = int(np.clip(minx, 0, width))
        xmax_px = int(np.clip(maxx, 0, width))
        ymin_px = int(np.clip(miny, 0, height))
        ymax_px = int(np.clip(maxy, 0, height))
        xmin_px, xmax_px = sorted([xmin_px, xmax_px])
        ymin_px, ymax_px = sorted([ymin_px, ymax_px])

        if xmin_px == xmax_px or ymin_px == ymax_px:
            print(f"  Skipping empty crop for sub-bucket {bucket_id}")
            continue

        try:
            zarr_img = np.asarray(z[time_idx, :3, ymin_px:ymax_px, xmin_px:xmax_px]).transpose(1, 2, 0)
            bucket_crowns_local = bucket_crowns_px.copy()
            bucket_crowns_local["geometry"] = bucket_crowns_px.geometry.apply(
                lambda g: translate(g, xoff=-xmin_px, yoff=-ymin_px)
            )
            boxes = crowns_to_boxes_local(bucket_crowns_local)
            if not boxes:
                print(f"  No valid boxes for sub-bucket {bucket_id}. Skipping.")
                continue
            boxes_np = np.asarray(boxes, dtype=np.float32)
            pil = Image.fromarray(zarr_img[..., :3]).convert("RGB")

            with torch.inference_mode():
                masks, scores = seg._predict_batch(pil, boxes_np)

            for idx, (thisscore, thiscrown) in enumerate(zip(scores, masks)):
                mask = np.squeeze(thiscrown)
                mask_np = (mask > 0).astype(np.uint8) * 255
                contours, _ = cv2.findContours(mask_np, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                crown_polygons = []
                for contour in contours:
                    contour = contour.reshape(-1, 2)
                    if cv2.contourArea(contour) < 4 or contour.shape[0] < 3:
                        continue
                    x = contour[:, 0] + xmin_px
                    y = contour[:, 1] + ymin_px
                    poly = Polygon(np.stack([x, y], axis=1)).buffer(0)
                    if poly.is_empty or not poly.is_valid or poly.area <= 0:
                        continue
                    crown_polygons.append(poly)
                if not crown_polygons:
                    continue
                best_poly = max(crown_polygons, key=lambda p: p.area)
                row = {
                    "geometry": best_poly,
                    "score": float(np.squeeze(thisscore)),
                    "time": time[time_idx],
                    "bucket_id": bucket_to_process,
                }
                if idx < len(bucket_crowns):
                    row["GlobalID"] = bucket_crowns.iloc[idx].get("GlobalID", None)
                    row["tag"] = bucket_crowns.iloc[idx].get("tag", None)
                    row["latin"] = bucket_crowns.iloc[idx].get("latin", None)
                time_rows.append(row)
        finally:
            for var in ['zarr_img', 'boxes', 'boxes_np', 'pil', 'masks', 'scores']:
                if var in locals():
                    del locals()[var]
            gc.collect()
    return time_rows

def process_direction_reseg(direction_indices):
    """Resegment only crowns below RESEG_THRESHOLD. Excellent-quality crowns are kept as-is
    and passed forward as the reference seed for the next time step."""
    # Reference stays untouched — read it as the starting seed.
    tile_crowns = _read_part(bucket_to_process, time[ref_idx]).copy()
    print(f"Loaded reference seed: {time[ref_idx]} ({len(tile_crowns)} crowns)")

    for i in direction_indices:
        current_time = time[i]

        if not _part_exists(bucket_to_process, current_time):
            print(f"⚠️  Partition missing for {current_time} — cannot resegment. Skipping.")
            continue

        current_part = _read_part(bucket_to_process, current_time)

        # --- classify crowns ---
        def _sim(r):
            s = r.get("similarity", None)
            return float(s) if pd.notna(s) else 0.0

        is_excellent = current_part.apply(
            lambda r: str(r.get("quality", "")).strip().lower() == "excellent", axis=1
        )
        is_good = (~is_excellent) & current_part.apply(lambda r: _sim(r) >= RESEG_THRESHOLD, axis=1)
        needs_reseg = ~is_excellent & ~is_good

        n_excellent = int(is_excellent.sum())
        n_good = int(is_good.sum())
        n_reseg = int(needs_reseg.sum())
        print(f"\nTime {current_time}: reseg={n_reseg}, already_good={n_good}, excellent(fixed)={n_excellent}")

        keep_as_is = current_part[~needs_reseg].copy()
        to_reseg_gids = set(current_part.loc[needs_reseg, "GlobalID"].astype(str))

        if not to_reseg_gids:
            print(f"  Nothing to resegment — updating seed and continuing.")
            # Excellent crowns propagate their geometry as seed for next step.
            tile_crowns = current_part.copy()
            continue

        # Use matching rows from tile_crowns (previous step's best geometry) as SAM prompts.
        seed_crowns = tile_crowns[tile_crowns["GlobalID"].astype(str).isin(to_reseg_gids)].copy()
        if seed_crowns.empty:
            print(f"  No seed crowns found for reseg candidates — keeping existing partition.")
            tile_crowns = current_part.copy()
            continue

        time_rows = _segment_crowns(seed_crowns, i)

        if not time_rows:
            print(f"  No new polygons produced — keeping existing partition for {current_time}.")
            tile_crowns = current_part.copy()
            continue

        gdf_reseg = gpd.GeoDataFrame(time_rows, crs=None)
        gdf_reseg["geometry"] = gdf_reseg["geometry"].apply(
            lambda geom: pixel_to_utm(geom, att_transform)
        )
        gdf_reseg = gdf_reseg.set_crs(crs)
        gdf_reseg = crownmap_metrics(original_crownmap=seed_crowns, segmented_crownmap=gdf_reseg)
        gdf_reseg = crown_avoid(gdf_reseg)
        gdf_reseg = crownmap_metrics(original_crownmap=seed_crowns, segmented_crownmap=gdf_reseg)

        ref_geoms = seed_crowns.set_index("GlobalID")["geometry"]
        gdf_reseg["geometry"] = gdf_reseg.apply(_pick_geometry_reseg, axis=1, ref_geoms=ref_geoms)
        fallback_count = int((
            gdf_reseg.apply(
                lambda r: (float(r["similarity"]) if pd.notna(r.get("similarity")) else 0.0) < RESEG_THRESHOLD
                          and r["GlobalID"] in ref_geoms.index,
                axis=1,
            )
        ).sum())
        print(f"  ↩ Fell back to seed for {fallback_count}/{len(gdf_reseg)} reseg crowns")
        gdf_reseg = gpd.GeoDataFrame(gdf_reseg, crs=crs)

        # Merge resegmented crowns back with the kept (Excellent + already-good) crowns.
        result = gpd.GeoDataFrame(
            pd.concat([keep_as_is, gdf_reseg], ignore_index=True), crs=crs
        )
        _write_part(result, bucket_to_process, current_time)
        print(f"  Saved partition: bucket={bucket_to_process}, time={current_time} ({len(result)} rows)")

        # Excellent crowns carry their fixed geometry forward into the next seed.
        tile_crowns = result.copy()

        if seg.device.type == "cuda":
            torch.cuda.empty_cache()
        del gdf_reseg, result
        print("  Cleared GPU cache.")

#####################################################################################
CROWNMAP_PATH = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_2022_2023_crownmap_raw.shp"
tiles_folder= r"D:\BCI_50ha_timeseries\tiles"
BUCKET_GRID_SHAPE = (2, 4)
BUCKET_BUFFER_SIZE = 5

crownmap_gdf, buck, bucket_attributes = build_bucket_attributes(
    crownmap_path=CROWNMAP_PATH,
    tiles_folder=tiles_folder,
    grid_shape=BUCKET_GRID_SHAPE,
    buffer_size=BUCKET_BUFFER_SIZE,
)

dir_address = r"D:\BCI_50ha_timeseries"
tiles_folder = os.path.join(dir_address, "tiles")

cfg = SegmenterConfig.from_yaml(r"C:\Users\vasquezvicente\repo\CanopyRS\canopyrs\config\segmenters\sam3_multi_selvamask_FT.yaml")
seg = Sam3PredictorWrapper(cfg)

bucket_to_process = "0_3"

z = zarr.open(
    os.path.join(tiles_folder, "aligned_local", bucket_to_process, "cube.zarr"),
    mode="r"
)

att_transform = Affine(*z.attrs['transform'])
att_inversed = ~att_transform
time = z.attrs['time']
crs = z.attrs['crs']

height, width = z.shape[2], z.shape[3]

ref = "2022-09-29T00:00:00"
ref_idx = time.index(ref)
backward_indices = list(range(ref_idx - 1, -1, -1))
forward_indices = list(range(ref_idx + 1, len(time)))

# ---------------------------------------------------
# MASTER GDF — PARTITIONED FOLDER
# Each (bucket_id, time) is a separate small parquet file.
# We never read the whole dataset to append a new slice.
# Layout: master_gdf_parts/{bucket_id}/{safe_time}.parquet
# ---------------------------------------------------
master_parts_dir = r"D:\BCI_50ha_timeseries\master_gdf_parts"
os.makedirs(master_parts_dir, exist_ok=True)

# Resegment pass: reference partition must already exist (written by crown_segment).
if not _part_exists(bucket_to_process, time[ref_idx]):
    raise RuntimeError(
        f"Reference partition missing for bucket={bucket_to_process}, time={time[ref_idx]}. "
        "Run 3.crown_segmentv2.py first."
    )
# ---------------------------------------------------
# RUN BACKWARD THEN FORWARD (resegment pass)
# ---------------------------------------------------
process_direction_reseg(backward_indices)
process_direction_reseg(forward_indices)
