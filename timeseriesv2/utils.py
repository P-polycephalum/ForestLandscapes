import numpy as np
from datetime import datetime 
# from arosics import COREG, COREG_LOCAL
import geopandas as gpd
import os
import rasterio
from shapely.geometry import box
import json
from shapely.geometry import box, MultiPolygon
from rasterio.warp import reproject, Resampling
import numpy as np
import tempfile
import shutil
from rasterio.enums import Resampling as ResamplingEnum
from tqdm import tqdm
import time
import numpy as np
from shapely.geometry import box
from matplotlib import pyplot as plt
import pandas as pd

# from sam2.build_sam import build_sam2
# from sam2.sam2_image_predictor import SAM2ImagePredictor
from rasterio.mask import mask
#import torch
import cv2
import gc
from shapely.geometry import GeometryCollection, Polygon, MultiPolygon

def split_gdf_buckets(gdf, x_size=250, y_size=250, n_tiles=None, grid_shape=None):
    """
    Split gdf containing polygons into buckets using centroid assignment.

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

def tile_ortho_with_grid(orthomosaic, grid_info, buffer, output_folder):
    """Split orthomosaic into tiles using a precomputed grid."""
    os.makedirs(output_folder, exist_ok=True)

    with rasterio.open(orthomosaic) as src:
        for idx, row in grid_info.iterrows():
            geom2 = box(row["xmin"] - buffer, row["ymin"] - buffer, row["xmax"] + buffer, row["ymax"] + buffer)
            out_image, out_transform = rasterio.mask.mask(src, [geom2], crop=True)
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            })

            output_filename = f"output_raster_{idx}.tif"
            filename = os.path.join(output_folder, output_filename)
            with rasterio.open(filename, "w", **out_meta) as dest:
                dest.write(out_image)

def crown_segment(tile_folder, shp, checkpoint, model_cfg, device, tile_to_shp=None):  
    sam2_model = build_sam2(model_cfg, checkpoint, device=device)
    predic = SAM2ImagePredictor(sam2_model)
    all=[]
    tiles = sorted(os.listdir(tile_folder))
    for tile in tiles:
        print("processing tile", tile)
        sub=os.path.join(tile_folder,tile)
        with rasterio.open(sub) as src:
            data=src.read()
            transposed_data=data.transpose(1,2,0)
            crs=src.crs
            affine_transform = src.transform 
            bounds=src.bounds
            main_box= box(bounds[0],bounds[1],bounds[2],bounds[3])
        crowns_src = tile_to_shp.get(tile, shp) if tile_to_shp else shp
        if crowns_src is None or crowns_src.empty:
            print(f"No polygon input for tile {tile}. Skipping.")
            continue
        crowns = crowns_src.to_crs(crs)
        mask = crowns['geometry'].within(main_box)
        test_crowns = crowns.loc[mask]

        print("starting box transformation from utm to xy")
        boxes=[]
        for index, row in test_crowns.iterrows():
            if isinstance(row.geometry, MultiPolygon):
                multi_polygon = row.geometry
                polygons = []
                for polygon in multi_polygon.geoms:
                    polygons.append(polygon)
                largest_polygon = max(polygons, key=lambda polygon: polygon.area)
                bounds = largest_polygon.bounds
                boxes.append(bounds)
            else:
                bounds = row.geometry.bounds
                boxes.append(bounds)
        box_mod=[]
        for box1 in boxes:
            xmin, ymin, xmax, ymax = box1
            x_pixel_min, y_pixel_min = ~affine_transform * (xmin, ymin)
            x_pixel_max, y_pixel_max = ~affine_transform * (xmax, ymax)
            trans_box=[x_pixel_min,y_pixel_max,x_pixel_max,y_pixel_min]
            box_mod.append(trans_box)
        if not box_mod:
            print(f"No valid boxes found for tile {tile}. Skipping.")
            continue
        print("The tile contains", len(box_mod), "polygons")

        input_boxes=torch.tensor(box_mod, device=device)
        print("about to set the image")
        predic.set_image(transposed_data[:,:,:3])
        with torch.inference_mode():
            masks, scores, logits = predic.predict(
                box=input_boxes,
                multimask_output=True,
            )
        predic = SAM2ImagePredictor(sam2_model)
        print("finish predicting now getting the utms for transformation")
        height, width, num_bands = transposed_data.shape
        utm_coordinates_and_values = np.empty((height, width, num_bands + 2))
        utm_transform = src.transform
       
        y_coords, x_coords = np.meshgrid(np.arange(height), np.arange(width), indexing='ij')
        utm_x, utm_y = rasterio.transform.xy(utm_transform, y_coords, x_coords)
        utm_coordinates_and_values[..., 0] = np.array(utm_x).reshape(height, width)
        utm_coordinates_and_values[..., 1] = np.array(utm_y).reshape(height, width)
        utm_coordinates_and_values[..., 2:] = transposed_data[..., :num_bands]

        all_polygons=[]
        for idx, (thisscore, thiscrown) in enumerate(zip(scores, masks)):
            maxidx=thisscore.tolist().index(max(thisscore.tolist()))
            thiscrown = thiscrown[maxidx]
            score=scores[1].tolist()[thisscore.tolist().index(max(thisscore.tolist()))]   
            mask = thiscrown.squeeze()
            utm_coordinates = utm_coordinates_and_values[:, :, :2]
            mask_np = mask.astype(np.uint8) if isinstance(mask, np.ndarray) else mask.cpu().numpy().astype(np.uint8)
            contours, _ = cv2.findContours(mask_np, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            polygons = []
            areas = []
            for contour in contours:
                contour_coords = contour.squeeze().reshape(-1, 2)
                contour_utm_coords = utm_coordinates[contour_coords[:, 1], contour_coords[:, 0]]
                if len(contour_utm_coords) >= 3:
                    polygon = Polygon(contour_utm_coords)
                    area = polygon.area
                    polygons.append(polygon)
                    areas.append(area)       
            if len(areas) == 0:
                print(f"No valid areas found for this crown. Skipping.")
                continue  
            largest_index = np.argmax(areas)
            gdf = gpd.GeoDataFrame(geometry=[polygons[largest_index]])
            gdf['area'] = areas[largest_index]
            gdf['score'] = score 
            gdf.crs = src.crs
            tag_value = test_crowns.iloc[idx]['tag']
            global_id= test_crowns.iloc[idx]['global_id']
            gdf['tag']=tag_value
            gdf['global_id']=global_id
            all_polygons.append(gdf)
        print("finish transforming back to utm")
        print(len(all_polygons),"crowns segmented")
        all.append(all_polygons)
        progress= len(all)/len(tiles)
        del input_boxes, masks, scores, logits
        torch.cuda.empty_cache()
        print(progress)
    final_gdfs = []
    for polygons_gdf_list in all:
        combined_gdf = gpd.GeoDataFrame(pd.concat(polygons_gdf_list, ignore_index=True), crs=src.crs)
        final_gdfs.append(combined_gdf)
    final_gdf = gpd.GeoDataFrame(pd.concat(final_gdfs, ignore_index=True), crs=src.crs)
    return final_gdf

def crownmap_QC(original_crownmap, segmented_crownmap, precision_threshold=0.5, recall_threshold=0.5, iou_threshold=0.5):
        """Calcalate quality metrics for segmented crowns and keep original polygons if quality is below threshold.

        Takes a GeoDataFrame of crown polygons.
        """
        crown_qc= crownmap_metrics(original_crownmap, segmented_crownmap)
        for crown in crown_qc['global_id'].unique():
            original_geom = original_crownmap.loc[original_crownmap['global_id'] == crown, 'geometry'].values[0]
            segmented_geom = crown_qc.loc[crown_qc['global_id'] == crown, 'geometry'].values[0]
            iou = crown_qc.loc[crown_qc['global_id'] == crown, 'IoU'].values[0]
            precision = crown_qc.loc[crown_qc['global_id'] == crown, 'Precision'].values[0]
            recall = crown_qc.loc[crown_qc['global_id'] == crown, 'Recall'].values[0]
        # if the crown is below any threshold, keep the original polygon
            if precision < precision_threshold or recall < recall_threshold or iou < iou_threshold:
                print(
                    f"Warning: Crown {crown} has low quality (IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}). "
                    "Keeping original polygon."
                )
                crown_qc.loc[segmented_crownmap['global_id'] == crown, 'geometry'] = original_geom

        return crown_qc

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
            print(f"Converted MultiPolygon to largest Polygon for index {index}.")

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
                    print(f"Polygon {idx} is more similar than {adj_idx}. Subtracting overlap from adjacent.")   
                elif polygon['similarity'] < adj_polygon['similarity']:
                    # Current polygon loses overlap
                    modifications[idx] = modifications.get(idx, polygon.geometry).difference(adj_polygon.geometry)
                    print(f"Polygon {idx} is less similar than {adj_idx}. Subtracting overlap from polygon.")

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

def crownmap_metrics(original_crownmap, segmented_crownmap):
    segmented_ids = set(segmented_crownmap['global_id'].unique())
    original_ids = set(original_crownmap['global_id'].unique())
    missing_in_segmented = original_ids - segmented_ids
    extra_in_segmented = segmented_ids - original_ids
    print(f"Missing IDs in segmented: {missing_in_segmented if missing_in_segmented else 'None'}")
    print(f"Extra IDs in segmented: {extra_in_segmented if extra_in_segmented else 'None'}")

    original_polys = {row['global_id']: row['geometry'] for _, row in original_crownmap.iterrows()}
    segmented_polys = {row['global_id']: row['geometry'] for _, row in segmented_crownmap.iterrows()}

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
        print(f"GlobalID: {crown}, IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}, F1: {f1:.2f}, Score: {score:.2f}")
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'IoU'] = iou
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'Precision'] = precision
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'Recall'] = recall
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'F1'] = f1
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'similarity'] = score
    
    return segmented_crownmap

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

def run_global_alignment_direction(files, start_idx, direction, reference_file, global_dir):
    """
    Run global alignment in one direction from a reference.
    direction: 'backward' or 'forward'
    """
    successful_alignments = [reference_file]
    current_reference = reference_file

    if direction == "backward":
        indices = range(start_idx - 1, -1, -1)
        dir_name = "BACKWARD"
    elif direction == "forward":
        indices = range(start_idx + 1, len(files))
        dir_name = "FORWARD"
    else:
        raise ValueError("direction must be 'backward' or 'forward'")

    for idx in indices:
        target = files[idx]
        orthomosaic_basename = os.path.basename(target)

        print(f"\n[{dir_name}] Processing: {orthomosaic_basename}")

        output_path = os.path.join(global_dir,orthomosaic_basename)

        if os.path.isfile(output_path):
            print(f"  Already processed. Skipping: {orthomosaic_basename}")
            successful_alignments.append(output_path)
            current_reference = output_path
            continue

        file_kwargs = {
            'path_out': output_path,
            'fmt_out': 'GTIFF',
            'r_b4match': 2,
            's_b4match': 2,
            'max_shift': 200,
            'max_iter': 20,
            'align_grids': True,
            'match_gsd': True,
            'binary_ws': False
        }

        target_date = None
        try:
            parts = orthomosaic_basename.split("_")
            target_date = datetime.strptime(f"{parts[2]}-{parts[3]}-{parts[4]}", "%Y-%m-%d")
        except (IndexError, ValueError):
            print(f"  Could not extract date from {orthomosaic_basename}. Will use unsorted candidates.")

        candidate_refs = []
        for ref_file in successful_alignments:
            try:
                fname = os.path.basename(ref_file)
                parts = fname.split("_")
                ref_date = datetime.strptime(f"{parts[2]}-{parts[3]}-{parts[4]}", "%Y-%m-%d")
                days_diff = abs((ref_date - target_date).days) if target_date else float('inf')
                candidate_refs.append((ref_file, days_diff))
            except (IndexError, ValueError):
                candidate_refs.append((ref_file, float('inf')))

        candidate_refs.sort(key=lambda x: x[1])

        alignment_successful = False
        for candidate_ref, days_diff in candidate_refs:
            print(f"  Trying reference: {os.path.basename(candidate_ref)} ({days_diff} days apart)")
            try:
                CR = COREG(candidate_ref, target, **file_kwargs, ws=(2048, 2048))
                CR.calculate_spatial_shifts()
                CR.correct_shifts()
                print(f"  ✓ Aligned successfully using: {os.path.basename(candidate_ref)}")
                successful_alignments.append(output_path)
                current_reference = output_path
                alignment_successful = True
                break
            except Exception as e:
                print(f"  ✗ Failed with {os.path.basename(candidate_ref)}: {str(e)[:120]}")

        if not alignment_successful:
            print(f"  ✗ All candidates failed for {orthomosaic_basename}. Skipping, continuing direction.")

    return successful_alignments, current_reference

def run_global_alignment(files, start_idx, reference_file, global_dir, categorical="both"):
    """
    categorical: 'both', 'backwards', or 'forward'
    """
    mode = categorical.lower().strip()
    if mode not in {"both", "backwards", "forward"}:
        raise ValueError("categorical must be one of: 'both', 'backwards', 'forward'")

    backward_alignments, backward_ref = [reference_file], reference_file
    forward_alignments, forward_ref = [reference_file], reference_file

    if mode in {"both", "backwards"}:
        print("=" * 70)
        print("Starting global alignment BACKWARD from reference...")
        print("=" * 70)
        backward_alignments, backward_ref = run_global_alignment_direction(
            files, start_idx, "backward", reference_file, global_dir
        )

    if mode in {"both", "forward"}:
        print("\n" + "=" * 70)
        print("Starting global alignment FORWARD from reference...")
        print("=" * 70)
        forward_alignments, forward_ref = run_global_alignment_direction(
            files, start_idx, "forward", reference_file, global_dir
        )

    return backward_alignments, forward_alignments, backward_ref, forward_ref



