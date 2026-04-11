# ALIGNMENT 
########################################################################################
# Libraries
import os
import shutil
import warnings
from matplotlib import dates
import rasterio
import numpy as np
from arosics import COREG, COREG_LOCAL
from datetime import datetime
from rasterio.warp import reproject, Resampling
import zarr
import time
import geopandas as gpd
from shapely.geometry import box
from tqdm import tqdm

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
        buffered_box = box(minx - buffer_size, miny - buffer_size, maxx + buffer_size,maxy + buffer_size)
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

#############################################################################################
## runs global alignment in both directions from reference, then local alignment in both directions from reference, with retries on failure
def _extract_date(filename):
    parts = filename.split("_")
    if len(parts) < 5:
        return None
    try:
        year = int(parts[2])
        month = int(parts[3])
        day = int(parts[4])
        return datetime(year, month, day)
    except ValueError:
        return None
def _make_file_kwargs(output_path):
    return {
        'path_out': output_path,
        'fmt_out': 'GTIFF',
        'r_b4match': 1,
        's_b4match': 1,
        'max_shift': 200,
        'max_iter': 20,
        'align_grids': True,
        'match_gsd': True,
        'binary_ws': False,
        'nodata': (0, 0)
    }
def _has_valid_pixels(path, band=1, nodata=0):
    try:
        with rasterio.open(path) as src:
            arr = src.read(band, masked=False)
            if arr.size == 0:
                return False
            return np.any(arr != nodata)
    except Exception:
        return False
def _try_align(target, output_path, successful_alignments_in_global_dir, target_date):
    """
    Try aligning target against all candidates in global_dir sorted by date proximity.
    All candidates are already resolved to global_dir paths.
    Returns True if any candidate succeeded.
    """
    # If target has no valid pixels in the matching band, keep it unchanged
    if not _has_valid_pixels(target, band=1, nodata=0):
        print(f"  Target has no valid pixels in match band. Copying without alignment: {os.path.basename(target)}")
        shutil.copy2(target, output_path)
        return True

    candidate_refs = []
    for ref_file in successful_alignments_in_global_dir:
        ref_date = _extract_date(os.path.basename(ref_file))
        days_diff = abs((ref_date - target_date).days) if (ref_date and target_date) else float('inf')
        candidate_refs.append((ref_file, days_diff))

    candidate_refs.sort(key=lambda x: x[1])

    file_kwargs = _make_file_kwargs(output_path)

    for candidate_ref, days_diff in candidate_refs:
        if not _has_valid_pixels(candidate_ref, band=1, nodata=0):
            print(f"  Skipping empty reference: {os.path.basename(candidate_ref)}")
            continue

        print(f"  Trying reference: {os.path.basename(candidate_ref)} ({days_diff} days apart)")
        try:
            CR = COREG(candidate_ref, target, **file_kwargs, ws=(2048, 2048))
            CR.calculate_spatial_shifts()
            CR.correct_shifts()
            print(f"  ✓ Aligned successfully using: {os.path.basename(candidate_ref)}")
            return True
        except Exception as e:
            print(f"  ✗ Failed with {os.path.basename(candidate_ref)}: {str(e)[:120]}")

    # Final fallback: keep original so file count is complete
    print(f"  No usable reference worked. Copying without alignment: {os.path.basename(target)}")
    shutil.copy2(target, output_path)
    return True
def run_global_alignment_direction(files, start_idx, direction, reference_file_in_global_dir, global_dir):
    """
    Run global alignment in one direction from a reference.
    successful_alignments tracks ONLY paths inside global_dir.
    """
    # Start pool with the reference already in global_dir
    successful_alignments = [reference_file_in_global_dir]
    failed_files = []

    if direction == "backward":
        indices = range(start_idx - 1, -1, -1)
        dir_name = "BACKWARD"
    elif direction == "forward":
        indices = range(start_idx + 1, len(files))
        dir_name = "FORWARD"
    else:
        raise ValueError("direction must be 'backward' or 'forward'")

    for idx in indices:
        target = files[idx]  # raw src — used only as INPUT to COREG, never as reference
        orthomosaic_basename = os.path.basename(target)
        output_path = os.path.join(global_dir, orthomosaic_basename)

        print(f"\n[{dir_name}] Processing: {orthomosaic_basename}")

        if os.path.isfile(output_path):
            print(f"  Already in global_dir. Skipping: {orthomosaic_basename}")
            successful_alignments.append(output_path)
            continue

        target_date = _extract_date(orthomosaic_basename)
        if target_date is None:
            print(f"  Could not extract date from {orthomosaic_basename}. Will use unsorted candidates.")

        ok = _try_align(target, output_path, successful_alignments, target_date)

        if ok:
            successful_alignments.append(output_path)
        else:
            print(f"  ✗ All candidates failed. Adding to retry queue.")
            failed_files.append((target, output_path, target_date))

    # --- RETRY PASS with full pool ---
    if failed_files:
        print(f"\n[{dir_name}] Retrying {len(failed_files)} failed files with full pool ({len(successful_alignments)} refs)...")
        still_failed = []

        for target, output_path, target_date in failed_files:
            orthomosaic_basename = os.path.basename(target)

            if os.path.isfile(output_path):
                successful_alignments.append(output_path)
                continue

            ok = _try_align(target, output_path, successful_alignments, target_date)

            if ok:
                successful_alignments.append(output_path)
                print(f"  ✓ Retry succeeded: {orthomosaic_basename}")
            else:
                still_failed.append(orthomosaic_basename)
                print(f"  ✗ Retry failed: {orthomosaic_basename}")

        if still_failed:
            print(f"\n[{dir_name}] Permanently failed ({len(still_failed)}):")
            for f in still_failed:
                print(f"    - {f}")

    return successful_alignments
def run_global_alignment(files, start_idx, reference_file, global_dir, categorical="both"):
    """
    categorical: 'both', 'backwards', or 'forward'
    """
    mode = categorical.lower().strip()
    if mode not in {"both", "backwards", "forward"}:
        raise ValueError("categorical must be one of: 'both', 'backwards', 'forward'")

    # Copy reference into global_dir FIRST — this is the seed for all alignment chains
    ref_basename = os.path.basename(reference_file)
    ref_in_global = os.path.join(global_dir, ref_basename)
    if not os.path.isfile(ref_in_global):
        shutil.copy2(reference_file, ref_in_global)
        print(f"  Copied reference to global_dir: {ref_basename}")
    else:
        print(f"  Reference already in global_dir: {ref_basename}")

    backward_alignments = [ref_in_global]
    forward_alignments = [ref_in_global]

    if mode in {"both", "backwards"}:
        print("=" * 70)
        print("Starting global alignment BACKWARD from reference...")
        print("=" * 70)
        backward_alignments = run_global_alignment_direction(
            files, start_idx, "backward", ref_in_global, global_dir
        )

    if mode in {"both", "forward"}:
        print("\n" + "=" * 70)
        print("Starting global alignment FORWARD from reference...")
        print("=" * 70)
        forward_alignments = run_global_alignment_direction(
            files, start_idx, "forward", ref_in_global, global_dir
        )

    total_out = len([f for f in os.listdir(global_dir) if f.endswith(".tif")])
    print("\n" + "=" * 70)
    print("Global alignment summary:")
    print(f"  Input files :  {len(files)}")
    print(f"  Output files:  {total_out}")
    print(f"  Missing:       {len(files) - total_out}")
    print("=" * 70)

    return backward_alignments, forward_alignments
#################################################################################################

def run_local_direction(indices, start_reference, source_files, local_dir):
    reference_local = start_reference
    indices_list = list(indices)
    with tqdm(indices_list, desc="Local alignment", unit="file") as pbar:
        for idx in pbar:
            source_file = source_files[idx]
            basename = os.path.basename(source_file)
            pbar.set_postfix(file=basename[:40])

            if not os.path.isfile(source_file):
                tqdm.write(f"  Missing input, skipping: {basename}")
                continue

            out_path = os.path.join(
                local_dir,
                basename.replace(".tif", "_local.tif"),
            )

            if os.path.isfile(out_path):
                tqdm.write(f"  Already processed, skipping: {basename}")
                reference_local = out_path
                continue

            try:
                kwargs_local = {
                    "grid_res": 200,
                    "window_size": (512, 512),
                    "path_out": out_path,
                    "fmt_out": "GTIFF",
                    "q": True,
                    "min_reliability": 30,
                    "r_b4match": 2,
                    "s_b4match": 2,
                    "max_shift": 100,
                    "nodata": (255, 255),
                    "ignore_errors": True,
                    "match_gsd": True,
                    "align_grids": True
                }
                CRL = COREG_LOCAL(reference_local, source_file, **kwargs_local)
                CRL.calculate_spatial_shifts()
                CRL.correct_shifts()
                reference_local = out_path
            except Exception as e:
                tqdm.write(f"  Failed ({basename}): {e}")
                tqdm.write(f"  Copying without local alignment: {basename}")
                shutil.copy2(source_file, out_path)
                reference_local = out_path

###################################################################################################
def _list_tifs(folder):
    return sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith(".tif")
    )

def _read_raster_template(ref_path):
    with rasterio.open(ref_path) as ref:
        return {
            "bands": ref.count,
            "height": ref.height,
            "width": ref.width,
            "dtype": ref.dtypes[0],
            "transform": ref.transform,
            "crs": ref.crs,
            "bounds": ref.bounds,
            "res_x": ref.transform.a,
            "res_y": abs(ref.transform.e),
        }

def _extract_dates_from_tiles(tiles):
    parsed_dates = []
    for tile in tiles:
        d = _extract_date(tile)  # uses your existing helper
        if d is None:
            raise ValueError(f"Could not parse date from filename: {tile}")
        parsed_dates.append(d)
    return parsed_dates

def build_zarr_cube_from_folder(tile_folder, zarr_name="cube.zarr", chunks=(1, -1, 256, 256)):
    tiles = _list_tifs(tile_folder)
    if not tiles:
        raise ValueError(f"No .tif files found in: {tile_folder}")

    ref_path = os.path.join(tile_folder, tiles[0])
    meta = _read_raster_template(ref_path)
    dates = _extract_dates_from_tiles(tiles)

    bands = meta["bands"]
    height = meta["height"]
    width = meta["width"]
    dtype = meta["dtype"]

    # allow -1 for "all bands"
    chunk_bands = bands if chunks[1] == -1 else chunks[1]
    zarr_chunks = (chunks[0], chunk_bands, chunks[2], chunks[3])

    zarr_path = os.path.join(tile_folder, zarr_name)
    cube = zarr.open(
        zarr_path,
        mode="w",
        shape=(len(tiles), bands, height, width),
        chunks=zarr_chunks,
        dtype=dtype,
    )

    cube.attrs["crs"] = str(meta["crs"])
    cube.attrs["transform"] = tuple(meta["transform"])
    cube.attrs["bounds"] = {
        "xmin": meta["bounds"].left,
        "ymin": meta["bounds"].bottom,
        "xmax": meta["bounds"].right,
        "ymax": meta["bounds"].top,
    }
    cube.attrs["resolution"] = (meta["res_x"], meta["res_y"])
    cube.attrs["shape"] = {
        "time": len(tiles),
        "bands": bands,
        "height": height,
        "width": width,
    }
    cube.attrs["time"] = [d.isoformat() for d in dates]
    cube.attrs["tile_files"] = tiles

    for i, tile in enumerate(tiles):
        tile_path = os.path.join(tile_folder, tile)
        with rasterio.open(tile_path) as src:
            cube[i] = src.read()  # (band, y, x)
        print(f"Written {i + 1}/{len(tiles)}: {tile}")

    print("\nDone.")
    print("Cube shape:", cube.shape)
    print("Saved at:", zarr_path)
    return zarr_path

def _get_cube_timesteps(cube_path):
    try:
        cube = zarr.open(cube_path, mode="r")
        if "time" in cube.attrs:
            return len(cube.attrs["time"])
        if hasattr(cube, "shape") and len(cube.shape) > 0:
            return int(cube.shape[0])
    except Exception as exc:
        print(f"Could not read cube timesteps for {cube_path}: {exc}")
    return None

dir_address = r"D:\BCI_50ha_timeseries"
tiles_folder = os.path.join(dir_address, "tiles")
#####################################################################################
CROWNMAP_PATH = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_2022_2023_crownmap_raw.shp"
BUCKET_GRID_SHAPE = (2, 4)
BUCKET_BUFFER_SIZE = 5

crownmap_gdf, buck, bucket_attributes = build_bucket_attributes(
    crownmap_path=CROWNMAP_PATH,
    tiles_folder=tiles_folder,
    grid_shape=BUCKET_GRID_SHAPE,
    buffer_size=BUCKET_BUFFER_SIZE,
)
#########################################################################
reference_token = "BCI_50ha_2022_09_29_orthomosaic" #and a given reference date
EXPECTED_TIMESTEPS = 176

for bucket_to_align in bucket_attributes:
    global_dir = os.path.join(tiles_folder, "aligned_global", bucket_to_align)
    local_dir = os.path.join(tiles_folder, "aligned_local", bucket_to_align)
    global_cube_path = os.path.join(global_dir, "cube.zarr")
    local_cube_path = os.path.join(local_dir, "cube.zarr")

    skip_global = False
    skip_local = False

    if os.path.isdir(global_cube_path): #if global cube exists, check timesteps
        global_steps = _get_cube_timesteps(global_cube_path)
        print(f"Bucket {bucket_to_align} global cube timesteps: {global_steps}")
        if global_steps == EXPECTED_TIMESTEPS:
            skip_global = True
            if os.path.isdir(local_cube_path): #if local cube exists, check timesteps
                local_steps = _get_cube_timesteps(local_cube_path) 
                print(f"Bucket {bucket_to_align} local cube timesteps: {local_steps}")
                if local_steps == EXPECTED_TIMESTEPS:
                    skip_local = True
                    print(
                        f"Bucket {bucket_to_align} has complete global and local cubes. Skipping alignment, keeping existing cubes.")
            if not skip_local:
                # Need global TIFs to run local alignment; check they still exist
                existing_global_tifs = len(_list_tifs(global_dir))
                if existing_global_tifs == EXPECTED_TIMESTEPS:
                    print(
                        f"Bucket {bucket_to_align} has complete global cube but local cube is missing/incomplete. "
                        "Skipping global alignment, reprocessing from global_dir."
                    )
                else:
                    skip_global = False
                    print(
                        f"Bucket {bucket_to_align} has complete global cube but global TIFs are missing "
                        f"({existing_global_tifs} found). Reprocessing global alignment."
                    )
        else:
            print(
                f"Bucket {bucket_to_align} global cube is incomplete ({global_steps}). "
                "Reprocessing from global_dir."
            )

    os.makedirs(global_dir, exist_ok=True)
    os.makedirs(local_dir, exist_ok=True)

    if bucket_to_align not in bucket_attributes:
        raise ValueError(f"Bucket {bucket_to_align} not found.")

    if not skip_global:
        files = bucket_attributes[bucket_to_align]["files"]
        if not files:
            raise ValueError(f"No files found for bucket {bucket_to_align}.")

        ref_index = next(
            ((i, f) for i, f in enumerate(files) if reference_token in os.path.basename(f)),
            None,
        )
        if ref_index is None:
            raise ValueError(
                f"Reference token '{reference_token}' not found in bucket {bucket_to_align}."
            )

        reference_file = files[ref_index[0]]
        print(f"Bucket {bucket_to_align}: {len(files)} files")
        print(f"Reference index: {ref_index} -> {os.path.basename(reference_file)}")

        backward_alignments, forward_alignments = run_global_alignment(
            files=files,
            start_idx=ref_index[0],
            reference_file=reference_file,
            global_dir=global_dir,
            categorical="both"
        )

    if not skip_local:    
        files_local = sorted([
            os.path.join(global_dir, f)
            for f in os.listdir(global_dir)
            if f.lower().endswith(".tif")
            ])

        if not files_local:
            raise ValueError(f"No global-aligned files found in {global_dir}")

        local_ref_index = next(
            (i for i, f in enumerate(files_local) if reference_token in os.path.basename(f)),
            None,
        )
        if local_ref_index is None:
            raise ValueError(f"Reference token '{reference_token}' not found in global_dir.")

        local_reference_file = files_local[local_ref_index]

        backward_indices = range(local_ref_index - 1, -1, -1)
        forward_indices = range(local_ref_index + 1, len(files_local))

        print("\nStarting local alignment backward from reference...")
        run_local_direction(backward_indices, local_reference_file, files_local, local_dir)

        print("\nStarting local alignment forward from reference...")
        run_local_direction(forward_indices, local_reference_file, files_local, local_dir)

        # copy reference into local output with _local suffix
        shutil.copy2(
            local_reference_file,
            os.path.join(
                local_dir,
                os.path.basename(local_reference_file).replace(".tif", "_local.tif")
            )
        )

    #test the length of the aligned tifs
    tile_folder = os.path.join(tiles_folder, "aligned_global", bucket_to_align)
    if len(_list_tifs(tile_folder)) == EXPECTED_TIMESTEPS:
        build_zarr_cube_from_folder(tile_folder=tile_folder)
        global_cube_path = os.path.join(global_dir, "cube.zarr") #after building the cube, check again if it has the expected timesteps, if so we can erase the tifs to save space, if not we keep them for debugging
        if os.path.isdir(global_cube_path):
            global_steps = _get_cube_timesteps(global_cube_path)
            print(f"Bucket {bucket_to_align} global cube timesteps: {global_steps}")
            if global_steps == EXPECTED_TIMESTEPS:
                tile_folder_aligned= os.path.join(tiles_folder, "aligned_local", bucket_to_align)
                if len(_list_tifs(tile_folder_aligned)) == EXPECTED_TIMESTEPS:
                    build_zarr_cube_from_folder(tile_folder=tile_folder_aligned)
                    local_cube_path = os.path.join(local_dir, "cube.zarr")
                    if os.path.isdir(local_cube_path):
                        local_steps = _get_cube_timesteps(local_cube_path)
                        print(f"Bucket {bucket_to_align} local cube timesteps: {local_steps}")
                        if local_steps == EXPECTED_TIMESTEPS:
                            print(f"Bucket {bucket_to_align} local cube complete. Erasing aligned_local tifs to save space.")
                            for f in os.listdir(tile_folder_aligned):
                                if f.lower().endswith(".tif"):
                                    os.remove(os.path.join(tile_folder_aligned, f))
                            print(f"Bucket {bucket_to_align} global cube complete. Erasing aligned_global tifs to save space.")
                            for f in os.listdir(tile_folder):
                                if f.lower().endswith(".tif"):
                                    os.remove(os.path.join(tile_folder, f))
                        else:
                            print(f"Bucket {bucket_to_align} local cube is incomplete ({local_steps}). Keeping local tifs for debugging.")
                    else:
                        print(f"Bucket {bucket_to_align} local cube was not saved correctly. Keeping local tifs for debugging.")
                else:
                    print(f"Bucket {bucket_to_align} local files are incomplete ({len(_list_tifs(tile_folder_aligned))} tifs). Keeping local tifs for debugging.")
            else:
                print(f"Bucket {bucket_to_align} global cube is incomplete ({global_steps}). Keeping tifs for debugging.")
        else:
            print(f"Bucket {bucket_to_align} global cube wasnt saved correctly. Keeping tifs for debugging.")
    else:
        print(f"Bucket {bucket_to_align} has incomplete global alignment ({len(_list_tifs(tile_folder))} tifs). Skipping cube build and keeping tifs for debugging.")

