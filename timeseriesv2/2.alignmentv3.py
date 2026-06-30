# ALIGNMENT
########################################################################################
import os
import shutil
import rasterio
import numpy as np
from arosics import COREG
from datetime import datetime
import zarr
import geopandas as gpd
from shapely.geometry import box

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
        'r_b4match': 2,
        's_b4match': 2,
        'max_iter': 20,      
        'align_grids': True,
        'match_gsd': True,
        'binary_ws': False,
        'nodata': (None, None)
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
        zarr_format=2,
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

CROWNMAP_PATH     = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_2022_2023_crownmap_raw.shp"
BUCKET_GRID_SHAPE = (2, 4)
BUCKET_BUFFER_SIZE = 5

crownmap_gdf, buck, bucket_attributes = build_bucket_attributes(
    crownmap_path=CROWNMAP_PATH,
    tiles_folder=tiles_folder,
    grid_shape=BUCKET_GRID_SHAPE,
    buffer_size=BUCKET_BUFFER_SIZE,
)

# ── Configuration ──────────────────────────────────────────────────────────────
FILE_SUBSET     = slice(None)   # all files
FIRST_REF_TOKEN = "2022_09_29"  # date token of the anchor image
MAX_FALLBACKS   = 5             # max number of already-aligned anchors to try per file
# ──────────────────────────────────────────────────────────────────────────────

def _align_global(ref_path, src_path, out_path, max_shift=100):
    """
    Align src to ref using COREG (global), writing to out_path.
    Returns True on success, False on failure.
    """
    try:
        CR = COREG(ref_path, src_path, **_make_file_kwargs(out_path),
                   ws=(2048, 2048), max_shift=max_shift)
        CR.calculate_spatial_shifts()

        reliability = getattr(CR, "shift_reliability", None)
        CR.correct_shifts()
        rel_str = f"{reliability:.1f}%" if reliability is not None else "n/a"
        print(f"    ✓ reliability={rel_str}  max_shift={max_shift}")
        return True

    except Exception as e:
        print(f"    ✗ COREG failed (max_shift={max_shift}): {str(e)[:120]}")
        if os.path.isfile(out_path):
            os.remove(out_path)
        return False


def _try_align(idx, source_files, out_folder, basenames, max_fallbacks=5):
    """
    Try to align source_files[idx] using up to max_fallbacks already-aligned
    files as reference, sorted by index proximity.  Tries max_shift 100 then 200
    for each candidate.  Returns True on success, False if all attempts exhausted.
    """
    basename = basenames[idx]
    out_path = os.path.join(out_folder, basename)

    if not _has_valid_pixels(source_files[idx], band=1, nodata=0):
        print(f"  ✗ No valid pixels: {basename}")
        return False

    aligned_by_dist = sorted(
        (i for i, b in enumerate(basenames)
         if os.path.isfile(os.path.join(out_folder, b))),
        key=lambda i: abs(i - idx)
    )[:max_fallbacks]

    for attempt, ref_cand in enumerate(aligned_by_dist):
        ref_path = os.path.join(out_folder, basenames[ref_cand])
        print(f"  [{attempt+1}/{len(aligned_by_dist)}] [{ref_cand}]→[{idx}] {basename}")
        for ms in [100, 200]:
            ok = _align_global(ref_path, source_files[idx], out_path, max_shift=ms)
            if ok:
                return True
            print(f"    ↓ max_shift={ms} failed")

    print(f"  ✗ Exhausted {len(aligned_by_dist)} anchors for {basename} — skipping (retry next run)")
    return False


for bucket_to_align in bucket_attributes:
    print(f"\n{'='*70}")
    print(f"=== Processing bucket: {bucket_to_align} ===")
    print(f"{'='*70}")

    out_folder = os.path.join(tiles_folder, "aligned_global", bucket_to_align)
    cube_path  = os.path.join(out_folder, "cube.zarr")

    if os.path.isdir(cube_path):
        print(f"  Cube already exists. Skipping bucket.")
        continue

    os.makedirs(out_folder, exist_ok=True)

    raw_files = bucket_attributes[bucket_to_align]["files"][FILE_SUBSET]
    if not raw_files:
        print(f"  No files in subset. Skipping.")
        continue

    source_files = sorted(
        raw_files,
        key=lambda f: _extract_date(os.path.basename(f)) or datetime.min
    )
    N         = len(source_files)
    basenames = [os.path.basename(f) for f in source_files]
    print(f"  Total files: {N}")

    # ── locate anchor ──────────────────────────────────────────────────────────────
    first_ref_token = f"BCI_50ha_{FIRST_REF_TOKEN}_orthomosaic_tile_{bucket_to_align}.tif"
    ref_idx = next((i for i, b in enumerate(basenames) if b == first_ref_token), None)
    if ref_idx is None:
        print(f"  Anchor '{first_ref_token}' not found. Skipping bucket.")
        continue
    if not _has_valid_pixels(source_files[ref_idx], band=1, nodata=0):
        print(f"  Anchor has no valid pixels. Skipping bucket.")
        continue
    ref_out = os.path.join(out_folder, first_ref_token)
    if not os.path.isfile(ref_out):
        shutil.copy2(source_files[ref_idx], ref_out)
    print(f"  Anchor: [{ref_idx}] {first_ref_token}")

    # ── backward chain: ref_idx-1 → 0 ────────────────────────────────────────────────
    print(f"\n  === Backward pass ===")
    for i in range(ref_idx - 1, -1, -1):
        out_path = os.path.join(out_folder, basenames[i])
        if os.path.isfile(out_path):
            continue
        aligned = sum(1 for b in basenames if os.path.isfile(os.path.join(out_folder, b)))
        print(f"\n  [{i}] {basenames[i]}  ({aligned}/{N} aligned)")
        _try_align(i, source_files, out_folder, basenames, MAX_FALLBACKS)

    # ── forward chain: ref_idx+1 → N-1 ──────────────────────────────────────────────
    print(f"\n  === Forward pass ===")
    for i in range(ref_idx + 1, N):
        out_path = os.path.join(out_folder, basenames[i])
        if os.path.isfile(out_path):
            continue
        aligned = sum(1 for b in basenames if os.path.isfile(os.path.join(out_folder, b)))
        print(f"\n  [{i}] {basenames[i]}  ({aligned}/{N} aligned)")
        _try_align(i, source_files, out_folder, basenames, MAX_FALLBACKS)

    final_aligned = sum(1 for b in basenames if os.path.isfile(os.path.join(out_folder, b)))
    print(f"\n  Done. {final_aligned}/{N} files aligned. Building cube...")
    if final_aligned > 0:
        build_zarr_cube_from_folder(tile_folder=out_folder)
        print(f"  Cube built: {_get_cube_timesteps(cube_path)} timesteps.")
    else:
        print(f"  No aligned files. Skipping cube build.")


