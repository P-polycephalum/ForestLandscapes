from arosics import COREG, COREG_LOCAL
from datetime import datetime
import rasterio
import random
import shutil
import zarr
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import os
import random
import pandas as pd
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")

# HELPERS
# -----------------------------------
def _extract_date(filename):
    parts = filename.split("_")
    try:
        year = int(parts[0][0:4])
        month = int(parts[0][4:6])
        day = int(parts[0][6:8])
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

# -----------------------------------
# CORE ALIGNMENT FUNCTION
# -----------------------------------
def _try_align(target, output_path, successful_alignments_in_global_dir, target_date):

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
            CR = COREG(candidate_ref, target, **file_kwargs, ws=(512, 512))
            CR.calculate_spatial_shifts()

            # -----------------------------------
            # EXTRACT SHIFTS
            # -----------------------------------
            dx_map = CR.x_shift_map
            dy_map = CR.y_shift_map
            dx_px = CR.x_shift_px
            dy_px = CR.y_shift_px

            # -----------------------------------
            # LOG SHIFTS
            # -----------------------------------
            shifts_log.append({
                "target": os.path.basename(target),
                "reference": os.path.basename(candidate_ref),
                "days_diff": days_diff,
                "dx_map": dx_map,
                "dy_map": dy_map,
                "dx_px": dx_px,
                "dy_px": dy_px
            })

            # -----------------------------------
            # APPLY SHIFT
            # -----------------------------------
            CR.correct_shifts()

            print(f"  ✓ Aligned successfully using: {os.path.basename(candidate_ref)}")
            print(f"    Shift (map): ({dx_map:.3f}, {dy_map:.3f})")

            return True

        except Exception as e:
            print(f"  ✗ Failed with {os.path.basename(candidate_ref)}: {str(e)[:120]}")

    return True

# -----------------------------------
# MAIN
# -----------------------------------

n_iterations = 4

base_path = r"C:\Users\vasquezvicente\Downloads\planet_examples\planet_examples"
out_root = r"C:\Users\vasquezvicente\Downloads\planet_examples\coreg_iterative"
os.makedirs(out_root, exist_ok=True)

all_shifts = []

current_path = base_path

for it in range(n_iterations):

    print(f"\n================ ITERATION {it+1} ================")

    out_folder = os.path.join(out_root, f"iter_{it+1}")
    os.makedirs(out_folder, exist_ok=True)

    list_files = [f for f in os.listdir(current_path) if f.endswith(".tif") and "rgb" in f]

    random_ref = random.choice(list_files)

    successful_alignments_in_global_dir = []
    successful_alignments_in_global_dir.append(os.path.join(current_path, random_ref))

    shifts_log = []

    for f in list_files:
        target = os.path.join(current_path, f)
        output_path = os.path.join(out_folder, f.replace(".tif", f"_iter{it+1}.tif"))

        try:
            _try_align(target, output_path, successful_alignments_in_global_dir, _extract_date(f))
            successful_alignments_in_global_dir.append(output_path)
        except Exception as e:
            print(f"Failed to align {f}: {str(e)[:120]}")

    df = pd.DataFrame(shifts_log)
    df["iteration"] = it + 1

    all_shifts.append(df)

    current_path = out_folder  # <-- key: feed next iteration

# -----------------------------------
# COMBINE ALL
# -----------------------------------
df_all = pd.concat(all_shifts, ignore_index=True)

csv_path = os.path.join(out_root, "shifts_all_iterations.csv")
df_all.to_csv(csv_path, index=False)

print(f"\nSaved: {csv_path}")

# -----------------------------------
# PLOT
# -----------------------------------
plt.figure(figsize=(8, 6))

for it in range(1, n_iterations + 1):
    dfi = df_all[df_all["iteration"] == it]
    plt.scatter(dfi["dx_map"], dfi["dy_map"], s=60, label=f"iter {it}")

plt.axhline(0)
plt.axvline(0)
plt.xlabel("dx_map (m)")
plt.ylabel("dy_map (m)")
plt.legend()
plt.title("Shift convergence across iterations")
plt.show()


#which iteration had the best convergence?
iteration_stats = df_all.groupby("iteration").agg(
    dx_mean=("dx_map", lambda x: x.abs().mean()),
    dx_std=("dx_map", lambda x: x.abs().std()),
    dy_mean=("dy_map", lambda x: x.abs().mean()),
    dy_std=("dy_map", lambda x: x.abs().std()),
).reset_index()

iteration_stats

#maybe align the grids if the dimensions are different, lets try without it.
from osgeo import gdal
def align_rasters_inplace(rasters):
    if not rasters:
        print("No rasters provided.")
        return False

    ref_path = rasters[0]
    ref_ds = gdal.Open(ref_path, gdal.GA_ReadOnly)
    if ref_ds is None:
        print(f"Could not open reference raster: {ref_path}")
        return False

    gt = ref_ds.GetGeoTransform()
    proj = ref_ds.GetProjection()
    width, height = ref_ds.RasterXSize, ref_ds.RasterYSize

    xmin, ymax = gt[0], gt[3]
    xres, yres = gt[1], abs(gt[5])
    xmax = xmin + width * xres
    ymin = ymax - height * yres

    print("Reference grid:")
    print(f"  File: {os.path.basename(ref_path)}")
    print(f"  Extent: {xmin}, {ymin}, {xmax}, {ymax}")
    print(f"  Resolution: {xres}, {yres}")
    print(f"  Size: {width}, {height}\n")

    ref_ds = None  # release lock on Windows

    for i, raster_path in enumerate(rasters):
        if i == 0:
            print(f"Skipping reference raster: {os.path.basename(raster_path)}")
            continue

        print(f"Aligning in place: {os.path.basename(raster_path)}")
        temp_path = raster_path.replace(".tif", "._tmp_align.tif")

        if os.path.exists(temp_path):
            os.remove(temp_path)

        warp_ds = gdal.Warp(
            destNameOrDestDS=temp_path,
            srcDSOrSrcDSTab=raster_path,
            format="GTiff",
            outputBounds=(xmin, ymin, xmax, ymax),
            xRes=xres,
            yRes=yres,
            width=width,
            height=height,
            dstSRS=proj,
            targetAlignedPixels=True,
            resampleAlg=gdal.GRA_Bilinear,
            multithread=True,
        )

        if warp_ds is None:
            raise RuntimeError(f"GDAL Warp failed for {raster_path}")

        warp_ds.FlushCache()
        warp_ds = None

        replaced = False
        for _ in range(5):
            try:
                os.replace(temp_path, raster_path)
                replaced = True
                break
            except PermissionError:
                time.sleep(0.5)

        if not replaced:
            raise PermissionError(f"Could not replace locked file: {raster_path}")

    print("\nAll rasters aligned in place.")
    return True

out_folder= r"C:\Users\vasquezvicente\Downloads\planet_examples\coreg_iterative\iter_2"
tiles = sorted([f for f in os.listdir(out_folder) if f.endswith(".tif")])
rasters = [os.path.join(out_folder, f) for f in tiles]
align_rasters_inplace(rasters)

#lets create the zarr cube from the aligned images
tiles= sorted([f for f in os.listdir(out_folder) if f.endswith(".tif")])
ref_path = os.path.join(out_folder, tiles[0])

with rasterio.open(ref_path) as ref:
    bands = ref.count
    height = ref.height
    width = ref.width
    dtype = ref.dtypes[0]

    transform = ref.transform
    crs = ref.crs
    bounds = ref.bounds
    res_x = transform.a
    res_y = abs(transform.e)


dates = []

for tile in tiles:
    date = _extract_date(tile)
    dates.append(date)

zarr_path = os.path.join(out_folder, "cube.zarr")
cube = zarr.open(
    zarr_path,
    mode="w",
    shape=(len(tiles), bands, height, width),
    chunks=(1, bands, 256, 256), 
    dtype=dtype
)

cube.attrs["crs"] = str(crs)
cube.attrs["transform"] = tuple(transform)

cube.attrs["bounds"] = {
    "xmin": bounds.left,
    "ymin": bounds.bottom,
    "xmax": bounds.right,
    "ymax": bounds.top
}

cube.attrs["resolution"] = (res_x, res_y)
cube.attrs["shape"] = {
    "time": len(tiles),
    "bands": bands,
    "height": height,
    "width": width
}
cube.attrs["time"] = [d.isoformat() for d in dates]

# --- Write each raster ---
for i, tile in enumerate(tiles):
    tile_path = os.path.join(out_folder, tile)

    with rasterio.open(tile_path) as src:
        data = src.read()  # (band, y, x)

    cube[i] = data  # ✅ correct

    print(f"Written {i+1}/{len(tiles)}: {tile}")
print("\nDone.")
print("Cube shape:", cube.shape)
print("Saved at:", zarr_path)


import dask.array as da
import napari
import zarr
zarr_path = r"C:\Users\vasquezvicente\Downloads\planet_examples\coreg_iterative\iter_2\cube.zarr"
z = zarr.open(zarr_path, mode="r")

rgb = da.from_zarr(z)[:, 0:3, :, :]
rgb = da.moveaxis(rgb, 1, -1)

times = list(z.attrs.get("time", []))

viewer = napari.Viewer()

viewer.add_image(
    rgb,
    rgb=True,
    name="cube",
    gamma=1.0
)