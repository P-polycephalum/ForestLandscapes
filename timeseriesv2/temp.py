import zarr
import dask.array as da
import napari
import numpy as np
import os
import geopandas as gpd
import pandas as pd
from affine import Affine
from shapely.geometry import MultiPolygon

master_parts_dir = r"D:\BCI_50ha_timeseries\master_gdf_parts"

# ── helpers ──────────────────────────────────────────────────────────────────

def _safe_time(time_str):
    """Convert time string to safe filename component."""
    return time_str.replace(":", "-").replace("/", "-")

def read_part_by_time(bucket_id, time_str):
    """Read a single parquet partition for one bucket and one time step."""
    part_path = os.path.join(master_parts_dir, str(bucket_id), f"{_safe_time(time_str)}.parquet")
    if not os.path.exists(part_path):
        return gpd.GeoDataFrame()
    try:
        part_gdf = gpd.read_parquet(part_path)
    except Exception as e:
        print(f"⚠️ Error reading partition {part_path}: {e}")
        return gpd.GeoDataFrame()
    if not part_gdf.empty:
        computed_crown_area = part_gdf.geometry.area
        if "crown_area" not in part_gdf.columns:
            part_gdf["crown_area"] = computed_crown_area
        else:
            missing = part_gdf["crown_area"].isna()
            if missing.any():
                part_gdf.loc[missing, "crown_area"] = computed_crown_area.loc[missing]
    return gpd.GeoDataFrame(part_gdf, crs="EPSG:32617")

def read_all_parts(bucket_id):
    """Read all parquet partitions for a tile (all time steps)."""
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    if not os.path.isdir(bucket_dir):
        print(f"⚠️ Bucket directory not found: {bucket_dir}")
        return gpd.GeoDataFrame()
    files = sorted(f for f in os.listdir(bucket_dir) if f.endswith(".parquet"))
    if not files:
        print(f"⚠️ No parquet files found in {bucket_dir}")
        return gpd.GeoDataFrame()
    parts = []
    for f in files:
        try:
            parts.append(gpd.read_parquet(os.path.join(bucket_dir, f)))
        except Exception as e:
            print(f"⚠️ Error reading {f}: {e}")
    if not parts:
        return gpd.GeoDataFrame()
    gdf = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs="EPSG:32617")
    if "crown_area" not in gdf.columns:
        gdf["crown_area"] = gdf.geometry.area
    return gdf

def _geom_to_pixel_parts(geom, inv_transform):
    """Convert a shapely geometry to a list of (row, col) pixel coord arrays."""
    if geom is None or geom.is_empty:
        return []
    polys = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    out = []
    for poly in polys:
        coords = np.array(poly.exterior.coords)
        px, py = inv_transform * (coords[:, 0], coords[:, 1])
        out.append(np.column_stack([py, px]))
    return out

def build_pixel_shapes_by_time(tile_gdf, inv_transform, img_h, img_w):
    """Return dict: time_str -> list of pixel-coord arrays for all polygons at that time."""
    pixel_shapes_by_time = {}
    tile_gdf.geometry = tile_gdf.geometry.simplify(0.1)
    for time_str, group in tile_gdf.groupby("time"):
        shapes = []
        for _, row in group.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue
            for coords_px in _geom_to_pixel_parts(geom, inv_transform):
                # Skip polygons entirely outside the image
                if (np.all(coords_px[:, 0] < 0) or np.all(coords_px[:, 0] > img_h) or
                        np.all(coords_px[:, 1] < 0) or np.all(coords_px[:, 1] > img_w)):
                    continue
                shapes.append(coords_px)
        pixel_shapes_by_time[str(time_str)] = shapes
    return pixel_shapes_by_time

# ── cube setup ────────────────────────────────────────────────────────────────

cube_root = r"D:\BCI_50ha_timeseries\tiles"
tile_type = "aligned_global"
tile_id = "0_0"
path = os.path.join(cube_root, tile_type, tile_id, "cube.zarr")
print(f"Opening {path}...")

if not os.path.exists(path):
    print(f"File not found: {path}")

z = zarr.open(path, mode="r")
times = list(z.attrs.get("time", []))
print(f"Shape: {z.shape}  —  {len(times)} time steps")

att_transform = Affine(*z.attrs["transform"])
inv_transform = ~att_transform

# (T, C, H, W) -> (T, H, W, C), keep RGB only
rgb = da.from_zarr(z)[:, :3, :, :]
rgb = da.moveaxis(rgb, 1, -1)
img_h, img_w = int(rgb.shape[1]), int(rgb.shape[2])

# ── load polygons ─────────────────────────────────────────────────────────────

print(f"Loading all polygons for tile {tile_id}...")
tile_gdf = read_all_parts(tile_id)
print(f"  {len(tile_gdf)} rows loaded across {tile_gdf['time'].nunique() if not tile_gdf.empty else 0} time steps.")

pixel_shapes_by_time = {}
if not tile_gdf.empty:
    pixel_shapes_by_time = build_pixel_shapes_by_time(tile_gdf, inv_transform, img_h, img_w)
    print(f"  Precomputed pixel shapes for {len(pixel_shapes_by_time)} time steps.")

# ── napari viewer ─────────────────────────────────────────────────────────────

viewer = napari.Viewer()
viewer.add_image(rgb, rgb=True, name=f"{tile_type}/{tile_id}")

shapes_layer = viewer.add_shapes(
    [],
    shape_type="polygon",
    edge_color="red",
    face_color="transparent",
    edge_width=6,
    name="polygons",
    opacity=1,
)

last_t_idx = {"value": None}

def _on_step(event=None):
    t = int(viewer.dims.current_step[0])
    label = times[t] if t < len(times) else str(t)
    viewer.title = f"t={t}  {label}"

    if t == last_t_idx["value"]:
        return
    last_t_idx["value"] = t

    time_str = str(times[t]) if t < len(times) else None
    shapes = pixel_shapes_by_time.get(time_str, []) if time_str else []
    shapes_layer.data = shapes
    if shapes:
        shapes_layer.edge_width = 5
        shapes_layer.edge_color = "red"
        shapes_layer.opacity = 1.0

viewer.dims.events.current_step.connect(_on_step)
_on_step(None)

napari.run()
