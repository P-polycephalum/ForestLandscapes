from time import time
import numpy as np
from torch import layout
import zarr
import dask.array as da
import napari
import geopandas as gpd
import pandas as pd
from affine import Affine
from shapely.geometry import Polygon, MultiPolygon
import os
from shapely.ops import transform as shp_transform
from qtpy.QtWidgets import (QPushButton, QVBoxLayout, QWidget, QInputDialog, 
                            QSlider, QLabel, QRadioButton, QButtonGroup, QHBoxLayout, QComboBox,
                            QLineEdit)
from qtpy.QtCore import QTimer, Qt

cube_root = r"D:\BCI_50ha_timeseries\tiles\aligned_local"
master_parts_dir = r"D:\BCI_50ha_timeseries\master_gdf_parts"  # new partitioned folder
# Camera control: 1.0 means "fit to view", >1 zooms in, <1 zooms out.
IMAGE_ZOOM_FACTOR = 5.0
# Zoom target in normalized image coordinates.
# (0.0, 0.0) = top-left, (1.0, 1.0) = bottom-right, (0.5, 0.5) = image center.
IMAGE_ZOOM_TARGET_X = 0.5
IMAGE_ZOOM_TARGET_Y = 0.5
# Dynamic zoom based on crown area: min/max zoom factors and area scaling.
MIN_ZOOM = 2.0
MAX_ZOOM = 15.0
# Area scaling factor: adjust to tune sensitivity. Higher = more zoom variation.
AREA_SCALING_FACTOR = 400 
######################################


def read_part_filtered(bucket_id, global_id=None):
    """Read all parquets for a bucket from the partitioned folder, optionally filtered by GlobalID.
    
    Args:
        bucket_id: Bucket identifier (e.g., tile name)
        global_id: Optional GlobalID to filter by. If None, returns all rows for the bucket.
    
    Returns:
        GeoDataFrame with all time slices for the bucket, optionally filtered by GlobalID.
    """
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    if not os.path.isdir(bucket_dir):
        print(f"⚠️ Bucket directory not found: {bucket_dir}")
        return gpd.GeoDataFrame()
    
    files = sorted(f for f in os.listdir(bucket_dir) if f.endswith(".parquet"))
    if not files:
        print(f"⚠️ No parquet files found in {bucket_dir}")
        return gpd.GeoDataFrame()
    
    # Read all time slices for this bucket
    if global_id is not None:
        parts = []
        for f in files:
            try:
                part = gpd.read_parquet(os.path.join(bucket_dir, f),filters=[("GlobalID", "==", global_id)]) 
                parts.append(part)
            except Exception as e:
                print(f"⚠️ Error reading {f}: {e}")
                continue
        if parts:
            bucket_gdf = gpd.GeoDataFrame(
                pd.concat(parts, ignore_index=True),
                crs="EPSG:32617"
            )
        else:
            bucket_gdf = gpd.GeoDataFrame()
    else:
        print(f"No global_id specified, reading entire bucket {bucket_id} for all trees...")
        parts = []
        for f in files:
            try:
                part = gpd.read_parquet(os.path.join(bucket_dir, f))
                parts.append(part)
            except Exception as e:
                print(f"⚠️ Error reading {f}: {e}")
                continue
        if parts:
            bucket_gdf = gpd.GeoDataFrame(
                pd.concat(parts, ignore_index=True),
                crs="EPSG:32617"
            )
        else:
            bucket_gdf = gpd.GeoDataFrame()

    if not bucket_gdf.empty:
        computed_crown_area = bucket_gdf.geometry.area
        if "crown_area" not in bucket_gdf.columns:
            bucket_gdf["crown_area"] = computed_crown_area
        else:
            missing_crown_area = bucket_gdf["crown_area"].isna()
            if missing_crown_area.any():
                bucket_gdf.loc[missing_crown_area, "crown_area"] = computed_crown_area.loc[missing_crown_area]

    return bucket_gdf

# Load species/tree metadata for this tile.

list_of_tiles = sorted(os.listdir(cube_root))
all_parquet_files = []
for tile in list_of_tiles:
    if not os.path.isdir(os.path.join(master_parts_dir, tile)):
        continue
    tile_gdf = gpd.read_parquet(os.path.join(master_parts_dir, tile, "2022-09-29T00-00-00.parquet"))
    all_parquet_files.append(tile_gdf)

species_df = pd.concat(all_parquet_files, ignore_index=True)

species_df["latin_str"] = species_df["latin"].astype(str).str.strip()
species_df.loc[species_df["latin"].isna(), "latin_str"] = ""
species_df["globalid_str"] = species_df["GlobalID"].astype(str).str.strip()
species_df.loc[species_df["GlobalID"].isna(), "globalid_str"] = ""
species_df["bucket_id_str"] = species_df["bucket_id"].astype(str).str.strip()
species_df.loc[species_df["bucket_id"].isna(), "bucket_id_str"] = ""

species = sorted(
    {
        s
        for s in species_df["latin_str"].unique()
        if s
    }
)
print(f"Loaded species metadata for {len(species)} unique species and {len(species_df)} trees.")

##########################################

def geom_to_pixels(geom, inv_transform):
    coords_world = np.array(geom.exterior.coords)
    coords_pixel = np.array([inv_transform * (x, y) for x, y in coords_world])
    return coords_pixel[:, ::-1]

# ---------------------------------------------------
# CUBE STATE (LAZY LOAD)
# ---------------------------------------------------
tile_id = ""
path = ""
z = None
rgb = None
times = []

#---------------------------------------------------

viewer = napari.Viewer()
image_layer = None
att_transform = None
inv_transform = None
shapes_layer = None
last_t_idx = {"value": None}
pixel_shapes_by_time = {}  # dict: time_str -> list of pixel coord arrays (precomputed at tile load)
pixel_shape_records_by_time = {}  # dict: time_str -> list[(globalid_str, pixel coords)]
tree_centroid_by_id = {}  # dict: GlobalID str -> (row, col) centroid in pixel coordinates
tree_crown_area_by_id = {}  # dict: GlobalID str -> crown_area (square meters)
current_selected_tree_id = ""
pending_polygon_edits = {}  # dict: (bucket_id, time_str, globalid_str) -> geometry
pending_tree_updates = {}  # dict: (bucket_id, time_str, globalid_str) -> metadata updates
original_geometry_by_key = {}  # dict: (bucket_id, time_str, globalid_str) -> original geometry
loaded_edited_by_key = {}  # dict: (time_str, globalid_str) -> bool edited flag from master_gdf
loaded_labels_by_key = {}  # dict: (time_str, globalid_str) -> {leafing, flowering, quality, ignore, notes}

def compute_dynamic_zoom(crown_area):
    """Compute zoom factor inversely proportional to crown area.
    Larger crowns get smaller zoom; smaller crowns get larger zoom.
    """
    if crown_area is None or crown_area <= 0:
        return IMAGE_ZOOM_FACTOR
    # Inverse scaling: zoom = MAX_ZOOM / (1 + area/scaling_factor)
    zoom = MAX_ZOOM / (1.0 + float(crown_area) / AREA_SCALING_FACTOR)
    zoom = max(MIN_ZOOM, min(MAX_ZOOM, zoom))  # Clamp to [MIN_ZOOM, MAX_ZOOM]
    return zoom


def set_camera_view(zoom_factor=IMAGE_ZOOM_FACTOR, center_rc=None):
    
    if rgb is None:
        return

    # Start from a reliable fit-to-data state, then apply an optional zoom factor.
    viewer.reset_view()
    if center_rc is None:
        h, w = int(rgb.shape[1]), int(rgb.shape[2])
        target_x = min(max(float(IMAGE_ZOOM_TARGET_X), 0.0), 1.0)
        target_y = min(max(float(IMAGE_ZOOM_TARGET_Y), 0.0), 1.0)
        center_row = target_y * (h - 1)
        center_col = target_x * (w - 1)
    else:
        center_row, center_col = center_rc
    viewer.camera.center = (center_row, center_col)
    viewer.camera.zoom = float(viewer.camera.zoom) * float(zoom_factor)


def _geom_to_pixel_parts(geom):
    if geom is None or geom.is_empty or inv_transform is None:
        return []
    polys = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    out = []
    for poly in polys:
        coords = np.array(poly.exterior.coords)
        px, py = inv_transform * (coords[:, 0], coords[:, 1])
        out.append(np.column_stack([py, px]))
    return out

def _current_time_str():
    if not times:
        return None
    t_idx = int(viewer.dims.current_step[0])
    if 0 <= t_idx < len(times):
        return str(times[t_idx])
    return None

def _update_pending_label():
    pending_count_label.setText(f"Pending updates: {len(pending_tree_updates)}")

def _selected_quality_value():
    btn = quality_value.checkedButton()
    return btn.text() if btn else "Not Checked"

def _selected_flowering_value():
    btn = flowering_group.checkedButton()
    return btn.text() if btn else "No"

def _selected_leafing_value():
    return int(leafing_slider.value())

def _selected_ignore_value():
    btn = ignore_group.checkedButton()
    if btn is None:
        return False
    return btn.text().strip().lower() == "yes"

def _selected_notes_value():
    return (notes_input.text() or "").strip()

def _normalize_edited_flag(value):
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)

def stage_current_polygon_edit():
    global pending_polygon_edits, pending_tree_updates
    if shapes_layer is None or att_transform is None:
        print("⚠️ No editable polygon layer loaded.")
        return

    time_str = _current_time_str()
    tree_id = (current_selected_tree_id or "").strip()
    if not tile_id or not time_str or not tree_id:
        print("⚠️ Select a tree and time step before staging edits.")
        return
    if len(shapes_layer.data) == 0:
        print("⚠️ No polygon geometry to stage.")
        return

    edited_polys = []
    for coords_pixel in shapes_layer.data:
        coords_world = np.array([att_transform * (col, row) for row, col in coords_pixel])
        if coords_world.shape[0] < 3:
            continue
        poly = Polygon(coords_world)
        if not poly.is_empty and poly.is_valid and poly.area > 0:
            edited_polys.append(poly)

    if not edited_polys:
        print("⚠️ Edited shape is invalid; nothing staged.")
        return

    key = (tile_id, time_str, tree_id)
    edited_geom = edited_polys[0] if len(edited_polys) == 1 else MultiPolygon(edited_polys)

    original_geom = original_geometry_by_key.get(key)
    if original_geom is not None and edited_geom.equals_exact(original_geom, tolerance=1e-6):
        is_edited = False
        pending_polygon_edits.pop(key, None)
    else:
        is_edited = True
        pending_polygon_edits[key] = edited_geom

    pending_tree_updates[key] = {
        "quality": _selected_quality_value(),
        "leafing": _selected_leafing_value(),
        "flowering": _selected_flowering_value(),
        "ignore": _selected_ignore_value(),
        "notes": _selected_notes_value(),
        "edited": is_edited,
    }

    # Keep camera behavior consistent with new edited geometry.
    centroid = edited_geom.centroid
    cx, cy = inv_transform * (centroid.x, centroid.y)
    tree_centroid_by_id[tree_id] = (float(cy), float(cx))
    tree_crown_area_by_id[tree_id] = float(edited_geom.area)

    _update_pending_label()
    print(
        f"Staged update for tree={tree_id}, time={time_str}, "
        f"quality={pending_tree_updates[key]['quality']}, edited={is_edited}, "
        f"ignore={pending_tree_updates[key]['ignore']}. "
        f"Pending updates: {len(pending_tree_updates)}"
    )

def clear_pending_polygon_edits():
    pending_polygon_edits.clear()
    pending_tree_updates.clear()
    _update_pending_label()
    print("Cleared all pending polygon edits.")

def _safe_time(time_str):
    """Convert time string to safe filename component."""
    return time_str.replace(":", "-").replace("/", "-")

def _write_part(gdf, bucket_id, time_str):
    """Write a GeoDataFrame partition to disk."""
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    os.makedirs(bucket_dir, exist_ok=True)
    part_path = os.path.join(bucket_dir, f"{_safe_time(time_str)}.parquet")
    gdf.to_parquet(part_path, index=False)

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
            missing_crown_area = part_gdf["crown_area"].isna()
            if missing_crown_area.any():
                part_gdf.loc[missing_crown_area, "crown_area"] = computed_crown_area.loc[missing_crown_area]

    return gpd.GeoDataFrame(part_gdf, crs="EPSG:32617")

def apply_pending_polygon_edits():
    if not pending_tree_updates and not pending_polygon_edits:
        print("No pending edits to apply.")
        return

    print(f"Applying {len(pending_tree_updates)} pending updates to partitioned folder {master_parts_dir}...")
    
    # Group edits by (bucket_id, time) to batch read/write per partition
    edits_by_partition = {}
    all_keys = set(pending_tree_updates.keys()) | set(pending_polygon_edits.keys())
    for (bucket_id_key, time_key, tree_id_key) in all_keys:
        partition_key = (bucket_id_key, time_key)
        if partition_key not in edits_by_partition:
            edits_by_partition[partition_key] = []
        edits_by_partition[partition_key].append((tree_id_key, bucket_id_key, time_key))

    updated_rows = 0
    unmatched = []
    updated_partitions = {}

    for (bucket_id_key, time_key), tree_edits in edits_by_partition.items():
        # Read only the target parquet for this bucket/time.
        part_gdf = read_part_by_time(bucket_id_key, time_key)
        if part_gdf.empty:
            print(f"⚠️ Partition not found: bucket={bucket_id_key}, time={time_key}")
            for tree_id_key, _, _ in tree_edits:
                unmatched.append((bucket_id_key, time_key, tree_id_key))
            continue

        # Ensure metadata columns exist
        for col in ["edited", "quality", "leafing", "flowering", "ignore", "notes"]:
            if col not in part_gdf.columns:
                part_gdf[col] = {"edited": False, "quality": "Not Checked", "leafing": 0, "flowering": "No", "ignore": False, "notes": ""}.get(col, "")

        # Apply all edits for this partition
        for tree_id_key, _, _ in tree_edits:
            mask = part_gdf["GlobalID"].astype(str).str.strip() == str(tree_id_key).strip()
            n_match = int(mask.sum())
            if n_match == 0:
                unmatched.append((bucket_id_key, time_key, tree_id_key))
                continue

            meta = pending_tree_updates.get((bucket_id_key, time_key, tree_id_key), {})
            part_gdf.loc[mask, "edited"] = bool(meta.get("edited", False))
            part_gdf.loc[mask, "quality"] = str(meta.get("quality", "Not Checked"))
            part_gdf.loc[mask, "leafing"] = int(meta.get("leafing", 0))
            part_gdf.loc[mask, "flowering"] = str(meta.get("flowering", "No"))
            part_gdf.loc[mask, "ignore"] = bool(meta.get("ignore", False))
            part_gdf.loc[mask, "notes"] = str(meta.get("notes", ""))

            edited_geom = pending_polygon_edits.get((bucket_id_key, time_key, tree_id_key))
            if edited_geom is not None:
                part_gdf.loc[mask, "geometry"] = edited_geom
                if "crown_area" in part_gdf.columns:
                    part_gdf.loc[mask, "crown_area"] = float(edited_geom.area)
            updated_rows += n_match

        # Cache the updated partition for writing
        updated_partitions[(bucket_id_key, time_key)] = part_gdf

    # Write all updated partitions back to disk
    for (bucket_id_key, time_key), updated_part_gdf in updated_partitions.items():
        _write_part(updated_part_gdf, bucket_id_key, time_key)
        print(f"  Saved partition: bucket={bucket_id_key}, time={time_key} ({len(updated_part_gdf)} rows)")

    print(f"Applied edits to {updated_rows} row(s).")
    if unmatched:
        print(f"⚠️ Unmatched staged edits: {len(unmatched)}")
    
    # Update in-memory dictionaries so viewer reflects changes immediately
    for (bucket_id_key, time_key, tree_id_key), meta in list(pending_tree_updates.items()):
        time_str_key = str(time_key)
        gid_key = str(tree_id_key).strip()
        loaded_edited_by_key[(time_str_key, gid_key)] = bool(meta.get("edited", False))
        if (time_str_key, gid_key) in loaded_labels_by_key:
            loaded_labels_by_key[(time_str_key, gid_key)].update({
                "leafing": meta.get("leafing", "—"),
                "flowering": meta.get("flowering", "—"),
                "quality": meta.get("quality", "—"),
                "ignore": meta.get("ignore", "—"),
            })
    
    pending_polygon_edits.clear()
    pending_tree_updates.clear()
    _update_pending_label()

def update_polygons(event=None):
    if rgb is None or len(times) == 0 or shapes_layer is None or not pixel_shapes_by_time:
        return

    t_idx = int(viewer.dims.current_step[0])
    if t_idx == last_t_idx["value"]:
        return
    last_t_idx["value"] = t_idx

    if t_idx >= len(times):
        shapes_layer.data = []
        return

    current_time_str = str(times[t_idx])
    shapes = pixel_shapes_by_time.get(current_time_str, []) #here it is 
    shapes_layer.data = shapes
    if shapes:
        shapes_layer.edge_width = 5
        shapes_layer.edge_color = "red"
        shapes_layer.opacity = 1.0

def update_title(event=None):
    if len(times) == 0:
        viewer.title = f"Tile: {tile_id} | No cube loaded"
        return

    t_idx = int(viewer.dims.current_step[0])
    if 0 <= t_idx < len(times):
        viewer.title = f"Tile: {tile_id} | Time: {times[t_idx]}"

def update_info_box(event=None):
    tree_id = (current_selected_tree_id or "").strip()
    time_str = _current_time_str()

    if not tree_id:
        info_box.setText("<i>No tree selected</i>")
        return

    # Species and tag from species_df
    tree_row = species_df.loc[species_df["globalid_str"] == tree_id]
    species_name = tree_row["latin_str"].iloc[0] if not tree_row.empty else ""
    tag_col = next((c for c in species_df.columns if c.lower() in {"tag", "tree_tag", "mnemonic", "sp_code"}), None)
    tag_val = str(tree_row[tag_col].iloc[0]) if (tag_col and not tree_row.empty) else tree_id

    # Area
    area = tree_crown_area_by_id.get(tree_id)
    area_str = f"{area:.2f} m²" if area is not None else "—"

    # Labels: staged takes priority over loaded
    labels = {}
    if time_str:
        loaded = loaded_labels_by_key.get((time_str, tree_id), {})
        staged = pending_tree_updates.get((tile_id, time_str, tree_id), {})
        labels = {**loaded, **staged}

    leafing  = labels.get("leafing",   "—")
    flowering = labels.get("flowering", "—")
    quality  = labels.get("quality",   "—")
    ignore   = labels.get("ignore",    "—")
    edited   = labels.get("edited",    "—")
    notes    = labels.get("notes",     "")

    date_str = (time_str or "").split("T")[0] if time_str else "—"

    lines = [
        f"<b>Date:</b> {date_str}",
        f"<b>Species:</b> {species_name}",
        f"<b>Tag:</b> {tag_val}",
        f"<b>Area:</b> {area_str}",
        "<hr/>",
        f"<b>Leafing:</b> {leafing}",
        f"<b>Flowering:</b> {flowering}",
        f"<b>Quality:</b> {quality}",
        f"<b>Edited:</b> {edited}",
        f"<b>Ignore:</b> {ignore}",
    ]
    if notes:
        lines.append(f"<b>Notes:</b> {notes}")
    info_box.setText("<br/>".join(lines))

def load_tree_data_for_tile(tile_id_value, selected_tree_id=None):
    global pixel_shapes_by_time, pixel_shape_records_by_time, tree_centroid_by_id
    global tree_crown_area_by_id, original_geometry_by_key, loaded_edited_by_key, loaded_labels_by_key

    selected_tree_str = (selected_tree_id or "").strip()
    if selected_tree_str:
        tile_gdf_target = read_part_filtered(tile_id_value, global_id=selected_tree_str)
        print(f"Loaded {len(tile_gdf_target)} rows for bucket={tile_id_value}, tree={selected_tree_str}")
    else:
        print("No tree selected, loading entire bucket...")
        tile_gdf_target = read_part_filtered(tile_id_value)

    if tile_gdf_target.empty:
        pixel_shapes_by_time = {}
        pixel_shape_records_by_time = {}
        tree_centroid_by_id = {}
        tree_crown_area_by_id = {}
        original_geometry_by_key = {}
        loaded_edited_by_key = {}
        loaded_labels_by_key = {}
        return

    # Only simplify the subset we are actually going to use.
    tile_gdf_target.geometry = tile_gdf_target.geometry.simplify(0.1)

    _h, _w = rgb.shape[1], rgb.shape[2]
    pixel_shapes_by_time = {}
    pixel_shape_records_by_time = {}
    tree_centroid_by_id = {}
    tree_crown_area_by_id = {}
    original_geometry_by_key = {}
    loaded_edited_by_key = {}
    loaded_labels_by_key = {}

    # Build per-tree centroid and crown area lookup for camera centering and dynamic zoom.
    for _gid, _group in tile_gdf_target.groupby("GlobalID"):
        if not _gid:
            continue
        _row = None
        _area = None
        for _, _record in _group.iterrows():
            _geom = _record.geometry
            if _geom is None or _geom.is_empty:
                continue
            if _row is None:
                _row = _geom.centroid
            if _area is None and "crown_area" in _record:
                _area = _record["crown_area"]
        if _row is None:
            continue
        _cx, _cy = inv_transform * (_row.x, _row.y)
        tree_centroid_by_id[_gid] = (float(_cy), float(_cx))
        if _area is not None:
            tree_crown_area_by_id[_gid] = float(_area)

    for _time_str, _group in tile_gdf_target.groupby("time"):
        _time_str_key = str(_time_str)
        _shapes = []
        _records = []
        for _, _record in _group.iterrows():
            _geom = _record.geometry
            if _geom is None or _geom.is_empty:
                continue
            _gid = str(_record["GlobalID"]).strip()
            original_geometry_by_key[(tile_id_value, _time_str_key, _gid)] = _geom
            loaded_edited_by_key[(_time_str_key, _gid)] = _normalize_edited_flag(_record.get("edited", False))
            loaded_labels_by_key[(_time_str_key, _gid)] = {
                "leafing":   _record.get("leafing", "—"),
                "flowering": _record.get("flowering", "—"),
                "quality":   _record.get("quality", "—"),
                "ignore":    _record.get("ignore", "—"),
                "edited":    _normalize_edited_flag(_record.get("edited", False)),
                "notes":     _record.get("notes", ""),
            }
            _polys = _geom.geoms if isinstance(_geom, MultiPolygon) else [_geom]
            for _poly in _polys:
                _coords = np.array(_poly.exterior.coords)
                _px, _py = inv_transform * (_coords[:, 0], _coords[:, 1])
                _cp = np.column_stack([_py, _px])
                if (np.all(_cp[:, 0] < 0) or np.all(_cp[:, 0] > _h) or
                        np.all(_cp[:, 1] < 0) or np.all(_cp[:, 1] > _w)):
                    continue
                _shapes.append(_cp)
                _records.append((_gid, _cp))
        pixel_shapes_by_time[_time_str_key] = _shapes
        pixel_shape_records_by_time[_time_str_key] = _records

def load_tile_cube(new_tile_id, selected_tree_id=None):
    global tile_id, path, z, rgb, times, att_transform, inv_transform, image_layer, shapes_layer, pixel_shapes_by_time, pixel_shape_records_by_time, tree_centroid_by_id, tree_crown_area_by_id, original_geometry_by_key, loaded_edited_by_key, loaded_labels_by_key

    new_path = os.path.join(cube_root, new_tile_id, "cube.zarr")
    if not os.path.exists(new_path):
        print(f"⚠️ Cube path not found for tile {new_tile_id}: {new_path}")
        return False

    print(f"Switching to tile: {new_tile_id}")
    path = new_path
    tile_id = new_tile_id

    z = zarr.open(path, mode="r")
    rgb = da.from_zarr(z)[:, 0:3, :, :]
    rgb = da.moveaxis(rgb, 1, -1)
    times = list(z.attrs.get("time", []))

    sample = rgb[0].compute()
    p2_new, p98_new = np.percentile(sample, [2, 98])

    if image_layer is None:
        image_layer = viewer.add_image(
            rgb,
            rgb=True,
            name="cube",
            contrast_limits=(float(p2_new), float(p98_new)),
            gamma=1.0
        )
    else:
        image_layer.data = rgb
        image_layer.contrast_limits = (float(p2_new), float(p98_new))

    # Create shapes layer on first cube load
    if shapes_layer is None:
        shapes_layer = viewer.add_shapes(
            [],
            shape_type="polygon",
            edge_color="red",
            face_color="transparent",
            edge_width=9,
            name="polygons",
            opacity=1
        )

    # Reconfigure time axis after replacing layer data so date swiping is available.
    if len(times) > 0:
        viewer.dims.set_range(0, (0, max(len(times) - 1, 0), 1))

    att_transform = Affine(*z.attrs["transform"])
    inv_transform = ~att_transform

    # Precompute pixel coordinates for all polygons at all time steps (paid once at tile load)
    print(f"Precomputing polygon pixel coords for tile {new_tile_id}...")
    start_time = time()
    load_tree_data_for_tile(new_tile_id, selected_tree_id=selected_tree_id)
    end_time = time()
    print(f"Precomputed {len(pixel_shapes_by_time)} time steps in {end_time - start_time:.2f} seconds.")

    last_t_idx["value"] = None
    time_step = 49

    viewer.dims.set_point(0, time_step)
    selected_tree_str = (selected_tree_id or "").strip()
    target_center = tree_centroid_by_id.get(selected_tree_str)
    target_area = tree_crown_area_by_id.get(selected_tree_str)
    dynamic_zoom = compute_dynamic_zoom(target_area)
    set_camera_view(zoom_factor=dynamic_zoom, center_rc=target_center)
    update_title()
    update_polygons()
    update_info_box()

    return True

control_widget = QWidget()
layout = QVBoxLayout()

info_box = QLabel("<i>No tree selected</i>")
info_box.setWordWrap(True)
info_box.setStyleSheet(
    "background:#1e1e1e; color:#d4d4d4; padding:6px; border-radius:4px; font-size:11px;"
)
layout.addWidget(info_box)
spp_selection = QLabel("Select species:")
layout.addWidget(spp_selection)
tree_selection = QLabel("Select tree:")
layout.addWidget(tree_selection)


spp_combo = QComboBox()
spp_combo.addItem("")
spp_combo.addItems(species)
layout.addWidget(spp_combo)

tree_combo = QComboBox()
tree_combo.addItem("")
layout.addWidget(tree_combo)

leafing_label = QLabel("Leafing: 0")
layout.addWidget(leafing_label)
leafing_slider = QSlider(Qt.Horizontal)
leafing_slider.setMinimum(0)
leafing_slider.setMaximum(100)
leafing_slider.setValue(0)
leafing_slider.setTickPosition(QSlider.TicksBelow)
leafing_slider.setTickInterval(10)
leafing_slider.valueChanged.connect(lambda v: leafing_label.setText(f"Leafing: {v}"))
layout.addWidget(leafing_slider)

# Flowering buttons (Yes, No, Maybe)
flowering_label = QLabel("Flowering:")
layout.addWidget(flowering_label)
flowering_layout = QHBoxLayout()
flowering_group = QButtonGroup()
flowering_no = QRadioButton("No")
flowering_yes = QRadioButton("Yes")
flowering_maybe = QRadioButton("Maybe")
flowering_no.setChecked(True)
flowering_group.addButton(flowering_no)
flowering_group.addButton(flowering_yes)
flowering_group.addButton(flowering_maybe)
flowering_layout.addWidget(flowering_no)
flowering_layout.addWidget(flowering_yes)
flowering_layout.addWidget(flowering_maybe)
layout.addLayout(flowering_layout)

ignore_label = QLabel("Ignore:")
layout.addWidget(ignore_label)
ignore_layout = QHBoxLayout()
ignore_group = QButtonGroup()
ignore_no = QRadioButton("No")
ignore_yes = QRadioButton("Yes")
ignore_no.setChecked(True)
ignore_group.addButton(ignore_no)
ignore_group.addButton(ignore_yes)
ignore_layout.addWidget(ignore_no)
ignore_layout.addWidget(ignore_yes)
layout.addLayout(ignore_layout)

notes_label = QLabel("Notes:")
layout.addWidget(notes_label)
notes_input = QLineEdit()
notes_input.setPlaceholderText("Why ignore/decommission this tree?")
layout.addWidget(notes_input)

quaility_label = QLabel("Quality:")
layout.addWidget(quaility_label)
quality_layout = QHBoxLayout()
quality_value = QButtonGroup()
quality_excellent = QRadioButton("Excellent")
quality_good = QRadioButton("Good")
alright_quality = QRadioButton("Alright")
quality_bad = QRadioButton("Bad")
very_bad_quality = QRadioButton("Very Bad")
quality_not_checked = QRadioButton("Not Checked")
quality_value.addButton(quality_excellent)
quality_value.addButton(quality_good)
quality_value.addButton(alright_quality)
quality_value.addButton(quality_bad)
quality_value.addButton(very_bad_quality)
quality_value.addButton(quality_not_checked)
quality_layout.addWidget(quality_excellent)
quality_layout.addWidget(quality_good)
quality_layout.addWidget(alright_quality)
quality_layout.addWidget(quality_bad)
quality_layout.addWidget(very_bad_quality)
quality_layout.addWidget(quality_not_checked)
layout.addLayout(quality_layout)

pending_count_label = QLabel("Pending geometry edits: 0")
layout.addWidget(pending_count_label)

stage_edit_button = QPushButton("Stage Polygon Edit")
stage_edit_button.clicked.connect(stage_current_polygon_edit)
layout.addWidget(stage_edit_button)

apply_edits_button = QPushButton("Apply All Pending Edits")
apply_edits_button.clicked.connect(apply_pending_polygon_edits)
layout.addWidget(apply_edits_button)

clear_edits_button = QPushButton("Clear Pending Edits")
clear_edits_button.clicked.connect(clear_pending_polygon_edits)
layout.addWidget(clear_edits_button)



def update_tree_combo(selected_species):
    selected_species = (selected_species or "").strip()
    tree_combo.clear()
    tree_combo.addItem("")

    if not selected_species or "GlobalID" not in species_df.columns:
        return

    matching_trees = sorted(
        {
            gid
            for gid in species_df.loc[species_df["latin_str"] == selected_species, "globalid_str"].unique()
            if gid
        }
    )

    tree_combo.addItems(matching_trees)

def update_cube_from_selection(event=None):
    global current_selected_tree_id
    selected_species = (spp_combo.currentText() or "").strip()
    selected_tree = (tree_combo.currentText() or "").strip()

    if not selected_species or not selected_tree:
        return

    row = species_df.loc[species_df["globalid_str"] == selected_tree, "bucket_id_str"]
    if row.empty:
        print(f"⚠️ No bucket found for tree={selected_tree}")
        return

    previous_tree_id = (current_selected_tree_id or "").strip()
    current_selected_tree_id = selected_tree

    target_tile = row.iloc[0]
    if image_layer is None or target_tile != tile_id:
        load_tile_cube(target_tile, selected_tree_id=selected_tree)
    else:
        if selected_tree != previous_tree_id:
            load_tree_data_for_tile(target_tile, selected_tree_id=selected_tree)
        last_t_idx["value"] = None
        viewer.dims.set_point(0, 49)
        target_center = tree_centroid_by_id.get(selected_tree)
        target_area = tree_crown_area_by_id.get(selected_tree)
        dynamic_zoom = compute_dynamic_zoom(target_area)
        set_camera_view(zoom_factor=dynamic_zoom, center_rc=target_center)
        update_polygons()
        update_polygons_for_selected_tree()
        update_title()
        update_info_box()

def update_polygons_for_selected_tree(event=None):
    if rgb is None or len(times) == 0 or shapes_layer is None:
        return

    time_str = _current_time_str()
    if time_str is None:
        shapes_layer.data = []
        return

    records = pixel_shape_records_by_time.get(time_str, [])
    tree_id = (current_selected_tree_id or "").strip()

    if not tree_id:
        shapes_layer.data = [coords for _, coords in records]
        shapes_layer.edge_color = "red"
        return

    pending_key = (tile_id, time_str, tree_id)
    print(f"[DEBUG] pending_key={pending_key!r}")
    print(f"[DEBUG] pending_polygon_edits keys={list(pending_polygon_edits.keys())}")
    if pending_key in pending_polygon_edits:
        shapes_layer.data = _geom_to_pixel_parts(pending_polygon_edits[pending_key])
        pending_meta = pending_tree_updates.get(pending_key, {})
        is_edited = bool(pending_meta.get("edited", True))
        shapes_layer.edge_color = "blue" if is_edited else "red"
        return

    shapes_layer.data = [coords for gid, coords in records if gid == tree_id]
    if pending_key in pending_tree_updates:
        is_edited = bool(pending_tree_updates[pending_key].get("edited", False))
    else:
        is_edited = loaded_edited_by_key.get((time_str, tree_id), False)
    shapes_layer.edge_color = "blue" if is_edited else "red"
    update_info_box()

spp_combo.currentTextChanged.connect(update_tree_combo)
tree_combo.currentTextChanged.connect(update_cube_from_selection)

control_widget.setLayout(layout)
viewer.window.add_dock_widget(control_widget, area="right", name="Species")
# ---------------------------------------------------
# CONNECT
# ---------------------------------------------------
viewer.dims.events.current_step.connect(update_polygons)
viewer.dims.events.current_step.connect(update_title)
viewer.dims.events.current_step.connect(update_polygons_for_selected_tree)
viewer.dims.events.current_step.connect(update_info_box)

update_polygons()
update_polygons_for_selected_tree()
update_title()

# ---------------------------------------------------
# RUN
# ---------------------------------------------------
napari.run()