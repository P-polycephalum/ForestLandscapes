import napari
import rasterio
import json
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, mapping
import os
import geopandas as gpd
from qtpy.QtWidgets import (QPushButton, QVBoxLayout, QWidget, QInputDialog, 
                            QSlider, QLabel, QRadioButton, QButtonGroup, QHBoxLayout, QComboBox)
from qtpy.QtCore import QTimer, Qt

plot="50ha"
spp="Cavanillesia platanifolia"
master_control_json=f"D:\\BCI_{plot}_timeseries\\master_control.json"
master_control_gdf= f"D:\\BCI_{plot}_timeseries\\master_control_polygons.gpkg"
print(f"Loading master control JSON from {master_control_json}...")
print(f"Loading master control GeoDataFrame from {master_control_gdf}...")

# Track current index
current_index = [0]  # Use list to make it mutable
current_transform = [None]  # Store transform for coordinate conversion
current_tree_id = [None]  # Store current tree ID (polygon_id)
saved_polygons = []  # Store all saved polygon edits
labeled_polygon_ids = set()  # Track which polygons are already labeled

# Load existing JSON if it exists
if os.path.exists(master_control_json):
    try:
        with open(master_control_json, 'r') as f:
            saved_polygons = json.load(f)
        
        if spp is not None:
            saved_polygons['crowns'] = {k: v for k, v in saved_polygons['crowns'].items() if v.get('latin') == spp}
        
        individuals = len(saved_polygons['crowns'].keys())
        unique_dates = tuple(sorted(set(crown['date'] for crown in saved_polygons['crowns'].values() if 'date' in crown)))
        # Populate labeled_polygon_ids with already-edited polygons
        for polygon_id, crown_data in saved_polygons['crowns'].items():
            if crown_data.get('edited', False):
                labeled_polygon_ids.add(polygon_id)
        print(f"Number of Individuals: {individuals}")
        print(f"Unique dates: {len(unique_dates)}")
        print(f"Already edited: {len(labeled_polygon_ids)}")
    except Exception as e:
        print(f"Could not load existing JSON: {e}")

geoms = gpd.read_file(master_control_gdf)
# filter geoms with polygon id in saved_polygons
if spp is not None:
    geoms = geoms[geoms['polygon_id'].isin(saved_polygons['crowns'].keys())]

viewer = napari.Viewer()

# Reference to image and shapes layers
image_layer = None
shapes_layer = None
reference_image_layer = None  # Reference orthomosaic layer
reference_shapes_layer = None  # Reference crown polygon layer
reference_transform = [None]  # Store reference transform for reuse
tree_rows = None  # Will be set when global_id is selected

# Function to load reference data (called once per global_id)
def load_reference(global_id):
    global reference_image_layer, reference_shapes_layer, reference_transform
    
    try:
        # Load reference orthomosaic
        reference_id = f"{global_id}_2022_09_29" # this is the reference
        ref_img_path = os.path.join(saved_polygons['paths']['crowns_dir'], global_id, f"{reference_id}.tif")
        if os.path.exists(ref_img_path):
            if reference_image_layer is not None:
                viewer.layers.remove(reference_image_layer)
            if reference_shapes_layer is not None:
                viewer.layers.remove(reference_shapes_layer)
            
            with rasterio.open(ref_img_path) as src:
                ref_img = src.read()
                reference_transform[0] = src.transform
            
            reference_image_layer = viewer.add_image(ref_img.transpose(1, 2, 0), name=f'Reference {global_id}', rgb=True, opacity=0.7, visible=False)
            
            # Load reference crown polygon
            try:
                ref_polygon_id = f"{global_id}_2022_09_29"
                ref_geom = geoms.loc[geoms['polygon_id'] == ref_polygon_id, 'geometry']
                if not ref_geom.empty:
                    ref_geom = ref_geom.values[0]
                    ref_geom = ref_geom.simplify(0.10, preserve_topology=True)
                    ref_shapes_to_add = []
                    
                    if isinstance(ref_geom, Polygon):
                        coords_world = np.array(ref_geom.exterior.coords)
                        coords_pixel = np.array([~reference_transform[0] * (x, y) for x, y in coords_world])
                        coords_pixel = coords_pixel[:, ::-1]
                        ref_shapes_to_add.append(coords_pixel)
                    elif isinstance(ref_geom, MultiPolygon):
                        for p in ref_geom.geoms:
                            coords_world = np.array(p.exterior.coords)
                            coords_pixel = np.array([~reference_transform[0] * (x, y) for x, y in coords_world])
                            coords_pixel = coords_pixel[:, ::-1]
                            ref_shapes_to_add.append(coords_pixel)
                    
                    if ref_shapes_to_add:
                        reference_shapes_layer = viewer.add_shapes(
                            ref_shapes_to_add,
                            shape_type='polygon',
                            edge_width=2,
                            edge_color='cyan',
                            face_color='transparent',
                            name='Reference Crown',
                            visible=False
                        )
            except Exception as e:
                print(f"Could not load reference crown polygon for {ref_polygon_id}: {e}")
        else:
            print(f"Reference image not found: {ref_img_path}")
    except Exception as e:
        print(f"Could not load reference data for {global_id}: {e}")

# Function to load a tree by index
def load_tree(polygon_id):
    global image_layer, shapes_layer
    geom = geoms.loc[geoms['polygon_id'] == polygon_id, 'geometry'].values[0]
    global_id = polygon_id.split("_")[0]
    try:
        img_path = os.path.join(saved_polygons['paths']['crowns_dir'],global_id, f"{polygon_id}.tif")
        with rasterio.open(img_path) as src:
            img = src.read()
            metadata = src.meta
            transform = src.transform
        
        # Store transform and tree_id globally
        current_transform[0] = transform
        current_tree_id[0] = polygon_id
        
        # Update or create image layer
        if image_layer is not None:
            viewer.layers.remove(image_layer)
            
        image_layer = viewer.add_image(img.transpose(1, 2, 0), name=f'Tree {polygon_id}', rgb=True)

        # Update or create shapes layer for current crown
        if shapes_layer is not None:
            viewer.layers.remove(shapes_layer)
        
        if geom is not None:
            
            #lets slightly simplify the polygon by 0.20
            geom = geom.simplify(0.10, preserve_topology=True)

            shapes_to_add = []
            
            if isinstance(geom, Polygon):
                coords_world = np.array(geom.exterior.coords)
                coords_pixel = np.array([~transform * (x, y) for x, y in coords_world])
                coords_pixel = coords_pixel[:, ::-1]
                shapes_to_add.append(coords_pixel)
            elif isinstance(geom, MultiPolygon):
                for p in geom.geoms:
                    coords_world = np.array(p.exterior.coords)
                    coords_pixel = np.array([~transform * (x, y) for x, y in coords_world])
                    coords_pixel = coords_pixel[:, ::-1]
                    shapes_to_add.append(coords_pixel)
            
            if shapes_to_add:
                # Check if already reviewed (either edited geometry or marked as excellent) - use green as visual indicator
                is_edited = saved_polygons['crowns'].get(str(polygon_id), {}).get('edited', False)
                is_excellent = saved_polygons['crowns'].get(str(polygon_id), {}).get('no_edits_needed', False)
                edge_color = 'green' if (is_edited or is_excellent) else 'yellow'
                
                shapes_layer = viewer.add_shapes(
                    shapes_to_add,
                    shape_type='polygon',
                    edge_width=3,
                    edge_color=edge_color,
                    face_color='transparent',
                    name='Crown Polygon'
                )
                # Set to select mode to enable editing
                shapes_layer.mode = 'direct'
        
        # Load reference crown polygon for same polygon_id (hidden by default)

    except Exception as e:
        print(f"Error loading tree {polygon_id}: {e}")

# Function to go to next tree (skip already labeled)
def next_tree():
    current_index[0] += 1
    # Skip already labeled polygons for this tree
    while current_index[0] < len(tree_rows):
        polygon_id = tree_rows.iloc[current_index[0]]['polygon_id']
        if polygon_id not in labeled_polygon_ids:
            load_tree(polygon_id)
            return
        current_index[0] += 1

    print("Reached the end of this tree's time series!")
    current_index[0] = len(tree_rows) - 1

def back_tree():
    current_index[0] -= 1
    if current_index[0] < 0:
        print("Already at the beginning of this tree's time series!")
        current_index[0] = 0
    else:
        load_tree(tree_rows.iloc[current_index[0]]['polygon_id'])

def save_polygon():
    if current_transform[0] is None or current_tree_id[0] is None:
        print("No tree loaded!")
        return
    
    transform = current_transform[0]
    tree_id = current_tree_id[0]
    
    # Get label values
    leafing_value = leafing_slider.value()
    flowering_value = flowering_group.checkedButton().text() if flowering_group.checkedButton() else "No"
    quality_choice = quality_value.checkedButton().text() if quality_value.checkedButton() else "not checked"
    
    # Check if polygon has been edited geometrically
    has_geometry_edits = shapes_layer is not None and len(shapes_layer.data) > 0
    
    # Update geometry if shapes were edited
    if has_geometry_edits:
        for i, coords_pixel in enumerate(shapes_layer.data):
            coords_world = np.array([transform * (col, row) for row, col in coords_pixel])
            poly = Polygon(coords_world)
            geoms.loc[geoms['polygon_id'] == tree_id, 'geometry'] = poly
    
    # Store the result with labels in JSON
    saved_polygons['crowns'][str(tree_id)]['labels'] = {
        "leafing": leafing_value,
        "flowering": flowering_value
    }
    saved_polygons['crowns'][str(tree_id)]['quality'] = quality_choice  # Store quality assessment
    
    # Mark as edited only if geometry was actually changed
    if has_geometry_edits:
        saved_polygons['crowns'][str(tree_id)]['edited'] = True
    
    # If marked as excellent, mark as no edits needed
    if quality_choice == "Excellent":
        saved_polygons['crowns'][str(tree_id)]['no_edits_needed'] = True
    if quality_choice == "Good":
        saved_polygons['crowns'][str(tree_id)]['no_edits_needed'] = True
    if quality_choice == "Alright":
        saved_polygons['crowns'][str(tree_id)]['no_edits_needed'] = True
    
    # Mark this polygon as labeled
    labeled_polygon_ids.add(tree_id)
    
    # Build status message
    status_msg = "Saved"
    if has_geometry_edits:
        status_msg += f" {len(shapes_layer.data)} polygon(s) with edits"
    else:
        status_msg += " labels only (no geometry changes)"
    if quality_choice == "Excellent":
        status_msg += " (Marked as Excellent - no edits needed)"
    
    print(f"{status_msg} for tree {tree_id}")

# Function to export all saved polygons to JSON and GeoDataFrame
def export_json():
    if len(saved_polygons) == 0:
        print("No polygons saved yet!")
        return
    
    # Load the full original JSON to preserve all crowns
    with open(master_control_json, 'r') as f:
        full_data = json.load(f)
    
    # Merge edited crowns back into the full dataset
    for crown_id, crown_data in saved_polygons['crowns'].items():
        full_data['crowns'][crown_id] = crown_data
    
    # Save the complete merged dataset
    with open(master_control_json, 'w') as f:
        json.dump(full_data, f, indent=2)
    print(f"Exported {len(saved_polygons['crowns'])} edited crowns (merged into full dataset with {len(full_data['crowns'])} total crowns)")
    
    # Load full original geoms and merge modified geometries
    full_geoms = gpd.read_file(master_control_gdf)
    modified_polygon_ids = geoms['polygon_id'].unique()
    
    # Update only the modified geometries in the full dataset
    for polygon_id in modified_polygon_ids:
        if polygon_id in geoms['polygon_id'].values:
            full_geoms.loc[full_geoms['polygon_id'] == polygon_id, 'geometry'] = geoms.loc[geoms['polygon_id'] == polygon_id, 'geometry'].values[0]
    
    # Save complete merged GeoDataFrame
    full_geoms.to_file(master_control_gdf, driver="GPKG", index=False)
    print(f"Exported modified geometries (merged into full dataset with {len(full_geoms)} total polygons)")
    
    viewer.close()

# Create control widget
control_widget = QWidget()
layout = QVBoxLayout()

# Global ID selector at top
global_id_label = QLabel("Select Tree (global_id):")
layout.addWidget(global_id_label)

global_ids = sorted(set(crown['global_id'] for crown in saved_polygons['crowns'].values() if 'global_id' in crown))
global_id_combo = QComboBox()
global_id_combo.addItems([str(gid) for gid in global_ids])
layout.addWidget(global_id_combo)

def on_global_id_changed(selected_id):
    global tree_rows, current_index, labeled_polygon_ids
    tree_rows = geoms[geoms['polygon_id'].astype(str).str.split("_").str[0] == selected_id]
    if tree_rows.empty:
        print(f"No rows found for global_id {selected_id}.")
        return
    # Load reference data once for this global_id
    load_reference(selected_id)
    # Filter labeled_polygon_ids to only include polygons for this global_id
    labeled_polygon_ids.intersection_update(set(tree_rows['polygon_id'].astype(str)))
    current_index[0] = min(48, len(tree_rows) - 1)  # Start at index 47 or latest if fewer rows
    load_tree(tree_rows.iloc[current_index[0]]['polygon_id'])
    print(f"\nLoaded tree {selected_id} with {len(tree_rows)} dates")

global_id_combo.currentTextChanged.connect(on_global_id_changed)

# Initialize with first global_id
first_global_id = global_id_combo.currentText()
on_global_id_changed(first_global_id)

next_btn = QPushButton("Next Tree (or press N)")
next_btn.clicked.connect(next_tree)
layout.addWidget(next_btn)

back_btn = QPushButton("Previous Tree (or press B)")
back_btn.clicked.connect(back_tree)
layout.addWidget(back_btn)

save_btn = QPushButton("Save Current Polygon (or press S)")
save_btn.clicked.connect(save_polygon)
layout.addWidget(save_btn)

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


export_btn = QPushButton("Export All & Close (or press Ctrl+S)")
export_btn.clicked.connect(export_json)
layout.addWidget(export_btn)

control_widget.setLayout(layout)
viewer.window.add_dock_widget(control_widget, area='right', name='Controls')

# Add keyboard shortcuts
@viewer.bind_key('s')
def save_on_key(viewer):
    save_polygon()

@viewer.bind_key('n')
def next_on_key(viewer):
     next_tree()

@viewer.bind_key('b')
def back_on_key(viewer):
     back_tree()

@viewer.bind_key('Control-S')
def export_on_key(viewer):
    export_json()

# Find the first unedited tree starting from index 128
def find_first_unedited(start_index):
    # Search forward from start_index
    for i in range(start_index, len(tree_rows)):
        polygon_id = tree_rows.iloc[i]['polygon_id']
        if polygon_id not in labeled_polygon_ids:
            return i
    # Search backward from start_index
    for i in range(start_index - 1, -1, -1):
        polygon_id = tree_rows.iloc[i]['polygon_id']
        if polygon_id not in labeled_polygon_ids:
            return i
    return -1  # All labeled

# Load the first unedited tree starting from index 47
first_index = min(48, len(tree_rows) - 1)
unedited_index = find_first_unedited(first_index)
if unedited_index >= 0:
    current_index[0] = unedited_index
    load_tree(tree_rows.iloc[current_index[0]]['polygon_id'])
    print(f"Starting at date {current_index[0] + 1}/{len(tree_rows)} (closest unedited from index 47)")
else:
    print("All polygons for this tree already edited!")
    current_index[0] = first_index
    load_tree(tree_rows.iloc[current_index[0]]['polygon_id'])

print("\n" + "="*60)
print("POLYGON EDITOR INSTRUCTIONS:")
print("="*60)
print(f"Already labeled: {len(labeled_polygon_ids)} dates")
print(f"Remaining: {len(tree_rows) - len(labeled_polygon_ids)} dates")
print("="*60)
print("1. Click on the polygon to select it")
print("2. Drag vertices to move them")
print("3. Set labels in the control panel:")
print("   - Leafing slider (0-100)")
print("   - Flowering buttons (Yes/No/Maybe)")
print("4. Press 'S' or click 'Save Current Polygon' to save edits")
print("5. Press 'N' to auto-skip to next unlabeled tree")
print("6. Press 'Ctrl+S' or click 'Export All & Close' when done")
print("7. App loads existing labels and skips already-labeled trees")
print("8. Results saved to D:\\edited_polygons.json")
print("="*60 + "\n")

napari.run()