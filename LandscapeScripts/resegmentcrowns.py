import os
import geopandas as gpd
import json
from LandscapeScripts.utils import tile_ortho, crown_segment, crown_avoid, crownmap_QC, crownmap_metrics

plot="50ha"
spp="Cavanillesia platanifolia"

master_control_json=f"D:\\BCI_{plot}_timeseries\\master_control.json"
master_control_gdf= f"D:\\BCI_{plot}_timeseries\\master_control_polygons.gpkg"

print(f"Loading master control JSON from {master_control_json}...")
with open(master_control_json, 'r') as f:
    master_control = json.load(f)
print(f"Loading master control GeoDataFrame from {master_control_gdf}...")
gdf = gpd.read_file(master_control_gdf)

print("Filtering for species and plot...")
json_filtered = {k: v for k, v in master_control['crowns'].items() if v.get('latin') == spp}
keys_to_keep = set(json_filtered.keys())
gdf_filtered = gdf[gdf['polygon_id'].isin(keys_to_keep)]

print(f"Found {len(gdf_filtered)} polygons for species {spp} in plot {plot}.")

# Extract unique global_ids from the filtered JSON
def calculate_metrics(reference_geom, target_geom, print_results=True):
    if reference_geom is None or target_geom is None:
            print(f"Warning: Missing geometry for reference or target. Skipping metrics calculation.")
            
    else:
        intersection_area = reference_geom.intersection(target_geom).area
        union_area = reference_geom.union(target_geom).area

        precision = intersection_area / target_geom.area if target_geom.area > 0 else 0
        recall = intersection_area / reference_geom.area if reference_geom.area > 0 else 0
        iou = intersection_area / union_area if union_area > 0 else 0
        f1= 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        score = 0.6 * iou + (0.4 * f1)
        if print_results:
            print(f"GlobalID:, IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}, F1: {f1:.2f}, Score: {score:.2f}")
        return iou, precision, recall, f1, score



unique_global_ids = sorted(set(crown.get('global_id') for crown in json_filtered.values() if 'global_id' in crown))
print(f"Unique global IDs in filtered JSON: {len(unique_global_ids)}")
# reference is index 49

for global_id in unique_global_ids:
    #print(f"Processing tree {global_id}...")
    tree_rows = gdf_filtered[gdf_filtered['polygon_id'].str.contains(global_id)]
    print(f"Tree {global_id} has {len(tree_rows)} polygons in the GeoDataFrame.")

    reference_geom = tree_rows.iloc[49]['geometry'] if len(tree_rows) > 49 else None
    for idx, row in tree_rows[50:].iterrows():
        #print(f"Processing polygon {row['polygon_id']}...")
        iou, precision, recall, f1, score = calculate_metrics(reference_geom, row['geometry'])

        if score < 0.5:
            print(f"Polygon {row['polygon_id']} has low similarity score ({score:.2f})")
            # we will transfer the reference geometry
            gdf_filtered.at[idx, 'geometry'] = reference_geom
            tree_rows.at[idx, 'geometry'] = reference_geom
            print(f"Replaced geometry of polygon {row['polygon_id']} with reference geometry.")
            reference_geom = reference_geom
        else: 
            reference_geom = row['geometry']  # Update reference to current for next iteration
    reference_geom = tree_rows.iloc[49]['geometry'] if len(tree_rows) > 49 else None
    for idx, row in tree_rows[48::-1].iterrows():
        #print(f"Processing polygon {row['polygon_id']}...")
        iou, precision, recall, f1, score = calculate_metrics(reference_geom, row['geometry'])

        if score < 0.5:
            print(f"Polygon {row['polygon_id']} has low similarity score ({score:.2f})")
            # we will transfer the reference geometry
            gdf_filtered.at[idx, 'geometry'] = reference_geom
            tree_rows.at[idx, 'geometry'] = reference_geom
            print(f"Replaced geometry of polygon {row['polygon_id']} with reference geometry.")
            reference_geom = reference_geom
        else: 
            reference_geom = row['geometry']  

# Update gdf with modified geometries from gdf_filtered
for idx, row in gdf_filtered.iterrows():
    polygon_id = row['polygon_id']
    gdf.loc[gdf['polygon_id'] == polygon_id, 'geometry'] = row['geometry']

# Save the updated GeoDataFrame to the master control GPKG
gdf.to_file(master_control_gdf, driver="GPKG", index=False)
print(f"Updated master control GeoDataFrame saved to {master_control_gdf}")
print(f"Total polygons updated: {len(gdf_filtered)}")
