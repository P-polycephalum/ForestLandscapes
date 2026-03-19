import os
import json
import geopandas as gpd
from shapely.geometry import mapping
# Paths

path = r"D:\BCI_50ha_timeseries"
path_polygons = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_crownmap_timeseries.gpkg"
master_control_json = r"D:\BCI_50ha_timeseries\master_control.json"

# Load crownmap
print("Loading crownmap...")
crownmap = gpd.read_file(path_polygons)

# Build master control JSON
master_control = {
    "paths": {
        "root": path,
        "crowns_dir": os.path.join(path, "crowns"),
        "crownmap_path": path_polygons,
        "orthomosaics_dir": os.path.join(path, "Product_local")  
    },
    "variables": {
        "global_id": {
            "description": "Unique tree identifier"
        },
        "date": {
            "description": "Acquisition date for crown polygon. Exception for reference crowns where date is 'reference'"
        },
        "latin": {
            "description": "Latin species name"
        },
        "labels": {
            "description": "Variable labels for the tree crown including leaf coverage, flowering, fruiting, and new leaves"
        },
        "area": {
            "description": "Area of the tree crown polygon in square meters"
        },
        "score": {
            "description": "Score from SAM2"
        },
        "tag": {
            "description": "Tag from ForestGEO census data" 
        },
        "IoU": {
            "description": "Intersection over Union similarity score between segmented crown and reference crown"
        },
        "Precision": {
            "description": "Precision score for the segmented crown compared to the reference crown"
        },
        "Recall": {
            "description": "Recall score for the segmented crown compared to the reference crown"
        },
        "F1": {
            "description": "F1 score for the segmented crown compared to the reference crown"
        },
        "similarity": {
            "description": "Overall similarity score between segmented crown and reference crown, calculated as the average of IoU, Precision, Recall, and F1 scores"
        }
    },
    "crowns": {}
}

# Process each row in crownmap
print("Building crown records...")
for idx, row in crownmap.iterrows():
    polygon_id = row['polygon_id']
    crown_record = {
        "global_id": row['global_id'],
        "date": row['date'],
        "latin": row['latin'],
        "labels": {},
        "area": row['area'],
        "score": row.get('score', None),
        "tag": row.get('tag', None),
        "IoU": row.get('IoU', None),
        "Precision": row.get('Precision', None),
        "Recall": row.get('Recall', None),
        "F1": row.get('F1', None),
        "similarity": row.get('similarity', None)
    }
    master_control["crowns"][str(polygon_id)] = crown_record
    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx + 1} rows...")


# Save to JSON
print(f"\nSaving master control to {master_control_json}...")
with open(master_control_json, 'w') as f:
    json.dump(master_control, f, indent=2)

master_control_gdf= crownmap[['polygon_id', 'geometry']].to_file(os.path.join(path, "master_control_polygons.gpkg"), driver="GPKG")    

print(f"✓ Saved master control polygons to {os.path.join(path, 'master_control_polygons.gpkg')}")
print(f"✓ Created master_control.json with {len(master_control['crowns'])} crowns")