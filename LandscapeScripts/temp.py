import os
import rasterio
import geopandas as gpd
import pandas as pd
from LandscapeScripts.utils import multi_to_polygon
import matplotlib.pyplot as plt
from shapely.geometry import box
from rasterio.mask import mask

wd=r"D:\BCI_ava_timeseries"
vis=os.path.join(wd,"visualization")
os.makedirs(vis, exist_ok=True)
segmented_crowns=os.path.join(wd,"segmented_crowns")
reference_crownmap=os.path.join(wd,"crownmap","BCI_ava_crownmap_2025.gpkg")

all_segmented=[]
for file in os.listdir(segmented_crowns):
    print(f"Processing {file}...")
    if file.endswith(".gpkg"):
        file_path=os.path.join(segmented_crowns,file)
        gdf=gpd.read_file(file_path)
        print(f"Number of rows in {file}: {len(gdf)}")
        gdf['date']= "_".join(file.split(".")[0].split("_")[3:6])
        all_segmented.append(gdf)

all_segmented_gdf=gpd.GeoDataFrame(pd.concat(all_segmented, ignore_index=True))
all_segmented_gdf['geometry'] = all_segmented_gdf.apply(lambda row: multi_to_polygon(row['geometry']) if row['geometry'].geom_type == 'MultiPolygon' else row['geometry'], axis=1)

all_segmented_gdf['date']= pd.to_datetime(all_segmented_gdf['date'], format="%Y_%m_%d")
all_segmented_gdf= all_segmented_gdf.merge(gpd.read_file(reference_crownmap)[['global_id','latin']], on='global_id', how='left')

#now we need the to add the reference crownmap
crownmap= gpd.read_file(reference_crownmap)
crownmap_subset = crownmap[['global_id', 'area_m2', 'latin', 'geometry','tag']].copy()
crownmap_subset = crownmap_subset.rename(columns={'area_m2': 'area'})

for index, row in crownmap_subset.iterrows():
    crownmap_subset.at[index, 'geometry'] = multi_to_polygon(row['geometry'])


crownmap_subset['date']= "reference"

all_segmented_gdf= pd.concat([all_segmented_gdf, crownmap_subset], ignore_index=True)

for idx, row in all_segmented_gdf.iterrows():
    if row['date'] != "reference":
        all_segmented_gdf.at[idx, 'date'] = row['date'].strftime("%Y_%m_%d")
    else:
        all_segmented_gdf.at[idx, 'date'] = "reference"


all_segmented_gdf['polygon_id'] = all_segmented_gdf['global_id'] + "_" + all_segmented_gdf['date']
all_segmented_gdf['similarity']= 1
all_segmented_gdf['F1']= 1
all_segmented_gdf['Precision']= 1
all_segmented_gdf['Recall']= 1
all_segmented_gdf['IoU']= 1

all_segmented_gdf.to_file(os.path.join(wd,'crownmap' ,"BCI_ava_crownmap_timeseries.gpkg"), driver="GPKG")

#pdf visualization
for global_id in global_ids_70_100[:20]:  # Process only the first 10 for demonstration
    subset = all_segmented_gdf[all_segmented_gdf['global_id'] == global_id]
    output_pdf = os.path.join(vis, f"{global_id}.pdf")
    orthomosaic_path = os.path.join(wd, "Product_local")
    generate_leafing_pdf(subset, output_pdf, orthomosaic_path, crowns_per_page=12, variables=['similarity','date'])

orthomosaics=os.path.join(wd,"Product_local")
orthomosaic_files= [os.path.join(orthomosaics, file) for file in os.listdir(orthomosaics) if file.endswith(".tif")]
orthomosaic_files.sort()  # Ensure files are in chronological order
orthomosaic_files[127]
#file visualization
vis_file=os.path.join(wd,'crowns')
os.makedirs(vis_file, exist_ok=True)
for tree in all_segmented_gdf.itertuples():
    global_id = tree.global_id
    date = tree.date if isinstance(tree.date, str) else tree.date.strftime("%Y_%m_%d")
    polygon_id = tree.polygon_id
    if not os.path.exists(os.path.join(vis_file, global_id, f"{polygon_id}.tif")):
        orthomosaic=os.path.join(wd, "Product_local", f"BCI_ava_{date}_orthomosaic.tif")
        with rasterio.open(orthomosaic) as src:
            bounds = tree.geometry.bounds
            box_crown_5 = box(bounds[0] - 5, bounds[1] - 5, bounds[2] + 5, bounds[3] + 5)

            out_image, out_transform = mask(src, [box_crown_5], crop=True)
            out_meta = src.meta.copy()
            x_min, y_min = out_transform * (0, 0)
            xres, yres = out_transform[0], out_transform[4]
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })
            tree_folder= os.path.join(vis_file, global_id)
            os.makedirs(tree_folder, exist_ok=True)
            with rasterio.open(os.path.join(tree_folder, f"{polygon_id}.tif"), "w", **out_meta) as dest:
                dest.write(out_image)
                print(f"✓ Saved visualization for tree {global_id} on {date} to {os.path.join(tree_folder, f'{polygon_id}.tif')}")
    else:
        print(f"✓ Visualization for tree {global_id} on {date} already exists at {os.path.join(tree_folder, f'{polygon_id}.tif')}")

all_segmented_gdf.to_file(os.path.join(wd,'crownmap' ,"BCI_ava_crownmap_timeseries.gpkg"), driver="GPKG")

date_target = .strftime("%Y_%m_%d")
path_orthomosaic = [os.path.join(orthomosaic_path, file) for file in os.listdir(orthomosaic_path) if date_target in file and file.endswith(".tif")]
plt.figure(figsize=(12, 6))
for global_id in all_segmented_gdf['global_id'].unique()[1:10]:
    subset = all_segmented_gdf[all_segmented_gdf['global_id'] == global_id].sort_values('date')
    plt.plot(subset['date'], subset['similarity'], marker='o', label=global_id, alpha=0.7)

plt.xlabel('Date')
plt.ylabel('Similarity')
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.show()


