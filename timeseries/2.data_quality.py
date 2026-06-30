import os
import pandas as pd
import geopandas as gpd
from shapely import box
import matplotlib.pyplot as plt
import shapely
import numpy as np
from statistics import mode
from PIL import Image

#load main dataset of labels 
labels_path=r"timeseries/dataset_raw/export_2026_04_23.csv"
labels=pd.read_csv(labels_path)

#Keep only good segmentation
crowns_labeled= labels[labels["segmentation"]=="good"]

#how many repeated ones
label_counts = crowns_labeled.groupby("polygon_id").size()
dist = label_counts.value_counts().sort_index()
dist.index.name = "n_labels"
dist.name = "n_crowns"
print(dist.to_frame())




#mnost labeled species
species_counts = crowns_labeled['latin'].value_counts()
species_counts_filtered = species_counts[species_counts >= 200]
fig, ax = plt.subplots(figsize=(14, max(6, len(species_counts_filtered) * 0.3)))
species_counts_filtered.plot.barh(ax=ax)
ax.set_xlabel('Count')
ax.set_ylabel('Species')
ax.set_title('Labeled Species Distribution')
ax.tick_params(axis='y', labelsize=8)
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.tight_layout()
plt.show()

#summary of flowering vs non flowering
crowns_labeled
flowering_counts = crowns_labeled['isFlowering'].value_counts()
flowering_summary = pd.DataFrame({'count': flowering_counts, 'percentage': flowering_counts / flowering_counts.sum() * 100})
print(flowering_summary)



flowering = crowns_labeled[(crowns_labeled['isFlowering']== "yes") | (crowns_labeled['isFlowering']== "maybe")]

species_counts_flowering = flowering['latin'].value_counts()
species_counts_flowering = species_counts_flowering[species_counts_flowering >= 20]
fig, ax = plt.subplots(figsize=(14, max(6, len(species_counts_flowering) * 0.3)))
species_counts_flowering.plot.barh(ax=ax)
ax.set_xlabel('Count')
ax.set_ylabel('Species')
ax.set_title('Flowering Labeled Species Distribution')
ax.tick_params(axis='y', labelsize=8)
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.tight_layout()
plt.show()



def customLeafing(leafing_values):
    values = list(leafing_values)
    sd_values = np.std(values)
    if len(values) == 1:
        return values[0]
    if len(values) >= 2 and sd_values <= 5:
        result= sum(values) / len(values)
        print("2 or more values with sd <=5",values )
        return result
    if len(values) >= 2 and sd_values > 5:
        try:
            reference_value = mode(values)
        except:
            reference_value = np.median(values)
        filtered_values = [v for v in values if abs(v - reference_value) <= 5]
        if filtered_values:
            result=sum(filtered_values) / len(filtered_values)
            print("2 or more values sd>5", values, "output: ", result)
            return result
        else:
            result= None
            print("2 or more values sd>5", values, "output: ", result)
            return result

def customFloweringNumeric(floweringN):
    values = list(floweringN)
    sd_values = np.std(values)
    if len(values) == 1:
        return values[0]
    if len(values) >= 2 and sd_values <= 5:
        result= sum(values) / len(values)
        return result
    if len(values) >= 2 and sd_values > 5:
        try:
            reference_value = mode(values)
        except:
            reference_value = np.median(values)
        filtered_values = [v for v in values if abs(v - reference_value) <= 5]
        if filtered_values:
            result=sum(filtered_values) / len(filtered_values)
            return result
        else:
            result= None
            return result

def customFlowering(floweringValues): 
    floweringValues=list(floweringValues)
    if len(floweringValues)==1:
        return floweringValues[0]
    if all(value == floweringValues[0] for value in floweringValues): 
        return floweringValues[0]         
    if "maybe" in floweringValues:
        return "maybe"
    if "yes" in floweringValues and "no" in floweringValues:
        return "maybe"
    else:
        return "no"
    
crowns_labeled_avg = crowns_labeled.groupby("polygon_id").agg({
    "leafing": customLeafing,
    "isFlowering":customFlowering,
    "floweringIntensity": customFloweringNumeric,
    "latin": "first",  
    "date":"first",
}).reset_index()

crowns_labeled_avg = crowns_labeled_avg[~crowns_labeled_avg['leafing'].isna()]
label_counts = crowns_labeled_avg.groupby("polygon_id").size()
dist = label_counts.value_counts().sort_index()
dist.index.name = "n_labels"
dist.name = "n_crowns"
print(dist.to_frame())

import geopandas as gpd
master_parts_dir = r"D:\BCI_50ha_timeseries\master_gdf_parts"  # new partitioned folder

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
            missing_crown_area = part_gdf["crown_area"].isna()
            if missing_crown_area.any():
                part_gdf.loc[missing_crown_area, "crown_area"] = computed_crown_area.loc[missing_crown_area]

    return gpd.GeoDataFrame(part_gdf, crs="EPSG:32617")

def _write_part(gdf, bucket_id, time_str):
    """Write a GeoDataFrame partition to disk."""
    bucket_dir = os.path.join(master_parts_dir, str(bucket_id))
    os.makedirs(bucket_dir, exist_ok=True)
    part_path = os.path.join(bucket_dir, f"{_safe_time(time_str)}.parquet")
    gdf.to_parquet(part_path, index=False)

buckets= os.listdir(master_parts_dir)
crowns_labeled_avg.columns

label_cols = ['leafing', 'isFlowering', 'floweringIntensity']
tmp_cols = {c: f"_new_{c}" for c in label_cols}

for bucket in buckets:
    bucket_path= os.path.join(master_parts_dir, bucket)
    times= os.listdir(bucket_path)
    times= [t.replace(".parquet","") for t in times]
    for t in times:
        data= read_part_by_time(bucket, t)
        data['polygon_id']= data['GlobalID'].astype(str)+"_"+data['time'].str[:10]

        if not data.empty:
            print(f"Successfully read partition for bucket {bucket} and time {t}")
            incoming = crowns_labeled_avg[['polygon_id'] + label_cols].rename(columns=tmp_cols)
            merged_data = data.merge(incoming, on='polygon_id', how='left')
            # Prefer master (data) values; only fill in from crowns_labeled_avg where master has no value
            for col, tmp_col in tmp_cols.items():
                if col not in merged_data.columns:
                    merged_data[col] = merged_data[tmp_col]
                else:
                    merged_data[col] = merged_data[col].combine_first(merged_data[tmp_col])
                merged_data.drop(columns=[tmp_col], inplace=True)
            

            n_with_labels = merged_data['leafing'].notna().sum()
            n_without_labels = merged_data['leafing'].isna().sum()
            if n_with_labels > 0:
                 print(f"Rows with labels: {n_with_labels}, Rows without labels: {n_without_labels}")
                 print(f"Sample of rows with labels for bucket {bucket} and time {t}:\n{merged_data[merged_data['leafing'].notna()].head()}")
            _write_part(merged_data, bucket, t)
        else:
            print(f"No data found for bucket {bucket} and time {t}")

tile= "0_0"
time= "2018-04-04T00-00-00"
crowns_labeled_avg['polygon_id']
data= read_part_by_time(tile, time)
data.columns

crowns_labeled_avg['polygon_id']
data['polygon_id']= data['GlobalID'].astype(str)+"_"+data['time'].str[:10]



#split the dataset to the flowers dataset
flower_dataset=crowns_labeled_avg[(crowns_labeled_avg["isFlowering"]=="yes")|(crowns_labeled_avg["isFlowering"]=="maybe")]
path_out= os.path.join(data_path,'flower_dataset')
#extract the features, crown based 
for i, (_, row) in enumerate(flower_dataset.iterrows()):
    print(f"Processing iteration {i + 1} of {len(flower_dataset)}")
    if not os.path.exists(os.path.join(path_out, row['polygon_id']+".png")):
        path_orthomosaic = os.path.join(data_path, 'orthomosaic_aligned_local', f"BCI_50ha_{row['date']}_local.tif")
        try:
            with rasterio.open(path_orthomosaic) as src:
                bounds = row.geometry.bounds
                box_crown_5 = box(bounds[0] - 5, bounds[1] - 5, bounds[2] + 5, bounds[3] + 5)
                print(box_crown_5)
                out_image, out_transform = mask(src, [box_crown_5], crop=True)
                x_min, y_min = out_transform * (0, 0)
                xres, yres = out_transform[0], out_transform[4]

                transformed_geom = shapely.ops.transform(
                        lambda x, y: ((x - x_min) / xres, (y - y_min) / yres),
                        row.geometry
                    )
                
                img_name = f"{row['polygon_id']}.png"
                img_path = os.path.join(path_out, img_name)

                fig, ax = plt.subplots(figsize=(10, 10))

                ax.imshow(out_image.transpose((1, 2, 0))[:, :, 0:3])

                ax.plot(*transformed_geom.exterior.xy, color='red')

                for interior in transformed_geom.interiors:
                    ax.plot(*interior.xy, color='red')

                ax.axis('off')

                fig.savefig(img_path, bbox_inches='tight', pad_inches=0)
                plt.close(fig)
                
                print(f"Saved: {img_path}")

        except Exception as e:
            print(f"Error processing {row['polygon_id']}: {e}")
    else:
        print("it already exists in dataset")


flower_csv= flower_dataset.drop(columns=['geometry'])
flower_csv= flower_csv[['isFlowering', 'leafing', 'floweringIntensity',
       'polygon_id', 'date','latin', 'area', 'score', 'tag', 'iou']]
flower_csv.to_csv(r'timeseries/dataset_corrections/flower.csv')

###break to call labelbox flowering


non_flower= crowns_labeled_avg[crowns_labeled_avg['isFlowering']=="no"]
#read the flower dataset back in 
flower_dataset=pd.read_csv(r"timeseries/dataset_corrections/flower_out.csv")
# get rid of the flowering lianas
flower_dataset= flower_dataset[flower_dataset['flowering_liana']!="yes"]

#merge no flower and flower datasets

all= pd.concat([non_flower,flower_dataset])
all.columns
all_crowns= all.merge(crowns[['latin','polygon_id','geometry']], left_on='polygon_id',right_on='polygon_id', how='left')
all_crowns['latin'] = all_crowns['latin_x'].combine_first(all_crowns['latin_y'])
all_crowns['geometry'] = all_crowns['geometry_x'].combine_first(all_crowns['geometry_y'])
all_crowns.drop(columns=['latin_x', 'latin_y','geometry_x','geometry_y'], inplace=True)

sgbt_dataset = gpd.GeoDataFrame(all_crowns , geometry='geometry')
sgbt_dataset.set_crs("EPSG:32617", allow_override=True, inplace=True)  

#make sure the combined dataset has sound estimates in leafing
print(sgbt_dataset['leafing'].describe())

#one comes out for the sgbt dataset
sgbt_dataset.to_file(r'timeseries/dataset_training/train.shp')


















##lets try to get prioria copaifera all flowering inividuals and the same amount of not leafing
prioria_all = crowns_labeled_avg[crowns_labeled_avg['latin'] == "Prioria copaifera"]
generate_leafing_pdf(prioria_all[0:12], r'plots/prico_normal.pdf',orthomosaic_path, crowns_per_page=12,variables=['floweringIntensity','isFlowering','date'] )



#i only want to keep dipteryx and jacaranda for my cnn of flowers
flower_cnn= crowns_labeled_avg[(crowns_labeled_avg['latin']=='Jacaranda copaia')|(crowns_labeled_avg['latin']=='Dipteryx oleifera')]

def classify(row):
    # If not flowering
    if row['isFlowering'] == 'no':
        return 0

    # If flowering (yes or maybe) and species is Dipteryx
    elif row['isFlowering'] in ['yes', 'maybe'] and row['latin'] == 'Dipteryx oleifera':
        return 1

    # If flowering (yes or maybe) and species is Jacaranda
    elif row['isFlowering'] in ['yes', 'maybe'] and row['latin'] == 'Jacaranda copaia':
        return 2

    else:
        print(row['isFlowering'], row['latin'])
        return None  

flower_cnn['class'] = flower_cnn.apply(classify, axis=1)
flower_cnn['date'] = flower_cnn['polygon_id'].str.split('_').str[1]
flower_cnn['date'] = flower_cnn['date'].str.replace('-', '_')

flower_cnn[flower_cnn['class'].isna()]

orthomosaic_path=os.path.join(data_path,"orthomosaic_aligned_local")
#list of orthomosaics
orthomosaic_list=os.listdir(orthomosaic_path)
import numpy as np
from PIL import Image
path_out=os.path.join(data_path,"flower_data")
os.makedirs(path_out, exist_ok=True)
for i, (_, row) in enumerate(flower_cnn.iterrows()):
    print(f"Processing iteration {i + 1} of {len(flower_cnn)}")
    if not os.path.exists(os.path.join(path_out, row['polygon_id']+".png")):
        path_orthomosaic = os.path.join(orthomosaic_path, f"BCI_50ha_{row['date']}_local.tif")
        try:
            with rasterio.open(path_orthomosaic) as src:
                out_image, out_transform = mask(src, [row.geometry], crop=True)
                img_array = np.moveaxis(out_image, 0, -1) 
                img_array = img_array.astype(np.uint8)
                img_name = f"{row['polygon_id']}.png"
                img_path = os.path.join(path_out, img_name)
                Image.fromarray(img_array).save(img_path)
                
                print(f"Saved: {img_path}")

        except Exception as e:
            print(f"Error processing {row['polygon_id']}: {e}")
    else:
        print("it already exists in dataset")

flower_cnn['file']= flower_cnn['polygon_id']+".png"
flower_cnn[['file','class']].to_csv(r'timeseries/dataset_training/train_cnn_flower.csv')


#verify you dont have empty values in leafing
crowns_labeled_avg=crowns_labeled_avg[~crowns_labeled_avg['leafing'].isna()]


cnn_dataset= crowns_labeled_avg[['polygon_id','leafing']]
cnn_dataset['polygon_id']=cnn_dataset["polygon_id"]+'.png'

cnn_dataset.to_csv(r'timeseries/dataset_training/train_cnn.csv')
path_out= os.path.join(data_path,"train_dataset")
for i, (_, row) in enumerate(crowns_labeled_avg.iterrows()):
    print(f"Processing iteration {i + 1} of {len(crowns_labeled_avg)}")
    if not os.path.exists(os.path.join(path_out, row['polygon_id']+".png")):
        path_orthomosaic = os.path.join(orthomosaic_path, f"BCI_50ha_{row['date']}_local.tif")
        try:
            with rasterio.open(path_orthomosaic) as src:
                out_image, out_transform = mask(src, [row.geometry], crop=True)
                img_array = np.moveaxis(out_image, 0, -1) 
                img_array = img_array.astype(np.uint8)
                img_name = f"{row['polygon_id']}.png"
                img_path = os.path.join(path_out, img_name)
                Image.fromarray(img_array).save(img_path)
                
                print(f"Saved: {img_path}")

        except Exception as e:
            print(f"Error processing {row['polygon_id']}: {e}")
    else:
        print("it already exists in dataset")


sgbt_dataset= crowns_labeled_avg[['polygon_id','leafing','isFlowering']].merge(crowns[['date','geometry','polygon_id']],
                              left_on="polygon_id",
                                right_on="polygon_id",
                                  how="left")





