import os
import pandas as pd
import json

folder_name = "export_2026_06_18"

location=os.path.join(r'timeseries', folder_name)
files= os.listdir(location)

# with open(os.path.join(location,files[1]), "r") as f:
#     observations = json.load(f)  

with open(os.path.join(location, files[0]), "r") as f:
    observations = [json.loads(line) for line in f]

with open(os.path.join(location, files[1]), "r") as f:
    plants= [json.loads(line) for line in f]

data = []
for obs in observations:
    path_items = obs.get("__key__", {}).get("path", "").replace('"', '').split(", ")
    second_item = path_items[1] if len(path_items) > 1 else None  
    data.append({
        "isFlowering": obs.get("isFlowering"),
        "leafing": obs.get("leafing"),
        "floweringIntensity": obs.get("floweringIntensity"),
        "segmentation": obs.get("segmentation"),
        "observation_id": obs.get("__key__", {}).get("name"),
        "polygon_id": second_item,
    })

df = pd.DataFrame(data)

data_plants= []

for obs in plants:
    data_plants.append({
        "date": obs.get("date"),
        "globalId": obs.get("globalId"),
        "latin": obs.get("latinName"),
        "polygon_id":  obs.get("__key__", {}).get("name"),
    })

df_plants= pd.DataFrame(data_plants)

df_merged= df.merge(df_plants, left_on='polygon_id', right_on='polygon_id', how='left')

#export as csv
df_merged.to_csv(f"timeseries/dataset_raw/{folder_name}.csv")



#lets examine df_merged
df_merged.columns
df_merged.value_counts("latin")
df_merged.value_counts("segmentation")

df_merged= df_merged[df_merged["segmentation"] == "good"]
df_merged.value_counts("isFlowering")

df_flowering= df_merged[df_merged["isFlowering"] == "yes"]
df_flowering.value_counts("latin")

path_crownmap_timeseries=r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_crownmap_timeseries.gpkg"
path_master_gdf_parts= r"C:\Users\vasquezvicente\repo\labeler\data\master_gdf_parts"






#replace yes with Yes and no with No in isFlowering'
df_merged["isFlowering"]= df_merged["isFlowering"].replace({"yes": "Yes", "no": "No"})
latin_counts= df_merged.value_counts("latin")
globalId_counts= df_merged.value_counts("globalId")[1:100]


for latin, count in latin_counts.items():
    print(f"{latin}: {count}")
for globalId, count in globalId_counts.items():
    which_species= df_merged[df_merged["globalId"] == globalId]["latin"].unique()
    print(f"{globalId}: {count}, Species: {', '.join(which_species)}")
    

import geopandas as gpd

# we will keep this self contain for safety of transfering the labels
target_globalid= "e35fe388-2ed7-43b6-bea8-e6f0a8f899c9"
#examine the master gdf tiles looking for the tile that contains the target globalId
for tile_folder in os.listdir(path_master_gdf_parts):
    parkeets= os.listdir(os.path.join(path_master_gdf_parts, tile_folder))
    #since all parkets have the same trees, we can just check the first parkeet
    first_parkeet= gpd.read_parquet(os.path.join(path_master_gdf_parts, tile_folder, parkeets[0]))
    if target_globalid in first_parkeet["GlobalID"].values:
        print(f"Found target globalId in tile: {tile_folder}")
        break

target_parkeets= os.listdir(os.path.join(path_master_gdf_parts, tile_folder))
label= df_merged[df_merged["globalId"] == target_globalid]

for _, row in label.iterrows():
    polygon_id = row["polygon_id"]
    date = row["date"].replace("_", "-")
    leafing = float(row["leafing"])
    parquet_path = os.path.join(
        path_master_gdf_parts,
        tile_folder,
        date + "T00-00-00.parquet"
    )
    parkeet_gdf = gpd.read_parquet(parquet_path)
    parkeet_gdf["polygon_id"] = (
        parkeet_gdf["GlobalID"].astype(str)
        + "_"
        + parkeet_gdf["time"].astype(str).str[:10]
    )
    mask = parkeet_gdf["polygon_id"] == polygon_id
    print(polygon_id, "matches:", mask.sum())
    parkeet_gdf.loc[mask, "leafing"] = leafing
    parkeet_gdf = parkeet_gdf.drop(columns="polygon_id")
    print(f"Updated leafing for polygon_id: {polygon_id} in parkeet: {date}T00-00-00.parquet")
    parkeet_gdf.to_parquet(parquet_path, index=False)

for parkeet in target_parkeets:
    parkeet_gdf= gpd.read_parquet(os.path.join(path_master_gdf_parts, tile_folder, parkeet))
    parkeet_gdf['polygon_id']= parkeet_gdf['GlobalID']+ "_" + parkeet_gdf["time"].astype(str).str[:10]
    
    #does label polygon_id exist in parkeet_gdf polygon_id?
    if not label["polygon_id"].isin(parkeet_gdf["polygon_id"]).any():
        print(f"polygon_id from label does not exist in parkeet_gdf for parkeet: {parkeet}")
        continue
    if "leafing" in parkeet_gdf.columns:
        print("leafing exist in parket gdf")
        parket_gdf_temp= parkeet_gdf.merge(label[["polygon_id", "leafing", "isFlowering"]], left_on="polygon_id", right_on="polygon_id", how="left")
        parket_gdf_temp['leafing']= parket_gdf_temp.apply(lambda row: row['leafing_x'] if pd.notna(row['leafing_x']) else row['leafing_y'], axis=1)
    else:
        parket_gdf_temp= parkeet_gdf.merge(label[["polygon_id", "leafing", "isFlowering"]], left_on="polygon_id", right_on="polygon_id", how="left")
    
    parket_gdf_temp["flowering"] = parket_gdf_temp.apply(
            lambda row: (
                row["flowering"]
                if row["flowering"] == "Yes" or pd.isna(row["flowering"])
                else row["isFlowering"]
            ),
            axis=1
        )
    parket_gdf_temp['polygon_id']
    #drop the columns that end with _x and _y
    parket_gdf_temp= parket_gdf_temp.drop(columns=[col for col in parket_gdf_temp.columns if col.endswith("_x") or col.endswith("_y")])
    #drop the polygon_id column
    parket_gdf_temp= parket_gdf_temp.drop(columns=["polygon_id", "isFlowering"])
    #ensure leafing and are both numericals
    parket_gdf_temp["leafing"]= pd.to_numeric(parket_gdf_temp["leafing"], errors="coerce")
    parket_gdf_temp.to_parquet(os.path.join(path_master_gdf_parts, tile_folder, parkeet), index=False)
    print(f"Updated parkeet: {parkeet} with flowering and leafing labels")


