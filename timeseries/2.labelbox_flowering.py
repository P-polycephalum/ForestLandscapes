import labelbox
import os
import pandas as pd
import numpy as np
from PIL import Image
import pickle
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
from shapely.geometry import box
import shapely
import matplotlib.pyplot as plt
from PIL import Image

## data paths
data_path=r"\\stri-sm01\ForestLandscapes\UAVSHARE\BCI_50ha_timeseries"
orthomosaic_path=os.path.join(data_path,"orthomosaic_aligned_local")
orthomosaic_list=os.listdir(orthomosaic_path)
data= pd.read_csv(r"timeseries/dataset_predictions/jacaranda_sgbt.csv")

#import the flower model
with open(r'timeseries/models/xgb_model_flower.pkl', 'rb') as file:
      model_flower = pickle.load(file)

X=data[['rccM', 'gccM', 'bccM', 'ExGM', 'gvM', 'npvM', 'shadowM','rSD', 'gSD', 'bSD',
       'ExGSD', 'gvSD', 'npvSD', 'gcorSD', 'gcorMD','entropy','elevSD']]

Y= data[['area', 'score', 'tag', 'GlobalID', 'iou',
       'date', 'latin', 'polygon_id']]

X_predict_flower= model_flower.predict(X)
df_final = Y.copy()  # Copy Y to keep the same structure
df_final['isFlowering_predicted'] = X_predict_flower

path_crowns=os.path.join(data_path,r"geodataframes\BCI_50ha_crownmap_timeseries.shp")
crowns=gpd.read_file(path_crowns)
crowns['polygon_id']= crowns['GlobalID']+"_"+crowns['date'].str.replace("_","-")
flower_cnn= df_final[(df_final['isFlowering_predicted']==1)|(df_final['isFlowering_predicted']==2)]
flower_cnn['polygon_id'] = flower_cnn['GlobalID'].astype(str) + "_" + flower_cnn['date'].str.replace("_", "-", regex=False)
                                                                                                     
flower_cnn_2= crowns[['polygon_id','geometry']].merge(flower_cnn, left_on='polygon_id', right_on='polygon_id',how='right')

path_out=os.path.join(data_path,"flower_dataset")
os.makedirs(path_out, exist_ok=True)
for i, (_, row) in enumerate(flower_cnn_2.iterrows()):
    print(f"Processing iteration {i + 1} of {len(flower_cnn_2)}")
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


#################################################
#call the labelbox API to pull the labels#


#get the dataset that is currently in labelbox
dataset = client.get_dataset("cm8bs9pgf00d40746btooascw")

#export the dataset that is currently in labelbox
export_task = dataset.export()
export_task.wait_till_done()

#function to parse the labelbox object to a json
def json_stream_handler(output: labelbox.BufferedJsonConverterOutput):
  print(output.json)

#parse the labelbox object

stream=export_task.get_buffered_stream(stream_type=labelbox.StreamType.RESULT).start(stream_handler=json_stream_handler)

#get the names that are already in labelbox
export_json = [data_row.json for data_row in export_task.get_buffered_stream()]
global_keys = [item["data_row"]["external_id"] for item in export_json]
print("Number of data rows in labelbox:", len(global_keys))

#path to the raw flower dataset hosted in ForestLandscapes server
data_path=r"\\stri-sm01\ForestLandscapes\UAVSHARE\BCI_50ha_timeseries"
path_out= os.path.join(data_path,'flower_dataset')
list_flower = [os.path.join(path_out, filename) for filename in os.listdir(path_out)]

#determine the ones the flowering trees that are in the server and not in labelbox
extra_files = [f for f in list_flower if f not in global_keys]
extra_file_ids = [os.path.basename(f) for f in extra_files]

#read in the raw flowering dataset
flowering_metadata=pd.read_csv(r'timeseries\dataset_corrections\flower.csv')  ##the predicted flowering does not have the metadata
flowering_metadata['polygon_id']=flowering_metadata['polygon_id']+".png"
extra_metadata_rows = flowering_metadata[flowering_metadata['polygon_id'].isin(extra_file_ids)]
extra_metadata_rows = extra_metadata_rows.drop_duplicates(subset='polygon_id', keep=False)


##extra_file_ids
extra_file_ids
extra_metadata_rows= pd.DataFrame({'polygon_id': extra_file_ids})
crowns['polygon_id']= crowns['GlobalID']+"_"+crowns['date'].str.replace("_","-")+".png"
extra_metadata_rows= extra_metadata_rows.merge(crowns[['polygon_id','latin','date']],how='left', right_on='polygon_id', left_on='polygon_id')

extra_metadata_rows['floweringIntensity']= 0
extra_metadata_rows['isFlowering']= 'maybe'
extra_metadata_rows['leafing']= 0
# create an assets package for labelbox
assets=[]
for index, row in extra_metadata_rows.iterrows():
    row_data_path = os.path.join(
        "\\\\stri-sm01\\ForestLandscapes\\UAVSHARE\\BCI_50ha_timeseries\\flower_dataset",
        f"{row['polygon_id']}"
    )
    dict_row = {
        'row_data': row_data_path,
        'global_key': row['polygon_id'],
        'media_type': 'IMAGE',
        'metadata_fields': [
            {"schema_id": "cm8bxtvk600q1074v0ad2cds7", "value": str(row['latin'])},
            {"schema_id": "cm8by8tfz08lh074lbwv5430o", "value": row['floweringIntensity']},
            {"schema_id": "cm8bydleg08dx0730cl2144z7", "value": row['isFlowering']},
            {"schema_id": "cm8bypte30abx075t61j4cvs7", "value": row['leafing']},
            {"schema_id": "cm8byrjeg08h8074mewqdf570", "value": row['date']}
        ]
    }
    assets.append(dict_row)

#push the data rows to labelbox
try:
    task = dataset.create_data_rows(assets)
    task.wait_till_done()
except Exception as err:
    print(f'Error while creating labelbox dataset -  Error: {err}')
    task.errors

import pandas as pd
import labelbox

def json_stream_handler(output: labelbox.BufferedJsonConverterOutput):
  print(output.json)
client = labelbox.Client(api_key="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJ1c2VySWQiOiJjbHRoYjZ4d2cwNXA3MDcxMGMzeTk1bTZjIiwib3JnYW5pemF0aW9uSWQiOiJjanF6bGVoNHcweGhxMDc3N2dqNDF0cHNrIiwiYXBpS2V5SWQiOiJjbW9zc2tpYTIwZWxlMDd0eTM5NmMwcTZoIiwic2VjcmV0IjoiMTFjYjFhY2E1Y2I0NmUzMWJjMjAzMjMyODI1MTMyNzUiLCJpYXQiOjE3Nzc5OTU0MjUsImV4cCI6MTc4MDQxNDYyNX0.I8vbnl7dBQEcxmgzWPZWa4sSqmanNKbDgHSYMTc2AXk")

#bring back the exports after labelling
project = client.get_project("cm8azfo2f037h074jgzid05f9")
export_task = project.export()
export_task.wait_till_done()
stream=export_task.get_buffered_stream(stream_type=labelbox.StreamType.RESULT).start(stream_handler=json_stream_handler)
export_json = [data_row.json for data_row in export_task.get_buffered_stream()]
global_keys = [item["data_row"]["external_id"] for item in export_json]


data = [] 
for row in export_json:
    polygon_id = row["data_row"]["external_id"]  # Extract polygon ID
    row_data = {"polygon_id": polygon_id}
    for project_id, project_data in row['projects'].items():  # every project is a columns
        for label in project_data['labels']:  # Access labels
            classifications = label['annotations']['classifications']
            
            for classification in classifications:
                if classification['name'] == 'newLeaves':  # Look for "newLeaves"
                    row_data['newLeaves']= classification['radio_answer']['value']
                elif classification['name'] == 'isFlowering':
                    row_data['isFlowering']= classification['radio_answer']['value']
                elif classification['name'] == 'isFruiting':
                    row_data['isFruiting']= classification['radio_answer']['value']
                elif classification['name']== 'segmentation':
                    row_data['segmentation']=classification['radio_answer']['value']
                elif classification['name'] == 'leafing':  # Fixed condition
                    row_data['leafing'] = classification['text_answer']['content']
                elif classification['name']== 'floweringIntensity':
                    row_data['floweringIntensity']=classification['text_answer']['content']
                elif classification['name'] == 'flowering_liana':
                    row_data['flowering_liana'] = [answer['name'] for answer in classification['checklist_answers']][0]
    data.append(row_data)
     
flowering_dataset= pd.DataFrame(data)

import os
import geopandas as gpd
flowering_dataset['polygon_id'] = flowering_dataset['polygon_id'].apply(os.path.basename)
flowering_dataset['polygon_id'] = flowering_dataset['polygon_id'].apply(lambda x: x.split(".")[0])
flowering_dataset=flowering_dataset[~flowering_dataset['leafing'].isna()]

flowering_dataset['isFlowering'].value_counts()
crownmap = r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_2022_2023_crownmap_raw.shp"
crownmap_gdf = gpd.read_file(crownmap)
flowering_dataset['GlobalID']=flowering_dataset['polygon_id'].apply(lambda x: "_".join(x.split("_")[:-1]))

crownmap_gdf['polygon_id']= crownmap_gdf['GlobalID']+"_"+crownmap_gdf['date'].str.replace("_","-")

flowering_dataset= flowering_dataset.merge(crownmap_gdf[['GlobalID','latin']], on='GlobalID', how='left')


flowering_dataset.to_csv('timeseries/dataset_corrections/flower_out.csv')

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
flowering_dataset[flowering_dataset['isFlowering'].isin(['maybe', 'yes'])]['latin'].value_counts().plot(kind='barh')
plt.title('Count of Flowering Trees by Latin Name')
plt.xlabel('Count')
plt.ylabel('Latin Name')
plt.show()