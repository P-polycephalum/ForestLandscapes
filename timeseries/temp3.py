import os
import geopandas as gpd
import json
import pandas as pd
import matplotlib.pyplot as plt


path= r"D:\BCI_50ha_timeseries\master_control.json"

with open(path) as f:
    data = json.load(f)

cavallinesia_crowns = [crown_id for crown_id in data['crowns'].keys() 
                       if data['crowns'][crown_id].get('latin') == "Cavanillesia platanifolia"]
print(f"Found {len(cavallinesia_crowns)} Cavanillesia platanifolia crowns:")
print(cavallinesia_crowns)
len(cavallinesia_crowns)

# lets create a dataframe with the crown ids and their latin names

crown_data = []
for crown_id, crown_info in data['crowns'].items():
    global_id = crown_info.get('global_id')
    date = crown_info.get('date')
    latin_name = crown_info.get('latin')
    leafing = crown_info.get('labels', {}).get('leafing')
    flowering = crown_info.get('labels', {}).get('flowering')
    area = crown_info.get('area')
    score = crown_info.get('score')
    tag = crown_info.get('tag')
    Iou = crown_info.get('IoU')
    precision = crown_info.get('Precision')
    recall = crown_info.get('Recall')
    F1 = crown_info.get('F1')
    score = crown_info.get('score')
    quality = crown_info.get('quality')
    edited = crown_info.get('edited')
    
    crown_data.append({
        'crown_id': crown_id,
        'global_id': global_id,
        'date': date,
        'latin': latin_name,
        'leafing': leafing,
        'flowering': flowering,
        'area': area,
        'score': score,
        'tag': tag,
        'IoU': Iou,
        'Precision': precision,
        'Recall': recall,
        'F1': F1,
        'quality': quality,
        'edited': edited
    })

crown_df = pd.DataFrame(crown_data)

target_species = "Cavanillesia platanifolia"
cavallinesia_df = crown_df[crown_df['latin'] == target_species]
print(cavallinesia_df.head())


cavallinesia_df['date'] = pd.to_datetime(cavallinesia_df['date'].str.replace('_', '-'))

# Find global_ids that have any "very bad" quality entries
bad_quality_ids = cavallinesia_df[cavallinesia_df['quality'] == 'Very Bad']['global_id'].unique()
print(f"Found {len(bad_quality_ids)} individuals with 'very bad' quality: {bad_quality_ids}")

cavallinesia_df = cavallinesia_df[cavallinesia_df['global_id'] != '87757645-6a6c-41c0-b862-6b192c82f4cb']

# Check for NA leafing values
na_leafing = cavallinesia_df[cavallinesia_df['leafing'].isna()]
na_leafing_ids = na_leafing['global_id'].unique()
print(f"Found {len(na_leafing_ids)} individuals with NA leafing values: {na_leafing_ids}")
print(f"Found {len(na_leafing)} entries with NA leafing values:")


cavallinesia_df['doy'] = cavallinesia_df['date'].dt.dayofyear

for global_id in cavallinesia_df['global_id'].unique():
    indv_data = cavallinesia_df[cavallinesia_df['global_id'] == global_id].sort_values('date')
    print(f"Global ID: {global_id}")
    print(indv_data[['date', 'leafing', 'quality']])
    print("\n")
    plt.figure(figsize=(10, 4))
    plt.plot(indv_data['date'], indv_data['leafing'], marker='o', label=f'ID {global_id}', linewidth=2)
    plt.title(f"Leafing of {target_species} Crown {global_id} Over Time")
    plt.xlabel("Date")
    plt.ylabel("Leafing")
    plt.xticks(rotation=45)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


