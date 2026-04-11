import os
import geopandas as gpd

path= r"D:\BCI_50ha_timeseries\master_gdf.gpkg"
if os.path.exists(path):
    master_gdf=gpd.read_file(path)


#unique dates
unique_dates = master_gdf['date'].unique()
print(f"Unique dates in master GDF:{len(unique_dates)} unique dates.")
unique_trees= master_gdf['global_id'].nunique()
print(f"Unique trees in master GDF: {unique_trees} unique trees.")

print(f"Total rows in master GDF: {len(master_gdf)}")
if len(master_gdf) != unique_trees * len(unique_dates):
    print("Warning: The number of rows in the master GDF does not match the expected number based on unique trees and dates.")



