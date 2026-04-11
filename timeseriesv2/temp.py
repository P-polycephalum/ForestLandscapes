import os
import geopandas as gpd
import pandas as pd

parquet_folder = r"D:\BCI_50ha_timeseries\master_gdf_parts"

# ── folder layout ──────────────────────────────────────────
# master_gdf_parts/
#   {bucket_id}/
#     {safe_time}.parquet   e.g. 2023-01-15T00-00-00.parquet

# 1. Read a single (bucket, time) slice
def read_part(bucket_id, time_str):
    safe_time = time_str.replace(":", "-").replace("/", "-")
    path = os.path.join(parquet_folder, str(bucket_id), f"{safe_time}.parquet")
    return gpd.read_parquet(path)

gdf = read_part("0_0", "2023-04-25T00:00:00")

# 2. Read all times for one bucket
def read_bucket(bucket_id):
    bucket_dir = os.path.join(parquet_folder, str(bucket_id))
    files = sorted(f for f in os.listdir(bucket_dir) if f.endswith(".parquet"))
    return gpd.GeoDataFrame(
        pd.concat([gpd.read_parquet(os.path.join(bucket_dir, f)) for f in files],
                  ignore_index=True),
        crs="EPSG:32617"
    )

bucket0 = read_bucket(0)

# 3. Read everything (all buckets, all times) — equivalent to the old monolithic file
def read_all():
    parts = []
    for bucket_id in sorted(os.listdir(parquet_folder)):
        bucket_dir = os.path.join(parquet_folder, bucket_id)
        if not os.path.isdir(bucket_dir):
            continue
        for f in sorted(os.listdir(bucket_dir)):
            if f.endswith(".parquet"):
                parts.append(gpd.read_parquet(os.path.join(bucket_dir, f)))
    return gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs="EPSG:32617")

master_gdf = read_all()

# 4. Inspect available buckets / times
buckets = sorted(os.listdir(parquet_folder))
print(buckets)

for b in buckets:
    times = sorted(os.listdir(os.path.join(parquet_folder, b)))
    print(f"bucket {b}: {len(times)} timesteps")