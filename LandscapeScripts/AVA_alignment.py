from arosics import COREG, COREG_LOCAL
import os
import shutil
import rasterio
import numpy as np
from shapely.geometry import box
from rasterio.mask import mask
import geopandas as gpd
from matplotlib import pyplot as plt
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


##We need to standarize the resolution of the orthomosaics and the uint8 format to be able to run the coregistration
# 0.05 is our target resolution, and uint8 is the target dtype
ava_path=r"D:\BCI_ava_timeseries\cropped"
ava_temp_path=r"D:\BCI_ava_timeseries\cropped_resampled"
reference_path= r"D:\BCI_ava_timeseries\cropped\AVA_plot_clipped.tif"
os.makedirs(ava_temp_path, exist_ok=True)

target_resolution = 0.05

for landscape in os.listdir(ava_path):
    input_path = os.path.join(ava_path, landscape)
    output_path = os.path.join(ava_temp_path, landscape)
    
    if os.path.exists(output_path):
        print(f"Already processed: {landscape}")
        continue
    print(f"Processing: {landscape}")
    with rasterio.open(input_path) as src:
        transform = src.transform
        current_res = (transform.a, -transform.e)
        
        scale_x = current_res[0] / target_resolution
        scale_y = current_res[1] / target_resolution
        
        new_width = int(src.width * scale_x)
        new_height = int(src.height * scale_y)
        
        # Create new transform
        new_transform = rasterio.transform.from_bounds(
            src.bounds.left, src.bounds.bottom,
            src.bounds.right, src.bounds.top,
            new_width, new_height
        )
        
        # Read and resample data
        from rasterio.enums import Resampling
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=Resampling.bilinear
        )
        
        # Convert to uint8 if not already
        if src.dtypes[0] != 'uint8':
            for band in range(data.shape[0]):
                # Check if this is the alpha band (last band for 4-band images)
                if band == data.shape[0] - 1 and data.shape[0] == 4:
                    # Set alpha to fully opaque
                    data[band] = np.full_like(data[band], 255, dtype='uint8')
                else:
                    band_data = data[band]
                    band_min, band_max = np.percentile(band_data, (1, 99))
                    band_data_clipped = np.clip(band_data, band_min, band_max)
                    data[band] = ((band_data_clipped - band_min) / (band_max - band_min) * 255).astype('uint8')
        
        # Update metadata
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": new_height,
            "width": new_width,
            "transform": new_transform,
            "dtype": 'uint8'
        })
        
        # Write output
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(data)
        
        print(f"  Resampled to {target_resolution}m, converted to uint8: {landscape}")
       
#this is to resample the reference.
with rasterio.open(reference_path) as src:
        output_path = reference_path.replace("cropped", "cropped_resampled")
        transform = src.transform
        current_res = (transform.a, -transform.e)
        
        scale_x = current_res[0] / target_resolution
        scale_y = current_res[1] / target_resolution
        
        new_width = int(src.width * scale_x)
        new_height = int(src.height * scale_y)
        
        # Create new transform
        new_transform = rasterio.transform.from_bounds(
            src.bounds.left, src.bounds.bottom,
            src.bounds.right, src.bounds.top,
            new_width, new_height
        )
        
        # Read and resample data
        from rasterio.enums import Resampling
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=Resampling.bilinear
        )
        data= np.where(np.isnan(data), 255, data) # Set NaN to 255 before converting to uint8
        
        # Convert to uint8 if not already
        if src.dtypes[0] != 'uint8':
            for band in range(data.shape[0]):
            # Check if this is the alpha band (last band for 4-band images)
                if band == data.shape[0] - 1 and data.shape[0] == 4:
                    # Set alpha to fully opaque
                    data[band] = np.full_like(data[band], 255, dtype='uint8')
                else:
                    band_data = data[band]
                    # Create mask for non-NaN values
                    valid_mask = ~np.isnan(band_data)
                    band_min, band_max = np.percentile(band_data[valid_mask], (1, 99))
                    band_data_clipped = np.clip(band_data, band_min, band_max)
                    data[band] = ((band_data_clipped - band_min) / (band_max - band_min) * 255).astype('uint8')
            
        # Update metadata
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": new_height,
            "width": new_width,
            "transform": new_transform,
            "dtype": 'uint8'
        })
        out_meta['nodata'] = 255  # Set nodata to 255 for uint8
        
        # Write output
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(data)
        
        print(f"  Resampled to {target_resolution}m, converted to uint8: {landscape}")
       
reference=r"D:\BCI_ava_timeseries\AVA_plot_clipped.tif"
reference_permanent=r"D:\BCI_ava_timeseries\AVA_plot_clipped.tif"
print("the reference is", reference)
global_dir= r"D:\BCI_ava_timeseries\Product_global"
os.makedirs(global_dir, exist_ok=True)
successful_alignments = [file for file in os.listdir(global_dir) if file.endswith(".tif")]
successful_alignments.append(os.path.basename(reference)) # Add the original reference to the list of successful alignments
orthomosaics= os.listdir(ava_temp_path)


#global alignment
for orthomosaic in orthomosaics[127::-1]:
    print(f"Processing {orthomosaic}...")
    target = os.path.join(ava_temp_path, orthomosaic)
    output_path = target.replace("cropped_resampled", "Product_global")
    if os.path.isfile(output_path):
        print(f"Global alignment for {orthomosaic} already processed. Skipping...")
        continue
    kwargs2 = { 'path_out': output_path,
                    'fmt_out': 'GTIFF',
                    'r_b4match': 2,
                    's_b4match': 2,
                    'max_shift': 200,
                    'max_iter': 20,
                    'align_grids':True,
                    'match_gsd': True,
                    'binary_ws': False
                }
    alignment_successful = False
    while not alignment_successful and successful_alignments:
        try:
            CR= COREG(reference, target, **kwargs2,ws=(2048,2048))
            CR.calculate_spatial_shifts()
            CR.correct_shifts()
            print("Global alignment successful")
            successful_alignments.append(output_path) # Add successful alignment to the list
            reference = output_path
            alignment_successful = True
        except Exception as e:
                print("Global alignment failed, retrying with closest successful alignment as reference")
                if successful_alignments: # if there are successful alignments to use as reference                    
                    # Extract target date
                    year = orthomosaic.split("_")[2]
                    month = orthomosaic.split("_")[3]
                    day = orthomosaic.split("_")[4]
                    target_date = datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
                    
                    # Extract dates from successful alignments and calculate proximity
                    alignment_dates = []
                    for aligned_file in successful_alignments:
                        try:
                            fname = os.path.basename(aligned_file)
                            y = fname.split("_")[2]
                            m = fname.split("_")[3]
                            d = fname.split("_")[4]
                            file_date = datetime.strptime(f"{y}-{m}-{d}", "%Y-%m-%d")
                            days_diff = abs((file_date - target_date).days)
                            alignment_dates.append((aligned_file, days_diff))
                        except (IndexError, ValueError):
                            continue
                    
                    # Sort by proximity (closest first)
                    alignment_dates.sort(key=lambda x: x[1])
                    
                    # Try each reference in order of proximity
                    for closest_ref, days_diff in alignment_dates:
                        print(f"Trying reference {os.path.basename(closest_ref)} ({days_diff} days difference)")
                        reference = os.path.join(r"D:\BCI_ava_timeseries\Product_global", closest_ref) if not os.path.isabs(closest_ref) else closest_ref
                        try:
                            CR = COREG(reference, target, **kwargs2, ws=(2048, 2048))
                            CR.calculate_spatial_shifts()
                            CR.correct_shifts()
                            print(f"Global alignment successful with reference {os.path.basename(closest_ref)}")
                            successful_alignments.append(output_path)
                            reference = output_path
                            alignment_successful = True
                            break
                        except Exception as e:
                            print(f"Failed with {os.path.basename(closest_ref)}: {e}")
                            continue
                    if not alignment_successful:
                        print("trying with the original")
                        try:
                            CR = COREG(reference_permanent, target, **kwargs2, ws=(2048, 2048))
                            CR.calculate_spatial_shifts()
                            CR.correct_shifts()
                            print("Global alignment successful with original reference")
                            successful_alignments.append(output_path)
                            reference = output_path
                            alignment_successful = True
                        except Exception as e:
                            print(f"Failed with original reference: {e}")
                            print("No successful alignments available, skipping this orthomosaic")
                            break
                                              
for orthomosaic in orthomosaics[128:]:
    print(f"Processing {orthomosaic}...")
    target = os.path.join(ava_temp_path, orthomosaic)
    output_path = target.replace("cropped_resampled", "Product_global")
    if os.path.isfile(output_path):
        print(f"Global alignment for {orthomosaic} already processed. Skipping...")
        continue
    kwargs2 = { 'path_out': output_path,
                    'fmt_out': 'GTIFF',
                    'r_b4match': 2,
                    's_b4match': 2,
                    'max_shift': 200,
                    'max_iter': 20,
                    'align_grids':True,
                    'match_gsd': True,
                    'binary_ws': False
                }
    alignment_successful = False
    while not alignment_successful and successful_alignments:
        try:
            CR= COREG(reference, target, **kwargs2,ws=(2048,2048))
            CR.calculate_spatial_shifts()
            CR.correct_shifts()
            print("Global alignment successful")
            successful_alignments.append(output_path) # Add successful alignment to the list
            reference = output_path
            alignment_successful = True
        except Exception as e:
                print("Global alignment failed, retrying with closest successful alignment as reference")
                if successful_alignments: # if there are successful alignments to use as reference                    
                    # Extract target date
                    year = orthomosaic.split("_")[2]
                    month = orthomosaic.split("_")[3]
                    day = orthomosaic.split("_")[4]
                    target_date = datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
                    
                    # Extract dates from successful alignments and calculate proximity
                    alignment_dates = []
                    for aligned_file in successful_alignments:
                        try:
                            fname = os.path.basename(aligned_file)
                            y = fname.split("_")[2]
                            m = fname.split("_")[3]
                            d = fname.split("_")[4]
                            file_date = datetime.strptime(f"{y}-{m}-{d}", "%Y-%m-%d")
                            days_diff = abs((file_date - target_date).days)
                            alignment_dates.append((aligned_file, days_diff))
                        except (IndexError, ValueError):
                            continue
                    
                    # Sort by proximity (closest first)
                    alignment_dates.sort(key=lambda x: x[1])
                    
                    # Try each reference in order of proximity
                    for closest_ref, days_diff in alignment_dates:
                        print(f"Trying reference {os.path.basename(closest_ref)} ({days_diff} days difference)")
                        reference = os.path.join(r"D:\BCI_ava_timeseries\Product_global", closest_ref) if not os.path.isabs(closest_ref) else closest_ref
                        try:
                            CR = COREG(reference, target, **kwargs2, ws=(2048, 2048))
                            CR.calculate_spatial_shifts()
                            CR.correct_shifts()
                            print(f"Global alignment successful with reference {os.path.basename(closest_ref)}")
                            successful_alignments.append(output_path)
                            reference = output_path
                            alignment_successful = True
                            break
                        except Exception as e:
                            print(f"Failed with {os.path.basename(closest_ref)}: {e}")
                            continue
                    if not alignment_successful:
                        print("trying with the original")
                        try:
                            CR = COREG(reference_permanent, target, **kwargs2, ws=(2048, 2048))
                            CR.calculate_spatial_shifts()
                            CR.correct_shifts()
                            print("Global alignment successful with original reference")
                            successful_alignments.append(output_path)
                            reference = output_path
                            alignment_successful = True
                        except Exception as e:
                            print(f"Failed with original reference: {e}")
                            print("No successful alignments available, skipping this orthomosaic")
                            break
                    

#local alignment attempt

ava_temp_local_dir= r"D:\BCI_ava_timeseries\Product_local"
os.makedirs(ava_temp_local_dir, exist_ok=True)

#reference = r"D:\BCI_ava_timeseries\reference.tif"
for orthomosaic in orthomosaics[127::-1]:
    print(f"Processing {orthomosaic}...")
    print(f"Reference for local alignment: {reference}")
    target = os.path.join(global_dir, orthomosaic)
    out_path = target.replace("Product_global", "Product_local")
    if os.path.isfile(out_path):
        print(f"Local alignment for {orthomosaic} already processed. Skipping...")
        # befor continuing, set the reference to the already aligned orthomosaic, so that the next one can be aligned to it if needed
        reference = out_path
        continue
    else:
        try:

            kwargs = {'grid_res': 200,
                    'window_size': (512, 512),
                    'path_out': out_path,
                    'fmt_out': 'GTIFF',
                    'q': False,
                    'min_reliability': 30,
                    'r_b4match': 2,
                    's_b4match': 2,
                    'max_shift': 100,
                    'nodata': (255,255),
                    'ignore_errors': True,
                    'match_gsd':False,
                            }
            CRL = COREG_LOCAL(reference, target, **kwargs)
            CRL.calculate_spatial_shifts()
            CRL.correct_shifts()
            CRL.CoRegPoints_table.to_csv(out_path.replace("orthomosaic.tif","aligned.csv"))
            reference = out_path
        except Exception as e:
            print(f"Local alignment failed for {orthomosaic}: {e}")
            #break the loop if local alignment fails, as it is likely that subsequent ones will also fail due to the same underlying issue
            break

reference=r"D:\BCI_ava_timeseries\AVA_plot_clipped.tif"
for orthomosaic in orthomosaics[128:]:
    print(f"Processing {orthomosaic}...")
    print(f"Reference for local alignment: {reference}")
    target = os.path.join(global_dir, orthomosaic)
    out_path = target.replace("Product_global", "Product_local")
    if os.path.isfile(out_path):
        print(f"Local alignment for {orthomosaic} already processed. Skipping...")
        # befor continuing, set the reference to the already aligned orthomosaic, so that the next one can be aligned to it if needed
        reference = out_path
        continue
    else:
        try:

            kwargs = {'grid_res': 200,
                    'window_size': (512, 512),
                    'path_out': out_path,
                    'fmt_out': 'GTIFF',
                    'q': False,
                    'min_reliability': 30,
                    'r_b4match': 2,
                    's_b4match': 2,
                    'max_shift': 100,
                    'nodata': (255,255),
                    'ignore_errors': True,
                    'match_gsd':False,
                            }
            CRL = COREG_LOCAL(reference, target, **kwargs)
            CRL.calculate_spatial_shifts()
            CRL.correct_shifts()
            CRL.CoRegPoints_table.to_csv(out_path.replace("orthomosaic.tif","aligned.csv"))
            reference = out_path
        except Exception as e:
            print(f"Local alignment failed for {orthomosaic}: {e}")
            #break the loop if local alignment fails, as it is likely that subsequent ones will also fail due to the same underlying issue
            break