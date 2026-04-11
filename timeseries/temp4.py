import os
import geopandas as gpd
import rasterio.mask as mask
import rasterio
from shapely.geometry import box

crownmap=gpd.read_file(r"D:\PNM\PNM_crownmap_2025.gpkg")
orthomosaic=gpd.read_file(r"D:\PNM\PNM_metrop_2024_12_20_orthomosaic.tif")

buffer=5
b_box= box(crownmap.total_bounds[0]-buffer, crownmap.total_bounds[1]-buffer, crownmap.total_bounds[2]+buffer, crownmap.total_bounds[3]+buffer)

def crop_raster(raster_path, output_path, bounds):
    '''
    Crops a raster to the specified bounding box and saves the output.
    Takes in shapely box bounds (minx, miny, maxx, maxy) and the path to the raster and output file.
    
    '''
    with rasterio.open(raster_path) as src:
        out_image, out_transform = mask.mask(src, [bounds], crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)

def resample_raster(input_path, output_path, target_resolution=0.05):
    '''
    Resamples a raster to a target resolution and saves the output.
    Takes in the path to the input raster, the path to the output raster, and the target resolution in meters (default: 0.05 m).
    '''
    with rasterio.open(input_path) as src:
        # Get current resolution
        current_resolution = abs(src.transform[0])  # pixel width in meters
        
        # Calculate scale factor
        scale_factor = current_resolution / target_resolution
        
        # Calculate new dimensions
        new_width = int(src.width * scale_factor)
        new_height = int(src.height * scale_factor)

        # Update metadata for the resampled raster
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": new_height,
            "width": new_width,
            "transform": src.transform * src.transform.scale(
                (src.width / new_width),
                (src.height / new_height)
            )
        })

        # Read and resample the data
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=rasterio.enums.Resampling.bilinear
        )

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(data)