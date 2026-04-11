
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.enums import Resampling as ResamplingEnum
import numpy as np

def resample_to_resolution_uint8(input_path, output_path, target_resolution, resampling_method=ResamplingEnum.bilinear):
    """
    Resample a raster to a target resolution and convert to uint8.
    
    Parameters:
    -----------
    input_path : str
        Path to input raster file
    output_path : str
        Path to output raster file
    target_resolution : float
        Target resolution in the same units as the raster CRS (e.g., meters)
    resampling_method : rasterio.enums.Resampling
        Resampling method (default: bilinear)
    """
    with rasterio.open(input_path) as src:
        transform = src.transform
        current_res = (transform.a, -transform.e)
        
        # Calculate scaling factors
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
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=resampling_method
        )
        
        # Convert to uint8 if not already
        if src.dtypes[0] != 'uint8':
            # Handle nodata values
            src_nodata = src.nodata
            for band in range(data.shape[0]):
                # Check if this is the alpha band (last band for 4-band images)
                if band == data.shape[0] - 1 and data.shape[0] == 4:
                    # Set alpha to fully opaque
                    data[band] = np.full_like(data[band], 255, dtype='uint8')
                else:
                    band_data = data[band]
                    # Mask nodata values
                    if src_nodata is not None:
                        valid_mask = band_data != src_nodata
                        if valid_mask.any():
                            band_min, band_max = np.percentile(band_data[valid_mask], (1, 99))
                            band_data_clipped = np.clip(band_data, band_min, band_max)
                            converted = ((band_data_clipped - band_min) / (band_max - band_min) * 255).astype('uint8')
                            # Set nodata pixels to 0
                            converted[~valid_mask] = 0
                            data[band] = converted
                        else:
                            data[band] = np.zeros_like(band_data, dtype='uint8')
                    else:
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
            "dtype": 'uint8',
        })
        # Remove nodata if converting to uint8 to avoid conflicts with valid data
        if 'nodata' in out_meta:
            del out_meta['nodata']
        
        # Write output
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(data)

def crop_raster(input_path, output_path, shapely_polygon):
    with rasterio.open(input_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, [shapely_polygon], crop=True)
        out_meta = src.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image)

def combine_ortho_dsm(ortho_path,dsm_path, output_path):
    with rasterio.open(ortho_path) as src:
        ortho_data = src.read()
        ortho_meta = src.meta.copy()
    with rasterio.open(dsm_path) as src:
        dem_data = src.read(1)
        dem_meta = src.meta
        dem_data=np.where(dem_data==dem_meta['nodata'],0,dem_data)
    resampled_dem = np.zeros((ortho_meta['height'], ortho_meta['width']), dtype=ortho_data.dtype)
    reproject(
    dem_data, resampled_dem,
    src_transform=dem_meta['transform'],
    src_crs=dem_meta['crs'],
    dst_transform=ortho_meta['transform'],
    dst_crs=ortho_meta['crs'],
    resampling=Resampling.nearest)
    ortho_data[4,:,:] = resampled_dem
    ortho_meta['count'] = 5
    with rasterio.open(output_path, 'w', **ortho_meta) as dst:
        dst.write(ortho_data)