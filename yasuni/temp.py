import os
from rasterio.mask import mask
from shapely.geometry import box
import rasterio
from rasterio.warp import transform_geom

path_shape= r"D:\BCI_ava_timeseries\cropped\BCI_ava_2020_10_26_orthomosaic.tif"

# Get bounds and CRS from reference raster
with rasterio.open(path_shape) as src:
    bounds = src.bounds
    ref_crs = src.crs
    print(f"Reference bounds: {bounds}")
    print(f"Reference CRS: {ref_crs}")

ha= r"D:\BCI_50ha_timeseries\BCI_50ha_2024_08_26_orthomosaic.tif"
with rasterio.open(ha) as src:
    print(f"Target CRS: {src.crs}")
    print(f"Target bounds: {src.bounds}")
    print(f"Number of bands: {src.count}")
    
    # Create box geometry in reference CRS
    box_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    
    # Transform geometry to target raster's CRS if needed
    if ref_crs != src.crs:
        print("CRS mismatch detected - reprojecting bounds")
        box_geom = transform_geom(ref_crs, src.crs, box_geom)
    
    # Mask with filled=False to preserve nodata
    out_image, out_transform = mask(src, [box_geom], crop=True, filled=False)
    
    # Select only RGB bands (0, 1, 2) and alpha band (7)
    selected_bands = out_image[[0, 1, 2, 7], :, :]
    
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": selected_bands.shape[1],
        "width": selected_bands.shape[2],
        "count": 4,
        "transform": out_transform
    })
    
    out_path= os.path.join(r"D:\BCI_ava_timeseries\cropped", "BCI_ava_2024_08_26_orthomosaic.tif")
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(selected_bands)