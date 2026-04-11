import numpy as np
from shapely.geometry import box

def tile_crownmap(crownmap_gdf, tile_size, overlap=0):
    """
    Tile a crownmap geodataframe with adjusted tile sizes to distribute residuals equally.
    
    Parameters:
    -----------
    crownmap_gdf : GeoDataFrame
        Input geodataframe to tile
    tile_size : float
        Target tile size (will be adjusted to fit bounds exactly)
    overlap : float
        Overlap between adjacent tiles (default: 0)
    
    Returns:
    --------
    list of shapely.geometry.box objects representing tile bounds
    """
    minx, miny, maxx, maxy = crownmap_gdf.total_bounds
    
    if tile_size <= 0:
        raise ValueError("tile_size must be greater than zero.")
    
    # Calculate ranges
    x_range = maxx - minx
    y_range = maxy - miny
    
    # Calculate number of tiles needed
    x_tiles = int(np.ceil(x_range / tile_size))
    y_tiles = int(np.ceil(y_range / tile_size))
    
    # Calculate residuals
    x_residual = x_range % tile_size
    y_residual = y_range % tile_size
    
    # Adjust tile size to distribute residuals equally
    if x_residual > 0:
        tile_size_x = tile_size + x_residual / x_tiles
    else:
        tile_size_x = tile_size
    
    if y_residual > 0:
        tile_size_y = tile_size + y_residual / y_tiles
    else:
        tile_size_y = tile_size
    
    if x_residual > 0 or y_residual > 0:
        print(f"Adjusted tile size - X: {tile_size_x:.2f}m, Y: {tile_size_y:.2f}m ({x_tiles}x{y_tiles} tiles)")
    
    # Generate tile boundaries with overlap
    effective_step_x = tile_size_x - overlap
    effective_step_y = tile_size_y - overlap
    
    tiles = []
    for x_idx in range(x_tiles):
        for y_idx in range(y_tiles):
            x_start = minx + x_idx * effective_step_x
            y_start = miny + y_idx * effective_step_y
            x_end = min(x_start + tile_size_x, maxx)
            y_end = min(y_start + tile_size_y, maxy)
            tile_bounds = box(x_start, y_start, x_end, y_end)
            tiles.append(tile_bounds)
    
    return tiles
