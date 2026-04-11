import os
import geopandas as gpd
from shapely.geometry import Point, Polygon
import pandas as pd
import matplotlib.pyplot as plt
# Create a 10-meter grid
from shapely.geometry import box
import numpy as np



path= r"crown_area_predictions.csv"
data= pd.read_csv(path)
gdf= gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(data.long, data.lat), crs="EPSG:4326")
gdf.to_crs("EPSG:32617", inplace=True)

ax = gdf.plot(figsize=(10, 10), alpha=0.1, markersize=0.5)
ax.set_title('Stem positons')
plt.show()

gdf['geometry'] = gdf.apply(lambda row: Point(row.geometry.x, row.geometry.y).buffer(row.radius), axis=1)
                                                                                                     


bounds = gdf.total_bounds  # (minx, miny, maxx, maxy)
minx, miny, maxx, maxy = bounds
cellsize = 20

# Generate grid cells
grid_cells = []
grid_ids = []
cell_id = 0

for x in np.arange(minx, maxx, cellsize):
    for y in np.arange(miny, maxy, cellsize):
        grid_cells.append(box(x, y, x + cellsize, y + cellsize))
        grid_ids.append(cell_id)
        cell_id += 1

grid = gpd.GeoDataFrame({'cell_id': grid_ids, 'geometry': grid_cells}, crs="EPSG:32617")

# Intersect polygons with grid cells
intersections = gpd.sjoin(gdf[['geometry', 'agb_Mg', 'radius']], grid, how='inner', predicate='intersects')

# Calculate overlap area and allocate AGB by area fraction
intersections['overlap_geom'] = intersections.apply(
    lambda row: gdf.loc[row.name, 'geometry'].intersection(grid.loc[row.index_right, 'geometry']),
    axis=1
)
intersections['overlap_area'] = intersections['overlap_geom'].area
intersections['polygon_area'] = gdf.loc[intersections.index, 'geometry'].area.values
intersections['agb_fraction'] = intersections['overlap_area'] / intersections['polygon_area']
intersections['agb_cell_contrib'] = intersections['agb_Mg'] * intersections['agb_fraction']

# Sum AGB by grid cell
grid_agb = intersections.groupby('cell_id')['agb_cell_contrib'].sum().reset_index()
grid_agb.columns = ['cell_id', 'agb_Mg']

# Merge back to grid
grid = grid.merge(grid_agb, on='cell_id', how='left')
grid['agb_Mg'] = grid['agb_Mg'].fillna(0)

print(f"Total AGB in polygons: {gdf['agb_Mg'].sum():.2f} Mg")
print(f"Total AGB in grid: {grid['agb_Mg'].sum():.2f} Mg")

# Plot grid with AGB values
ax = grid.plot(column='agb_Mg', cmap='YlGn', figsize=(12, 10), alpha=0.8)
ax.set_title('AGB per 10m Grid Cell')
plt.colorbar(ax.collections[0], label='AGB (Mg)')
plt.show()


quadrats=r"C:\Users\vasquezvicente\repo\quad_poly.gpkg"
quad_gdf=gpd.read_file(quadrats)
quad_gdf.columns

quad=quad_gdf[quad_gdf['plot_id'] == 'BCI 50 ha plot'].copy()
quad.to_crs("EPSG:32617", inplace=True)


# so for every quad polygon we extract the sum of the AGB in the grid that intersects with it. we can also calculate the area of the quad and get the AGB per unit area.

for idx, row in quad.iterrows():
    quad_polygon = row.geometry
    intersecting_cells = grid[grid.intersects(quad_polygon)]
    total_agb = intersecting_cells['agb_Mg'].sum()
    quad_area = quad_polygon.area / 10000  # convert to hectares
    agb_per_ha = total_agb / quad_area if quad_area > 0 else 0
    print(f"Plot {row['plot_id']}: Total AGB = {total_agb:.2f} Mg, Area = {quad_area:.2f} ha, AGB/ha = {agb_per_ha:.2f} Mg/ha")
    quad.at[idx, 'total_agb'] = total_agb
    quad.at[idx, 'quad_area_ha'] = quad_area
    quad.at[idx, 'agb_per_ha'] = agb_per_ha


ax=quad.plot(column='agb_per_ha', cmap='YlOrRd', figsize=(12, 10), alpha=0.8)
ax.set_title('AGB per Hectare by Quadrant')
plt.colorbar(ax.collections[0], label='AGB per Hectare (Mg/ha)')
plt.show()

# 