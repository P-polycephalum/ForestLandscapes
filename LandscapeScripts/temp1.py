import os
import gc
import cv2
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.geometry import box as box1
from shapely import ops as shapely_ops



def generate_leafing_pdf(geodataframe, output_pdf, orthomosaic_path, crowns_per_page=12, variables=[]):
    """
    Generates a PDF with deciduous crowns plotted.

    Parameters:
        geodataframe (GeoDataFrame): DataFrame containing crown geometries and metadata.
        output_pdf (str): Path to the output PDF file.
        orthomosaic_path (str): Path to the orthomosaic folder containing image files.
        crowns_per_page (int): Number of crowns to plot per page (default: 12).
        variables(tupple): must be numeric variables. 
    """
    crowns_plotted = 0

    with PdfPages(output_pdf) as pdf_pages:
        fig, axes = plt.subplots(4, 3, figsize=(15, 20))
        axes = axes.flatten()

        for i, (_, row) in enumerate(geodataframe.iterrows()):
            date_target = row['date']
            path_orthomosaic = [
                os.path.join(orthomosaic_path, file)
                for file in os.listdir(orthomosaic_path)
                if date_target in file and file.endswith(".tif")
            ]
            print(path_orthomosaic)
            try:
                with rasterio.open(path_orthomosaic[0]) as src:
                    bounds = row.geometry.bounds
                    box_crown_5 = box1(bounds[0] - 5, bounds[1] - 5, bounds[2] + 5, bounds[3] + 5)

                    out_image, out_transform = mask(src, [box_crown_5], crop=True)
                    x_min, y_min = out_transform * (0, 0)
                    xres, yres = out_transform[0], out_transform[4]

                    # Transform geometry
                    transformed_geom = shapely_ops.transform(
                        lambda x, y: ((x - x_min) / xres, (y - y_min) / yres),
                        row.geometry
                    )

                    ax = axes[crowns_plotted % crowns_per_page]
                    ax.imshow(out_image.transpose((1, 2, 0))[:, :, 0:3])
                    ax.plot(*transformed_geom.exterior.xy, color='red', linewidth=2)
                    ax.axis('off')

                    # Add text label
                    annotation_text = f"{row['latin']}\n"
                    for var in variables:
                        if var in row:
                            try:
                                val = float(row[var])
                                annotation_text += f"{var}: {val:.2f}\n"
                            except (ValueError, TypeError):
                                annotation_text += f"{var}: {row[var]}\n"
                    # Add text label
                    ax.text(5, 5, annotation_text.strip(),
                            fontsize=12, color='white', backgroundcolor='black', verticalalignment='top')
                    crowns_plotted += 1

            except Exception as e:
                print(f"Error processing {path_orthomosaic}: {e}")
                continue  # Skip the current iteration if an error occurs

            # Save PDF and start a new page every `crowns_per_page` crowns
            if crowns_plotted % crowns_per_page == 0 or i == len(geodataframe) - 1:
                plt.tight_layout()
                pdf_pages.savefig(fig)
                plt.close(fig)

                # Create new figure for the next batch
                if i != len(geodataframe) - 1:  # Prevent unnecessary re-creation at end
                    fig, axes = plt.subplots(4, 3, figsize=(15, 20))
                    axes = axes.flatten()
    print(f"PDF saved: {output_pdf}")


timeseries_path=r"D:\BCI_50ha_timeseries\crownmap\BCI_50ha_crownmap_timeseries.gpkg"
gdf = gpd.read_file(timeseries_path)

#lets filter the latin= "Pseudobombax septenatum"

orthomosaic_folder = r"E:\backups\Product_local"

pseudoseptenatum_gdf = gdf[gdf['latin'] == "Pseudobombax septenatum"]
pseudoseptenatum_gdf['date']

for global_id in pseudoseptenatum_gdf['global_id'].unique():
    print(f"Processing GlobalID: {global_id}")
    gdf_tile = pseudoseptenatum_gdf[pseudoseptenatum_gdf['global_id'] == global_id].copy()

    generate_leafing_pdf(
        geodataframe=gdf_tile,
        output_pdf=f"{global_id}_leafing.pdf",
        orthomosaic_path=orthomosaic_folder,
        crowns_per_page=12,
        variables=['date'])
    print(f"Finished processing GlobalID: {global_id}\n")