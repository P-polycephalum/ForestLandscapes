import os
import gc
import cv2
import numpy as np
import pandas as pd
import geopandas as gpd
import torch
import rasterio
from rasterio.mask import mask
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, box
from shapely import ops as shapely_ops
from sam2.build_sam import build_sam2
from sam2.sam2_image_predictor import SAM2ImagePredictor

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
            date_target= row['date'].strftime("%Y_%m_%d")
            path_orthomosaic = [os.path.join(orthomosaic_path, file) for file in os.listdir(orthomosaic_path) if date_target in file and file.endswith(".tif")]
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

def multi_to_polygon(geom):
    """Resolve MultiPolygon geometries to their largest Polygon.
    Takes a geometry object, checks if it's a MultiPolygon, and if so, extracts the largest Polygon based on area.
    If it's not a MultiPolygon, returns the original geometry.
    """
    if isinstance(geom, MultiPolygon):
        multi_polygon = geom
        polygons = []
        for polygon in multi_polygon.geoms:
            polygons.append(polygon)
        largest_polygon = max(polygons, key=lambda polygon: polygon.area)
        return largest_polygon
    else:
        print("Geometry is not a MultiPolygon. Returning original geometry.")
        return geom

def tile_ortho(orthomosaic, tile_size, buffer, output_folder):
    """
    Tiles an orthomosaic image into smaller raster tiles with optional buffer.
    
    :param orthomosaic: Path to the input orthomosaic raster file.
    :param tile_size: Size of each tile in pixels.
    :param buffer: Buffer to add around each tile.
    :param output_folder: Folder to save the output tiles.
    """
    with rasterio.open(orthomosaic) as src:
        bounds = src.bounds
        xmin, ymin, xmax, ymax = bounds
        if tile_size <= 0:
            raise ValueError("tile_size must be greater than zero.")      
        x_range = xmax - xmin
        y_range = ymax - ymin
        x_tiles = int(np.ceil(x_range / tile_size))
        y_tiles = int(np.ceil(y_range / tile_size))
        x_residual = x_range % tile_size
        y_residual = y_range % tile_size
        if x_residual > 0:
            tile_size_x = tile_size + x_residual / x_tiles
        else:
            tile_size_x = tile_size
        if y_residual > 0:
            tile_size_y = tile_size + y_residual / y_tiles
        else:
            tile_size_y = tile_size
        if x_residual > 0 or y_residual > 0:
            print(f"Warning: Adjusted tile size used for residual coverage - X: {tile_size_x}, Y: {tile_size_y}")
        xmins = np.arange(xmin, (xmax - tile_size_x + 1), tile_size_x)
        xmaxs = np.arange((xmin + tile_size_x), xmax + 1, tile_size_x)
        ymins = np.arange(ymin, (ymax - tile_size_y + 1), tile_size_y)
        ymaxs = np.arange((ymin + tile_size_y), ymax + 1, tile_size_y)
        X, Y = np.meshgrid(xmins, ymins)
        Xmax, Ymax = np.meshgrid(xmaxs, ymaxs)
        gridInfo = pd.DataFrame({
            'xmin': X.flatten(),
            'ymin': Y.flatten(),
            'xmax': Xmax.flatten(),
            'ymax': Ymax.flatten(),
        })
        print(gridInfo)
    with rasterio.open(orthomosaic) as src:
        for idx, row in gridInfo.iterrows():
            geom2 = box(row['xmin']-buffer, row['ymin']-buffer, row['xmax']+buffer, row['ymax']+buffer)
            out_image, out_transform = rasterio.mask.mask(src, [geom2], crop=True)
            # Update metadata for the output raster
            out_meta = src.meta
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })
            output_filename = f"output_raster_{idx}.tif"
            filename=os.path.join(output_folder,output_filename)
            with rasterio.open(filename, "w", **out_meta) as dest:
                dest.write(out_image)

def crown_segment(tile_folder,shp,checkpoint, model_cfg, device):
    
    sam2_model = build_sam2(model_cfg, checkpoint, device=device)
    predic = SAM2ImagePredictor(sam2_model)
    all=[]
    tiles= os.listdir(tile_folder)
    for tile in tiles:
        print("processing tile", tile)
        sub=os.path.join(tile_folder,tile)
        with rasterio.open(sub) as src:
            data=src.read()
            transposed_data=data.transpose(1,2,0)
            crs=src.crs
            affine_transform = src.transform 
            bounds=src.bounds
            main_box= box(bounds[0],bounds[1],bounds[2],bounds[3])
        crowns=shp
        crowns= crowns.to_crs(crs)
        mask = crowns['geometry'].within(main_box)
        test_crowns = crowns.loc[mask]

        print("starting box transformation from utm to xy")
        boxes=[]
        for index, row in test_crowns.iterrows():
            if isinstance(row.geometry, MultiPolygon):
                multi_polygon = row.geometry
                polygons = []
                for polygon in multi_polygon.geoms:
                    polygons.append(polygon)
                largest_polygon = max(polygons, key=lambda polygon: polygon.area)
                bounds = largest_polygon.bounds
                boxes.append(bounds)
            else:
                bounds = row.geometry.bounds
                boxes.append(bounds)
        box_mod=[]
        for box in boxes:
            xmin, ymin, xmax, ymax = box
            x_pixel_min, y_pixel_min = ~affine_transform * (xmin, ymin)
            x_pixel_max, y_pixel_max = ~affine_transform * (xmax, ymax)
            trans_box=[x_pixel_min,y_pixel_max,x_pixel_max,y_pixel_min]
            box_mod.append(trans_box)
        if not box_mod:
            print(f"No valid boxes found for tile {tile}. Skipping.")
            continue
        print("The tile contains", len(box_mod), "polygons")

        input_boxes=torch.tensor(box_mod, device=device)
        print("about to set the image")
        predic.set_image(transposed_data[:,:,:3])
        with torch.inference_mode():
            masks, scores, logits = predic.predict(
                box=input_boxes,
                multimask_output=True,
            )
        predic = SAM2ImagePredictor(sam2_model)
        print("finish predicting now getting the utms for transformation")
        height, width, num_bands = transposed_data.shape
        utm_coordinates_and_values = np.empty((height, width, num_bands + 2))
        utm_transform = src.transform
       
        y_coords, x_coords = np.meshgrid(np.arange(height), np.arange(width), indexing='ij')
        utm_x, utm_y = rasterio.transform.xy(utm_transform, y_coords, x_coords)
        utm_coordinates_and_values[..., 0] = np.array(utm_x).reshape(height, width)
        utm_coordinates_and_values[..., 1] = np.array(utm_y).reshape(height, width)
        utm_coordinates_and_values[..., 2:] = transposed_data[..., :num_bands]

        all_polygons=[]
        for idx, (thisscore, thiscrown) in enumerate(zip(scores, masks)):
            maxidx=thisscore.tolist().index(max(thisscore.tolist()))
            thiscrown = thiscrown[maxidx]
            score=scores[1].tolist()[thisscore.tolist().index(max(thisscore.tolist()))]   
            mask = thiscrown.squeeze()
            utm_coordinates = utm_coordinates_and_values[:, :, :2]
            mask_np = mask.astype(np.uint8) if isinstance(mask, np.ndarray) else mask.cpu().numpy().astype(np.uint8)
            contours, _ = cv2.findContours(mask_np, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            polygons = []
            areas = []
            for contour in contours:
                contour_coords = contour.squeeze().reshape(-1, 2)
                contour_utm_coords = utm_coordinates[contour_coords[:, 1], contour_coords[:, 0]]
                if len(contour_utm_coords) >= 3:
                    polygon = Polygon(contour_utm_coords)
                    area = polygon.area
                    polygons.append(polygon)
                    areas.append(area)       
            if len(areas) == 0:
                print(f"No valid areas found for this crown. Skipping.")
                continue  
            largest_index = np.argmax(areas)
            gdf = gpd.GeoDataFrame(geometry=[polygons[largest_index]])
            gdf['area'] = areas[largest_index]
            gdf['score'] = score 
            gdf.crs = src.crs
            tag_value = test_crowns.iloc[idx]['tag']
            global_id= test_crowns.iloc[idx]['global_id']
            gdf['tag']=tag_value
            gdf['global_id']=global_id
            all_polygons.append(gdf)
        print("finish transforming back to utm")
        print(len(all_polygons),"crowns segmented")
        all.append(all_polygons)
        progress= len(all)/len(tiles)
        del input_boxes, masks, scores, logits
        torch.cuda.empty_cache()
        print(progress)
    final_gdfs = []
    for polygons_gdf_list in all:
        combined_gdf = gpd.GeoDataFrame(pd.concat(polygons_gdf_list, ignore_index=True), crs=src.crs)
        final_gdfs.append(combined_gdf)
    final_gdf = gpd.GeoDataFrame(pd.concat(final_gdfs, ignore_index=True), crs=src.crs)
    return final_gdf

def crown_avoid(crown_gdf):
    """Resolve overlapping crown polygons and keep valid polygon geometries.

    Takes a GeoDataFrame of crown polygons, converts multipolygons to their
    largest part, subtracts overlaps based on relative area, and drops any
    non-polygon geometries that remain.
    """
    
    crown_avoidance = crown_gdf.copy()
    crown_avoidance['geometry'] = crown_avoidance.geometry.buffer(0)

    for index, row in crown_avoidance.iterrows():
        if isinstance(row.geometry, MultiPolygon):
            multi_polygon = row.geometry
            polygons = [polygon for polygon in multi_polygon.geoms]
            largest_polygon = max(polygons, key=lambda polygon: polygon.area)
            crown_avoidance.at[index, 'geometry'] = largest_polygon
            print(f"Converted MultiPolygon to largest Polygon for index {index}.")

    sindex = crown_avoidance.sindex
    modifications = {}  # Dictionary to collect modifications
    for idx, polygon in crown_avoidance.iterrows():
        possible_matches_index = list(sindex.intersection(polygon['geometry'].bounds))
        possible_matches = crown_avoidance.iloc[possible_matches_index]
        adjacents = possible_matches[possible_matches.geometry.intersects(polygon['geometry']) & (possible_matches.index != idx)]
        if adjacents.empty:
            continue
        else:
            for adj_idx, adj_polygon in adjacents.iterrows():
                if polygon['similarity'] > adj_polygon['similarity']:
                    # Adjacent loses overlap
                    modifications[adj_idx] = modifications.get(adj_idx, adj_polygon.geometry).difference(polygon.geometry)
                    print(f"Polygon {idx} is more similar than {adj_idx}. Subtracting overlap from adjacent.")   
                elif polygon['similarity'] < adj_polygon['similarity']:
                    # Current polygon loses overlap
                    modifications[idx] = modifications.get(idx, polygon.geometry).difference(adj_polygon.geometry)
                    print(f"Polygon {idx} is less similar than {adj_idx}. Subtracting overlap from polygon.")

    for idx, new_geom in modifications.items():
        crown_avoidance.at[idx, 'geometry'] = new_geom
    for index, row in crown_avoidance.iterrows():
        if isinstance(row.geometry, MultiPolygon):
            multi_polygon = row.geometry
            polygons = [polygon for polygon in multi_polygon.geoms]
            largest_polygon = max(polygons, key=lambda polygon: polygon.area)
            crown_avoidance.at[index, 'geometry'] = largest_polygon

    for index, row in crown_avoidance.iterrows():
        geom = row["geometry"]
        if isinstance(geom, GeometryCollection):
            polygons = [g for g in geom.geoms if isinstance(g, Polygon)]
            if polygons:
                crown_avoidance.at[index, "geometry"] = polygons[0]
            else:
                crown_avoidance.at[index, "geometry"] = pd.NA
        elif not isinstance(geom, Polygon):
            crown_avoidance.at[index, "geometry"] = pd.NA
    return crown_avoidance

def crownmap_metrics(original_crownmap, segmented_crownmap):
    segmented_ids = set(segmented_crownmap['global_id'].unique())
    original_ids = set(original_crownmap['global_id'].unique())
    missing_in_segmented = original_ids - segmented_ids
    extra_in_segmented = segmented_ids - original_ids
    print(f"Missing IDs in segmented: {missing_in_segmented if missing_in_segmented else 'None'}")
    print(f"Extra IDs in segmented: {extra_in_segmented if extra_in_segmented else 'None'}")

    original_polys = {row['global_id']: row['geometry'] for _, row in original_crownmap.iterrows()}
    segmented_polys = {row['global_id']: row['geometry'] for _, row in segmented_crownmap.iterrows()}

    common_ids = original_ids & segmented_ids
    for crown in common_ids:
        original_geom = original_polys.get(crown)
        segmented_geom = segmented_polys.get(crown)
        original_geom= original_geom.buffer(0)
        segmented_geom= segmented_geom.buffer(0)
        if original_geom is None or segmented_geom is None:
            print(f"Warning: Crown {crown} missing geometry in one of the inputs. Skipping.")
            continue

        intersection_area = original_geom.intersection(segmented_geom).area
        union_area = original_geom.union(segmented_geom).area
        precision = intersection_area / segmented_geom.area if segmented_geom.area > 0 else 0
        recall = intersection_area / original_geom.area if original_geom.area > 0 else 0
        iou = intersection_area / union_area if union_area > 0 else 0
        f1= 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        score = 0.6 * iou + (0.4 * f1)
        print(f"GlobalID: {crown}, IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}, F1: {f1:.2f}, Score: {score:.2f}")
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'IoU'] = iou
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'Precision'] = precision
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'Recall'] = recall
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'F1'] = f1
        segmented_crownmap.loc[segmented_crownmap['global_id'] == crown, 'similarity'] = score
    
    return segmented_crownmap

def crownmap_QC(original_crownmap, segmented_crownmap, precision_threshold=0.5, recall_threshold=0.5, iou_threshold=0.5):
        """Calcalate quality metrics for segmented crowns and keep original polygons if quality is below threshold.

        Takes a GeoDataFrame of crown polygons.
        """
        crown_qc= crownmap_metrics(original_crownmap, segmented_crownmap)
        for crown in crown_qc['global_id'].unique():
            original_geom = original_crownmap.loc[original_crownmap['global_id'] == crown, 'geometry'].values[0]
            segmented_geom = crown_qc.loc[crown_qc['global_id'] == crown, 'geometry'].values[0]
            iou = crown_qc.loc[crown_qc['global_id'] == crown, 'IoU'].values[0]
            precision = crown_qc.loc[crown_qc['global_id'] == crown, 'Precision'].values[0]
            recall = crown_qc.loc[crown_qc['global_id'] == crown, 'Recall'].values[0]
        # if the crown is below any threshold, keep the original polygon
            if precision < precision_threshold or recall < recall_threshold or iou < iou_threshold:
                print(
                    f"Warning: Crown {crown} has low quality (IoU: {iou:.2f}, Precision: {precision:.2f}, Recall: {recall:.2f}). "
                    "Keeping original polygon."
                )
                crown_qc.loc[segmented_crownmap['global_id'] == crown, 'geometry'] = original_geom

        return crown_qc

