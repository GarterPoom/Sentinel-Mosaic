import geopandas as gpd # working with geospatial data in vector format.
import rasterio # reading and writing raster data.
from rasterio.mask import mask  # mask or clip raster data using vector geometries.
from rasterio.enums import Resampling # resampling methods for raster data.
from rasterio.features import shapes # extracting shapes from raster data.
from rasterio.shutil import copy as rio_copy # copying raster files with metadata and overviews.
import os # interacting with the operating system, such as file paths and directories.
import tempfile # creating temporary files.
import glob # finding files matching a specified pattern.
import logging # logging messages for debugging and information.
import numpy as np # numerical operations, especially with arrays.
from shapely.geometry import shape, box # working with geometric objects.
from fiona.crs import from_epsg # handling coordinate reference systems.

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s') # Set up logging configuration

def polygonize_clipped_raster(clipped_raster_path, output_polygon_path, nodata_fill_value):
    """
    Polygonizes a clipped raster, creating polygons for valid data areas (pixels set to 1).

    Args:
        clipped_raster_path (str): Path to the clipped raster file.
        output_polygon_path (str): Path to save the output polygon GeoJSON file.
        nodata_fill_value (float or int): The value in the clipped raster that represents nodata
                                          (i.e., the fill value used by rasterio.mask).
        # Note: This function assumes the input raster is already clipped to the desired area.
    """
    try:
        with rasterio.open(clipped_raster_path) as src:
            image = src.read(1)  # Read the first band
            transform = src.transform
            crs = src.crs
            
            # Create a binary mask: 1 for valid data, 0 for nodata_fill_value
            # Pixels equal to nodata_fill_value become 0, others (valid data) become 1.
            binary_mask = np.where(image == nodata_fill_value, 0, 1).astype(np.uint8)
            
            # Extract shapes (polygons) from the binary_mask where pixel value is 1
            # The mask=(binary_mask == 1) ensures we only process regions of 1s.
            results = [
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v) in enumerate(
                    shapes(binary_mask, mask=(binary_mask == 1), transform=transform)
                ) if v == 1 # Filter for polygons derived from pixels that were set to 1
            ]

            if not results:
                logging.info(f"No valid data (value 1) polygons found in {os.path.basename(clipped_raster_path)} to polygonize.")
                return

            geometries = [shape(result['geometry']) for result in results]
            gdf_polygons = gpd.GeoDataFrame(geometry=geometries, crs=crs)
            gdf_polygons['value'] = 1 # Add a column indicating the (binary) value of the polygonized area

            os.makedirs(os.path.dirname(output_polygon_path), exist_ok=True)
            gdf_polygons.to_file(output_polygon_path, driver='ESRI ShapeFile') # Save the polygons to a shapefile
            logging.info(f"Polygonized raster saved to {output_polygon_path}")
    except Exception as e:
        logging.error(f"Error polygonizing {os.path.basename(clipped_raster_path)}: {e}")

def clip_rasters_by_subdistrict_hierarchy(
    folder_path, shapefile_path,
    province_attr, district_attr, subdistrict_attr,
    output_dir, polygon_dir, 
    clip_level_attr: str, # New parameter to specify the clipping level attribute
    use_district_for_path: bool = True,
    use_subdistrict_for_path: bool = True):

    """
    Clip rasters by subdistrict hierarchy using a shapefile and save the clipped
    rasters to a directory. Also, polygonize the clipped rasters and save the
    resulting polygons to a separate directory.

    Args:
        folder_path (str): The folder path containing the raster files to be clipped.
        shapefile_path (str): The path to the shapefile containing the subdistrict boundaries.
        province_attr (str): The attribute name for the province in the shapefile.
        district_attr (str): The attribute name for the district in the shapefile.
        subdistrict_attr (str): The attribute name for the subdistrict in the shapefile.
        output_dir (str): The directory where the clipped rasters will be saved.
        clip_level_attr (str): The attribute name ('PV_TN', 'AP_TN', or 'TB_TN') to use for the primary clipping unit.
        polygon_dir (str): The directory where the polygonized clipped rasters will be saved.
        use_district_for_path (bool): Whether to include district in the output path hierarchy.
        use_subdistrict_for_path (bool): Whether to include subdistrict in the output path hierarchy.

    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True) # Ensure the output directory exists
    os.makedirs(polygon_dir, exist_ok=True) # Ensure the polygon directory exists

    gdf = gpd.read_file(shapefile_path) # Read the shapefile into a GeoDataFrame

    for attr in [province_attr, district_attr, subdistrict_attr]: # Check if the required attributes are present in the GeoDataFrame
        if attr not in gdf.columns: # Check if the required attributes are present in the GeoDataFrame
            logging.error(f"Attribute '{attr}' not found in shapefile.")
            return

    if clip_level_attr not in [province_attr, district_attr, subdistrict_attr]: # Validate clip_level_attr
        logging.error(f"Invalid clip_level_attr '{clip_level_attr}'. Must be one of '{province_attr}', '{district_attr}', or '{subdistrict_attr}'.")
        return

    # The following line is kept for informational logging but not directly used for iteration if clip_level_attr is different
    unique_subdistricts_for_info = gdf[subdistrict_attr].dropna().unique() 
    logging.info(f"Found {len(unique_subdistricts_for_info)} unique subdistricts in attribute '{subdistrict_attr}' (for informational purposes).")

    raster_files = glob.glob(os.path.join(folder_path, '*.tif')) # Find all .tif files in the specified folder path
    if not raster_files:
        logging.warning("No .tif files found in the folder.")
        return

    logging.info(f"Found {len(raster_files)} raster files. Starting batch clipping...")

    unique_clip_units = gdf[clip_level_attr].dropna().unique() # Get unique values for the specified clipping level
    logging.info(f"Clipping by attribute '{clip_level_attr}'. Found {len(unique_clip_units)} unique units.")

    clipped_raster_info_list = [] # To store info for polygonization
    
    for raster_path in raster_files: # Iterate through each raster file found
        raster_name = os.path.basename(raster_path)
        raster_basename_no_ext = os.path.splitext(raster_name)[0]
        try:
            with rasterio.open(raster_path) as src: # Open the raster file using rasterio
                raster_crs = src.crs # Get the coordinate reference system of the raster
                raster_bounds = src.bounds # Get the bounds of the raster
                original_src_nodata = src.nodata # Get nodata value from original raster
                raster_bbox = box(*raster_bounds) # Create a bounding box from the raster bounds
            
            nodata_fill_value_for_clipped = original_src_nodata if original_src_nodata is not None else 0
            gdf_proj = gdf.to_crs(raster_crs) # Reproject the GeoDataFrame to the raster's CRS

            for current_clip_unit_name in unique_clip_units: # Iterate through each unique unit at the specified clipping level
                # Filter the GeoDataFrame for the current clipping unit
                filtered_gdf = gdf_proj[gdf_proj[clip_level_attr] == current_clip_unit_name] 
                if filtered_gdf.empty: 
                    continue
                if not filtered_gdf.intersects(raster_bbox).any():
                    continue

                # Extract the province, district, and subdistrict names associated with this clip unit
                # Assume the first row of the filtered GDF is representative of the hierarchy
                # Ensure there's at least one row before trying to access .iloc[0]
                if not filtered_gdf.empty:
                    province = str(filtered_gdf[province_attr].iloc[0])
                    district = str(filtered_gdf[district_attr].iloc[0])
                    subdistrict = str(filtered_gdf[subdistrict_attr].iloc[0])
                else: # Should not happen if filtered_gdf.empty check passed, but as a safeguard
                    province = "UnknownProvince"
                    district = "UnknownDistrict"
                    subdistrict = "UnknownSubdistrict"
                
                # --- Path and Filename Construction ---
                path_parts = [output_dir]
                poly_path_parts = [polygon_dir]

                # Sanitize names for use in filenames and potentially paths if needed
                province_clean_fn = province.replace(" ", "_").replace(os.sep, "_")
                district_clean_fn = district.replace(" ", "_").replace(os.sep, "_") 
                subdistrict_clean_fn = subdistrict.replace(" ", "_").replace(os.sep, "_") 
                clip_unit_clean_fn = str(current_clip_unit_name).replace(" ", "_").replace(os.sep, "_") # Sanitize the current clip unit name

                path_parts.append(province_clean_fn) # Province is always in path (use cleaned name for path)
                poly_path_parts.append(province_clean_fn)

                # Determine which parts of the hierarchy to include in the path
                current_path_level_name_for_district = district_clean_fn
                current_path_level_name_for_subdistrict = subdistrict_clean_fn

                if use_district_for_path:
                    path_parts.append(current_path_level_name_for_district)
                    poly_path_parts.append(current_path_level_name_for_district)
                
                if use_subdistrict_for_path:
                    # Only add subdistrict to path if district is also in path, or if it's the primary clip level
                    if use_district_for_path or clip_level_attr == subdistrict_attr:
                         path_parts.append(current_path_level_name_for_subdistrict)
                         poly_path_parts.append(current_path_level_name_for_subdistrict)


                clipped_folder = os.path.join(*path_parts)
                polygon_output_folder = os.path.join(*poly_path_parts)

                geom_json = [geom.__geo_interface__ for geom in filtered_gdf.geometry] # Convert geometries to GeoJSON format

                # Filename includes the clip unit name and the raster name for uniqueness
                unique_file_identifier = f"{clip_unit_clean_fn}_{raster_basename_no_ext}"
                
                with rasterio.Env(GDAL_CACHEMAX=512): # Set GDAL cache size to 512 MB for better performance
                    with rasterio.open(raster_path) as src:
                        out_image, out_transform = mask(src, geom_json, crop=True, nodata=nodata_fill_value_for_clipped, filled=True)
                        out_meta = src.meta.copy()

                out_meta.update({ # Update metadata for the output raster
                    "driver": "GTiff", # Output format
                    "height": out_image.shape[1], # Height of the output raster
                    "width": out_image.shape[2], # Width of the output raster
                    "transform": out_transform, # Affine transform for the output raster
                    "nodata": nodata_fill_value_for_clipped, # Ensure nodata is set in meta
                    "tiled": True, # Enable tiling for the output raster
                    "blockxsize": 256, # Block size in x direction
                    "blockysize": 256, # Block size in y direction
                    "compress": "lzw" # Compression method for the output raster
                })

                os.makedirs(clipped_folder, exist_ok=True) # Ensure the clipped folder exists
                os.makedirs(polygon_output_folder, exist_ok=True) # Ensure polygon folder exists

                clipped_filename = f"clipped_{unique_file_identifier}.tif" 
                clipped_path = os.path.join(clipped_folder, clipped_filename)
                
                with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp: 
                    temp_output = tmp.name

                with rasterio.open(temp_output, "w", **out_meta) as dest: 
                    dest.write(out_image)

                with rasterio.open(temp_output, 'r+') as dataset: 
                    dataset.build_overviews([2, 4, 8, 16, 32], Resampling.nearest)
                    dataset.update_tags(ns='rio_overview', resampling='nearest')

                rio_copy(temp_output, clipped_path, copy_src_overviews=True) 
                os.remove(temp_output) 

                logging.info(f"Saved clipped raster: {clipped_path}")

                polygon_filename = f"polygon_{unique_file_identifier}.shp" 
                polygon_final_path = os.path.join(polygon_output_folder, polygon_filename)

                clipped_raster_info_list.append({
                    'clipped_path': clipped_path,
                    'nodata_fill_value': nodata_fill_value_for_clipped,
                    'polygon_final_path': polygon_final_path
                })

        except Exception as e: 
            logging.error(f"Error processing {raster_name} for {clip_level_attr} '{current_clip_unit_name if 'current_clip_unit_name' in locals() else 'N/A'}': {e}")
    
    logging.info("Batch clipping finished.")

    if not clipped_raster_info_list:
        logging.info("No rasters were successfully clipped, skipping polygonization.")
    else:
        logging.info(f"Starting polygonization for {len(clipped_raster_info_list)} clipped rasters...")
        for info in clipped_raster_info_list:
            clipped_path = info['clipped_path']
            nodata_value = info['nodata_fill_value']
            polygon_final_path_item = info['polygon_final_path']
            polygonize_clipped_raster(clipped_path, polygon_final_path_item, nodata_value)
        logging.info("Batch polygonizing finished.")

    logging.info("Overall batch clipping and polygonizing process complete.")

def main():
    clip_rasters_by_subdistrict_hierarchy(
        folder_path="LANDSAT_9/", # Path to the folder containing raster files
        shapefile_path="Thailand/Thailand - Subnational Administrative Boundaries.shp", # Path to the shapefile with subdistrict boundaries
        province_attr="PV_TN", # Attribute for province in the shapefile
        district_attr="AP_TN", # Attribute for district in the shapefile
        subdistrict_attr="TB_TN", # Attribute for subdistrict in the shapefile
        output_dir="Clipped_Rasters", # Directory to save clipped rasters
        polygon_dir="Polygonized_Rasters", # Directory to save polygonized rasters
        
        clip_level_attr="PV_TN",  # Attribute to use for clipping level (can be 'PV_TN', 'AP_TN', or 'TB_TN')
        
        use_district_for_path=False, # Whether to include district in the output path hierarchy
        use_subdistrict_for_path=False # Whether to include subdistrict in the output path hierarchy
    )

if __name__ == "__main__":
    main()
