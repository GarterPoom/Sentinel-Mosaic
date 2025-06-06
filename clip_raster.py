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

def clip_rasters_by_subdistrict_hierarchy(
    folder_path, shapefile_path,
    province_attr, district_attr, subdistrict_attr,
    output_dir, polygon_dir):

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
        polygon_dir (str): The directory where the polygonized clipped rasters will be saved.

    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True) # Ensure the output directory exists
    os.makedirs(polygon_dir, exist_ok=True) # Ensure the polygon directory exists

    gdf = gpd.read_file(shapefile_path) # Read the shapefile into a GeoDataFrame

    for attr in [province_attr, district_attr, subdistrict_attr]: # Check if the required attributes are present in the GeoDataFrame
        if attr not in gdf.columns:
            logging.error(f"Attribute '{attr}' not found in shapefile.")
            return

    unique_subdistricts = gdf[subdistrict_attr].dropna().unique() # Get unique subdistricts from the specified attribute
    logging.info(f"Found {len(unique_subdistricts)} unique subdistricts in attribute '{subdistrict_attr}'.") # Log the number of unique subdistricts found

    raster_files = glob.glob(os.path.join(folder_path, '*.tif')) # Find all .tif files in the specified folder path
    if not raster_files:
        logging.warning("No .tif files found in the folder.")
        return

    logging.info(f"Found {len(raster_files)} raster files. Starting batch clipping...")

    clipped_raster_info_list = [] # To store info for polygonization

    for raster_path in raster_files: # Iterate through each raster file found
        raster_name = os.path.basename(raster_path)
        try:
            with rasterio.open(raster_path) as src: # Open the raster file using rasterio
                raster_crs = src.crs # Get the coordinate reference system of the raster
                raster_bounds = src.bounds # Get the bounds of the raster
                original_src_nodata = src.nodata # Get nodata value from original raster
                raster_bbox = box(*raster_bounds) # Create a bounding box from the raster bounds
            
            nodata_fill_value_for_clipped = original_src_nodata if original_src_nodata is not None else 0
            gdf_proj = gdf.to_crs(raster_crs) # Reproject the GeoDataFrame to the raster's CRS

            for subdistrict in unique_subdistricts: # Iterate through each unique subdistrict
                filtered_gdf = gdf_proj[gdf_proj[subdistrict_attr] == subdistrict] # Filter the GeoDataFrame for the current subdistrict
                if filtered_gdf.empty: 
                    continue
                if not filtered_gdf.intersects(raster_bbox).any():
                    continue

                provinces = filtered_gdf[province_attr].unique() # Get unique provinces in the filtered GeoDataFrame
                districts = filtered_gdf[district_attr].unique() # Get unique districts in the filtered GeoDataFrame
                if len(provinces) != 1 or len(districts) != 1: # Ensure only one province and one district per subdistrict
                    continue

                province = str(provinces[0]) # Get the province name
                district = str(districts[0]) # Get the district name

                geom_json = [geom.__geo_interface__ for geom in filtered_gdf.geometry] # Convert geometries to GeoJSON format

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

                clipped_folder = os.path.join(output_dir, province, district, subdistrict) # Create output directory structure for clipped rasters
                os.makedirs(clipped_folder, exist_ok=True) # Ensure the clipped folder exists
                clipped_path = os.path.join(clipped_folder, f"clipped_{raster_name}") # Save path for the clipped raster

                with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp: # Create a temporary file for the clipped raster
                    temp_output = tmp.name

                with rasterio.open(temp_output, "w", **out_meta) as dest: # Open the temporary file for writing
                    dest.write(out_image)

                with rasterio.open(temp_output, 'r+') as dataset: # Open the temporary raster file for reading and writing
                    dataset.build_overviews([2, 4, 8, 16, 32], Resampling.nearest)
                    dataset.update_tags(ns='rio_overview', resampling='nearest')

                rio_copy(temp_output, clipped_path, copy_src_overviews=True) # Copy the temporary raster to the final output path with overviews
                os.remove(temp_output) # Remove the temporary file after copying

                logging.info(f"Saved clipped raster: {clipped_path}")

                # Store info for later polygonization
                clipped_raster_info_list.append({
                    'clipped_path': clipped_path,
                    'nodata_fill_value': nodata_fill_value_for_clipped,
                    'output_structure': (province, district, subdistrict) # To recreate polygon output path
                })

        except Exception as e: # Handle any exceptions that occur during processing
            logging.error(f"Error processing {raster_name} for subdistrict {subdistrict if 'subdistrict' in locals() else 'N/A'}: {e}")
    
    logging.info("Batch clipping finished.")

def main():
    clip_rasters_by_subdistrict_hierarchy(
        folder_path="LANDSAT_9/", 
        shapefile_path="Thailand/Thailand - Subnational Administrative Boundaries.shp",
        province_attr="PV_TN",
        district_attr="AP_TN",
        subdistrict_attr="TB_TN",
        output_dir="Clipped_Rasters",
    )

if __name__ == "__main__":
    main()
