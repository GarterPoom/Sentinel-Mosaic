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

def polygonize_clipped_raster(clipped_raster_path, output_polygon_path, nodata_fill_value, output_geojson_path=None):
    """
    Polygonizes raster files found in the directory specified by `clipped_raster_path`,
    creating polygons for valid data areas. Output polygons are saved as Shapefiles
    in the directory specified by `output_polygon_path`.

    Args:
        clipped_raster_path (str): Path to the DIRECTORY containing input raster files (e.g., .tif).
                                   The function will search for '*.tif' files within this directory.
        output_polygon_path (str): Path to the DIRECTORY where output polygon Shapefiles will be saved.
                                   Each Shapefile will be named after its corresponding raster file.
        nodata_fill_value (float or int): The value in the input rasters that represents nodata.
                                          (i.e., the fill value used by rasterio.mask).
        output_geojson_path (str, optional): Path to the DIRECTORY where output polygon GeoJSON files will be saved.
                                             If None, GeoJSON export is skipped. Defaults to None.
        # Note: This function assumes each input raster is already clipped to the desired area of interest.
    """
    # Arguments are now expected to be directories
    input_dir = clipped_raster_path
    output_dir = output_polygon_path
    geojson_output_dir = output_geojson_path

    try:
        os.makedirs(output_dir, exist_ok=True) # Ensure output directory exists

        # Use glob to find all .tif files recursively in the input_dir
        search_pattern = os.path.join(input_dir, '**', '*.tif')
        raster_files = glob.glob(search_pattern, recursive=True)

        if geojson_output_dir:
            os.makedirs(geojson_output_dir, exist_ok=True) # Ensure GeoJSON output directory exists

        if not raster_files:
            logging.info(f"No .tif files found in {input_dir} to polygonize.")
            return

        logging.info(f"Found {len(raster_files)} raster files in {input_dir} for polygonization.")

        for individual_raster_file in raster_files:
            base_name = os.path.basename(individual_raster_file)
            
            # Determine the relative path of the raster file from the input_dir
            relative_path_to_file = os.path.relpath(individual_raster_file, input_dir)
            # Determine the relative directory structure
            relative_dir_structure = os.path.dirname(relative_path_to_file)
            
            # Create corresponding output subdirectory
            current_output_sub_dir = os.path.join(output_dir, relative_dir_structure)
            os.makedirs(current_output_sub_dir, exist_ok=True)
            
            # Construct the output polygon file paths
            output_filename_base = os.path.splitext(base_name)[0]
            current_output_shapefile_path = os.path.join(current_output_sub_dir, output_filename_base)

            if geojson_output_dir:
                current_geojson_output_sub_dir = os.path.join(geojson_output_dir, relative_dir_structure)
                os.makedirs(current_geojson_output_sub_dir, exist_ok=True)
                current_output_geojson_path = os.path.join(current_geojson_output_sub_dir, output_filename_base + '.geojson')

            try: # Inner try-except for individual file processing
                with rasterio.open(individual_raster_file) as src:
                    image = src.read(1)  # Read the first band
                    transform = src.transform
                    crs = src.crs
                    
                    # Create a binary mask: 1 for valid data, 0 for nodata_fill_value
                    binary_mask = np.where(image == nodata_fill_value, 0, 1).astype(np.uint8)
                    
                    # Extract shapes (polygons) from the binary_mask where pixel value is 1
                    results = [
                        {'properties': {'raster_val': v}, 'geometry': s}
                        for i, (s, v) in enumerate(
                            shapes(binary_mask, mask=(binary_mask == 1), transform=transform)
                        ) if v == 1 # Filter for polygons derived from pixels that were set to 1
                    ]

                    if not results:
                        logging.info(f"No valid data (value 1) polygons found in {base_name} to polygonize.")
                        continue # Skip to the next file

                    geometries = [shape(result['geometry']) for result in results]
                    gdf_polygons = gpd.GeoDataFrame(geometry=geometries, crs=crs)
                    gdf_polygons['value'] = 1 # Add a column indicating the (binary) value of the polygonized area

                    gdf_polygons.to_file(current_output_shapefile_path, driver='ESRI ShapeFile')
                    logging.info(f"Polygonized {base_name} and saved to {current_output_shapefile_path}")

                    if geojson_output_dir and current_output_geojson_path:
                        gdf_polygons.to_file(current_output_geojson_path, driver='GeoJSON')
                        logging.info(f"Exported GeoJSON for {base_name} to {current_output_geojson_path}")

            except Exception as e_file:
                logging.error(f"Error polygonizing {base_name}: {e_file}")
                # Continue to the next file even if one fails

    except Exception as e_main: # Outer try-except for issues like directory access
        logging.error(f"Error in polygonization process for directory {input_dir}: {e_main}")

def main():
    """
    Main function to demonstrate and run the raster polygonization process.
    Sets up input/output directories and nodata value, then calls the polygonization function.
    """
    # --- Configuration ---
    # TODO: Replace these placeholder paths with actual paths to your data.
    # Example using raw strings for Windows paths:
    # input_raster_directory = r"D:\Landsat_Mosaic_Data\clipped_rasters"
    # output_polygon_directory = r"D:\Landsat_Mosaic_Data\output_polygons"
    
    # Using relative paths for easier example execution:
    input_raster_directory = "Clipped_Rasters" 
    output_polygon_directory = "Polygonized_Rasters"
    output_geojson_directory = "GeoJSON_Polygons" # New directory for GeoJSON outputs
    
    # This is the value in your rasters that signifies "no data".
    # Adjust this to your specific nodata value (e.g., -9999, 0, 255).
    nodata_value_in_raster = 0 

    # --- Setup (Optional: for demonstration, create dummy directories if they don't exist) ---
    if not os.path.exists(input_raster_directory):
        os.makedirs(input_raster_directory)
        logging.info(f"Created dummy input directory: {os.path.abspath(input_raster_directory)}")
        logging.info(f"Please place your .tif raster files in this directory for processing.")
    
    if not os.path.exists(output_polygon_directory):
        os.makedirs(output_polygon_directory)
        logging.info(f"Created dummy output directory: {os.path.abspath(output_polygon_directory)}")
        logging.info(f"Output .shp files will be saved in this directory.")
    
    if not os.path.exists(output_geojson_directory): # Create GeoJSON output directory if it doesn't exist
        os.makedirs(output_geojson_directory)
        logging.info(f"Created dummy GeoJSON output directory: {os.path.abspath(output_geojson_directory)}")
        logging.info(f"Output .geojson files will be saved in this directory.")

    # --- Execution ---
    logging.info(f"Starting polygonization process...")
    logging.info(f"Input raster directory: {os.path.abspath(input_raster_directory)}")
    logging.info(f"Output polygon directory: {os.path.abspath(output_polygon_directory)}")
    logging.info(f"Output GeoJSON directory: {os.path.abspath(output_geojson_directory)}")
    logging.info(f"Nodata value for masking: {nodata_value_in_raster}")

    polygonize_clipped_raster(
        clipped_raster_path=input_raster_directory,   # This is now a directory path
        output_polygon_path=output_polygon_directory, # This is now a directory path
        nodata_fill_value=nodata_value_in_raster,
        output_geojson_path=output_geojson_directory  # Pass the new GeoJSON directory
    )
    
    logging.info("Polygonization process finished.")

if __name__ == "__main__":
    main()
