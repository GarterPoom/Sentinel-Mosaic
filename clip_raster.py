import glob
import os
import logging
from osgeo import gdal #, osr # osr might not be strictly needed here but often imported with gdal


# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def clip_raster_to_shapefile(raster_path, shapefile_path, output_path):
    """
    Clips a raster file to the extent of a shapefile using gdal.Warp and saves the output.

    Args:
        raster_path (str): Path to the input raster file.
        shapefile_path (str): Path to the shapefile defining the clipping area.
        output_path (str): Path where the clipped raster will be saved.
    """
    logger.info(f"Clipping raster {raster_path} to shapefile {shapefile_path} using gdal.Warp...")

    # Ensure output directory exists
    output_dir_for_file = os.path.dirname(output_path)
    if output_dir_for_file: # Check if dirname is not empty (e.g. for files in current dir)
        os.makedirs(output_dir_for_file, exist_ok=True)

    # Check if input raster exists
    if not os.path.exists(raster_path):
        logger.error(f"Input raster {raster_path} not found. Skipping.")
        return

    # Check if shapefile exists
    if not os.path.exists(shapefile_path):
        logger.error(f"Shapefile {shapefile_path} not found. Skipping clipping for {raster_path}.")
        return

    # Attempt to get NoData value from source raster
    src_nodata = None
    try:
        src_ds_check = gdal.Open(raster_path)
        if not src_ds_check:
            logger.warning(f"Could not open {raster_path} with GDAL to check for NoData. Proceeding without explicit NoData.")
        else:
            band = src_ds_check.GetRasterBand(1) # Check first band
            nodata_val = band.GetNoDataValue()
            if nodata_val is not None:
                src_nodata = nodata_val
                logger.info(f"Source NoData value for {raster_path}: {src_nodata}")
            else:
                logger.info(f"No explicit NoData value found in {raster_path}.")
            src_ds_check = None # Close dataset
    except Exception as e:
        logger.warning(f"Error reading NoData value from {raster_path}: {e}. Proceeding without explicit NoData.")

    warp_options_dict = {
        'format': 'GTiff',
        'cutlineDSName': shapefile_path,
        'cropToCutline': True,
        'dstAlpha': True,  # Add an alpha band for transparency outside cutline
        'creationOptions': ['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES', 'PREDICTOR=2'],
        'resampleAlg': gdal.GRIORA_NearestNeighbour, # Explicitly set resample algorithm
        'multithread': True, # Optional, if GDAL is built with it enabled
    }

    if src_nodata is not None:
        warp_options_dict['srcNodata'] = src_nodata
        warp_options_dict['dstNodata'] = src_nodata # Propagate NoData value

    try:
        ds = gdal.Warp(output_path, raster_path, **warp_options_dict)
        if ds is None:
            logger.error(f"gdal.Warp failed for {raster_path} (returned None). Output may be invalid or not created.")
            if os.path.exists(output_path) and os.path.getsize(output_path) == 0:
                 logger.warning(f"Output file {output_path} was created but is empty.")
            elif not os.path.exists(output_path):
                 logger.warning(f"Output file {output_path} was not created.")
            return
        ds = None  # Explicitly close the dataset
        logger.info(f"Successfully clipped raster and saved to {output_path}")

    except RuntimeError as e:
        logger.error(f"gdal.Warp failed for {raster_path}. GDAL RuntimeError: {e}")
        if os.path.exists(output_path): # Clean up failed output
            try:
                os.remove(output_path)
                logger.info(f"Removed partially created/failed output file: {output_path}")
            except OSError as ose:
                logger.warning(f"Could not remove partially created/failed output file {output_path}: {ose}")
    except Exception as e: # Catch any other unexpected errors
        logger.error(f"An unexpected error occurred during gdal.Warp for {raster_path}: {e}")
        if os.path.exists(output_path): # Clean up failed output
            try:
                os.remove(output_path)
                logger.info(f"Removed partially created/failed output file due to unexpected error: {output_path}")
            except OSError as ose:
                logger.warning(f"Could not remove partially created/failed output file {output_path} after unexpected error: {ose}")

def build_overviews(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest'):
    """
    Builds raster pyramid overviews for a given GeoTIFF file.

    Args:
        filepath (str): Path to the GeoTIFF file.
        overview_levels (list): List of integers representing the downsampling factors
                                for each overview level. Default is [2, 4, 8, 16, 32].
        resampling_method (str): Resampling method to use for overview creation.
                                 Common options include 'average', 'nearest', 'cubic', 'mode',
                                 'lanczos'. Default is 'nearest'.
    """
    logger.info(f"Building overviews for {filepath} using levels {overview_levels} with {resampling_method} resampling...")
    try:
        ds = gdal.Open(filepath, gdal.GA_Update)
        if ds is None:
            logger.error(f"Cannot open {filepath} to build overviews.")
            return
        
        # Build overviews
        ds.BuildOverviews(resampling_method, overview_levels)
        ds = None # Close the dataset to flush changes

        logger.info(f"Successfully built overviews for {filepath}")
    except Exception as e:
        logger.error(f"Failed to build overviews for {filepath}: {e}")

def main_process(): # Renamed original main to avoid conflict if this script is imported
    raster_path = r"Sentinel-2"
    shapefile_path = r"Thailand\Thailand - Subnational Administrative Boundaries.shp" # Ensure this path is correct
    output_dir = r"output_clipped_rasters" # Define an output directory
    os.makedirs(output_dir, exist_ok=True)

    raster_files = glob.glob(os.path.join(raster_path, '*.tif')) + \
                       glob.glob(os.path.join(raster_path, '*.tiff'))
    
    if not raster_files:
        logger.error(f"No raster files found in {raster_path}. Exiting.")
        return

    if not os.path.exists(shapefile_path):
        logger.error(f"Shapefile not found at {shapefile_path}. Exiting.")
        return

    for raster_file_path in raster_files:
        base_name = os.path.basename(raster_file_path)
        output_filename = os.path.join(output_dir, f"clipped_{base_name}")
        clip_raster_to_shapefile(raster_file_path, shapefile_path, output_filename)
        if os.path.exists(output_filename) and os.path.getsize(output_filename) > 0: # Check if clipping was successful
            build_overviews(output_filename)

if __name__ == "__main__":
    # Enable GDAL exceptions for more informative error messages from GDAL operations
    gdal.UseExceptions()
    main_process()