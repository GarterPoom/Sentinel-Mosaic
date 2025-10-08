import os
import logging
import shutil
from osgeo import gdal

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def build_overviews(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest', output_dir=None):
    """
    Builds raster pyramid overviews for a given GeoTIFF file.
    """
    logger.info(f"Building overviews for {filepath} using levels {overview_levels} with {resampling_method} resampling...")
    try:
        if output_dir:
            os.makedirs(output_dir, exist_ok=True) # Ensure output directory exists
            filename = os.path.basename(filepath) # Get the filename
            output_filepath = os.path.join(output_dir, filename) # Define output path
            if not os.path.exists(output_filepath): # Avoid overwriting existing files
                shutil.copy2(filepath, output_filepath) # Copy file to output directory
                logger.info(f"Copied {filepath} to {output_filepath}") # Log the copy action
            target_filepath = output_filepath # Use the copied file for building overviews
        else:
            target_filepath = filepath # Use the original file if no output directory is specified

        ds = gdal.Open(target_filepath, gdal.GA_ReadOnly) # Open the dataset in read-only mode
        if ds is None: # Check if the dataset was opened successfully
            logger.error(f"Cannot open {target_filepath} to build overviews.") # Log error if opening fails
            return  # Return early if the dataset cannot be opened

        ds.BuildOverviews(resampling_method, overview_levels) # Build the overviews
        ds = None # Close the dataset
        logger.info(f"Successfully built overviews for {target_filepath}") # Log success
    except Exception as e: # Catch any exceptions that occur
        logger.error(f"Failed to build overviews for {filepath}: {e}") # Log the exception

def main():
    """
    Recursively build raster pyramid overviews for all GeoTIFFs in a directory tree.
    """
    search_directory = r'Root directory' # Root directory to search for GeoTIFF files
    output_directory = r'Raster_Pyramid_Output' # Directory to save processed files

    logger.info(f"Searching recursively for GeoTIFF files under: {search_directory}") # Log the search directory
    logger.info(f"Output directory: {output_directory}") # Log the output directory

    # Walk through all subdirectories recursively
    for root, dirs, files in os.walk(search_directory): # Iterate through files in the current directory     
        for file in files: # Check each file
            if file.lower().endswith('.tif'): # Process only .tif files
                filepath = os.path.join(root, file) # Get the full file path
                logger.info(f"Found GeoTIFF file: {filepath}") # Log the found file
                build_overviews(filepath, output_dir=output_directory) # Build overviews and save to output directory
                logger.info(f"Finished processing {filepath}") # Log completion of processing

    logger.info("All GeoTIFF files processed successfully.") # Log overall completion

if __name__ == "__main__":
    main()
