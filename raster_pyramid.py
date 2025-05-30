import os
import glob
import logging
from collections import defaultdict
from osgeo import gdal, osr

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def build_overviews(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest'):
    """
    Builds raster pyramid overviews for a given GeoTIFF file.

    Args:
        filepath (str): Path to the GeoTIFF file.
        overview_levels (list): List of integers representing the downsampling factors
                                for each overview level. Default is [2, 4, 8, 16, 32].
        resampling_method (str): Resampling method to use for overview creation.
                                 Common options include 'average', 'nearest', 'cubic', 'mode',
                                 'lanczos'. Default is 'average'.
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


def main():
    # Define the root directory to search for GeoTIFF files
    search_directory = r'Raster_Mosaic'
    # Define the pattern for GeoTIFF files (e.g., .tif, .tiff)
    # The '**' will match any files and directories recursively
    pattern = os.path.join(search_directory, '**', '*.tif')

    logger.info(f"Searching for GeoTIFF files in {search_directory} using pattern: {pattern}")

    # Use glob to find all files matching the pattern recursively
    for filepath in glob.glob(pattern, recursive=True):
        logger.info(f"Found GeoTIFF file: {filepath}")
        build_overviews(filepath)


if __name__ == "__main__":
    main()