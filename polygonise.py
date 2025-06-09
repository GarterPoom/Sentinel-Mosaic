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

