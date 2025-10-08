import os # Library with methods for interacting with the operating system, like creating files and directories.
import glob # Library for finding files and directories that match a specified pattern
import logging # Library for logging messages for debugging and information
import shutil # Library for copying files and directories, including metadata and overviews
from collections import defaultdict # Library for creating dictionaries with default values
from osgeo import gdal, osr, ogr # Library for working with geospatial data, like GDAL (Geospatial Data Abstraction Library) and OGR (OpenGIS)

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s") # Set up logging configuration
logger = logging.getLogger(__name__) # Create a logger object

# --- HELPER FUNCTIONS (No changes here) ---
def create_clean_cutline(source_shp_path):
    """
    Reads a shapefile, dissolves all features into a single geometry,
    fixes the geometry using a buffer(0) operation, and saves it
    to a new in-memory shapefile.

    Args:
        source_shp_path (str): Path to the source shapefile.

    Returns:
        str: The in-memory path to the new, clean shapefile, or None on failure.
    """
    memory_shp_path = f"/vsimem/{os.path.basename(source_shp_path)}_clean.shp" # Create a temporary in-memory path for the shapefile
    logger.info(f"Creating a clean, dissolved cutline from {source_shp_path}...") # Log the operation

    try:
        source_ds = ogr.Open(source_shp_path) # Open the source shapefile using OGR
        if not source_ds: # Check if the shapefile could be opened
            logger.error(f"Unable to open shapefile: {source_shp_path}") # Log the error message
            return None # Return None to indicate failure
        source_layer = source_ds.GetLayer() # Get the layer from the shapefile
        srs = source_layer.GetSpatialRef() # Get the spatial reference system from the layer

        dissolved_geom = ogr.Geometry(ogr.wkbMultiPolygon) # Create an empty multi-polygon geometry
        for feature in source_layer: # Loop through all features in the layer
            geom = feature.GetGeometryRef() # Get the geometry of the feature
            if geom: # Check if the feature has a geometry
                dissolved_geom = dissolved_geom.Union(geom) # Add the geometry to the dissolved geometry

        logger.info("Applying buffer(0) to fix potential geometry issues...") # Log the operation
        clean_geom = dissolved_geom.Buffer(0) # Apply a buffer(0) operation to the dissolved geometry to fix potential geometry issues

        if clean_geom.IsEmpty(): # Check if the geometry is empty after cleaning 
            logger.error("Geometry is empty after cleaning. Cannot create cutline.") # Log the error message
            return None # Return None to indicate failure

        driver = ogr.GetDriverByName('ESRI Shapefile') # Get the ESRI Shapefile driver to create the in-memory shapefile
        if os.path.exists(memory_shp_path): # Check if the in-memory shapefile already exists
              driver.DeleteDataSource(memory_shp_path) # Delete the existing in-memory shapefile

        mem_ds = driver.CreateDataSource(memory_shp_path) # Create the in-memory shapefile
        mem_layer = mem_ds.CreateLayer('clean', srs, geom_type=ogr.wkbPolygon) # Create a layer in the in-memory shapefile

        feature_defn = mem_layer.GetLayerDefn() # Get the definition of the layer
        new_feature = ogr.Feature(feature_defn) # Create a new feature
        new_feature.SetGeometry(clean_geom) # Set the geometry of the feature
        mem_layer.CreateFeature(new_feature) # Create the feature in the layer

        new_feature = None # Clean up the feature
        mem_ds = None   # Clean up the in-memory shapefile
        source_ds = None # Clean up the source shapefile

        logger.info(f"Clean cutline created at in-memory path: {memory_shp_path}") # Log the operation success
        return memory_shp_path # Return the in-memory path to the new, clean shapefile
    
    except Exception as e: # Catch any exceptions that may occur
        logger.error(f"Failed to create clean cutline: {e}") # Log the error message
        return None # Return None to indicate failure

def analyze_rasters(files):
    """
    Analyze rasters to determine the most common projection and average resolution.
    """
    logger.info("Analyzing rasters to determine optimal mosaic parameters...")
    srs_counts = defaultdict(int) # Create a dictionary to store projection counts
    x_res_list = [] # Create a list to store x-resolution values
    y_res_list = [] # Create a list to store y-resolution values

    for f in files: # Loop through all rasters
        ds = gdal.Open(f) # Open the raster
        if ds is None: # Check if the raster could be opened
            logger.warning(f"Cannot open {f} for analysis") # Log the warning
            continue # Skip to the next raster

        srs = osr.SpatialReference() # Create a SpatialReference object
        proj_wkt = ds.GetProjection() # Get the projection WKT
        if proj_wkt: # Check if the projection WKT is not empty
            srs.ImportFromWkt(proj_wkt) # Import the projection WKT into the SpatialReference object
            authority_name = srs.GetAuthorityName(None) # Get the authority name
            authority_code = srs.GetAuthorityCode(None) # Get the authority code
            
            if authority_name and authority_name.upper() == "EPSG" and authority_code: # Check if the projection is EPSG
                srs_counts[f"EPSG:{authority_code}"] += 1 # Increment the count for the EPSG code

        gt = ds.GetGeoTransform() # Get the geotransform values
        x_res_list.append(abs(gt[1])) # Append the absolute value of the x-resolution
        y_res_list.append(abs(gt[5])) # Append the absolute value of the y-resolution
        ds = None

    # Determine target EPSG, prioritizing the local projected CRS for Thailand (32647)
    if "EPSG:32647" in srs_counts: # If 32647 is found
        target_epsg = "EPSG:32647" # Set target EPSG to 32647
        logger.info("Detected EPSG:32647 (UTM Zone 47N). This will be the target.") # Log the operation
    elif "EPSG:4326" in srs_counts: # If 4326 is found
        target_epsg = "EPSG:4326" # Set target EPSG to 4326
        logger.info("Detected EPSG:4326. This will be the target.") # Log the operation
    elif srs_counts:  # Other projections found, but not 32647 or 4326
        target_epsg = "EPSG:32647"
        logger.info(f"Detected other projected CRS. Defaulting to {target_epsg} for Thailand.") # Log the operation
    else:  # No projections found
        target_epsg = "EPSG:4326" # Set target EPSG to 4326 as a default
        logger.warning("Could not determine a common EPSG code. Defaulting to EPSG:4326.") # Log the warning message

    if not x_res_list or not y_res_list: # If there are no resolution values 
        logger.error("Could not determine average resolution from any input files.") # Log the error
        return target_epsg, None, None # Return the target EPSG and None for resolution
    
    avg_x_res = sum(x_res_list) / len(x_res_list) # Calculate the average x-resolution
    avg_y_res = sum(y_res_list) / len(y_res_list) # Calculate the average y-resolution

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}") # Log the analysis results
    return target_epsg, avg_x_res, avg_y_res # Return the target EPSG, average x-resolution, and average y-resolution

def raster_pyramid(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest'):
    """
    Build internal raster pyramid for a raster file. For efficient to display in GIS Software like QGIS.

    Args:
        filepath (str): Path to the raster file.
    """
    logger.info(f"Building raster pyramid for {filepath} using levels {overview_levels} with {resampling_method} resampling...") # Log the operation message

    try:
        ds = gdal.Open(filepath, gdal.GA_Update) # Open the raster file in internal raster pyramid 
        if ds is None:
            logger.error(f"Cannot open {filepath} to build raster pyramid.") # Log the error message
            return # Return

        ds.BuildOverviews(resampling_method.upper(), overview_levels) # Build the internal raster pyramid using the specified resampling method and overview levels

        logger.info(f"Successfully built raster pyramid for {filepath}") # Log the success message for the operation

        ds = None # Close the dataset to flush changes 
    except Exception as e: # Catch any exceptions that may occur
        logger.error(f"Failed to build overviews for {filepath}: {e}") # Log the error message with the exception details

def raster_overlays_shapefile(raster_path, shapefile_path):
    """
    Return True if raster spatially overlaps with any feature in the shapefile
    """
    try:
        shapefile = ogr.Open(shapefile_path) # Open the shapefile using OGR
        if not shapefile: # Check if the shapefile could be opened 
            logger.error(f"Unable to open shapefile: {shapefile_path}") # Log the error message
            return False # Return False to indicate failure
        layer = shapefile.GetLayer() # Get the layer from the shapefile
 
        ds = gdal.Open(raster_path) # Open the raster using GDAL
        if not ds: # Check if the raster could be opened
            logger.warning(f"Cannot open raster {raster_path} to check overlay") # Log the warning
            return False # Return False to indicate failure
        gt = ds.GetGeoTransform() # Get the geotransform
        cols = ds.RasterXSize # Get the number of columns
        rows = ds.RasterYSize # Get the number of rows
        minx, maxx = gt[0], gt[0] + cols * gt[1] # Calculate the minimum and maximum x-coordinates
        maxy, miny = gt[3], gt[3] + rows * gt[5] # Calculate the maximum and minimum y-coordinates
        ds = None # Close the dataset

        ring = ogr.Geometry(ogr.wkbLinearRing) # Create a linear ring
        ring.AddPoint(minx, miny) # Add the minimum x- and y-coordinates 
        ring.AddPoint(minx, maxy) # Add the maximum x- and y-coordinates
        ring.AddPoint(maxx, maxy) # Add the maximum x- and y-coordinates
        ring.AddPoint(maxx, miny) # Add the minimum x- and y-coordinates
        ring.AddPoint(minx, miny) # Add the minimum x- and y-coordinates

        raster_geom = ogr.Geometry(ogr.wkbPolygon) # Create a polygon
        raster_geom.AddGeometry(ring) # Add the linear ring

        for feature in layer: # Loop through all features in the layer
            geom = feature.GetGeometryRef() # Get the geometry of the feature
            if geom and raster_geom.Intersects(geom): # Check if the raster intersects with the feature
                return True # Return True to indicate success

        return False # Return False to indicate failure

    except Exception as e: # Catch any exceptions that may occur
        logger.error(f"Overlay test failed for {raster_path}: {e}") # Log the error message
        return False # Return False to indicate failure

def main():
    """
    Script entry point.

    This script takes a directory of GeoTIFFs with different resolutions and projections
    and creates a mosaic for each subfolder within. The mosaic is reprojected to the
    average resolution of the input rasters and aligned to the Thailand subnational
    administrative boundaries shapefile.

    The script will:

    1. Create a clean version of the shapefile by dissolving all features into a single
       geometry and fixing any geometry issues with a buffer(0) operation.
    2. Loop through all subfolders in the input directory and check which rasters
       overlap with the shapefile.
    3. For each subfolder, determine the average resolution of the overlapping rasters.
    4. Reproject each overlapping raster to the average resolution and aligned to the
       shapefile.
    5. Build a VRT file from the reprojected rasters.
    6. Create a final mosaic from the VRT file using gdal.Translate.
    7. Build overviews for the final mosaic.
    8. Clean up temporary files.

    The script will output one mosaic per subfolder in the input directory, with the
    filename format <subfolder_name>_Canopy_Mosaic.tif.

    Example usage:

        python canopy_raster_mosiac.py

    """
    root_dir = r'Canopy' # Root directory containing subfolders of GeoTIFFs
    output_dir = r'Raster_Mosaic' # Directory to store the final mosaics
    shapefile_path = r'Thailand\L05_Province_ESRI_2559.shp' # Path to the shapefile
    os.makedirs(output_dir, exist_ok=True) # Create the output directory

    nodata_value = 0 # Value to use for nodata

    clean_shapefile_path = create_clean_cutline(shapefile_path) # Create a clean version of the shapefile
    if not clean_shapefile_path: # Check if the clean shapefile could be created
        logger.error("Could not create a clean cutline from the source shapefile. Exiting.") # Log the error message
        return # Exit the script

    subfolders = [f.path for f in os.scandir(root_dir) if f.is_dir()] # Get a list of subfolders in the root directory

    if not subfolders: # Check if any subfolders were found 
        logger.warning(f"No subfolders found in {root_dir}. Exiting.") # Log the warning
        return # Exit the script

    for subfolder_path in subfolders: # Loop through all subfolders
        subfolder_name = os.path.basename(subfolder_path) # Get the name of the subfolder
        logger.info(f"\n--- Processing subfolder: {subfolder_name} ---") # Log the subfolder name

        final_output_filename = f"{subfolder_name}_Canopy_Mosaic.tif" # Create the output filename
        final_output_path = os.path.join(output_dir, final_output_filename) # Create the output path

        # Find all raster files within the current subfolder
        raster_files = glob.glob(os.path.join(subfolder_path, '*.tif')) + \
                       glob.glob(os.path.join(subfolder_path, '*.tiff'))

        if not raster_files: # Check if any raster files were found
            logger.warning(f"No raster files found in {subfolder_path}. Skipping.") # Log the warning
            continue # Skip to the next subfolder

        overlaying_rasters = [] # List to store the rasters that overlap the shapefile
        for raster_file in raster_files: # Loop through all raster files
            logger.info(f"Checking overlap for: {raster_file}") # Log the raster file name
            if raster_overlays_shapefile(raster_file, clean_shapefile_path): # Check if the raster overlaps the shapefile
                logger.info(" -> Overlaps shapefile") # Log the result
                overlaying_rasters.append(raster_file) # Add the raster to the list
            else: # The raster does not overlap the shapefile
                logger.info(" -> No overlap") # Log the result

        if not overlaying_rasters: # Check if any rasters overlap the shapefile
            logger.warning(f"No rasters overlap the shapefile in {subfolder_path}. Skipping.") # Log the warning
            continue # Skip to the next subfolder

        logger.info(f"Selected {len(overlaying_rasters)} raster(s) for mosaicking") # Log the number of selected rasters

        target_epsg, x_res, y_res = analyze_rasters(overlaying_rasters) # Analyze the rasters
        if x_res is None or y_res is None: # Check if the average resolution could be determined
            logger.error(f"Failed to determine average resolution for {subfolder_name}. Skipping this subfolder.") # Log the error
            continue # Skip to the next subfolder

        all_reprojected = [] # List to store the reprojected rasters
        reprojected_temp_dir = os.path.join(output_dir, f"temp_{subfolder_name}_reprojected") # Create a temporary directory
        os.makedirs(reprojected_temp_dir, exist_ok=True) # Create the temporary directory

        for i, raster_file in enumerate(overlaying_rasters): # Loop through all selected rasters
            try: 
                # Check if reprojection is needed 
                ds_check = gdal.Open(raster_file) # Open the dataset
                if ds_check is None: # Check if the dataset could be opened
                    logger.warning(f"Cannot open {raster_file} to check CRS. Skipping.") # Log the warning
                    continue # Skip to the next raster

                srs_check = osr.SpatialReference() # Create a SpatialReference object
                source_epsg_str = None # String to store the EPSG code 
                proj_wkt_check = ds_check.GetProjection() # Get the projection
                if proj_wkt_check: # Check if the projection is not empty
                    srs_check.ImportFromWkt(proj_wkt_check) # Import the projection
                    authority_name = srs_check.GetAuthorityName(None) # Get the authority name
                    authority_code = srs_check.GetAuthorityCode(None) # Get the authority code
                    if authority_name and authority_name.upper() == "EPSG" and authority_code: # Check if the projection is EPSG
                        source_epsg_str = f"EPSG:{authority_code}" # Get the EPSG code
                ds_check = None  # Close the dataset

                if source_epsg_str == target_epsg: # Check if the source and target EPSG codes are the same
                    logger.info(f"'{os.path.basename(raster_file)}' is already in target CRS ({target_epsg}). Skipping reprojection.") # Log the message
                    all_reprojected.append(raster_file) # Add the raster to the list
                    continue # Skip to the next raster

                base_name = os.path.basename(raster_file) # Get the base name of the raster
                reprojected_path = os.path.join(reprojected_temp_dir, f"reproj_{i}_{base_name}") # Create the reprojected path

                logger.info(f"Reprojecting '{base_name}' from {source_epsg_str or 'Unknown CRS'} to {target_epsg}...") # Log the message
                warp_options = gdal.WarpOptions( # Create the warp options
                    dstSRS=target_epsg, # Set the target CRS
                    xRes=x_res, # Set the x resolution
                    yRes=y_res, # Set the y resolution
                    targetAlignedPixels=True, # Align to the shapefile
                    srcNodata=0, # Assuming 0 is nodata in source
                    dstNodata=0, # Set nodata in output
                    resampleAlg='near', # Set the resampling algorithm
                    dstNodata=nodata_value, # Set the nodata value
                    outputType=gdal.GDT_UInt16, # Set the output type
                    creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'], # Set the creation options
                    errorThreshold=0.0 # Set the error threshold
                )
                ds = gdal.Warp(reprojected_path, raster_file, options=warp_options) # Reproject the raster
                if ds is not None: # Check if the dataset was successfully reprojected
                    all_reprojected.append(reprojected_path) # Add the reprojected path to the list
                ds = None # Close the dataset

            except Exception as e: # Catch any exceptions
                logger.error(f"Failed to process or reproject {raster_file}: {e}") # Log the error
                continue # Skip to the next raster

        if not all_reprojected: # Check if any rasters were successfully reprojected
            logger.error(f"No rasters were successfully reprojected for subfolder {subfolder_name}!") # Log the error
            shutil.rmtree(reprojected_temp_dir) # Remove the temporary directory
            continue # Skip to the next subfolder

        vrt_path = os.path.join(reprojected_temp_dir, f'aligned_mosaic_{subfolder_name}.vrt') # Create the VRT path
        logger.info(f"Building VRT for {subfolder_name} from reprojected files...") # Log the message

        vrt = gdal.BuildVRT( # Build the VRT file for the reprojected rasters
            vrt_path, # Set the VRT path
            all_reprojected, # Set the list of reprojected rasters
            options=gdal.BuildVRTOptions( # Set the VRT options
                resampleAlg='nearest', # Set the resampling algorithm
                addAlpha=False, # Set to add alpha channel
                separate=False, # Set to separate bands
                srcNodata=nodata_value, # Set the source nodata value
                VRTNodata=nodata_value # Set the VRT nodata value
            )
        )
        if vrt is None: # Check if the VRT was successfully built
            logger.error(f"Failed to build VRT for {subfolder_name}") # Log the error
            shutil.rmtree(reprojected_temp_dir) # Remove the temporary directory
            continue # Skip to the next subfolder
        vrt = None # Close the VRT
 
        logger.info(f"Creating final mosaic for {subfolder_name}...") # Log the message

        gdal.Translate( # Translate the VRT to a TIFF file
            final_output_path, # Set the output path
            vrt_path, # Set the VRT path
            options=gdal.TranslateOptions( # Set the translation options
                format='GTiff', # Set the output format
                creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES', 'PREDICTOR=2'], # Set the creation options
                noData=nodata_value # Set the nodata value
            )
        )
        logger.info(f"Final mosaic for {subfolder_name} saved to: {final_output_path}") # Log the message

        raster_pyramid(final_output_path) # Build the internal raster pyramid

        if os.path.exists(final_output_path): # Check if the final output file exists
            logger.info(f"Cleaning up temporary files for {subfolder_name}...") # Log the message
            try:
                shutil.rmtree(reprojected_temp_dir) # Remove the temporary directory
            except OSError as e: # Catch any exceptions
                logger.warning(f"Could not remove temporary directory {reprojected_temp_dir}: {e}.") # Log the warning

        logger.info(f"Mosaic creation for {subfolder_name} complete!") # Log the message

    logger.info("All processed!") # Log the message

if __name__ == "__main__":
    main()
