import os
import glob
import logging
import shutil # Import shutil for rmtree
from collections import defaultdict
from osgeo import gdal, osr, ogr

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# --- NEW HELPER FUNCTION ---
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
    memory_shp_path = f"/vsimem/{os.path.basename(source_shp_path)}_clean.shp"
    logger.info(f"Creating a clean, dissolved cutline from {source_shp_path}...")

    # Use a try-except block for robustness
    try:
        # Open source shapefile
        source_ds = ogr.Open(source_shp_path)
        if not source_ds:
            logger.error(f"Unable to open shapefile: {source_shp_path}")
            return None
        source_layer = source_ds.GetLayer()
        srs = source_layer.GetSpatialRef()

        # Union all geometries into a single one
        dissolved_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        for feature in source_layer:
            geom = feature.GetGeometryRef()
            if geom:
                dissolved_geom = dissolved_geom.Union(geom)

        # Fix the geometry using the buffer(0) trick
        logger.info("Applying buffer(0) to fix potential geometry issues...")
        clean_geom = dissolved_geom.Buffer(0)

        # If the geometry is empty after cleaning, it's an error
        if clean_geom.IsEmpty():
            logger.error("Geometry is empty after cleaning. Cannot create cutline.")
            return None

        # Create an in-memory shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        # Delete if it already exists in memory from a previous failed run
        if os.path.exists(memory_shp_path):
             driver.DeleteDataSource(memory_shp_path)

        mem_ds = driver.CreateDataSource(memory_shp_path)
        # Ensure the layer is created with the correct geometry type and SRS
        mem_layer = mem_ds.CreateLayer('clean', srs, geom_type=ogr.wkbPolygon)

        # Add the clean geometry as a single feature
        feature_defn = mem_layer.GetLayerDefn()
        new_feature = ogr.Feature(feature_defn)
        new_feature.SetGeometry(clean_geom)
        mem_layer.CreateFeature(new_feature)

        # Flush to /vsimem/ and close datasources
        new_feature = None
        mem_ds = None
        source_ds = None

        logger.info(f"Clean cutline created at in-memory path: {memory_shp_path}")
        return memory_shp_path

    except Exception as e:
        logger.error(f"Failed to create clean cutline: {e}")
        return None

def analyze_rasters(files):
    """
    Analyze rasters to determine the most common projection and average resolution.
    This version finds the most common EPSG code regardless of whether it is
    projected or geographic to minimize unnecessary reprojections.
    """
    logger.info("Analyzing rasters to determine optimal mosaic parameters...")
    srs_counts = defaultdict(int)
    x_res_list = []
    y_res_list = []

    for f in files:
        ds = gdal.Open(f)
        if ds is None:
            logger.warning(f"Cannot open {f} for analysis")
            continue

        # Get projection and resolution
        srs = osr.SpatialReference()
        proj_wkt = ds.GetProjection()
        if proj_wkt:
            srs.ImportFromWkt(proj_wkt)
            # Use a generic authority lookup that works for both PROJCS and GEOGCS
            authority_name = srs.GetAuthorityName(None)
            authority_code = srs.GetAuthorityCode(None)
            
            # Ensure we only count valid EPSG codes
            if authority_name and authority_name.upper() == "EPSG" and authority_code:
                srs_counts[f"EPSG:{authority_code}"] += 1

        gt = ds.GetGeoTransform()
        x_res_list.append(abs(gt[1]))
        y_res_list.append(abs(gt[5]))
        ds = None

    # --- SIMPLIFIED DECISION LOGIC ---
    # If any EPSG code was found, use the most common one as the target.
    # This prevents reprojecting files that are already in the most common CRS.
    if not srs_counts:
        target_epsg = "EPSG:4326"
        logger.warning("Could not determine a common EPSG code. Defaulting to EPSG:4326.")

    elif "EPSG:4326" in srs_counts:
        target_epsg = "EPSG:4326"
        logger.info("Detected EPSG:4326. This will be the target.")

    else:
        target_epsg = "EPSG:32647"
        logger.info(f"Determined CRS to be: {target_epsg}. UTM Zone 47 North (UTM Zone 47N) Coverage Thailand.")

    # Average resolution
    if not x_res_list or not y_res_list:
        logger.error("Could not determine average resolution from any input files.")
        return target_epsg, None, None
    
    avg_x_res = sum(x_res_list) / len(x_res_list)
    avg_y_res = sum(y_res_list) / len(y_res_list)

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}")
    return target_epsg, avg_x_res, avg_y_res

def build_overviews(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest'):
    """
    Builds raster pyramid overviews for a given GeoTIFF file.
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

def raster_overlays_shapefile(raster_path, shapefile_path):
    """
    Return True if raster spatially overlaps (fully or partially) with any feature in the shapefile
    """
    try:
        shapefile = ogr.Open(shapefile_path)
        if not shapefile:
            logger.error(f"Unable to open shapefile: {shapefile_path}")
            return False
        layer = shapefile.GetLayer()
        # Since our clean shapefile has only one feature, we can simplify this
        shapefile_feature = layer.GetNextFeature()
        shapefile_geom = shapefile_feature.GetGeometryRef()

        ds = gdal.Open(raster_path)
        if not ds:
            logger.warning(f"Cannot open raster {raster_path} to check overlay")
            return False
        gt = ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        minx, maxx, miny, maxy = gt[0], gt[0] + cols * gt[1], gt[3] + rows * gt[5], gt[3]
        ds = None

        # Ensure correct extent calculation for rotated rasters
        if gt[2] != 0 or gt[4] != 0:
             # More complex calculation needed for rotated rasters, but for now we assume they are not.
             pass

        # Create raster bounding box geometry
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(minx, miny)
        ring.AddPoint(minx, maxy)
        ring.AddPoint(maxx, maxy)
        ring.AddPoint(maxx, miny)
        ring.AddPoint(minx, miny) # Close the ring

        raster_geom = ogr.Geometry(ogr.wkbPolygon)
        raster_geom.AddGeometry(ring)

        return shapefile_geom.Intersects(raster_geom)
    except Exception as e:
        logger.error(f"Overlay test failed for {raster_path}: {e}")
        return False

def main():
    # Configure input and output directories/paths
    root_dir = r'Canopy'
    output_dir = r'Raster_Mosaic'
    shapefile_path = r'Thailand\L05_Province_ESRI_2559.shp'
    os.makedirs(output_dir, exist_ok=True)

    # --- MODIFICATION ---
    # Create a clean, dissolved version of the cutline shapefile in memory first
    clean_shapefile_path = create_clean_cutline(shapefile_path)
    if not clean_shapefile_path:
        logger.error("Could not create a clean cutline from the source shapefile. Exiting.")
        return

    # Find all direct subdirectories within root_dir
    subfolders = [f.path for f in os.scandir(root_dir) if f.is_dir()]

    if not subfolders:
        logger.warning(f"No subfolders found in {root_dir}. Exiting.")
        return

    for subfolder_path in subfolders:
        subfolder_name = os.path.basename(subfolder_path)
        logger.info(f"\n--- Processing subfolder: {subfolder_name} ---")

        final_output_filename = f"{subfolder_name}_Canopy_Mosaic.tif"
        final_output_path = os.path.join(output_dir, final_output_filename)

        raster_files = glob.glob(os.path.join(subfolder_path, '*.tif')) + \
                       glob.glob(os.path.join(subfolder_path, '*.tiff'))

        if not raster_files:
            logger.warning(f"No raster files found in {subfolder_path}. Skipping.")
            continue

        logger.info(f"Found {len(raster_files)} raster files in {subfolder_name} to process")

        # --- MODIFICATION ---
        # Use the clean shapefile for the overlay check
        overlaying_rasters = [f for f in raster_files if raster_overlays_shapefile(f, clean_shapefile_path)]
        if not overlaying_rasters:
            logger.warning(f"No rasters overlay the shapefile in {subfolder_path}. Skipping.")
            continue

        target_epsg, x_res, y_res = analyze_rasters(overlaying_rasters)
        if x_res is None or y_res is None:
            logger.error(f"Failed to determine average resolution for {subfolder_name}. Skipping this subfolder.")
            continue

        all_reprojected = []
        reprojected_temp_dir = os.path.join(output_dir, f"temp_{subfolder_name}_reprojected")
        os.makedirs(reprojected_temp_dir, exist_ok=True)

        for i, raster_file in enumerate(overlaying_rasters):
            try:
                base_name = os.path.basename(raster_file)
                reprojected_path = os.path.join(reprojected_temp_dir, f"reproj_{i}_{base_name}")

                # --- MODIFICATION ---
                # Use the clean shapefile for the gdal.Warp cutline
                warp_options = gdal.WarpOptions(
                    dstSRS=target_epsg,
                    xRes=x_res,
                    yRes=y_res,
                    targetAlignedPixels=True,
                    resampleAlg='near',
                    srcNodata=0,
                    dstNodata=0,
                    cutlineDSName=clean_shapefile_path, # Use the clean shapefile
                    cropToCutline=True,
                    outputType=gdal.GDT_UInt16,
                    creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'],
                    errorThreshold=0.0
                )
                ds = gdal.Warp(reprojected_path, raster_file, options=warp_options)
                if ds is not None:
                    all_reprojected.append(reprojected_path)
                ds = None

            except Exception as e:
                logger.error(f"Failed to reproject {raster_file}: {e}")
                continue

        if not all_reprojected:
            logger.error(f"No rasters were successfully reprojected for subfolder {subfolder_name}!")
            # Clean up the empty temp directory
            shutil.rmtree(reprojected_temp_dir) # Use shutil.rmtree
            continue

        vrt_path = os.path.join(reprojected_temp_dir, f'aligned_mosaic_{subfolder_name}.vrt')
        logger.info(f"Building VRT for {subfolder_name} from reprojected files...")
        vrt = gdal.BuildVRT(
            vrt_path,
            all_reprojected,
            options=gdal.BuildVRTOptions(
                resampleAlg='nearest',
                addAlpha=False,
                separate=False,
                srcNodata=0,
                VRTNodata=0
            )
        )
        if vrt is None:
            logger.error(f"Failed to build VRT for {subfolder_name}")
            shutil.rmtree(reprojected_temp_dir) # Clean up on failure
            continue
        vrt = None

        logger.info(f"Creating final mosaic for {subfolder_name}...")
        gdal.Translate(
            final_output_path,
            vrt_path,
            options=gdal.TranslateOptions(
                format='GTiff',
                creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES', 'PREDICTOR=2']
            )
        )
        logger.info(f"Final mosaic for {subfolder_name} saved to: {final_output_path}")

        build_overviews(final_output_path)

        if os.path.exists(final_output_path):
            logger.info(f"Cleaning up temporary files for {subfolder_name}...")
            try:
                # --- MODIFICATION ---
                # Use shutil.rmtree to remove the directory and all its contents
                shutil.rmtree(reprojected_temp_dir)
            except OSError as e:
                logger.warning(f"Could not remove temporary directory {reprojected_temp_dir}: {e}.")

        logger.info(f"Mosaic creation for {subfolder_name} complete!")

    logger.info("All subfolders processed!")

if __name__ == "__main__":
    main()