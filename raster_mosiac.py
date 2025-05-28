import os
import glob
import logging
from collections import defaultdict
from osgeo import gdal, osr

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def analyze_rasters(files):
    """Analyze rasters to determine projection and average resolution"""
    logger.info("Analyzing rasters to determine optimal mosaic parameters...")
    proj_counts = defaultdict(int)
    x_res_list = []
    y_res_list = []

    for f in files:
        ds = gdal.Open(f)
        if ds is None:
            logger.warning(f"Cannot open {f} for analysis")
            continue

        # Projection
        srs = osr.SpatialReference()
        srs.ImportFromWkt(ds.GetProjection())
        if srs.IsProjected():
            epsg = srs.GetAuthorityCode(None)
            proj_counts[epsg] += 1

        # Resolution
        gt = ds.GetGeoTransform()
        x_res_list.append(abs(gt[1]))
        y_res_list.append(abs(gt[5]))
        ds = None

    # Most common projection
    if not proj_counts:
        logger.warning("No projections found. Defaulting to EPSG:4326")
        target_epsg = "EPSG:4326"
    else:
        target_epsg = "EPSG:32647" # Use EPSG:32647 as default which is WGS84 UTM Zone 47 North (UTM Zone 47N) Coverage Thailand.

    # Average resolution - check if lists are empty to avoid ZeroDivisionError
    if not x_res_list or not y_res_list:
        logger.error("Could not determine average resolution from any input files.")
        return target_epsg, None, None # Indicate failure to calculate resolution
    avg_x_res = sum(x_res_list) / float(len(x_res_list)) # Use float for division
    avg_y_res = sum(y_res_list) / float(len(y_res_list)) # Use float for division

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}")
    return target_epsg, avg_x_res, avg_y_res

def build_overviews(filepath, overview_levels=[2, 4, 8, 16, 32], resampling_method='nearest'):
    """
    Builds raster pyramid overviews for a given GeoTIFF file.

    Args:
        filepath (str): Path to the GeoTIFF file.
        overview_levels (list): List of integers representing the downsampling factors
                                for each overview level. Default is [2, 4, 8, 16, 32].
        resampling_method (str): Resampling method to use for overview creation.
                                 Common options include 'average', 'nearest', 'cubic',
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
    # Configure input and output directories/paths
    root_dir = r'LANDSAT_9' # This is now the parent directory containing subfolders
    output_dir = r'Raster_Mosaic'

    os.makedirs(output_dir, exist_ok=True)

    # Find all direct subdirectories within root_dir
    subfolders = [f.path for f in os.scandir(root_dir) if f.is_dir()]

    if not subfolders:
        logger.warning(f"No subfolders found in {root_dir}. Exiting.")
        return

    for subfolder_path in subfolders:
        subfolder_name = os.path.basename(subfolder_path)
        logger.info(f"\n--- Processing subfolder: {subfolder_name} ---")

        # Define output path for the current subfolder's mosaic
        final_output_filename = f"{subfolder_name}_Mosaic.tif"
        final_output_path = os.path.join(output_dir, final_output_filename)

        # Find all raster files within the current subfolder
        # Note: We're only looking directly in the subfolder now, not recursively deeper
        raster_files = glob.glob(os.path.join(subfolder_path, '*.tif')) + \
                       glob.glob(os.path.join(subfolder_path, '*.tiff'))

        if not raster_files:
            logger.warning(f"No raster files found in {subfolder_path}. Skipping.")
            continue

        logger.info(f"Found {len(raster_files)} raster files in {subfolder_name} to process")

        # Analyze rasters for the current subfolder
        target_epsg, x_res, y_res = analyze_rasters(raster_files)
        if x_res is None or y_res is None:
            logger.error(f"Failed to determine average resolution for {subfolder_name}. Skipping this subfolder.")
            continue

        # Reproject all rasters for the current subfolder
        all_reprojected = []
        for i, raster_file in enumerate(raster_files):
            base_name = os.path.basename(raster_file)
            # Create a temporary reprojected path within the output_dir, specific to this subfolder
            reprojected_temp_dir = os.path.join(output_dir, f"temp_{subfolder_name}_reprojected")
            os.makedirs(reprojected_temp_dir, exist_ok=True)
            reprojected_path = os.path.join(reprojected_temp_dir, f"reproj_{i}_{base_name}")

            try:
                logger.info(f"Reprojecting {base_name} to {target_epsg} with resolution {x_res}, {y_res} and aligned pixels")
                warp_options = gdal.WarpOptions(
                                dstSRS=target_epsg,
                                xRes=x_res,
                                yRes=y_res,
                                targetAlignedPixels=True, # Align pixels
                                resampleAlg='near',
                                srcNodata=0, # Consider making this configurable or auto-detected if possible
                                dstNodata=0, # Consider making this configurable
                                outputType=gdal.GDT_UInt16,
                                creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'],
                                errorThreshold=0.0 # Default is 0.125, 0.0 means exact reprojection
                            )
                ds = gdal.Warp(
                    reprojected_path,
                    raster_file,
                    options=warp_options
                )
                if ds is None:
                    logger.error(f"gdal.Warp failed for {raster_file} and returned None, but did not raise an exception.")
                else:
                    all_reprojected.append(reprojected_path)
                ds = None
            except Exception as e:
                logger.error(f"Failed to reproject {raster_file}: {e}")
                continue

        if not all_reprojected:
            logger.error(f"No rasters were successfully reprojected for subfolder {subfolder_name}!")
            # Clean up empty temp directory if it was created
            if os.path.exists(reprojected_temp_dir) and not os.listdir(reprojected_temp_dir):
                os.rmdir(reprojected_temp_dir)
            continue

        # Build VRT for the current subfolder's reprojected files
        vrt_path = os.path.join(reprojected_temp_dir, f'aligned_mosaic_{subfolder_name}.vrt') # Unique VRT name
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
            continue
        vrt = None

        # Translate VRT to final GeoTIFF for the current subfolder
        logger.info(f"Creating final mosaic for {subfolder_name}...")
        gdal.Translate(
            final_output_path,
            vrt_path,
            options=gdal.TranslateOptions(
                format='GTiff',
                creationOptions=[
                    'TILED=YES',
                    'COMPRESS=LZW',
                    'BIGTIFF=YES',
                    'PREDICTOR=2'
                ]
            )
        )
        logger.info(f"Final mosaic for {subfolder_name} saved to: {final_output_path}")

        # Build overviews for the final mosaic
        build_overviews(final_output_path)

        # Clean up temporary files for the current subfolder
        if os.path.exists(final_output_path):
            logger.info(f"Cleaning up temporary files for {subfolder_name}...")
            try:
                os.remove(vrt_path)
            except Exception as e:
                logger.warning(f"Failed to remove VRT {vrt_path}: {e}")

            for file in all_reprojected:
                try:
                    os.remove(file)
                except Exception as e:
                    logger.warning(f"Failed to remove {file}: {e}")
            try:
                # Remove the temporary reprojected files directory
                os.rmdir(reprojected_temp_dir)
            except OSError as e:
                logger.warning(f"Could not remove temporary directory {reprojected_temp_dir}: {e}. It might not be empty.")

        logger.info(f"Mosaic creation for {subfolder_name} complete!")

    logger.info("All subfolders processed!")

if __name__ == "__main__":
    main()