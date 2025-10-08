import os
import glob
import logging
import shutil
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
                                 'lanczos'. Default is 'nearest'.
    """
    logger.info(f"Building overviews for {filepath} using levels {overview_levels} with {resampling_method} resampling...")
    try:
        ds = gdal.Open(filepath, gdal.GA_ReadOnly) # Open the dataset in read-only mode
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
    root_dir = r'Raster_Pyramid_Output' # This is now the parent directory containing subfolders

    output_dir = r'Raster_Mosaic'

    os.makedirs(output_dir, exist_ok=True)

    processing_jobs = []
    # Find all direct subdirectories within root_dir
    subfolders = [f.path for f in os.scandir(root_dir) if f.is_dir()]

    if subfolders:
        logger.info(f"Found {len(subfolders)} sub-directories to process.")
        for subfolder_path in subfolders:
            processing_jobs.append({
                'name': os.path.basename(subfolder_path),
                'path': subfolder_path
            })
    else:
        # If no subfolders, check for rasters directly in the root directory
        raster_files_in_root = glob.glob(os.path.join(root_dir, '*.tif')) + \
                               glob.glob(os.path.join(root_dir, '*.tiff'))
        if raster_files_in_root:
            logger.info("No sub-directories found. Processing raster files directly in the root directory.")
            processing_jobs.append({
                'name': os.path.basename(root_dir), # Use root directory name for the mosaic
                'path': root_dir
            })

    if not processing_jobs:
        logger.warning(f"No sub-directories or raster files found to process in '{root_dir}'. Exiting.")
        return

    for job in processing_jobs:
        job_name = job['name']
        job_path = job['path']
        logger.info(f"\n--- Processing job: {job_name} ---")

        # Define output path for the current subfolder's mosaic
        final_output_filename = f"{job_name}_THEOS-2_Mosaic.tif"
        final_output_path = os.path.join(output_dir, final_output_filename)

        # Find all raster files within the current subfolder
        raster_files = glob.glob(os.path.join(job_path, '*.tif')) + \
                       glob.glob(os.path.join(job_path, '*.tiff'))

        if not raster_files:
            logger.warning(f"No raster files found in {job_path}. Skipping.")
            continue

        logger.info(f"Found {len(raster_files)} raster files in {job_name} to process")

        # Analyze rasters for the current subfolder
        target_epsg, x_res, y_res = analyze_rasters(raster_files)
        if x_res is None or y_res is None:
            logger.error(f"Failed to determine average resolution for {job_name}. Skipping this job.")
            continue

        # Define the temporary directory for this subfolder's reprojected files
        reprojected_temp_dir = os.path.join(output_dir, f"temp_{job_name}_reprojected")
        os.makedirs(reprojected_temp_dir, exist_ok=True)

        # Reproject all rasters for the current subfolder
        all_reprojected = []
        for i, raster_file in enumerate(raster_files):
            try:
                ds = gdal.Open(raster_file)
                if ds is None:
                    logger.warning(f"Cannot open {raster_file}. Skipping.")
                    continue

                srs = osr.SpatialReference()
                srs.ImportFromWkt(ds.GetProjection())
                source_epsg = srs.GetAuthorityCode(None)
                ds = None

                # If projection of raster is EPSG:4326 or EPSG:32647 just skipping reprojection
                if source_epsg in ("4326", "32647"):
                    logger.info(f"{os.path.basename(raster_file)} is in a supported CRS (EPSG:{source_epsg}), skipping reprojection.")
                    all_reprojected.append(raster_file)
                    continue

                base_name = os.path.basename(raster_file)
                reprojected_path = os.path.join(reprojected_temp_dir, f"reproj_{i}_{base_name}")

                logger.info(f"Reprojecting {base_name} to {target_epsg} with resolution {x_res}, {y_res} and aligned pixels")
                warp_options = gdal.WarpOptions(
                    dstSRS=target_epsg,
                    xRes=x_res,
                    yRes=y_res,
                    targetAlignedPixels=True, # Ensure pixels are aligned to the resolution grid for reprojection
                    resampleAlg='near',
                    srcNodata=0, # Assuming 0 is nodata in source
                    dstNodata=0, # Set nodata in output
                    outputType=gdal.GDT_UInt16,
                    creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'],
                    errorThreshold=0.0 # Strict reprojection, no tolerance
                )
                ds = gdal.Warp(reprojected_path, raster_file, options=warp_options)
                if ds is None:
                    logger.error(f"gdal.Warp failed for {raster_file} and returned None.")
                else:
                    all_reprojected.append(reprojected_path)
                ds = None

            except Exception as e:
                logger.error(f"Failed to reproject {raster_file}: {e}")
                continue

        if not all_reprojected:
            logger.error(f"No rasters were successfully reprojected for job '{job_name}'!")
            # Clean up empty temp directory if it was created
            if os.path.exists(reprojected_temp_dir):
                shutil.rmtree(reprojected_temp_dir, ignore_errors=True)
            continue

        # Build VRT for the current subfolder's reprojected files
        vrt_path = os.path.join(reprojected_temp_dir, f'aligned_mosaic_{job_name}.vrt') # Unique VRT name
        logger.info(f"Building VRT for {job_name} from reprojected files...")
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
            logger.error(f"Failed to build VRT for {job_name}")
            continue
        vrt = None

        # Translate VRT to final GeoTIFF for the current subfolder
        logger.info(f"Creating final mosaic for {job_name}...")
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
        logger.info(f"Final mosaic for {job_name} saved to: {final_output_path}")

        # Build overviews for the final mosaic
        build_overviews(final_output_path)

        # Clean up temporary files for the current subfolder
        if os.path.exists(final_output_path):
            logger.info(f"Cleaning up temporary files for {job_name}...")
            # The reprojected_temp_dir contains the VRT and any reprojected files.
            # We can remove the whole directory.
            try:
                # Remove the temporary reprojected files directory
                shutil.rmtree(reprojected_temp_dir)
            except OSError as e:
                logger.warning(f"Could not remove temporary directory {reprojected_temp_dir}: {e}. It might not be empty.")

        logger.info(f"Mosaic creation for {job_name} complete!")

    logger.info("All processing jobs complete!")

if __name__ == "__main__":
    main()