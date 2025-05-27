import os
import re
import glob
import logging
from collections import defaultdict
from osgeo import gdal, osr

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def analyze_rasters(files):
    logger.info("Analyzing rasters to determine optimal mosaic parameters...")
    proj_counts = defaultdict(int)
    x_res_list = []
    y_res_list = []

    for f in files:
        try:
            ds = gdal.Open(f)
            if ds is None:
                logger.warning(f"Cannot open {f} for analysis. Skipping this file.")
                continue

            projection_wkt = ds.GetProjection()
            if not projection_wkt:
                logger.debug(f"No projection information found in {f}.")
            else:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(projection_wkt)
                if srs.IsProjected():
                    epsg_code = srs.GetAuthorityCode(None)
                    if epsg_code:
                        proj_counts[epsg_code] += 1
                    else:
                        logger.debug(f"Could not get Authority Code for projection in {f}")

            gt = ds.GetGeoTransform()
            if gt and len(gt) == 6 and gt[1] != 0 and gt[5] != 0: # Ensure valid transform and non-zero resolution
                x_res_list.append(abs(gt[1]))
                y_res_list.append(abs(gt[5]))
            else:
                logger.warning(f"Invalid, missing, or zero-resolution geotransform for {f}. Skipping resolution analysis for this file.")
            ds = None # Close dataset
        except RuntimeError as e:
            logger.warning(f"GDAL RuntimeError while analyzing {f}: {e}. Skipping this file.")
            continue
        except Exception as e:
            logger.warning(f"Unexpected error while analyzing {f}: {e}. Skipping this file.")
            continue

    if not x_res_list or not y_res_list: # If no resolutions were collected
        logger.error("No valid raster data with geotransform found in the input files for analysis. Cannot determine mosaic parameters.")
        return None, None, None

    # Determine target_epsg based on original logic
    if not proj_counts: # No projected CRS identified in valid rasters
        logger.warning("No projected CRS identified in valid rasters. Defaulting to EPSG:4326 (WGS84).")
        determined_epsg = "EPSG:4326" # Use EPSG:4326 as default which is WGS84 (World Geodetic System 1984).
    else:
        # Original logic: always use EPSG:32647 if any projection is found.
        logger.info(f"Found projected CRS codes: {dict(proj_counts)}. Using predefined EPSG:32647 for Thailand.")
        determined_epsg = "EPSG:32647" # Use EPSG:32647 as default which is WGS84 UTM Zone 47 North (UTM Zone 47N) Coverage Thailand.

    avg_x_res = sum(x_res_list) / len(x_res_list)
    avg_y_res = sum(y_res_list) / len(y_res_list)

    logger.info(f"Determined Target EPSG: {determined_epsg}, Avg XRes: {avg_x_res:.6f}, Avg YRes: {avg_y_res:.6f}")
    return determined_epsg, avg_x_res, avg_y_res

def main():
    root_dir = 'Raster_Resample'
    output_dir = 'Raster_Mosaic'
    os.makedirs(output_dir, exist_ok=True)

    raster_files = glob.glob(os.path.join(root_dir, '*.tif')) + glob.glob(os.path.join(root_dir, '*.tiff'))
    logger.info(f"Found {len(raster_files)} raster files to process")

    # Regex to extract date and orbit ID
    orbit_date_pattern = re.compile(r'_(\d{8})_R(\d{3})_')

    # Group files by (date, orbit ID)
    orbit_groups = defaultdict(list)

    for f in raster_files:
        basename = os.path.basename(f)
        match = orbit_date_pattern.search(os.path.basename(f))
        if match:
            date_str = match.group(1)  # e.g., 20250316
            orbit_id = f"R{match.group(2)}" # e.g., R061
            key = (date_str, orbit_id)
            orbit_groups[key].append(f)

    if not orbit_groups:
        logger.error("No valid orbit IDs and dates found in filenames.")
        return

    for (date_str, orbit_id), files in orbit_groups.items():
        logger.info(f"\nProcessing Orbit {orbit_id} on {date_str} with {len(files)} files...")

        orbit_dir = os.path.join(output_dir, orbit_id)
        os.makedirs(orbit_dir, exist_ok=True)
        
        # Output file: yyyymmdd_orbitID_Mosaic.tif
        final_output_path = os.path.join(output_dir, f"{date_str}_{orbit_id}_Mosaic.tif")

        target_epsg_str, x_res, y_res = analyze_rasters(files)

        if target_epsg_str is None:
            logger.error(f"Skipping Orbit {orbit_id} on {date_str} as mosaic parameters could not be determined from input rasters.")
            continue

        all_reprojected = []
        for i, raster_file in enumerate(files):
            base_name = os.path.basename(raster_file)
            reprojected_path = os.path.join(orbit_dir, f"reproj_{i}_{base_name}")

            logger.info(f"Attempting to reproject {base_name} to {target_epsg_str} with resolution {x_res:.6f}, {y_res:.6f}")

            # Pre-check if the source raster can be opened by GDAL
            src_ds_check = None
            try:
                src_ds_check = gdal.Open(raster_file)
                if src_ds_check is None:
                    logger.warning(f"Source raster {raster_file} cannot be opened by GDAL. Skipping reprojection for this file.")
                    continue # Skip to the next file
            except RuntimeError as e:
                logger.warning(f"GDAL RuntimeError while trying to open source raster {raster_file} for pre-check: {e}. Skipping.")
                continue
            finally:
                if src_ds_check:
                    src_ds_check = None # Dereference/close

            try:
                gdal.Warp(
                    reprojected_path,
                    raster_file,
                    options=gdal.WarpOptions(
                        dstSRS=target_epsg_str,
                        xRes=x_res,
                        yRes=y_res,
                        targetAlignedPixels=True,
                        resampleAlg='near',
                        srcNodata=0,
                        dstNodata=0,
                        creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES']
                    )
                )
                logger.info(f"Successfully reprojected {base_name} to {reprojected_path}")
                all_reprojected.append(reprojected_path)
            except RuntimeError as e:
                logger.warning(f"Failed to reproject {raster_file} to {reprojected_path}. GDAL Error: {e}. Skipping this file.")
                if os.path.exists(reprojected_path):
                    try:
                        os.remove(reprojected_path)
                        logger.debug(f"Removed partially created/failed reprojected file: {reprojected_path}")
                    except OSError as ose:
                        logger.warning(f"Could not remove partially created/failed reprojected file {reprojected_path}: {ose}")
                continue
            except Exception as e:
                logger.error(f"An unexpected error occurred while reprojecting {raster_file} to {reprojected_path}: {e}. Skipping this file.")
                if os.path.exists(reprojected_path):
                    try:
                        os.remove(reprojected_path)
                        logger.debug(f"Removed partially created/failed reprojected file due to unexpected error: {reprojected_path}")
                    except OSError as ose:
                        logger.warning(f"Could not remove partially created/failed reprojected file {reprojected_path} after unexpected error: {ose}")
                continue

        if not all_reprojected:
            logger.error(f"No rasters were successfully reprojected for Orbit {orbit_id} on {date_str}. Skipping mosaic creation for this group.")
            continue

        vrt_path = os.path.join(orbit_dir, 'aligned_mosaic.vrt')
        logger.info("Building VRT from reprojected files...")
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
        if vrt is None: # gdal.BuildVRT returns a Dataset object on success, or None on failure.
            logger.error(f"Failed to build VRT for {orbit_id}")
            continue
        vrt = None

        logger.info("Creating final mosaic...")
        gdal.Translate(
            final_output_path,
            vrt_path,
            options=gdal.TranslateOptions(
                format='GTiff',
                creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES', 'PREDICTOR=2']
            )
        )
        logger.info(f"Final mosaic saved to: {final_output_path}")

        if os.path.exists(final_output_path):
            logger.info("Cleaning up temporary files...")
            try:
                if os.path.exists(vrt_path):
                    os.remove(vrt_path)
                    logger.debug(f"Removed VRT: {vrt_path}")
                else:
                    logger.debug(f"VRT file {vrt_path} not found for cleanup.")
            except OSError as e:
                logger.warning(f"Failed to remove VRT: {vrt_path}. Error: {e}")

            for file_to_remove in all_reprojected:
                if os.path.exists(file_to_remove):
                    try:
                        os.remove(file_to_remove)
                        logger.debug(f"Removed temporary reprojected file: {file_to_remove}")
                    except OSError as e:
                        logger.warning(f"Failed to remove temporary reprojected file {file_to_remove}. Error: {e}")

        logger.info(f"Mosaic creation complete for orbit {orbit_id}!")

if __name__ == "__main__":
    main()
