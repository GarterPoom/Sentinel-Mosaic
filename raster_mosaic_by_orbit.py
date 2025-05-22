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
        ds = gdal.Open(f)
        if ds is None:
            logger.warning(f"Cannot open {f} for analysis")
            continue

        srs = osr.SpatialReference()
        srs.ImportFromWkt(ds.GetProjection())
        if srs.IsProjected():
            epsg = srs.GetAuthorityCode(None)
            proj_counts[epsg] += 1

        gt = ds.GetGeoTransform()
        x_res_list.append(abs(gt[1]))
        y_res_list.append(abs(gt[5]))
        ds = None

    if not proj_counts:
        logger.warning("No projections found. Defaulting to EPSG:4326")
        target_epsg = "EPSG:4326" # Use EPSG:4326 as default which is WGS84 (World Geodetic System 1984).
    else:
        target_epsg = "EPSG:32647" # Use EPSG:32647 as default which is WGS84 UTM Zone 47 North (UTM Zone 47N) Coverage Thailand.

    avg_x_res = sum(x_res_list) / len(x_res_list)
    avg_y_res = sum(y_res_list) / len(y_res_list)

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}")
    return target_epsg, avg_x_res, avg_y_res

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

        target_epsg, x_res, y_res = analyze_rasters(files)

        all_reprojected = []
        for i, raster_file in enumerate(files):
            base_name = os.path.basename(raster_file)
            reprojected_path = os.path.join(orbit_dir, f"reproj_{i}_{base_name}")

            logger.info(f"Reprojecting {base_name} to {target_epsg} with aligned pixels")
            gdal.Warp(
                reprojected_path,
                raster_file,
                options=gdal.WarpOptions(
                    dstSRS=target_epsg,
                    xRes=x_res,
                    yRes=y_res,
                    targetAlignedPixels=True,
                    resampleAlg='near',
                    srcNodata=0,
                    dstNodata=0,
                    creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES']
                )
            )
            all_reprojected.append(reprojected_path)

        if not all_reprojected:
            logger.error(f"No rasters were successfully reprojected for {orbit_id}!")
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
        if vrt is None:
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
                os.remove(vrt_path)
            except Exception as e:
                logger.warning(f"Failed to remove VRT: {e}")

            for file in all_reprojected:
                try:
                    os.remove(file)
                except Exception as e:
                    logger.warning(f"Failed to remove {file}: {e}")

        logger.info(f"Mosaic creation complete for orbit {orbit_id}!")

if __name__ == "__main__":
    main()
