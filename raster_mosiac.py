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
        target_epsg = f"EPSG:{max(proj_counts.items(), key=lambda x: x[1])[0]}"

    # Average resolution - check if lists are empty to avoid ZeroDivisionError
    if not x_res_list or not y_res_list:
        logger.error("Could not determine average resolution from any input files.")
        return target_epsg, None, None # Indicate failure to calculate resolution
    avg_x_res = sum(x_res_list) / float(len(x_res_list)) # Use float for division
    avg_y_res = sum(y_res_list) / float(len(y_res_list)) # Use float for division

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}")
    return target_epsg, avg_x_res, avg_y_res

def main():
    # Configure input and output directories/paths
    root_dir = 'Raster_Resample'
    output_dir = 'Raster_Mosaic'
    final_output_path = os.path.join(output_dir, '20250221_25_Mosaic.tif')

    os.makedirs(output_dir, exist_ok=True)

    # Find all raster files
    raster_files = glob.glob(os.path.join(root_dir, '*.tif')) + glob.glob(os.path.join(root_dir, '*.tiff'))
    logger.info(f"Found {len(raster_files)} raster files to process")

    # Analyze rasters
    target_epsg, x_res, y_res = analyze_rasters(raster_files)
    if x_res is None or y_res is None:
        logger.error("Failed to determine average resolution from input rasters. Cannot proceed.")
        # Attempt to clean up output directory if it was created and is empty
        if not os.listdir(output_dir): os.rmdir(output_dir)
        exit(1)

    # Reproject all rasters
    all_reprojected = []
    for i, raster_file in enumerate(raster_files):
        base_name = os.path.basename(raster_file)
        reprojected_path = os.path.join(output_dir, f"reproj_{i}_{base_name}")

        try:
            logger.info(f"Reprojecting {base_name} to {target_epsg} with resolution {x_res}, {y_res} and aligned pixels")
            warp_options = gdal.WarpOptions(
                    dstSRS=target_epsg,
                    xRes=x_res,
                    yRes=y_res,
                    targetAlignedPixels=True,
                    resampleAlg='near',
                    srcNodata=0, # Consider making this configurable or auto-detected if possible
                    dstNodata=0, # Consider making this configurable
                    outputType=gdal.GDT_UInt16,
                    creationOptions=['TILED=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'],
                    # Adding error handling for GDAL internal errors
                    errorThreshold=0.0 # Default is 0.125, 0.0 means exact reprojection
                )
            ds = gdal.Warp(
                reprojected_path,
                raster_file,
                options=warp_options
            )
            if ds is None: # Should not happen if an exception is raised, but good for robustness
                logger.error(f"gdal.Warp failed for {raster_file} and returned None, but did not raise an exception.")
            else:
                all_reprojected.append(reprojected_path)
            ds = None # Release dataset
        except Exception as e:
            logger.error(f"Failed to reproject {raster_file}: {e}")
            # Optionally, include traceback for detailed debugging:
            # import traceback
            # logger.error(f"Traceback: {traceback.format_exc()}")
            continue # Skip to the next file

    if not all_reprojected:
        logger.error("No rasters were successfully reprojected!")
        exit(1)

    # Build VRT
    vrt_path = os.path.join(output_dir, 'aligned_mosaic.vrt')
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
        logger.error("Failed to build VRT")
        exit(1)
    vrt = None

    # Translate VRT to final GeoTIFF
    logger.info("Creating final mosaic...")
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
    logger.info(f"Final mosaic saved to: {final_output_path}")

    # Clean up
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

    logger.info("Mosaic creation complete!")

if __name__ == "__main__":
    main()
