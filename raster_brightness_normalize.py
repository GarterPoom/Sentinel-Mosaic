import os
import glob
import logging
import numpy as np
from collections import defaultdict
from osgeo import gdal, osr

# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Configure input and output directories/paths
root_dir = 'Raster_แม่น้ำโขง'
output_dir = 'Raster_Mosaic'
temp_dir = os.path.join('Raster_Normalized')
final_output_path = os.path.join(output_dir, '20250323_27_Mosaic.tif')

os.makedirs(output_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)

# Find all raster files
raster_files = glob.glob(os.path.join(root_dir, '*.tif')) + glob.glob(os.path.join(root_dir, '*.tiff'))
logger.info(f"Found {len(raster_files)} raster files to process")

def normalize_brightness(input_path, output_path):
    """Normalize raster pixel values to 0–255 range per band"""
    ds = gdal.Open(input_path)
    if ds is None:
        logger.warning(f"Could not open {input_path} for normalization.")
        return None

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(
        output_path,
        ds.RasterXSize,
        ds.RasterYSize,
        ds.RasterCount,
        gdal.GDT_Byte,
        options=['TILED=YES', 'COMPRESS=LZW']
    )
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())

    for band_index in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_index)
        arr = band.ReadAsArray().astype(np.float32)

        min_val = np.nanmin(arr)
        max_val = np.nanmax(arr)
        if max_val - min_val > 0:
            scaled = ((arr - min_val) / (max_val - min_val) * 255).astype(np.uint8)
        else:
            scaled = np.zeros_like(arr, dtype=np.uint8)

        out_ds.GetRasterBand(band_index).WriteArray(scaled)

    ds = None
    out_ds = None
    return output_path

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

    # Average resolution
    avg_x_res = sum(x_res_list) / len(x_res_list)
    avg_y_res = sum(y_res_list) / len(y_res_list)

    logger.info(f"Target EPSG: {target_epsg}, Avg Res: {avg_x_res}, {avg_y_res}")
    return target_epsg, avg_x_res, avg_y_res

# Analyze rasters
target_epsg, x_res, y_res = analyze_rasters(raster_files)

# Normalize and reproject all rasters
all_reprojected = []
for i, raster_file in enumerate(raster_files):
    base_name = os.path.basename(raster_file)
    normalized_path = os.path.join(temp_dir, f"norm_{i}_{base_name}")
    reprojected_path = os.path.join(output_dir, f"reproj_{i}_{base_name}")

    logger.info(f"Normalizing brightness of {base_name}")
    norm_result = normalize_brightness(raster_file, normalized_path)
    if not norm_result:
        continue

    logger.info(f"Reprojecting {base_name} to {target_epsg}")
    gdal.Warp(
        reprojected_path,
        norm_result,
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

# Cleanup
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