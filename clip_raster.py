import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.shutil import copy as rio_copy
import os
import tempfile
import glob
import logging
from shapely.geometry import box

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def clip_rasters_by_subdistrict_hierarchy(
    folder_path, shapefile_path,
    province_attr, district_attr, subdistrict_attr,
    output_dir):
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load shapefile
    gdf = gpd.read_file(shapefile_path)
    
    # Ensure all needed attributes exist
    for attr in [province_attr, district_attr, subdistrict_attr]:
        if attr not in gdf.columns:
            logging.error(f"Attribute '{attr}' not found in shapefile.")
            return
    
    # Get unique subdistricts
    unique_subdistricts = gdf[subdistrict_attr].dropna().unique()
    logging.info(f"Found {len(unique_subdistricts)} unique subdistricts in attribute '{subdistrict_attr}'.")
    
    raster_files = glob.glob(os.path.join(folder_path, '*.tif'))
    if not raster_files:
        logging.warning("No .tif files found in the folder.")
        return
    
    logging.info(f"Found {len(raster_files)} raster files. Starting batch clipping...")

    for raster_path in raster_files:
        raster_name = os.path.basename(raster_path)
        try:
            with rasterio.open(raster_path) as src:
                raster_crs = src.crs
                raster_bounds = src.bounds
                raster_bbox = box(*raster_bounds)
            
            gdf_proj = gdf.to_crs(raster_crs)
            
            for subdistrict in unique_subdistricts:
                filtered_gdf = gdf_proj[gdf_proj[subdistrict_attr] == subdistrict]
                if filtered_gdf.empty:
                    logging.warning(f"No matching features for subdistrict '{subdistrict}' in {raster_name}. Skipping.")
                    continue
                
                if not filtered_gdf.intersects(raster_bbox).any():
                    logging.warning(f"Features for subdistrict '{subdistrict}' do not intersect raster {raster_name}. Skipping.")
                    continue
                
                # Extract province and district values (assuming single unique values in filtered_gdf)
                provinces = filtered_gdf[province_attr].unique()
                districts = filtered_gdf[district_attr].unique()
                
                if len(provinces) != 1 or len(districts) != 1:
                    logging.warning(f"Multiple provinces or districts found for subdistrict '{subdistrict}', skipping.")
                    continue
                
                province = str(provinces[0])
                district = str(districts[0])
                
                geom_json = [geom.__geo_interface__ for geom in filtered_gdf.geometry]
                
                with rasterio.Env(GDAL_CACHEMAX=512):
                    with rasterio.open(raster_path) as src:
                        out_image, out_transform = mask(src, geom_json, crop=True)
                        out_meta = src.meta.copy()
                
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                    "tiled": True,
                    "blockxsize": 256,
                    "blockysize": 256,
                    "compress": "lzw"
                })
                
                # Construct nested output directory
                value_dir = os.path.join(output_dir, province, district, subdistrict)
                os.makedirs(value_dir, exist_ok=True)
                
                output_path = os.path.join(value_dir, f"clipped_{raster_name}")
                
                with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
                    temp_output = tmp.name
                
                with rasterio.open(temp_output, "w", **out_meta) as dest:
                    dest.write(out_image)
                
                with rasterio.open(temp_output, 'r+') as dataset:
                    dataset.build_overviews([2, 4, 8, 16, 32], Resampling.average)
                    dataset.update_tags(ns='rio_overview', resampling='average')
                
                rio_copy(temp_output, output_path, copy_src_overviews=True)
                os.remove(temp_output)
                
                logging.info(f"Saved: {output_path}")
        
        except Exception as e:
            logging.error(f"Error processing {raster_name}: {e}")
    
    logging.info("Batch clipping by subdistrict hierarchy complete.")


def main():
    clip_rasters_by_subdistrict_hierarchy(
        folder_path="LANDSAT_9/", # Replace with your raster folder path
        shapefile_path="Thailand/Thailand - Subnational Administrative Boundaries.shp", # Replace with your shapefile path
        province_attr="PV_TN",       # Replace with actual province attribute name
        district_attr="AP_TN", # Replace with actual district attribute name
        subdistrict_attr="TB_TN", # Replace with actual subdistrict attribute name
        output_dir="Clipped_Rasters"
    )

if __name__ == "__main__":
    main()
