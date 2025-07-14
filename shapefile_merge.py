import geopandas as gpd
import pandas as pd
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

# Define the directory containing your shapefiles
shapefile_dir = r'/path/to/shapefiles' 

# Create an empty list to store GeoDataFrames
gdfs = []

# Iterate through the directory and read shapefiles
# make code to handle with subfolder in the shapefile_dir
for root, dirs, files in os.walk(shapefile_dir):
    for filename in files:
        if filename.endswith('.shp'):
            filepath = os.path.join(root, filename)
            gdf = gpd.read_file(filepath)
            logging.info(f"Loaded shapefile: {filepath}")

            # Reproject the GeoDataFrame to WGS1984 (EPSG:4326)
            gdf = gdf.to_crs(epsg=4326)
            gdfs.append(gdf)
            logging.info(f"Reprojected shapefile to WGS1984: {filepath}")


logging.info(f"Found {len(gdfs)} shapefiles to merge.")

# Concatenate all GeoDataFrames into a single GeoDataFrame
merged_gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))

# Define the output path for the merged shapefile
output_path = os.path.join(shapefile_dir, 'Wildfire_merged_shapefile.shp')

# Export the merged GeoDataFrame to a new shapefile
merged_gdf.to_file(output_path)

print(f"Shapefiles merged successfully to: {output_path}")