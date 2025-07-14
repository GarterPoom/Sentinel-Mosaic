import geopandas as gpd
import os
import glob

# --- Configuration ---
# Root directory containing the first set of shapefiles to be processed.
# The script will search this directory recursively for all .shp files.
INPUT_ROOT_DIR = r"path/to/input/shapefiles"

# Path to the second shapefile that will be unioned with each file from the input directory.
# This file is loaded only once for efficiency.
BASE_SHAPEFILE_PATH = r"path/to/base/shapefile.shp"

# Root directory where the output shapefiles will be saved.
# The original directory structure from INPUT_ROOT_DIR will be preserved.
OUTPUT_ROOT_DIR = r"path/to/output/shapefiles"

# Suffix to add to the output filenames to distinguish them from the originals.
OUTPUT_FILENAME_SUFFIX = "_unioned"

# --- Main Script ---

def main():
    """
    Finds all shapefiles in an input directory, performs a union operation
    with a base shapefile, and saves the results to an output directory,
    preserving the original file structure.
    """
    # 1. Find all input shapefiles using glob
    search_pattern = os.path.join(INPUT_ROOT_DIR, '**', '*.shp')
    shp1_paths = glob.glob(search_pattern, recursive=True)

    if not shp1_paths:
        print(f"No .shp files found in '{INPUT_ROOT_DIR}'. Exiting.")
        return

    print(f"Found {len(shp1_paths)} shapefiles to process.")
    

    # 2. Load the base shapefile once
    print(f"Loading base shapefile: {BASE_SHAPEFILE_PATH}")
    try:
        gdf2 = gpd.read_file(BASE_SHAPEFILE_PATH)
    except Exception as e:
        print(f"Error loading base shapefile '{BASE_SHAPEFILE_PATH}': {e}. Exiting.")
        return

    # 3. Loop through each input shapefile and process it
    for shp1_path in shp1_paths:
        print(f"\n--- Processing: {shp1_path} ---")
        try:
            # Load the current shapefile
            gdf1 = gpd.read_file(shp1_path)

            # --- Ensure CRS match ---
            # We reproject the current shapefile (gdf1) to match the base shapefile's CRS (gdf2).
            if gdf1.crs != gdf2.crs:
                print(f"CRS mismatch. Reprojecting '{os.path.basename(shp1_path)}' to match base CRS (EPSG:{gdf2.crs.to_epsg()}).")
                gdf1 = gdf1.to_crs(gdf2.crs)

            # --- Perform Overlay Union ---                                                                       
            print("Performing overlay union...")
            union_gdf = gpd.overlay(gdf1, gdf2, how='union', keep_geom_type=True)

            # --- Prepare dynamic output path ---
            relative_path = os.path.relpath(shp1_path, INPUT_ROOT_DIR)
            relative_dir = os.path.dirname(relative_path)
            filename_base, ext = os.path.splitext(os.path.basename(relative_path))

            output_filename = f"{filename_base}{OUTPUT_FILENAME_SUFFIX}{ext}"
            output_dir = os.path.join(OUTPUT_ROOT_DIR, relative_dir)
            output_path = os.path.join(output_dir, output_filename)

            os.makedirs(output_dir, exist_ok=True)

            # --- Save the result ---
            print(f"Saving unioned shapefile to: {output_path}")
            union_gdf.to_file(output_path)

        except Exception as e:
            print(f"!!! ERROR processing {shp1_path}: {e}")
            continue

    print("\nDone. All files processed.")

if __name__ == "__main__":
    main()