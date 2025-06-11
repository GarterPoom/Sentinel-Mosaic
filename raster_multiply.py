import os
import rasterio
import numpy as np
import glob # Import the glob module
from rasterio.enums import Resampling # For overviews

def multiply_raster_pixels_in_directory(directory, factor, output_suffix="_multiplied", output_directory=''):
    """
    Multiplies pixel values of all raster files in a given directory by a factor.

    Args:
        directory (str): The path to the directory containing raster files.
        factor (float): The multiplication factor.
        output_suffix (str): Suffix to append to the processed filenames.
        output_directory (str): The directory where processed files will be saved.
    Returns:
        None
    """
    # Store extensions in lowercase for easier and case-insensitive comparison
    raster_extensions_lower = {ext.lower() for ext in ('.tif', '.tiff', '.img', '.jp2', '.envi', '.hdr', '.dat')}
    processed_files_count = 0
    failed_files_count = 0

    print(f"Scanning directory: {directory}")
    print(f"Multiplying pixel values by: {factor}")

    # Use glob to find all entries in the directory
    all_entries_pattern = os.path.join(directory, '*')
    all_entries = glob.glob(all_entries_pattern)

    files_to_process = []
    for filepath in all_entries:
        # Check if it's a file and if its extension matches one of the raster extensions
        if os.path.isfile(filepath):
            _, ext = os.path.splitext(filepath)
            if ext.lower() in raster_extensions_lower:
                files_to_process.append(filepath)

    for filepath in files_to_process: # Iterate over the filtered list of raster files
        filename = os.path.basename(filepath) # Get the filename for logging
        try:
            with rasterio.open(filepath) as src:
                print(f"Processing {filename}...")
                    
                # Read raster data
                data = src.read() # Reads all bands
                profile = src.profile # Get metadata

                # Ensure data is float for multiplication to avoid type issues
                # and to correctly store the result of multiplication by a float.
                if data.dtype != np.float32 and data.dtype != np.float64:
                    # Promote to float32 if not already float
                    # Using float32 is often a good balance of precision and file size
                    data_float = data.astype(np.float32)
                else:
                    data_float = data

                # This block should be outside the if/else for data type conversion
                multiplied_data = data_float * factor

                # Update profile for the output file
                profile.update(dtype=multiplied_data.dtype.name) # e.g., 'float32'

                # Construct output filename
                name, original_ext = os.path.splitext(filename) # Use the original file extension
                output_filename = f"{name}{output_suffix}{original_ext}"
                output_filepath = os.path.join(output_directory, output_filename)

                # Write the modified data to a new file
                with rasterio.open(output_filepath, 'w', **profile) as dst:
                    dst.write(multiplied_data)

                # Build overviews for the newly created raster
                build_overviews_for_raster(output_filepath)

                print(f"Successfully processed. Output saved as: {output_filename}")
                processed_files_count += 1

        except Exception as e:
            print(f"Error processing {filename}: {e}")
            failed_files_count += 1

    print("\n--- Processing Summary ---")
    print(f"Successfully processed files: {processed_files_count}")
    print(f"Failed/skipped files: {failed_files_count}")

def build_overviews_for_raster(filepath, levels=[2, 4, 8, 16, 32], resampling_method=Resampling.nearest):
    """
    Builds raster pyramid overviews for a given raster file using rasterio.

    Args:
        filepath (str): Path to the raster file.
        levels (list): List of integers representing the downsampling factors
                       for each overview level.
        resampling_method (rasterio.enums.Resampling): Resampling method to use.
    """
    try:
        print(f"Building overviews for {os.path.basename(filepath)}...")
        with rasterio.open(filepath, 'r+') as dst:
            dst.build_overviews(levels, resampling_method)
            # Optionally, update tags if needed, though build_overviews often handles this.
            # dst.update_tags(ns='rio_overview', resampling=resampling_method.name)
        print(f"Successfully built overviews for {os.path.basename(filepath)}.")
    except Exception as e:
        print(f"Error building overviews for {os.path.basename(filepath)}: {e}")


def main():

    # The script will process files in the directory where it is located.
    # You can change this to a specific path if needed.
    root_directory = r'LANDSAT_9'

    # Define the path for the output directory
    output_dir_path = 'LANDSAT_9_Multiplied'

    # Ensure the output directory exists
    # os.makedirs will create the directory if it doesn't exist
    # and do nothing if it already exists (due to exist_ok=True).
    os.makedirs(output_dir_path, exist_ok=True)

    multiplication_factor = 0.47

    # Call the function with the correct arguments
    multiply_raster_pixels_in_directory(
        directory=root_directory,
        factor=multiplication_factor,
        output_directory=output_dir_path
    )

if __name__ == "__main__":
    main()