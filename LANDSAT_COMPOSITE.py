import os
import glob
import rasterio
import numpy as np

# --- CONFIG ---
input_dir      = r'LANDSAT_9'
output_dir     = r'Raster_Composite'
band_indices   = [6, 4, 3]          # Band 7 (R), Band 5 (G), Band 4 (B)
output_suffix  = '_Com754.tif'  # Output file suffix
overwrite_ok   = False              # Set to True to overwrite existing files

# -------------------------------------------------------------
def composite_by_index(src_path, band_indices):
    with rasterio.open(src_path) as src:
        if src.count < max(band_indices):
            raise ValueError(f"{os.path.basename(src_path)} has only {src.count} bands, "
                             f"but requested bands {band_indices}")

        arrays = [src.read(i) for i in band_indices]
        meta   = src.meta.copy()
        meta.update(count=len(arrays), dtype=arrays[0].dtype)

    return np.stack(arrays), meta

def export_composite(array3d, meta, dst_path, overwrite_ok=False):
    if (not overwrite_ok) and os.path.exists(dst_path):
        raise FileExistsError(f'{dst_path} already exists (set overwrite_ok=True to replace).')

    os.makedirs(os.path.dirname(dst_path), exist_ok=True)

    with rasterio.open(dst_path, 'w', **meta) as dst:
        for i in range(array3d.shape[0]):
            dst.write(array3d[i], i + 1)

# -------------------------------------------------------------
def main():
    tif_paths = sorted(glob.glob(os.path.join(input_dir, '*.tif')))
    if not tif_paths:
        raise FileNotFoundError(f'No *.tif files found in {input_dir}')

    for tif in tif_paths:
        try:
            composite, meta = composite_by_index(tif, band_indices)

            base        = os.path.splitext(os.path.basename(tif))[0]
            out_name    = f'{base}{output_suffix}'
            out_path    = os.path.join(output_dir, out_name)

            export_composite(composite, meta, out_path, overwrite_ok)
            print(f'✓ Exported {out_path}')

        except ValueError as ve:
            print(f'⚠️ Skipped {os.path.basename(tif)}: {ve}')
        except Exception as e:
            print(f'❌ Error with {os.path.basename(tif)}: {e}')

    print('Done.')

if __name__ == '__main__':
    main()
