from raster_multiply import main as multiply_raster_pixels
from cliping_raster import main as clip_raster
from polygonise_GeoJSON import main as polygonize_raster

def main():
    multiply_raster_pixels()
    clip_raster()
    polygonize_raster()

if __name__ == "__main__":
    main()                                 