import cv2
from PIL import Image, ImageDraw
import skimage.morphology, skimage.measure
import imageio
import numpy as np
import matplotlib.pyplot as plt
from FractalDimension.MBdimension import fractal_dimension
import cell_data
from density_calculation import get_cell_data_to_draw_image_bw, get_image_dimensions, draw_image_bw


def save_fig(cells, fig_export_path, fig_name):
    """
    @param cells cellStates dict
    @param hull object convex hull
    @param points numpy array bacteria end points
    @param  fig_export_path directory to export images
    @param fig_name image name

    """
    # vertexes
    cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations = \
        get_cell_data_to_draw_image_bw(cells, um_pixel_ratio=0.144)

    # Create image dimensions
    img_dimensions = get_image_dimensions(cell_centers_x, cell_centers_y)

    # Draw black and white image
    bw_img = draw_image_bw(img_dimensions, cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations)

    # part of `calculate_colony_density` function

    # Perform morphological closing to get rid of any gaps
    img_close = skimage.morphology.closing(bw_img, footprint=np.ones((7, 7), dtype='uint8'))

    # Get contours
    contours, hierarchy = cv2.findContours(np.copy(img_close), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Fill contours so that the image is a filled envelope of the microcolony
    img_close = cv2.drawContours(img_close, contours, -1, 255, cv2.FILLED)

    img_close_fromarray = Image.fromarray(img_close)
    img_close_fromarray.save(fig_export_path + fig_name + "_fill.png")


def calc_fractal_dimension(cells, fig_export_path, fig_name):
    """
    goal: calculation of fractal dimension (Minkowski–Bouligand dimension)

    @param  fig_export_path directory to export images
    @param fig_name image name
    @return fractal dimension
    """
    save_fig(cells, fig_export_path, fig_name)
    # Import the image in greyscale
    greyscale_img = imageio.v2.imread(fig_export_path + fig_name + '_fill.png', as_gray="True") / 255.0
    return fractal_dimension(greyscale_img)
