"""
draw_image_bw and calculate_colony_density are based on code from Sohaib Nadeem and Jonathan Chalaturnik:
https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/
Modifications by Aaron Yip:
-Repurposed to work with CellModeller data
-Fixed bugs with colony area calculations
-Added documentation
Instructions:
Call the main() function, giving a cellStates dict as an argument.
Optional: include path to a directory to export plots
"""
# Installed modules
import cv2
# Standard modules
import numpy as np
import pandas as pd
import skimage.measure
import skimage.morphology
from PIL import Image, ImageDraw


def get_cell_data_to_draw_image_bw(cells, um_pixel_ratio):
    """
    Convert CellModeller data into a format suitable for drawing black/white images

    @param  cells           cellStates dict
    @return ...             see variable names in return statement
    """
    # Add pixels to img border to avoid cutting off cells
    add_pixels_to_img_border = 5

    # Initialize storage variables
    cell_centers_x = cells['x_center'].to_numpy()
    cell_centers_y = cells['y_center'].to_numpy()
    cell_lengths = cells['major'].to_numpy()
    cell_radii = cells['minor'].to_numpy()
    cell_orientations = cells['orientation'].to_numpy()

    # convert um to pixel
    cell_centers_x = cell_centers_x / um_pixel_ratio
    cell_centers_y = cell_centers_y / um_pixel_ratio
    cell_lengths = cell_lengths / um_pixel_ratio
    cell_radii = cell_radii / um_pixel_ratio

    # Ensure all cells have positive coordinates for drawing
    min_x = min(cell_centers_x) - add_pixels_to_img_border
    min_y = min(cell_centers_y) - add_pixels_to_img_border
    cell_centers_x = cell_centers_x + abs(min_x)
    cell_centers_y = cell_centers_y + abs(min_y)

    return cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations


def get_image_dimensions(cell_centers_x, cell_centers_y):
    """
    Obtains dimensions for an image in pixels based on colony size

    @param  cell_centers_x  list or nparray of x-coordinates for cell centroids
    @param  cell_centers_y  list or nparray of y-coordinates for cell centroids
    @return img_dimensions  (width, height)
    """
    min_x = min(cell_centers_x) - 5
    max_x = max(cell_centers_x) + 5
    min_y = min(cell_centers_y) - 5
    max_y = max(cell_centers_y) + 5
    width = abs(max_x) + abs(min_x)
    height = abs(max_y) + abs(min_y)
    img_dimensions = (round(width), round(height))

    return img_dimensions


def draw_cell(draw, center_x, center_y, length, radius, orientation, fill, outline):
    """
    Draws a cell on a black/white canvas
    Based on code of Sohaib Nadeem:
    https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/image_drawing.py

    @param  draw        pillow Image
    @param  center_x    cell centroid x-coordinate
    @param  center_y    cell centroid y-coordinate
    @param  length      cell length (does not include ends of cell)
    @param  radius      cell radius
    @param  orientation cell angle in radians
    @param  fill        intensity of the cell fill (ranges from 0 to 255 for bw image)
    @param  outline     intensity of the cell outline (ranges from 0 to 255 for bw image)
    """

    # Adjust orientations to image coordinate system
    img_orientation = -orientation  # flip y-coordinates because (0,0) is in top left corner
    img_orientation_norm = img_orientation + np.pi / 2

    # Calculate lengths
    half_length_along_axis_x = length / 2 * np.cos(img_orientation)
    half_length_along_axis_y = length / 2 * np.sin(img_orientation)
    radius_perpendicular_to_axis_x = radius * np.cos(img_orientation_norm)
    radius_perpendicular_to_axis_y = radius * np.sin(img_orientation_norm)

    # Draw rectangle
    p1 = (center_x + half_length_along_axis_x + radius_perpendicular_to_axis_x,
          center_y + half_length_along_axis_y + radius_perpendicular_to_axis_y)
    p2 = (center_x + half_length_along_axis_x - radius_perpendicular_to_axis_x,
          center_y + half_length_along_axis_y - radius_perpendicular_to_axis_y)
    p3 = (center_x - half_length_along_axis_x - radius_perpendicular_to_axis_x,
          center_y - half_length_along_axis_y - radius_perpendicular_to_axis_y)
    p4 = (center_x - half_length_along_axis_x + radius_perpendicular_to_axis_x,
          center_y - half_length_along_axis_y + radius_perpendicular_to_axis_y)
    draw.polygon([p1, p2, p3, p4], fill=fill, outline=outline)

    # Draw ends of cell
    p5 = (center_x + half_length_along_axis_x - radius, center_y + half_length_along_axis_y - radius)
    p6 = (center_x + half_length_along_axis_x + radius, center_y + half_length_along_axis_y + radius)
    p7 = (center_x - half_length_along_axis_x - radius, center_y - half_length_along_axis_y - radius)
    p8 = (center_x - half_length_along_axis_x + radius, center_y - half_length_along_axis_y + radius)
    end_1 = img_orientation_norm * 180 / np.pi
    start_1 = end_1 - 180
    draw.pieslice([p5, p6], start=start_1, end=end_1, fill=fill, outline=outline)  # start and end angle in degrees
    draw.pieslice([p7, p8], start=end_1, end=start_1, fill=fill, outline=outline)


def draw_image_bw(img_dimensions, cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations):
    """
    Draw a black and white image from cell data (cells are white, background is black)

    Any cell dimension and position must be given in pixels

    @param  img_dimensions      dimensions of the world in pixels (width, height)
    @param  cell_centers_x      list or nparray of x-coordinates for cell centroids
    @param  cell_centers_y      list or nparray of y-coordinates for cell centroids
    @param  cell_lengths        list or nparray of cell lengths
    @param  cell_radii          list or nparray of cell radii
    @param  cell_orientations   list or nparray of cell orientations
    @return img                 black and white image of cells
    """
    cell_count = len(cell_centers_x)
    img = Image.new(mode='L', size=img_dimensions, color=0)
    draw = ImageDraw.Draw(img)
    fill = 255
    outline = None

    for i in range(cell_count):
        draw_cell(draw, cell_centers_x[i], cell_centers_y[i], cell_lengths[i], cell_radii[i], cell_orientations[i],
                  fill, outline)

    return img


def calculate_colony_density(img):
    """
    Based on code of J. Chalaturnik and S. Nadeem:
    https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/img_processing.py
    Calculates the fraction of the colony occupied by cells.
    @param  img             black/white pillow image of a colony
    @return density         colony density parameter (1 = colony is completely filled in)
    """
    cell_area = np.count_nonzero(np.array(img) == 255)

    # Perform morphological closing to get rid of any gaps
    img_close = skimage.morphology.closing(img, footprint=np.ones((7, 7), dtype='uint8'))

    # Get contours
    contours, hierarchy = cv2.findContours(np.copy(img_close), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Fill contours so that the image is a filled envelope of the micro-colony
    img_close = cv2.drawContours(img_close, contours, -1, 255, cv2.FILLED)

    # Calculate filled area of colony
    img_props = pd.DataFrame(skimage.measure.regionprops_table(img_close, properties=["area"]))
    filled_contour_area = img_props["area"][0]

    # Calculate fraction of area occupied by cells
    density = cell_area / filled_contour_area

    return density


def calc_density(cells_df, um_pixel_ratio):
    """
    The main function for calculating colony density

    @param  cells_df     dataframe     bacteria features value
    @return density         colony density parameter (1 = colony is completely filled in)
    """
    # Obtain cell data in format suitable for creating black/white images
    cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations = \
        get_cell_data_to_draw_image_bw(cells_df, um_pixel_ratio)

    # Create image dimensions
    img_dimensions = get_image_dimensions(cell_centers_x, cell_centers_y)

    # Draw black and white image
    bw_img = draw_image_bw(img_dimensions, cell_centers_x, cell_centers_y, cell_lengths, cell_radii, cell_orientations)

    # compute colony density
    density = calculate_colony_density(bw_img)

    return density
