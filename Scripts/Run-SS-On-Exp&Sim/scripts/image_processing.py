import cv2
import numpy as np
from PIL import Image
import BacteriaFeatures
from matplotlib import pyplot as plt
import alphashape
from shapely.geometry import Polygon, MultiPolygon
from matplotlib.patches import Polygon as plg


def save_fig(contours, fig_export_path, fig_name):

    """
    This function takes in a 2D NumPy array, saves it as an image file with .png format.
    @param contours numpy representing the image to be saved.
    @param fig_export_path string a string representing the file path where the saved image should be exported to.
    @param fig_name  string  representing the name of the saved image.
    """

    contours_from_array = Image.fromarray(contours)
    contours_from_array.save(fig_export_path + fig_name + ".png")


def interpolate_points(start, end, num_points=100):
    x_values = np.linspace(start[0], end[0], num_points)
    y_values = np.linspace(start[1], end[1], num_points)
    return list(zip(x_values, y_values))


def fill_contours(cs, um_pixel_ratio=0.144, margin=10, shape=0.048, fig_export_path=None, fig_name=None):
    """
    This function takes in a dictionary of cell states (cs),
    and generates a 2D NumPy array representing the filled envelope of the micro-colony.
    The envelope is filled using a morphological closing operation on a grayscale image generated from
    the input cell state dictionary.
    @param cs cellStates dict representing the cell states of a colony / micro-colony.
    @param um_pixel_ratio float representing the ratio of micrometers to pixels in the generated image.
    @param margin int pixel representing the margin to add to the generated image.
    @param shape float shape parameter
    """

    # bacteria endpoints
    bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y = \
        BacteriaFeatures.find_bacteria_endpoints(cs)

    bacteria_endpoint_x = bacteria_end_point1_x + bacteria_end_point2_x
    bacteria_endpoint_y = bacteria_end_point1_y + bacteria_end_point2_y
    # convert to pixel
    bacteria_endpoint_x = BacteriaFeatures.convert_um_pixel(bacteria_endpoint_x, um_pixel_ratio)
    bacteria_endpoint_y = BacteriaFeatures.convert_um_pixel(bacteria_endpoint_y, um_pixel_ratio)

    points = [(x, y) for x, y in zip(bacteria_endpoint_x, bacteria_endpoint_y)]

    fig, ax = plt.subplots()

    alpha_shape = alphashape.alphashape(points, shape)
    if isinstance(alpha_shape, MultiPolygon) or alpha_shape.boundary is None:
        shape = round(shape - 0.05, 3)
        alpha_shape = alphashape.alphashape(points, shape)

    coords = alpha_shape.boundary.coords.xy

    vert = [(x, y) for x, y in zip(coords[0], coords[1])]

    # Create a Polygon object
    # print(points)
    for point in points:
        plt.scatter(point[0], point[1])
    polygon_plt = plg(vert, facecolor='none', edgecolor='red')
    ax.add_patch(polygon_plt)

    # Set the limits of the axis
    ax.set_xlim([min(coords[0]) - margin, max(coords[0]) + margin])
    ax.set_ylim([min(coords[1]) - margin, max(coords[1]) + margin])

    polygon_shapely = Polygon(vert)
    exterior_coords = polygon_shapely.exterior.coords

    # Generate all boundary points
    boundary_points = []
    for i in range(len(exterior_coords) - 1):
        start = exterior_coords[i]
        end = exterior_coords[i + 1]
        points = interpolate_points(start, end, num_points=100)
        boundary_points.extend(points)

    x, y = zip(*boundary_points)
    plt.scatter(x, y, color='blue', s=2, alpha=0.5)

    plt.gca().invert_yaxis()
    # Remove the x-axis and y-axis
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    if fig_export_path is not None:
        plt.savefig(fig_export_path + fig_name + ".png")

    return boundary_points


def find_external_contours(gray, fig_export_path, fig_name):

    """
    This function takes in a 2D NumPy array representing the filled envelope of a colony / micro-colony,
    and returns a 2D NumPy array representing the boundary of the colony / micro-colony.
    The boundary is found by performing a contour finding operation on the input image.
    @param filled_contours numpy representing the filled envelope of a micro-colony.
    @param fig_export_path str representing the file path where the saved image should be exported to (optional)
    @param fig_name str representing the name of the saved image.
"""

    # Threshold the grayscale image
    ret, thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY)

    # Find the contours in the binary image
    external_contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Draw the contours on a black image
    boundary = np.zeros(gray.shape, dtype=np.uint8)
    cv2.drawContours(boundary, external_contours, -1, 255, 1)

    if fig_export_path:
        save_fig(boundary, fig_export_path, fig_name)

    return boundary


def pixel_to_coordinate(contours):

    """
    This function takes in a list of 2D NumPy arrays representing contours,
    and returns a tuple of two lists representing the x and y coordinates of the white pixels in the input contours.
    @param contours numpy a list of 2D NumPy arrays representing contours.
    """

    x_coordinates = []
    y_coordinates = []
    for x_index, contour in enumerate(contours):
        white_pixel = np.where(contour == 255)
        if white_pixel[0].size > 0:
            x_coordinates.extend([x_index] * white_pixel[0].size)
            y_coordinates.extend(list(white_pixel[0]))

    return x_coordinates, y_coordinates
    
