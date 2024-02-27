from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
import alphashape
from shapely.geometry import MultiPolygon
#from BacteriaFeatures import convert_um_pixel, find_bacteria_endpoints
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.BacteriaFeatures import convert_um_pixel, find_bacteria_endpoints


def cal_perimeter(ordered_points):
    """
    Calculates the perimeter of some ordered list of points
    @param ordered_points numpy array
    """
    perimeter_length = 0

    for it in range(len(ordered_points)):
        if it == 0:
            continue
        pre = ordered_points[it - 1]
        now = ordered_points[it]
        perimeter_length = perimeter_length + np.sqrt((pre[0] - now[0]) ** 2 + (pre[1] - now[1]) ** 2)

    perimeter_length = perimeter_length + np.sqrt((ordered_points[0][0] - ordered_points[-1][0]) ** 2 +
                                                  (ordered_points[0][1] - ordered_points[-1][1]) ** 2)

    return perimeter_length


def fit_alpha_shape(cs, um_pixel_ratio=0.144, shape=0.048):
    """
    This function takes in a dictionary of cell states (cs),
    and generates a 2D NumPy array representing the filled envelope of the micro-colony.
    The envelope is filled using a morphological closing operation on a grayscale image generated from
    the input cell state dictionary.
    @param cs cellStates dict representing the cell states of a colony / micro-colony.
    @param um_pixel_ratio float representing the ratio of micrometers to pixels in the generated image.
    @param shape float shape parameter
    """
    # bacteria endpoints
    bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y = \
        find_bacteria_endpoints(cs)

    bacteria_endpoint_x = bacteria_end_point1_x + bacteria_end_point2_x
    bacteria_endpoint_y = bacteria_end_point1_y + bacteria_end_point2_y
    # convert to pixel
    bacteria_endpoint_x = convert_um_pixel(bacteria_endpoint_x, um_pixel_ratio)
    bacteria_endpoint_y = convert_um_pixel(bacteria_endpoint_y, um_pixel_ratio)

    points = [(x, y) for x, y in zip(bacteria_endpoint_x, bacteria_endpoint_y)]

    alpha_shape = alphashape.alphashape(points, shape)

    return alpha_shape


def optimal_alpha(cs, um_pixel_ratio):
    min_area = None
    best_shape = None
    best_boundary = None
    shape_values = np.arange(0.001, 0.1, 0.001)
    area_list = []
    for shape_val in shape_values:
        fitted_alpha_shape = fit_alpha_shape(cs, um_pixel_ratio=um_pixel_ratio, shape=shape_val)
        if not isinstance(fitted_alpha_shape, MultiPolygon) and fitted_alpha_shape.boundary is not None:
            coords = fitted_alpha_shape.boundary.coords.xy
            boundary_points = np.array([(x, y) for x, y in zip(coords[0], coords[1])])
            area = polygon_area(boundary_points)
            area_list.append(area)
            if not min_area or area < min_area:
                min_area = area
                best_shape = shape_val
                best_boundary = boundary_points
        else:
            area_list.append(-1)

    return best_boundary, best_shape


def polygon_area(coords):
    """
    This function calculate the area formed by a set of coordinates
    """
    n = len(coords)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += coords[i][0] * coords[j][1]
        area -= coords[j][0] * coords[i][1]
    area = abs(area) / 2.0
    return area


def fit_convex_hull(cs, um_pixel_ratio):
    """
    This function fits a convex hull around the bacterial colony
    """
    # bacteria endpoints
    bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y = \
        find_bacteria_endpoints(cs)

    bacteria_endpoint_x = bacteria_end_point1_x + bacteria_end_point2_x
    bacteria_endpoint_y = bacteria_end_point1_y + bacteria_end_point2_y
    # convert to pixel
    bacteria_endpoint_x = convert_um_pixel(bacteria_endpoint_x, um_pixel_ratio)
    bacteria_endpoint_y = convert_um_pixel(bacteria_endpoint_y, um_pixel_ratio)

    points = list(zip(bacteria_endpoint_x, bacteria_endpoint_y))
    points = np.array(points)
    hull = ConvexHull(points)

    return hull, points


def cal_convexity(cs, um_pixel_ratio=0.144, fig_export_path=None, fig_name=None):
    """
    Calculates the ratio of the convex hull perimeter and perimeter of the fitted colony boundary
    @param cs: dictionary, cell status
    @param fig_export_path: string
    @param fig_name: string, the name of the scene
    @param um_pixel_ratio: float, um to pixel conversion constant
    """

    # fit convex hull
    fitted_hull, coordinates = fit_convex_hull(cs, um_pixel_ratio)
    # calculated the length of the convex perimeter
    convex_hull_perimeter = cal_perimeter(coordinates[fitted_hull.vertices])

    # find the best alphashape value
    boundary, best_shape = optimal_alpha(cs, um_pixel_ratio)
    contour_perimeter = cal_perimeter(boundary)

    # plot
    if fig_export_path:
        # plot all endpoints of the colony
        plt.plot(coordinates[:, 0], coordinates[:, 1], 'o')

        # Plot convex hull:
        plt.plot(coordinates[fitted_hull.vertices, 0], coordinates[fitted_hull.vertices, 1], 'r--', lw=2)
        plt.plot(coordinates[fitted_hull.vertices[0], 0], coordinates[fitted_hull.vertices[0], 1], 'ro')

        # plot contour
        plt.plot(boundary[:, 0], boundary[:, 1], 'o')
        plt.plot(boundary[:, 0], boundary[:, 1], 'b--', lw=2)
        plt.title("convex hull + boundary.png")
        plt.savefig(fig_export_path + fig_name + "_" + "hull_n_boundary.png")
        plt.close()

    return convex_hull_perimeter / contour_perimeter
