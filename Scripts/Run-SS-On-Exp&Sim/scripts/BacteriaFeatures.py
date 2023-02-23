import pandas as pd
import numpy as np
from lownerJohnEllipseMaster.src.lownerJohnEllipse import welzl


def min_max_normalization(feature_val, magnitude=400):

    if type(feature_val).__module__ != 'numpy':
        feature_val = np.array(feature_val)

    normalized_val = magnitude * ((feature_val - min(feature_val)) / (max(feature_val) - min(feature_val)))

    return normalized_val


def convert_um_pixel(feature_val, um_pixel_ratio=0.144):
    """
    goal: Convert um to pixels
    @param feature_val list bacteria feature values
    """

    return np.array(feature_val) / um_pixel_ratio


def convert_radian_degree(radian_val):

    return np.array(radian_val) * 180 / np.pi

# find co-vertices of an ellipse (bacteria):
# references:
# https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
# https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio


# find co vertexes
def find_co_vertex(center_x_list, center_y_list, minor_list, angle_rotation_list):

    """
    goal: find co-vertex of bacteria in a specif micro colony(local mode) or in specific time step (global mode)
    @param center_x_list list x coordinate of center of bacteria
    @param center_y_list list y coordinate of center of bacteria
    @param minor_list list  length of minor axis of bacteria
    @param angle_rotation_list list orientation of bacteria

    Returns:
        x_co_vertex list x coordinate of co-vertex of bacteria
        y_co_vertex list y coordinate of co-vertex of bacteria

    """
    # list of co-vertices points
    x_co_vertex = []
    y_co_vertex = []

    for i in range(len(center_x_list)):
        center_x = center_x_list[i]
        center_y = center_y_list[i]
        minor = minor_list[i] / 2
        angle_rotation = angle_rotation_list[i]

        # equations
        # np.power((x- center_x) * np.sin(angle_rotation) - (y-center_y) * np.cos(angle_rotation), 2) =
        # np.power(minor, 2)
        # np.power((x - center_x) * np.cos(angle_rotation) + (y - center_y) * np.sin(angle_rotation), 2) = 0
        # so: (x - center_x) * np.cos(angle_rotation) + (y - center_y) * np.sin(angle_rotation) = 0

        co_vertex_1_y = minor / (np.tan(angle_rotation) * np.sin(angle_rotation) + np.cos(angle_rotation)) + center_y
        co_vertex_1_x = -(co_vertex_1_y - center_y) * np.tan(angle_rotation) + center_x
        co_vertex_2_y = -minor / (np.tan(angle_rotation) * np.sin(angle_rotation) + np.cos(angle_rotation)) + center_y
        co_vertex_2_x = -(co_vertex_2_y - center_y) * np.tan(angle_rotation) + center_x

        x_co_vertex.extend([float(co_vertex_1_x), float(co_vertex_2_x)])
        y_co_vertex.extend([float(co_vertex_1_y), float(co_vertex_2_y)])

    return x_co_vertex, y_co_vertex


def find_bacteria_endpoints(cs):
    # end points of bacteria
    bacteria_end_point1_x = [cs[it].ends[0][0] for it in cs]
    bacteria_end_point1_y = [cs[it].ends[0][1] for it in cs]
    bacteria_end_point2_x = [cs[it].ends[1][0] for it in cs]
    bacteria_end_point2_y = [cs[it].ends[1][1] for it in cs]

    return bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y


def find_bacteria_major(cs):
    bacteria_major = [cs[it].length for it in cs]
    return bacteria_major


def find_bacteria_minor(cs):
    bacteria_minor = [cs[it].radius for it in cs]
    return bacteria_minor


def find_bacteria_id(cs):
    bacteria_id = [cs[it].id for it in cs]
    return bacteria_id


def find_bacteria_label(cs):
    bacteria_label = [cs[it].label for it in cs]
    return bacteria_label


def find_bacteria_orientation(cs):
    # orientation
    # direction vector: [x, y] --> orientation: arctan (y / x)
    bacteria_orientation = [np.arctan2(cs[it].dir[1], cs[it].dir[0]) for it in cs if cs]
    return bacteria_orientation


def find_bacteria_center_position(cs):
    # center position
    bacteria_x_center = [cs[it].pos[0] for it in cs]
    bacteria_y_center = [cs[it].pos[1] for it in cs]
    return bacteria_x_center, bacteria_y_center


def find_bacteria_cell_age(cs):
    # cell age
    bacteria_cell_age = [cs[it].cellAge for it in cs]
    return bacteria_cell_age


def find_bacteria_growth_rate(cs, dt):
    # growth rate
    strainRate_rolling = [cs[it].strainRate_rolling for it in cs]
    strainRate_rolling = [np.nan if str(x) == 'nan' else x for x in strainRate_rolling]
    bacteria_growth_rate = [element / dt for element in strainRate_rolling]
    return bacteria_growth_rate


def func_current_bacteria_features(cs, dt):
    """
    Goal: This function retrieves and/or calculates important features of bacteria ( it can be from sim. or exp.)
          in a specific time step and stores them in a dataframe.
          It can be used to check whether bacteria micro colonies have merged, fit an ellipse around micro colonies,
          and calculate summary statistics.
          These features includes:
          'id', 'label', 'minor', 'major', 'x_end_point1', 'y_end_point1', 'x_end_point2', 'y_end_point2',
                    'orientation', 'x_center', 'y_center', 'cell_age', 'growth_rate'
    @param cs: dictionary  bacteria features value
    @param  dt: float  interval time
    Return: df_current_time_step  dataframe   Values of important features of bacteria in a dataframe
    """
    # get important features of bacteria to draw them

    bacteria_id = find_bacteria_id(cs)
    bacteria_label = find_bacteria_label(cs)
    bacteria_minor = find_bacteria_minor(cs)
    bacteria_major = find_bacteria_major(cs)

    bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y = \
        find_bacteria_endpoints(cs)

    bacteria_orientation = find_bacteria_orientation(cs)
    bacteria_x_center, bacteria_y_center = find_bacteria_center_position(cs)
    bacteria_cell_age = find_bacteria_cell_age(cs)
    bacteria_growth_rate = find_bacteria_growth_rate(cs, dt)

    # convert to dataframe
    bacteria_features = list(zip(bacteria_id, bacteria_label, bacteria_minor, bacteria_major,
                             bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y,
                             bacteria_orientation, bacteria_x_center, bacteria_y_center, bacteria_cell_age,
                             bacteria_growth_rate))

    columns_name = ['id', 'label', 'minor', 'major', 'x_end_point1', 'y_end_point1', 'x_end_point2', 'y_end_point2',
                    'orientation', 'x_center', 'y_center', 'cell_age', 'growth_rate']
    df_current_time_step = pd.DataFrame(bacteria_features, columns=columns_name)

    return df_current_time_step


def fit_enclosing_ellipse(points):
    """
    goal: An ellipse is fitted around bacteria in a specific micro colony (local mode) or
    in specific time step (global mode).

    @param points dataframe endpoints (vertex & co-vertex) of bacteria in specific micro colony(local mode) or
    in specific time step (global mode)

    Returns: enclosing_ellipse tuple a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    """

    # convert dataframe to numpy array
    points = points.to_numpy()
    # print(points)
    # finds the smallest ellipse covering a finite set of points
    # https://github.com/dorshaviv/lowner-john-ellipse
    enclosing_ellipse = welzl(interior=points)

    return enclosing_ellipse


def fit_ellipse(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or in specific time step (global mode)

    Returns: ellipse_params  tuple a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    """
    # fit ellipse to micro colony
    # endpoints
    vertex1_x, vertex1_y, vertex2_x, vertex2_y = find_bacteria_endpoints(cs)

    bacteria_center_x, bacteria_center_y = find_bacteria_center_position(cs)
    bacteria_minor = find_bacteria_minor(cs)
    bacteria_orientation = find_bacteria_orientation(cs)

    # now I calculate co-vertexes of bacteria ellipse
    co_vertex_x, co_vertex_y = find_co_vertex(bacteria_center_x, bacteria_center_y, bacteria_minor,
                                              bacteria_orientation)

    endpoints_x = vertex1_x + vertex2_x + co_vertex_x
    endpoints_y = vertex1_y + vertex2_y + co_vertex_y
    endpoints = pd.DataFrame(zip(endpoints_x, endpoints_y))

    # fit ellipse
    ellipse_params = fit_enclosing_ellipse(endpoints)

    return ellipse_params

