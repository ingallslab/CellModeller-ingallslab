# Extracting bacteria features
import numpy as np


# find co-vertices of an ellipse (bacteria):
# references:
# https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
# https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio


# find co vertexes
def find_co_vertex(cs):
    """
    goal: find co-vertex of bacteria in a specif micro colony(local mode) or in specific time step (global mode)
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)

    Returns:
        x_co_vertex list x coordinate of co-vertex of bacteria
        y_co_vertex list y coordinate of co-vertex of bacteria

    """

    # get cells features
    minor_list = find_bacteria_minor(cs)
    angle_rotation_list = find_bacteria_orientation(cs)
    center_x_list, center_y_list = find_bacteria_center(cs)

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


def find_vertex(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)
    @return vertex1_x, vertex1_y, vertex2_x, vertex2_y  x,y vertices
    """
    vertex1_x = [cs[it].ends[0][0] for it in cs]
    vertex1_y = [cs[it].ends[0][1] for it in cs]
    vertex2_x = [cs[it].ends[1][0] for it in cs]
    vertex2_y = [cs[it].ends[1][1] for it in cs]

    return vertex1_x, vertex1_y, vertex2_x, vertex2_y


def find_bacteria_center(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)
    @return bacteria_center_x, bacteria_center_y  bacteria center coordinate
    """
    bacteria_center_x = [cs[it].pos[0] for it in cs]
    bacteria_center_y = [cs[it].pos[1] for it in cs]

    return bacteria_center_x, bacteria_center_y


def find_bacteria_orientation(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)
    """
    # direction vector: [x, y] --> orientation: arctan (y / x)
    bacteria_orientation = [np.arctan2(cs[it].dir[1], cs[it].dir[0]) for it in cs if cs]

    return bacteria_orientation


def find_bacteria_minor(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)
    """

    bacteria_minor = [cs[it].radius for it in cs]

    return bacteria_minor


def bacteria_end_points(cs):
    """
    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)
    @return endpoints_x, endpoints_y
    """

    # vertexes
    vertex1_x, vertex1_y, vertex2_x, vertex2_y = find_vertex(cs)

    # now I calculate co-vertexes of bacteria ellipse
    co_vertex_x, co_vertex_y = find_co_vertex(cs)

    endpoints_x = vertex1_x + vertex2_x + co_vertex_x
    endpoints_y = vertex1_y + vertex2_y + co_vertex_y

    return endpoints_x, endpoints_y
