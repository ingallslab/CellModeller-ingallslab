from lownerJohnEllipse import welzl
import numpy as np
import pandas as pd
import scipy.linalg as la


def fit_enclosing_ellipse(points):
    # convert dataframe to numpy array
    points = points.to_numpy()
    # print(points)
    # finds the smallest ellipse covering a finite set of points
    # https://github.com/dorshaviv/lowner-john-ellipse
    enclosing_ellipse = welzl(points)
    return enclosing_ellipse


# find co-vertices of an ellipse (bacteria):
# references:
# https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
# https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio


# find co vertexes
def find_co_vertex(center_x_list, center_y_list, minor_list, angle_rotation_list):
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


def fit_ellipse(bac_in_micro_colony):
    # fit ellipse to micro colony
    # endpoints
    vertex1_x = bac_in_micro_colony['x_end_point1'].values.tolist()
    vertex1_y = bac_in_micro_colony['y_end_point1'].values.tolist()
    vertex2_x = bac_in_micro_colony['x_end_point2'].values.tolist()
    vertex2_y = bac_in_micro_colony['y_end_point2'].values.tolist()

    bacteria_center_x = bac_in_micro_colony['x_center'].values.tolist()
    bacteria_center_y = bac_in_micro_colony['y_center'].values.tolist()
    bacteria_minor = bac_in_micro_colony['minor'].values.tolist()
    bacteria_orientation = bac_in_micro_colony['orientation'].values.tolist()

    # now I calculate co-vertexes of bacteria ellipse
    co_vertex_x, co_vertex_y = find_co_vertex(bacteria_center_x, bacteria_center_y, bacteria_minor,
                                              bacteria_orientation)

    endpoints_x = vertex1_x + vertex2_x + co_vertex_x
    endpoints_y = vertex1_y + vertex2_y + co_vertex_y
    endpoints = pd.DataFrame(zip(endpoints_x, endpoints_y))

    # fit ellipse
    ellipse_params = fit_enclosing_ellipse(endpoints)

    return ellipse_params


def aspect_ratio_calc(ellipse_params):
    # calculate aspect ratio
    center_pos, major, minor, theta = ellipse_params
    aspect_ratio = minor / major

    return aspect_ratio


def anisotropy_calc(bac_in_micro_colony, neighbour_distance):
    # main idea: https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/cell_data_processing.py#L184

    local_anisotropies = []

    # orientation of bacteria
    bacteria_orientation = bac_in_micro_colony["orientation"]
    bacteria_x_end_point1 = bac_in_micro_colony['x_end_point1']
    bacteria_y_end_point1 = bac_in_micro_colony['y_end_point1']
    bacteria_x_end_point2 = bac_in_micro_colony['x_end_point2']
    bacteria_y_end_point2 = bac_in_micro_colony['y_end_point2']

    for bacterium_index in range(bac_in_micro_colony.shape[0]):
        num_neighbours = 0
        # Projection matrix
        projection_matrix = np.zeros(shape=(2, 2))
        for other_bacterium_index in range(bac_in_micro_colony.shape[0]):
            if other_bacterium_index != bacterium_index:
                distance_between_bacteria_end_point1 = np.hypot(bacteria_x_end_point1.iloc[other_bacterium_index] -
                                                                bacteria_x_end_point1.iloc[bacterium_index],
                                                                bacteria_y_end_point1.iloc[other_bacterium_index]
                                                                - bacteria_y_end_point1.iloc[bacterium_index])

                distance_between_bacteria_end_point2 = np.hypot(bacteria_x_end_point2.iloc[other_bacterium_index] -
                                                                bacteria_x_end_point2.iloc[bacterium_index],
                                                                bacteria_y_end_point2.iloc[other_bacterium_index]
                                                                - bacteria_y_end_point2.iloc[bacterium_index])

                distance_between_bacteria_end_point1_2 = np.hypot(bacteria_x_end_point1.iloc[other_bacterium_index] -
                                                                  bacteria_x_end_point2.iloc[bacterium_index],
                                                                  bacteria_y_end_point1.iloc[other_bacterium_index]
                                                                  - bacteria_y_end_point2.iloc[bacterium_index])

                distance_between_bacteria_end_point2_1 = np.hypot(bacteria_x_end_point2.iloc[other_bacterium_index] -
                                                                  bacteria_x_end_point1.iloc[bacterium_index],
                                                                  bacteria_y_end_point2.iloc[other_bacterium_index]
                                                                  - bacteria_y_end_point1.iloc[bacterium_index])

                distance_list = [distance_between_bacteria_end_point1, distance_between_bacteria_end_point2,
                                 distance_between_bacteria_end_point1_2, distance_between_bacteria_end_point2_1]

                if min(distance_list) <= neighbour_distance:
                    # Compute the sum of the projection matrices on the orientation vectors of the neighbouring bacteria
                    # projection matrix
                    """
                    cos(angle)                  cos(angle)*sin(angle)
                    cos(angle)*sin(angle)       sin(angle)
                    """
                    projection_matrix += np.matrix([[np.cos(bacteria_orientation.iloc[other_bacterium_index]) ** 2,
                                                     np.cos(bacteria_orientation.iloc[other_bacterium_index]) * np.sin(
                                                         bacteria_orientation.iloc[other_bacterium_index])],
                                                    [np.cos(bacteria_orientation.iloc[other_bacterium_index]) * np.sin(
                                                        bacteria_orientation.iloc[other_bacterium_index]),
                                                     np.sin(bacteria_orientation.iloc[other_bacterium_index]) ** 2]])

                    num_neighbours += 1

        if num_neighbours > 0:
            # Compute the mean of the projection matrices on the orientation vectors of the neighbouring bacteria
            projection_matrix = projection_matrix / num_neighbours
            # Get the max real eigenvalues of the mean projection matrix; this is the local anisotropy
            local_anisotropies.append(max(la.eigvals(projection_matrix).real))

    if local_anisotropies:
        # calculate mean anisotropy
        mean_anisotropy = np.mean(local_anisotropies)
    else:
        mean_anisotropy = 'Nan'

    return mean_anisotropy
