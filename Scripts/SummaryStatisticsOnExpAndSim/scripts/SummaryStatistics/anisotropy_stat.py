import numpy as np
import scipy.linalg as la
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.process import cell_data


def anisotropy_calc(cs, max_neighbour_distance):
    """
    goal: calculation of Anisotropy
    @param cs dict bacteria features value in specific micro colony(local mode) or in specific time step (global mode)
    @param max_neighbour_distance float There is a maximum distance between bacteria that can be neighbours.

    Returns: mean_anisotropy float average value of anisotropy of bacteria in specific micro colony(local mode) or
    in specific time step (global mode)

    """
    # main idea: https://github.com/ingallslab/bsim-related/blob/main/bsim_related/data_processing/cell_data_processing.py#L184

    local_anisotropies = []

    if len(cs) >= 2:
        # orientation of bacteria
        bacteria_orientation = cell_data.find_bacteria_orientation(cs)
        bacteria_x_end_point1, bacteria_y_end_point1, bacteria_x_end_point2, bacteria_y_end_point2 = \
            cell_data.find_vertex(cs)

        for bacterium_index in range(len(cs)):
            num_neighbours = 0
            # Projection matrix
            projection_matrix = np.zeros(shape=(2, 2))
            for other_bacterium_index in range(len(cs)):
                if other_bacterium_index != bacterium_index:
                    distance_between_bacteria_end_point1 = np.hypot(bacteria_x_end_point1[other_bacterium_index] -
                                                                    bacteria_x_end_point1[bacterium_index],
                                                                    bacteria_y_end_point1[other_bacterium_index]
                                                                    - bacteria_y_end_point1[bacterium_index])

                    distance_between_bacteria_end_point2 = np.hypot(bacteria_x_end_point2[other_bacterium_index] -
                                                                    bacteria_x_end_point2[bacterium_index],
                                                                    bacteria_y_end_point2[other_bacterium_index]
                                                                    - bacteria_y_end_point2[bacterium_index])

                    distance_between_bacteria_end_point1_2 = np.hypot(bacteria_x_end_point1[other_bacterium_index] -
                                                                      bacteria_x_end_point2[bacterium_index],
                                                                      bacteria_y_end_point1[other_bacterium_index]
                                                                      - bacteria_y_end_point2[bacterium_index])

                    distance_between_bacteria_end_point2_1 = np.hypot(bacteria_x_end_point2[other_bacterium_index] -
                                                                      bacteria_x_end_point1[bacterium_index],
                                                                      bacteria_y_end_point2[other_bacterium_index]
                                                                      - bacteria_y_end_point1[bacterium_index])

                    distance_list = [distance_between_bacteria_end_point1, distance_between_bacteria_end_point2,
                                     distance_between_bacteria_end_point1_2, distance_between_bacteria_end_point2_1]

                    if min(distance_list) <= max_neighbour_distance:
                        # Compute the sum of the projection matrices on the orientation vectors of
                        # the neighbouring bacteria
                        # projection matrix
                        """
                        cos(angle)                  cos(angle)*sin(angle)
                        cos(angle)*sin(angle)       sin(angle)
                        """
                        projection_matrix += np.matrix([[np.cos(bacteria_orientation[other_bacterium_index]) ** 2,
                                                         np.cos(bacteria_orientation[other_bacterium_index]) * np.sin(
                                                             bacteria_orientation[other_bacterium_index])],
                                                        [np.cos(bacteria_orientation[other_bacterium_index]) * np.sin(
                                                            bacteria_orientation[other_bacterium_index]),
                                                         np.sin(bacteria_orientation[other_bacterium_index]) ** 2]])

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
            mean_anisotropy = np.nan
    else:
        mean_anisotropy = np.nan
    return mean_anisotropy
