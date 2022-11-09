import sys
import glob
import networkx as nx
import pandas as pd
import numpy as np
import pickle
from scipy.spatial import distance_matrix
import CellModeller
import ProcessCellProfilerData
from summary_statistics import fit_ellipse, aspect_ratio_calc, anisotropy_calc
from fractal_dimension import calc_fractal_dimension
import density_calculation
import growth_rate_stats

# set number of recursion
sys.setrecursionlimit(10000)


def finding_sub_graphs(graph):
    """
    Goal: Finding the largest disconnected subparagraphs
    @param graph graph neighboring bacteria graph in current time-step
    return sub_graphs list  List elements including bacteria ids for each subgraph
    """
    sub_graphs = sorted(nx.connected_components(graph), key=len, reverse=True)

    return sub_graphs


def make_adjacency_matrix(cs, max_distance_between_cells):
    """
    goal: make adjacency matrix
    @param cs dictionary  bacteria information
    @return: adjacency_matrix_df   dataframe  adjacency matrix dataframe
    """

    # bacteria id
    cell_id = [cs[it].id for it in cs]

    # end points of bacteria
    bacteria_first_end_points_x_in_current_timestep = [cs[it].ends[0][0] for it in cs]
    bacteria_first_end_points_y_in_current_timestep = [cs[it].ends[0][1] for it in cs]
    bacteria_second_end_points_x_in_current_timestep = [cs[it].ends[1][0] for it in cs]
    bacteria_second_end_points_y_in_current_timestep = [cs[it].ends[1][1] for it in cs]

    bacteria_end_point_1 = list(zip(bacteria_first_end_points_x_in_current_timestep,
                                    bacteria_first_end_points_y_in_current_timestep))

    bacteria_end_point_2 = list(zip(bacteria_second_end_points_x_in_current_timestep,
                                    bacteria_second_end_points_y_in_current_timestep))

    # each row shows distance of a bacterium  to the other bacteria
    dist_endpoints_1 = pd.DataFrame(distance_matrix(bacteria_end_point_1, bacteria_end_point_1))
    dist_endpoints_2 = pd.DataFrame(distance_matrix(bacteria_end_point_2, bacteria_end_point_2))
    dist_endpoints_1_2 = pd.DataFrame(distance_matrix(bacteria_end_point_1, bacteria_end_point_2))
    dist_endpoints_2_1 = pd.DataFrame(distance_matrix(bacteria_end_point_2, bacteria_end_point_1))

    # find the min distance of each two bacterium
    df_min_distance = pd.concat([dist_endpoints_1, dist_endpoints_2, dist_endpoints_1_2,
                                 dist_endpoints_2_1]).groupby(level=0).min()

    # ignore diagonal values
    # because it's bacteria distance of itself (zero value)
    # I change them to `nans`
    np.fill_diagonal(df_min_distance.values, np.nan)

    # extract i,j pairs where distance < threshold
    paires = pd.DataFrame(np.argwhere(df_min_distance.to_numpy() <= max_distance_between_cells),
                          columns=['cell', 'neighbour cell'])

    # create dataframe
    adjacency_matrix = pd.crosstab(paires['cell'], paires['neighbour cell'])
    # change row index and column name to bacteria id
    dictionary_index_to_id = dict(zip(list(range(len(cs))), cell_id))
    # 1. rename
    adjacency_matrix = adjacency_matrix.rename(columns=dictionary_index_to_id, index=dictionary_index_to_id)
    # 2. reindex
    adjacency_matrix = adjacency_matrix.reindex(index=cell_id, columns=cell_id, fill_value=0)

    return adjacency_matrix


def micro_colony_analysis(pickle_files_directory, summary_statistic_method, dt, min_size_of_micro_colony,
                          max_distance_between_cells, um_pixel_ratio, fig_export_path):
    """
    Goal: this is the main function; the function calls other function that are described above and
    calculates summary statistics for each micro colony in each time-step and store the outputs in a dictionary.

    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply on the micro-colonies
    @param dt float interval time
    @param min_size_of_micro_colony int minimum size of micro colony
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    @param fig_export_path directory to export images
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic
    
    """

    # store summary statistics of micro colonies
    aspect_ratio_list = []
    anisotropy_list = []
    density_list = []
    dist_vs_growth_rate_list = []
    fractal_dimension_list = []

    # read pickle files
    path = pickle_files_directory + "/*.pickle"
    filename_list = [filename for filename in sorted(glob.glob(path))]

    # In order to identify merged micro colonies, micro colony labels are stored.
    micro_colonies = []

    for cnt, filename in enumerate(filename_list):

        # store ellipses
        ellipses = []
        # store anisotropies of micro colonies of this time step
        local_anisotropy_list = []
        # store anisotropies of micro colonies of this time step
        local_aspect_ratio_list = []
        # store density of micro colonies of this time step
        local_density_list = []
        # store Correlate growth penalty based on location in micro colony of this time step
        local_dist_vs_growth_rate_list = []
        # store fractal dimension of this time step
        local_fractal_dimension_list = []

        micro_colonies_in_current_time_step = []

        timestep = cnt + 1
        print('time step:' + str(timestep))

        # read current pickle file
        current_bacteria_info = pickle.load(open(filename_list[cnt], 'rb'))
        cs = current_bacteria_info['cellStates']

        # find neighbours
        adjacency_matrix = make_adjacency_matrix(cs, max_distance_between_cells)

        # create graph
        graph = nx.from_pandas_adjacency(adjacency_matrix)

        # finding sub-graphs
        sub_graphs = finding_sub_graphs(graph)

        for sub_graph_index, subgraph in enumerate(sub_graphs):
            if len(subgraph) >= min_size_of_micro_colony:
                subgraph_nodes = list(subgraph)
                # find unique bacteria
                subgraph_bacteria_label = list(set([cs[it].label for it in subgraph_nodes]))

                intersection_list = [(mico_colony_index, np.intersect1d(subgraph_bacteria_label, micro_colony)) for
                                     mico_colony_index, micro_colony in enumerate(micro_colonies)]
                # remove empty array (without intersection)
                intersection_list = [list(x) for x in intersection_list if np.size(x[1]) > 0]

                if len(intersection_list) <= 1:  # As a result, the micro_colonies were not merged
                    # Add new labels to its micro_colony list
                    if len(intersection_list) == 1:
                        micro_colonies[intersection_list[0][0]].extend(subgraph_bacteria_label)
                        micro_colonies[intersection_list[0][0]] = \
                            list(set(micro_colonies[intersection_list[0][0]]))
                    else:
                        micro_colonies.append(subgraph_bacteria_label)

                    # store micro colony bacteria id
                    micro_colonies_in_current_time_step.append(subgraph_nodes)

                    # features value of bacteria in this micro colony
                    bacteria_in_this_micro_colony = {k: v for k, v in cs.items() if k in subgraph_nodes}

                    # fit ellipse
                    ellipse_params = fit_ellipse(bacteria_in_this_micro_colony)
                    # append ellipse
                    ellipses.append(ellipse_params)

                    """
                                calculation of aspect ratio
                    """
                    if "Aspect Ratio" in summary_statistic_method:
                        aspect_ratio = aspect_ratio_calc(ellipse_params)
                        # store aspect ratio
                        aspect_ratio_list.append(aspect_ratio)
                        local_aspect_ratio_list.append(aspect_ratio)
                    """
                                Finish
                    """

                    """
                        calculation of Anisotropy
                    """
                    if "Anisotropy" in summary_statistic_method:
                        mean_anisotropy = anisotropy_calc(bacteria_in_this_micro_colony, max_distance_between_cells)
                        # store anisotropy
                        anisotropy_list.append(mean_anisotropy)
                        local_anisotropy_list.append(mean_anisotropy)
                    """
                        Finish
                    """

                    """
                        calculation of Density
                    """
                    if "Density" in summary_statistic_method:
                        density = density_calculation.main(bacteria_in_this_micro_colony)
                        # store density
                        density_list.append(density)
                        local_density_list.append(density)
                    """
                        Finish
                    """

                    """
                        calculation of Correlate growth penalty based on location in microcolony
                    """
                    # Note that this summary stat currently is calc. only for micro colonies containing more than 2
                    # bacteria because the center of the bacteria is used to fit ellipse but the ellipse can be
                    # when at least three points are available
                    if "dist_vs_growth_rate" in summary_statistic_method and len(bacteria_in_this_micro_colony) > 2:
                        dist_vs_growth_rate = growth_rate_stats.main(bacteria_in_this_micro_colony, dt)
                        # store dist_vs_growth_rate
                        dist_vs_growth_rate_list.append(dist_vs_growth_rate)
                        local_dist_vs_growth_rate_list.append(dist_vs_growth_rate)
                    """
                        Finish
                    """

                    """
                        calculation of fractal dimension in microcolony
                    """
                    if "fractal_dimension" in summary_statistic_method:
                        fig_name = str(timestep) + '_' + str(sub_graph_index)
                        fractal_dimension = calc_fractal_dimension(bacteria_in_this_micro_colony, fig_export_path,
                                                                   fig_name)
                        # store fractal dimension
                        fractal_dimension_list.append(fractal_dimension)
                        local_fractal_dimension_list.append(fractal_dimension)
                    """
                        Finish
                    """

    report_mean_summary_statistics = {}
    if "Aspect Ratio" in summary_statistic_method:
        report_mean_summary_statistics["Aspect Ratio"] = np.mean(aspect_ratio_list)
    if "Anisotropy" in summary_statistic_method:
        # remove nan values
        anisotropy_list = [x for x in anisotropy_list if str(x) != 'nan']
        report_mean_summary_statistics["Anisotropy"] = np.mean(anisotropy_list)
    if "Density" in summary_statistic_method:
        report_mean_summary_statistics["Density"] = np.mean(density_list)
    if "dist_vs_growth_rate" in summary_statistic_method:
        # remove nan values
        dist_vs_growth_rate_list = [x for x in dist_vs_growth_rate_list if str(x) != 'nan']
        report_mean_summary_statistics["dist_vs_growth_rate"] = np.mean(dist_vs_growth_rate_list)
    if "fractal dimension" in summary_statistic_method:
        report_mean_summary_statistics["fractal dimension"] = np.mean(fractal_dimension_list)

    return report_mean_summary_statistics


def global_analysis(pickle_files_directory, summary_statistic_method, dt, max_distance_between_cells, um_pixel_ratio,
                    fig_export_path):
    """
    Goal: this is the main function; the function calls other functions described above and
    calculates summary statistics for each time-step and store the outputs in a dictionary.

    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply over time-lapse
    @param dt float interval time
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    @param  fig_export_path directory to export images
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic

    """

    # store summary statistics over time-lapse (for each individual time-step)
    aspect_ratio_list = []
    anisotropy_list = []
    density_list = []
    dist_vs_growth_rate_list = []
    fractal_dimension_list = []

    # read pickle files
    path = pickle_files_directory + "/*.pickle"
    filename_list = [filename for filename in sorted(glob.glob(path))]

    for cnt, filename in enumerate(filename_list):

        timestep = cnt + 1
        print('time step:' + str(timestep))

        # read current pickle file
        current_bacteria_info = pickle.load(open(filename_list[cnt], 'rb'))
        cs = current_bacteria_info['cellStates']

        # fit ellipse
        # if  number of bacteria < 2 : ellipse_params = `nan`
        ellipse_params = fit_ellipse(cs)

        """
           calculation of aspect ratio
        """
        # if  number of bacteria < 2 : ellipse_params = `nan`
        # so: the aspect ratio can not be calculated
        if "Aspect Ratio" in summary_statistic_method and len(cs) > 1:
            aspect_ratio = aspect_ratio_calc(ellipse_params)
            # store aspect ratio
            aspect_ratio_list.append(aspect_ratio)
        """
           Finish
        """

        """
           calculation of Anisotropy
        """
        if "Anisotropy" in summary_statistic_method and len(cs) > 1:
            mean_anisotropy = anisotropy_calc(cs, max_distance_between_cells)
            # store anisotropy
            anisotropy_list.append(mean_anisotropy)
        """
           Finish
        """

        """
           calculation of Density
        """
        if "Density" in summary_statistic_method:
            density = density_calculation.main(cs)
            # store anisotropy
            density_list.append(density)
        """
           Finish
        """

        """
           calculation of Correlate growth penalty based on location 
        """
        # Note that this summary stat currently is calc. only for time steps containing more than 2 bacteria
        # because the center of the bacteria is used to fit ellipse but the ellipse can be fitted when at
        # least three points are available
        if "dist_vs_growth_rate" in summary_statistic_method and len(cs) > 2:
            dist_vs_growth_rate = growth_rate_stats.main(cs, dt)
            # store dist_vs_growth_rate
            dist_vs_growth_rate_list.append(dist_vs_growth_rate)
        """
           Finish
        """

        """
           calculation of fractal dimension
        """
        if "fractal dimension" in summary_statistic_method:
            # image name
            fig_name = str(timestep)
            fractal_dimension = calc_fractal_dimension(cs, fig_export_path, fig_name)
            # store fractal dimension
            fractal_dimension_list.append(fractal_dimension)
        """
           Finish
        """

    report_mean_summary_statistics = {}
    if "Aspect Ratio" in summary_statistic_method:
        report_mean_summary_statistics["Aspect Ratio"] = np.mean(aspect_ratio_list)
    if "Anisotropy" in summary_statistic_method:
        # remove nan values
        anisotropy_list = [x for x in anisotropy_list if str(x) != 'nan']
        report_mean_summary_statistics["Anisotropy"] = np.mean(anisotropy_list)
    if "Density" in summary_statistic_method:
        report_mean_summary_statistics["Density"] = np.mean(density_list)
    if "dist_vs_growth_rate" in summary_statistic_method:
        # remove nan values
        dist_vs_growth_rate_list = [x for x in dist_vs_growth_rate_list if str(x) != 'nan']
        report_mean_summary_statistics["dist_vs_growth_rate"] = np.mean(dist_vs_growth_rate_list)
    if "fractal dimension" in summary_statistic_method:
        report_mean_summary_statistics["fractal dimension"] = np.mean(fractal_dimension_list)

    return report_mean_summary_statistics


def data_analysis(pickle_files_directory, summary_statistic_method, dt, mode='global', max_distance_between_cells=3.4,
                  um_pixel_ratio=0.144, min_size_of_micro_colony=2, fig_export_path=''):
    """
    goal: calculation of summary statistics in global or local mode
    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply on the micro-colonies
    @param dt float interval time
    @param mode str it shows we want to calculate global summary statistics or local summary statistics
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    @param min_size_of_micro_colony int minimum size of micro colony (only for local mode)
    @param  fig_export_path directory to export images
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic

    Returns:

    """

    if mode == 'local':
        return micro_colony_analysis(pickle_files_directory, summary_statistic_method, dt, min_size_of_micro_colony,
                                     max_distance_between_cells, um_pixel_ratio, fig_export_path)
    elif mode == 'global':
        return global_analysis(pickle_files_directory, summary_statistic_method, dt, max_distance_between_cells,
                               um_pixel_ratio, fig_export_path)
