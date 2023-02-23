import sys
import glob
import networkx as nx
import pandas as pd
import numpy as np
import pickle
from scipy.spatial import distance_matrix
import CellModeller
import CellProfilerAnalysis
from summary_statistics import fit_ellipse, aspect_ratio_calc, anisotropy_calc
from density_calculation import calc_density
from growth_rate_stats import calc_dist_vs_growth_rate
from fourier_descriptor import calc_fourier_descriptor

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


def func_current_bacteria_info(cs, dt):
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
    bacteria_id_in_current_timestep = [cs[it].id for it in cs]
    bacteria_label_in_current_timestep = [cs[it].label for it in cs]
    bacteria_minor_in_current_timestep = [cs[it].radius for it in cs]
    bacteria_major_in_current_timestep = [cs[it].length for it in cs]
    # end points of bacteria
    bacteria_first_end_points_x_in_current_timestep = [cs[it].ends[0][0] for it in cs]
    bacteria_first_end_points_y_in_current_timestep = [cs[it].ends[0][1] for it in cs]
    bacteria_second_end_points_x_in_current_timestep = [cs[it].ends[1][0] for it in cs]
    bacteria_second_end_points_y_in_current_timestep = [cs[it].ends[1][1] for it in cs]
    # orientation
    # direction vector: [x, y] --> orientation: arctan (y / x)
    bacteria_orientation_in_current_time_step = [np.arctan2(cs[it].dir[1], cs[it].dir[0]) for it in cs if cs]

    # center position
    bacteria_x_center_in_current_timestep = [cs[it].pos[0] for it in cs]
    bacteria_y_center_in_current_timestep = [cs[it].pos[1] for it in cs]

    # cell age
    bacteria_cell_age_in_current_timestep = [cs[it].cellAge for it in cs]

    # growth rate
    strainRate_rolling = [cs[it].strainRate_rolling for it in cs]
    strainRate_rolling = [np.nan if x == 'nan' else x for x in strainRate_rolling]
    bacteria_growth_rate_in_current_timestep = [element / dt for element in strainRate_rolling]

    # convert to dataframe
    bacteria_info = list(zip(bacteria_id_in_current_timestep,
                             bacteria_label_in_current_timestep,
                             bacteria_minor_in_current_timestep,
                             bacteria_major_in_current_timestep,
                             bacteria_first_end_points_x_in_current_timestep,
                             bacteria_first_end_points_y_in_current_timestep,
                             bacteria_second_end_points_x_in_current_timestep,
                             bacteria_second_end_points_y_in_current_timestep,
                             bacteria_orientation_in_current_time_step,
                             bacteria_x_center_in_current_timestep,
                             bacteria_y_center_in_current_timestep,
                             bacteria_cell_age_in_current_timestep,
                             bacteria_growth_rate_in_current_timestep))

    columns_name = ['id', 'label', 'minor', 'major', 'x_end_point1', 'y_end_point1', 'x_end_point2', 'y_end_point2',
                    'orientation', 'x_center', 'y_center', 'cell_age', 'growth_rate']
    df_current_time_step = pd.DataFrame(bacteria_info, columns=columns_name)

    return df_current_time_step


def micro_colony_analysis(pickle_files_directory, summary_statistic_method, dt, min_size_of_micro_colony,
                          max_distance_between_cells, um_pixel_ratio, fig_path):
    """
    Goal: this is the main function; the function calls other function that are described above and
    calculates summary statistics for each micro colony in each time-step and store the outputs in a dictionary.
    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply on the micro-colonies
    @param dt float interval time
    @param min_size_of_micro_colony int minimum size of micro colony
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic
    
    """

    # store summary statistics of micro colonies
    aspect_ratio_list = []
    anisotropy_list = []
    density_list = []
    dist_vs_growth_rate_list = []
    fourier_descriptor_list = []

    # read pickle files
    path = pickle_files_directory + "/*.pickle"
    filename_list = [filename for filename in sorted(glob.glob(path))]

    # In order to identify merged micro colonies, micro colony labels are stored.
    micro_colonies = []
    micro_colony_num = 0

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
        # store fourier descriptor of micro colonies of this time step
        local_fourier_descriptor_list = []
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

        # get important features of bacteria
        df_current_time_step = func_current_bacteria_info(cs, dt)

        for sub_graph_index, subgraph in enumerate(sub_graphs):
            if len(subgraph) >= min_size_of_micro_colony:
                subgraph_nodes = list(subgraph)
                # find unique bacteria labels
                subgraph_bacteria_label = list(set(df_current_time_step[df_current_time_step["id"].isin(subgraph_nodes)]
                                                   ["label"].values.tolist()))

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
                    micro_colony_num += 1

                    # fit ellipse
                    ellipse_params = fit_ellipse(bacteria_in_this_micro_colony)
                    # append ellipse
                    ellipses.append(ellipse_params)

                    # calculation of aspect ratio
                    if "Aspect Ratio" in summary_statistic_method:
                        aspect_ratio = aspect_ratio_calc(ellipse_params)
                        # store aspect ratio
                        aspect_ratio_list.append(aspect_ratio)
                        local_aspect_ratio_list.append(aspect_ratio)

                    # calculation of Anisotropy
                    if "Anisotropy" in summary_statistic_method:
                        mean_anisotropy = anisotropy_calc(bacteria_in_this_micro_colony, max_distance_between_cells)
                        # store anisotropy
                        anisotropy_list.append(mean_anisotropy)
                        local_anisotropy_list.append(mean_anisotropy)

                    # calculation of Density
                    if "Density" in summary_statistic_method:
                        density = calc_density(bacteria_in_this_micro_colony)
                        # store density
                        density_list.append(density)
                        local_density_list.append(density)

                    # calculation of Correlate growth penalty based on location in microcolony
                    if "dist_vs_growth_rate" in summary_statistic_method:
                        dist_vs_growth_rate = calc_dist_vs_growth_rate(bacteria_in_this_micro_colony, dt)
                        # store dist vs groeth rate
                        dist_vs_growth_rate_list.append(dist_vs_growth_rate)
                        local_dist_vs_growth_rate_list.append(dist_vs_growth_rate)

                    # calculation of fourier_descriptor in microcolony
                    if "fourier_descriptor" in summary_statistic_method:
                        fourier_descriptor = calc_fourier_descriptor(bacteria_in_this_micro_colony, fig_path,
                                                                     'micro_colony_img' + str(micro_colony_num+1) +
                                                                     '_fill')
                        # store fourier_descriptor
                        fourier_descriptor_list.append(fourier_descriptor)
                        local_fourier_descriptor_list.append(fourier_descriptor)

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
    if "fourier_descriptor" in summary_statistic_method:
        # remove nan values
        fourier_descriptor_list = [x for x in fourier_descriptor_list if str(x) != 'nan']
        report_mean_summary_statistics["fourier_descriptor"] = np.mean(fourier_descriptor_list)

    return report_mean_summary_statistics


def global_analysis(pickle_files_directory, summary_statistic_method, dt, max_distance_between_cells, um_pixel_ratio,
                    fig_path):
    """
    Goal: this is the main function; the function calls other functions described above and
    calculates summary statistics for each time-step and store the outputs in a dictionary.
    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply over time-lapse
    @param dt float interval time
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic
    """

    # store summary statistics over time-lapse (for each individual time-step)
    aspect_ratio_list = []
    anisotropy_list = []
    density_list = []
    dist_vs_growth_rate_list = []
    fourier_descriptor_list = []

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

        # calculation of aspect ratio
        # if  number of bacteria < 2 : ellipse_params = `nan`
        # so: the aspect ratio can not be calculated
        if "Aspect Ratio" in summary_statistic_method and len(cs) > 1:
            aspect_ratio = aspect_ratio_calc(ellipse_params)
            # store aspect ratio
            aspect_ratio_list.append(aspect_ratio)

        # calculation of Anisotropy
        if "Anisotropy" in summary_statistic_method and len(cs) > 1:
            mean_anisotropy = anisotropy_calc(cs, max_distance_between_cells)
            # store anisotropy
            anisotropy_list.append(mean_anisotropy)

        # calculation of Density
        if "Density" in summary_statistic_method:
            density = calc_density(cs)
            # store density
            density_list.append(density)

        # calculation of Correlate growth penalty based on location in microcolony
        if "dist_vs_growth_rate" in summary_statistic_method:
            dist_vs_growth_rate = calc_dist_vs_growth_rate(cs, dt)
            # store dist vs growth rate
            dist_vs_growth_rate_list.append(dist_vs_growth_rate)
        
        # calculation of fourier_descriptor
        if "fourier_descriptor" in summary_statistic_method:
            fourier_descriptor = calc_fourier_descriptor(cs, fig_path, 'timestep_' + str(timestep+1) + '_fill')
            # store fourier descriptor
            fourier_descriptor_list.append(fourier_descriptor)

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
    if "fourier_descriptor" in summary_statistic_method:
        # remove nan values
        fourier_descriptor_list = [x for x in fourier_descriptor_list if str(x) != 'nan']
        report_mean_summary_statistics["fourier_descriptor"] = np.mean(fourier_descriptor_list)

    return report_mean_summary_statistics


def data_analysis(pickle_files_directory, summary_statistic_method, dt, mode='global', max_distance_between_cells=3.4,
                  um_pixel_ratio=0.144, min_size_of_micro_colony=2, fig_path=None):
    """
    goal: calculation of summary statistics in global or local mode
    @param pickle_files_directory str directory of simulation / experimental pickle files
    @param summary_statistic_method   str     the method that we apply on the micro-colonies
    @param dt float interval time
    @param mode str it shows we want to calculate global summary statistics or local summary statistics
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param um_pixel_ratio float  convert um to pixel (requires for calculation of density summary statistic)
    @param min_size_of_micro_colony int minimum size of micro colony (only for local mode)
    @param fig_path str Directory of saved contours as images (optional)
    Return report_mean_summary_statistics dictionary  According to the summary statics calculated for micro colonies,
    the average value of each summary statistic is reported as follows:
    Summary statistic name: average value of summary statistic
    Returns:
    """

    if mode == 'local':
        return micro_colony_analysis(pickle_files_directory, summary_statistic_method, dt, min_size_of_micro_colony,
                                     max_distance_between_cells, um_pixel_ratio, fig_path)
    elif mode == 'global':
        return global_analysis(pickle_files_directory, summary_statistic_method, dt, max_distance_between_cells,
                               um_pixel_ratio, fig_path)
