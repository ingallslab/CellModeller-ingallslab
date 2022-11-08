import networkx as nx
from scipy.spatial import distance_matrix
import numpy as np
import pandas as pd


def finding_sub_graphs(cs, max_distance_between_cells):
    """
    Goal: Finding the largest disconnected subparagraphs
    @param cs dict bacteria features value in specific time step (global mode)
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    return sub_graphs list  List elements including bacteria ids for each subgraph
    """
    graph = make_graph(cs, max_distance_between_cells)
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


def make_graph(cs, max_distance_between_cells):
    """
    @param cs dict bacteria features value in specific time step (global mode)
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @return graph graph
    """
    # find neighbours
    adjacency_matrix = make_adjacency_matrix(cs, max_distance_between_cells)

    # create graph
    graph = nx.from_pandas_adjacency(adjacency_matrix)

    return graph
