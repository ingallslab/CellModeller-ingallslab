"""
This script summarizes the total anisotropy in a list of cells
using the orientational order parameter. 
The script can be used for both local and global order parameters.
Order parameter values for certain configurations:
order_parameter is 1/sqrt(N) if all cells are randomly oriented
order_parameter = 1 if all cells are oriented the same way
Order parameter formula from liquid crystal theory applied to the xy plane.
See refs https://www.nature.com/articles/s41467-018-06370-3 and https://arxiv.org/pdf/2003.10509.pdf
"""

import numpy as np

from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.density_calculation import convert_cellmodeller_orientation_to_radians


def get_cell_orientations(cells):
    """
    Create a list of cell angles from CellState.cell_dir in radians
    @param  cells               cellStates dict
    @return cell_orientations   nparray of cell orientations in radians
    """
    # Initialize storage variables
    n_cells = len(cells)
    cell_orientations = np.zeros(n_cells)
    
    # Get cell data
    for i, cell in enumerate(cells.values()):
        cell_orientations[i] = convert_cellmodeller_orientation_to_radians(cell.dir)
        
    return cell_orientations

def get_order_parameter(cell_angles):
    """
    Calculates the order parameter for a list of cells
    
    @param cell_angles      list of cell orientation angles in radians
    @return order_parameter Scalar order parameter (data_type = double)
    """
    # Calculate director vector (average direction of cells in colony)
    director_angle = np.mean(cell_angles)
    director_vector = np.array([np.cos(director_angle), np.sin(director_angle)])
    
    # Calculate angle (theta) between individual cell and director
    # Order parameter is the average of all thetas
    sum_count = 0
    for angle in cell_angles:
        cell_vector = np.array([np.cos(angle), np.sin(angle)])
        theta = np.arccos(np.dot(director_vector, cell_vector)) #no need to divide by mags of each vector since they are unit vectors
        sum_count += np.absolute(np.cos(2*theta))
    order_parameter = np.sqrt(sum_count/len(cell_angles))
    
    return order_parameter
    
def get_global_order_parameter(cells):
    """
    The main function for calculating anisotropy
    @param  cells           cellStates dict
    @return order_parameter Scalar order parameter
    """
    # Run calculation
    angles = get_cell_orientations(cells)
    order_parameter = get_order_parameter(angles)
    
    return order_parameter
