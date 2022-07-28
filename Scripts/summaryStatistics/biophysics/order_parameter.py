"""
This script summarizes the total anisotropy in a list of cells
using the orientational order parameter. 
The script can be used for both local and global order parameters.

Order parameter values for certain configurations:
order_parameter \approx 0 if all cells are randomly oriented
order_parameter = 1 if all cells are oriented the same way

Order parameter formula from liquid crystal theory: https://en.wikipedia.org/wiki/Liquid_crystal#Order_parameter
"""

import pickle
import pandas as pd
import numpy as np
import math

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
        sum_count += (3*np.cos(theta)**2 - 1)/2
    order_parameter = sum_count/len(cell_angles)
    
    return order_parameter
    
if __name__ == "__main__":
    """
    For testing purposes
    """
    angles = [math.pi/2, math.pi/2] #list of angles for each cell in radians
    order_parameter = get_order_parameter(angles)
    print('order_parameter = ', order_parameter)