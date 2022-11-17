# Standard modules
import sys
import numpy as np
import pickle #only needed for testing

# Installed modules 
import CellModeller

# Local/custom modules
from lje.src.lownerJohnEllipse import plot_ellipse, welzl

def get_centroids(cells):    
    """
    Obtain cell controids
    @param  cells           cellStates dict
    @return centroids       nparray of centroids (n_cells, 2)
    """                
    # Initialize storage variables
    n_cells = len(cells.keys())
    centroids = np.zeros((n_cells, 2))
    
    # Collect centroid values from cellStates dict
    for i, cell in enumerate(cells.values()):
        centroids[i][0] = cell.pos[0]
        centroids[i][1] = cell.pos[1]
    
    return centroids    
     
def calc_aspect_ratio(a, b):
    """
    Compute the aspect ratio of an ellipse

    @param  a               (float) major axis of the ellipse
    @param  b               (float) minor axis of the ellipse
    @return aspect_ratio    (float) minor/major axis
    """
    aspect_ratio = b/a

    return aspect_ratio

def get_aspect_ratio(cells):
    """
    The main function that finds the aspect ratio of a colony.
    To be called by pyabc.
    
    @param  cells           cellStates dict
    @return aspect_ratio    aspect_ratio of a colony (float) 
    """
    if len(cells.keys()) >= 1000:
        sys.setrecursionlimit(10000) # necessary to run welzl for colonies with n_cells > 1000   
        
    # Find cell centroids
    centroids = get_centroids(cells)
    
    # Fit ellipse around colony
    center, major, minor, theta = welzl(centroids)
    
    # Calculate aspect ratio
    aspect_ratio = calc_aspect_ratio(major, minor)
    
    return aspect_ratio   
