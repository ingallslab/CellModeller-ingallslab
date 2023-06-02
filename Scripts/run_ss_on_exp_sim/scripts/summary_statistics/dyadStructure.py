import numpy as np
from numpy.linalg import norm

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates

def find_bacteria_center(cs):
    bacteria_center_x = [cs[it].pos[0] for it in cs]
    bacteria_center_y = [cs[it].pos[1] for it in cs]
    return bacteria_center_x, bacteria_center_y


def calcDyadStructure(cs):
    '''
    Function used to measure cell dyad structure (d2 quantity). Method of calculation
    adapted from Storck et. al.

    Inputs:
        cell states must be of colony of 2 cells
    Outputs:
        numpy array (dtype = float) -> dyad structure (d2 quantity)
    '''

    # get cell center positions 
    cell_center_x, cell_center_y = find_bacteria_center(cs) # lists

    cell_1 = np.array([cell_center_x[0],cell_center_y[0]])
    cell_2 = np.array([cell_center_x[1],cell_center_y[1]])

    # get cell_1 pole position
    cell_end_point1_x = [cs[it].ends[0][0] for it in cs] 
    cell_end_point1_y = [cs[it].ends[0][1] for it in cs] 
    #print(cell_end_point1_y)

    x_pole = np.array([cell_end_point1_x[0],cell_end_point1_y[0]])

    cells_vecter = cell_1 - cell_2
    x_orientation = cell_1 - x_pole
    d2 = abs(np.dot(cells_vecter, x_orientation)/(norm(cells_vecter)*norm(x_orientation)))
    return d2


"""
if __name__ == '__main__':
    picklefile = 'step-000011.pickle'
    cs = load_cellStates("", picklefile)
    #print(vars(cs[4]))
    #print(cs)
    if len(cs) == 2:
        print(calcDyadStructure(cs))
"""
