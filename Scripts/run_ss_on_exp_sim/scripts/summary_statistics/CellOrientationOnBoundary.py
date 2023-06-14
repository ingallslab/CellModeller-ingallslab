import alphashape as ashape
import numpy as np
from scipy import interpolate
from numpy.linalg import norm
import matplotlib.pyplot as plt
import math

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates


"""
idea from https://github.com/ingallslab/bsim-hpc-package/blob/master/examples/approximate-bayesian-computing/abc_helpers.py line 160
Ati fixed essential bug
"""

# alterd from BacteriaFeatures.find_bacteria_endpoints
def find_bacteria_endpoints(cs):
    """
    find the end points of the cells in cs
    @param cs dict bacteria features value

    Returns: 4 dict:   end point 1 x value
                        end point 1 y value
                        end point 2 x value
                        end point 2 y value
    """
    # end points of bacteria
    bacteria_end_point1_x = {}
    bacteria_end_point1_y = {}
    bacteria_end_point2_x = {}
    bacteria_end_point2_y = {}
    for it in cs:
        bacteria_end_point1_x[it] = cs[it].ends[0][0]
        bacteria_end_point1_y[it] = cs[it].ends[0][1]
        bacteria_end_point2_x[it] = cs[it].ends[1][0]
        bacteria_end_point2_y[it] = cs[it].ends[1][1]

    return bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y


def find_bacteria_center(cs):
    """
    find the center of the cells in cs
    @param cs dict bacteria features value

    Returns: 2 lists:   cell center x value
                        cell center y value

    """
    # end points of bacteria
    bacteria_center_x = [cs[it].pos[0] for it in cs]
    bacteria_center_y = [cs[it].pos[1] for it in cs]

    return bacteria_center_x, bacteria_center_y


def calc_cell_orientation_on_boundary(cs, fig_path=None, scene_name=None, display = False):
    '''
    Function that calculates an alphashape around one cell colony to
    obtain a set of cells (points) that defines the boundary and fits a cubic
    spline to generate continuous representation of that boundary.

    Distribution of cell angle with respect to colony boundary is then
    calculated and returned.
    
    goal: calculation of CellOrientationOnBoundary
    @param cs dict bacteria features value  in specific time step (global mode)
    @param fig_path: string, plot export path
    @param scene_name string the name of the scene

    Returns: distribution of cell orientation on boundary in specific time step (global mode)
    Effects: plots a barplot of distribution of the cell orientation (should not if running in unit test)

    '''

    # empty list to save cell angles later on
    angles = []

    # endpoints
    bacteria_x_end_point1, bacteria_y_end_point1, bacteria_x_end_point2, bacteria_y_end_point2 = find_bacteria_endpoints(cs)

    # get cell center positions
    bacteria_center_x, bacteria_center_y = find_bacteria_center(cs)

    # convert list to numpy array
    bacteria_center_x = np.array(bacteria_center_x)
    bacteria_center_y = np.array(bacteria_center_y)

    # create 2 column array containing x and y positions
    cell_center_coordinates = np.hstack((bacteria_center_x.reshape(-1, 1), bacteria_center_y.reshape(-1, 1)))

    # alpha value that determines alphashape fitting
    # larger alpha value than 0.01 will likely overfit boundary
    alpha = 0.01

    # calculate alphashape and obtain (x,y) coordinates that define shape
    shape = ashape.alphashape(cell_center_coordinates, alpha)
    coords = shape.boundary.coords.xy
    #print("coords: ", coords[0][:-1])

    # convert alphashape data to numpy arrays
    alpha_x = np.array(coords[0][:-1])
    alpha_y = np.array(coords[1][:-1])

    # find actual cell id where alphashape coordinates match the cell coordinates
    indices = []
    for i in range(len(alpha_x)):
        for ind in cs:
            if (cs[ind].pos[0] == alpha_x[i]) & (cs[ind].pos[1] == alpha_y[i]):
                indices.append(ind)
                continue

    # remove duplicates, do we need this?
    # indices = list(set(indices))
    #print("indices: ", len(indices), "\n", indices)

    # using cell id, calculate vector along scene_nameection of each cell that
    # makes up alphashape boundary
    diff_x = []
    diff_y = []
    for id in indices:
        diff_x.append(bacteria_x_end_point2[id] - bacteria_x_end_point1[id])
        diff_y.append(bacteria_y_end_point2[id] - bacteria_y_end_point1[id])
    diff_x = np.array(diff_x)
    diff_y = np.array(diff_y)
    cells_on_boundary = np.hstack((diff_x.reshape(-1, 1), diff_y.reshape(-1, 1)))
    #print("cells_on_boundary: ",cells_on_boundary)

    # fit splines to x=f(u) and y=g(u), treating both as periodic. also note
    # that s=0 is needed in order to force the spline fit to pass through all
    # the input points.
    tck, u = interpolate.splprep([alpha_x, alpha_y], s=0, per=1)

    # small step size for calculating vector tangent to boundary
    du = 1e-6

    # loop over all cells on the boundary (positions are saved in u)
    for i in range(len(u)):
        # find point u + du on spline to calculate tangent vector
        xd, yd = interpolate.splev(u[i] + du, tck)

        boundary_vec = np.array([alpha_x[i] - xd, alpha_y[i] - yd])
        cell_vec = cells_on_boundary[i, :]
        #print("cell_vec: ", cell_vec)

        # calculate angle of cell on boundary with respect to boundary
        angleNumerator = cell_vec.dot(boundary_vec)
        angleDenominator = norm(cell_vec) * norm(boundary_vec)
        angle = math.degrees(np.arccos(angleNumerator / angleDenominator))
        angles.append(angle)

    # turn the domian from [0,360] to [0,180]
    for cnt in range(len(angles)):
        if angles[cnt] > 90:
            angles[cnt] = abs(angles[cnt]-180) 

    if fig_path:
        # get distribution of angles with binsize = 360 (360 refers to angle in radian?)
        #print("angles: ", angles)
        y, x = np.histogram(angles, bins=180)

        # plot the distribution
        plt.hist(angles, bins=180)
        # Add labels and title
        plt.xlabel('angles')
        plt.ylabel('Frequency')
        plt.title('distribution of angles of cells on boundary')
        # save the plot
        plt.savefig(fig_path+scene_name+'_cell_orientation.png')
        if display:
            # Show the plot
            plt.show()
        plt.close()

    # omit final bin edge in x because np.histogram returns N+1 edges
    # and N frequencies
    # return [y, x[:-1]]
    return sum(angles) / len(angles)


"""
if __name__ == '__main__':
    picklefile = 'step-00800.pickle'
    cs = load_cellStates("", picklefile)
    # print(vars(cs[31]))
    # print(cs)
    print(calc_cell_orientation_on_boundary(cs, True))
"""
