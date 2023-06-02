import numpy as np
from numpy.linalg import norm
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates

"""
new age distance returning the slope of the scatter plot
"""


def calc_distances_from_centroid(cs):
    '''
    Some auxillary calculations are:
        - distance of cells to colony centroid
    Inputs:
        cs (dictionary of cell status objects)
    Outputs:
        panda dataframe with 'distance' column and cell id as indecies
    '''

    # obtain cell centers
    cell_centers_x = [cs[it].pos[0] for it in cs]
    cell_centers_y = [cs[it].pos[1] for it in cs] # a list of centers
    cell_centers_x = np.array(cell_centers_x)
    cell_centers_y = np.array(cell_centers_y)
    # colony center
    centroid_x = np.sum(cell_centers_x)/np.size(cell_centers_x)
    centroid_y = np.sum(cell_centers_y)/np.size(cell_centers_y)

    # calculate distance of all cells to colony centroid
    xdist = (cell_centers_x - centroid_x).reshape(-1,1)
    ydist = (cell_centers_y - centroid_y).reshape(-1,1)
    coord = np.hstack((xdist,ydist))
    #print(coord)
    distances_from_centroid = norm(coord, axis = 1)
    #print(type(distances_from_centroid))

        
    return distances_from_centroid


def calcAgeDistanceDistribution(cs, fig_path=None, dir=None, display=False):
    '''
    Function that groups cells into their respective age group and calculates
    the distance of each cell in an age group to the cell colony's centroid.
    The average distance to colony centroid is calculated for each age group.

    Inputs:
        cs (dictionary of cell status objects)
        fig_path: string, plot export path
        dir: string, the name of the scene
    Outputs:
        dataframe -> [age, distance]
    Effects:
        display a line plot of the age group and their mean distance (should not if running in unit test)
    '''


    # obtain distances_from_centroid
    distances_from_centroid = calc_distances_from_centroid(cs) # numpy array

    # obtain all ages in cs
    cell_ages = [cs[it].cellAge for it in cs]
    cell_ages = np.array(cell_ages)

    # fit a linear regression line
    model = LinearRegression().fit(cell_ages.reshape(-1, 1), distances_from_centroid)
    slope = model.coef_[0]
    intercept = model.intercept_

    if fig_path:
        # plot the plot
        plt.scatter(cell_ages, distances_from_centroid)
        # plot the regression line
        plt.plot(cell_ages, slope*cell_ages + intercept, color='red')

        plt.xlabel('cell ages')
        plt.ylabel('distance to centroid')
        plt.title('mean distance to centroid of each cell age groups in '+ dir)
        # save plot
        #print('save')
        plt.savefig(fig_path+dir+'_age_distance.png')
        if display:
            plt.show()
        plt.close()

    return slope


"""
if __name__ == '__main__':
    picklefile = 'jan3_step-000097.pickle'
    cs = load_cellStates("", picklefile)
    #print(vars(cs[23]))
    #print(len(cs))
    print(calcAgeDistanceDistribution(cs, "hi" ,True))
"""
