import os
import sys
sys.path.append(os.path.join(os.path.expanduser('~'), 'CellModeller-ingallslab/Scripts/summaryStatistics'))
import numpy as np
import CellModeller
from scipy.optimize import curve_fit
from countPopulations import count_cell_types
from helperFunctions import create_pickle_list, read_time_step, load_cellStates, get_max_cell_type
import matplotlib.pyplot as plt

'''
The main summary statistic is calculated in get_exp_deviation.
'''
def get_exp_deviation(file_dir_path, dt):
    """
    Summary metric for deviation from an exponential fit.
    Larger standard deviations indicate deviation from exponential growth.
    @param  file_dir_path   string containing path to directory containing .pickle files
    @param  dt              time step
    @return std_residuals   standard deviation of residuals of exponential fit
    """
    time_list, population_list = get_population_curve(file_dir_path, dt)
    std_residuals = exp_error_std(time_list, population_list)
    
    return std_residuals

def get_population_curve(file_dir_path, dt):
    """
    Obtain total number of cells at each time point.
    @param  pickle_list     list of pickle file names
    @param  dt              time step
    @return time_list       list of time points
    @return population_list list of total cell count at each time point
    """
    # Preliminary actions
    pickle_list = create_pickle_list(file_dir_path)
    max_cell_type = get_max_cell_type(load_cellStates(file_dir_path, pickle_list[0]))

    # Read time-step from pickle file and count populations; store in a list
    time_list = []
    population_list = []
    for file in pickle_list:
        # Get time
        time = read_time_step(file)*dt
        time_list.append(time)
        # Read all populations
        cells = load_cellStates(file_dir_path, file)
        populations = count_cell_types(cells, max_cell_type)
        total_population = np.sum(populations)
        population_list.append(total_population)
        
    return time_list, population_list  
    
def exp_error_std(xdata, ydata):
    """
    Calculate standard deviation of sum of squared residuals of an exponential fit.
    @param  xdata           list of data for x-axis
    @param  ydata           list of data for y-axis
    @return std_residual    standard deviation of residuals
    """
    opt_parms, parm_cov = curve_fit(exp_func, xdata, ydata, maxfev=10000)
    f = 1*np.exp(opt_parms[0]*np.array(xdata)) #exponential growth equation
    residuals = np.array(ydata) - f
    std_residuals = np.std(residuals)
    
    return std_residuals
    
def exp_func(t, K):
    """
    Standard exponential function
    """
    return np.exp(K * t)
    
"""
For plotting purposes only
"""
def get_growth_rate(t, y):
    """
    Fit exponential function to data
    """
    opt_parms, parm_cov = curve_fit(exp_func, t, y, maxfev=10000)
    K = opt_parms[0]
    
    return K

def plot_population_curve(time_list, population_list):
    """
    Add a population curve to a plot.
    @param  time_list       list of time points
    @param  population_list list of total cell count at each time point
    """
    plt.plot(time_list, population_list, 'o')
    
def plot_fit(time_list, growth_rate):
    t = np.asarray(time_list)
    pop = np.asarray(population_list)
    start = min(time_list)
    end = max(time_list)
    t = np.linspace(start, end, 100)
    y = 1*np.exp(growth_rate*t)
    plt.plot(t, y, '--')    
