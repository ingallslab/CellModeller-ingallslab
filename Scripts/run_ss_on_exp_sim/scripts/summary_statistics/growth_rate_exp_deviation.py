import os
import sys
sys.path.append(os.path.join(os.path.expanduser('~'), 'CellModeller-ingallslab/Scripts/summaryStatistics'))
import numpy as np
import CellModeller
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt


from Scripts.run_ss_on_exp_sim.scripts.helper_functions.countPopulations import count_cell_types
from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import create_pickle_list, get_max_cell_type,\
    load_cellStates, read_time_step


'''
The main summary statistic is calculated in get_exp_deviation.
'''
def get_exp_deviation(file_dir_path, dt):
    """
    Summary metric for deviation from an exponential fit.
    Larger standard deviations indicate deviation from exponential growth.
    @param  file_dir_path   string containing path to directory containing .pickle files
    @param  dt              time step
    @return r_squared       correlation coefficient (R**2 value) of regression to ln(population_curve) vs. time
    """
    time_list, population_list = get_population_curve(file_dir_path, dt)
    
    x = np.array(time_list)
    y = np.log(np.array(population_list))
    slope, intercept, r_squared = perform_linear_regression(x, y)
    
    return r_squared

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
        
def perform_linear_regression(x, y):
    """
    Fit a line to [x, y] data. Get regression coefficients and R**2.
    
    @param  x           (1D nparray) X data
    @param  y           (1D nparray) Y data
    @return slope       Slope of linear fit
    @return intercept   Intercept of linear fit
    @return r_squared   R**2 of linear fit
    """
    # Necessary to reshape data for LinearRegression model
    x_reshape = x.reshape((-1,1))
    
    # Linear regression
    model = LinearRegression().fit(x_reshape, y)
    r_squared = model.score(x_reshape, y)
    model_prediction = model.predict(x_reshape)
    slope = model.coef_[0] # slope
    intercept = model.intercept_ # intercept
    r_squared = model.score(x_reshape, y)
    
    return slope, intercept, r_squared
    
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
