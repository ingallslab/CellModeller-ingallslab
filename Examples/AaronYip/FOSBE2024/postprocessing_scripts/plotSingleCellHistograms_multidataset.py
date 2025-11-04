"""
This script plots single-cell measurements exported from
ProcessCellProfilerdata.ProcessData() 

Use to get statistics for multiple data sets.

Instructions:
1. Optional: adjust plotting parameters for min/max of x-axis
2. Optional: adjust filtering range for growth rate
3. Run the script and select the .csv files to process
4. Save figures manually
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm, weibull_min
from scipy.stats import lognorm
from tkinter import Tk
from tkinter.filedialog import askopenfilenames

plt.rcParams.update({'font.size': 12})

def plot_radius_pdf(data, fig, distribution='norm', xmin=0, xmax=2.0):
    """
    Plots a histogram of the cell radii as probability density function
    
    @param data     Analyzed data from CelProfiler as dataframe
    @param fig      pyplot figure
    @param xmin     Start plotting the distribution at xmin
    @param xmax     Stop plotting the distribution at xmax
    """
    # Fetch mean readius of each cell
    radius_list = []
    for data_set in data:
        radius_list.append(data_set.groupby(['id'])['radius'].mean())
    radius = pd.concat(radius_list) 
    
    # Create histogram
    plt.hist(radius, bins = 'auto', density = True, facecolor = 'b', alpha = 1)
    fit_dist(radius, xmin, xmax, distribution)
        
    # labels
    plt.xlabel('Radius (um)')
    plt.ylabel('Density')
    #plt.xlim(xmin, xmax)
    
def plot_division_length_pdf(data, fig, distribution='norm', xmin=0, xmax=20, filter_min=2, filter_max=20):
    """
    Plots a histogram of division length as probability density function.
    Currently, the function assumes that 'length' is the total length of the cell;
    CellModeller takes the length as L_total - 2*r
    
    @param data     Analyzed data from CelProfiler as dataframe
    @param fig      pyplot figure
    @param xmin     Start plotting the distribution at xmin
    @param xmax     Stop plotting the distribution at xmax
    """
    # Fetch division length of each cell
    div_length_list = [] # The list of dataframes to fill
    for data_set in data:
        # Filter between min and max values
        div_length_filter = (
            data_set.
            groupby("id", dropna=True).
            filter(lambda x: x['targetVol'].mean() - 2*x['radius'].mean() >= filter_min and x['targetVol'].mean() - 2*x['radius'].mean() <= filter_max)
        )
        div_length_list.append(div_length_filter.groupby(['id'])['targetVol'].mean() - 2*div_length_filter.groupby(['id'])['radius'].mean())   
    division_length = pd.concat(div_length_list)
    
    # plot historgram
    plt.hist(division_length, bins = 'auto', density = True, facecolor = 'b', alpha = 1)
    fit_dist(division_length, xmin, xmax, distribution)
    
    # labels
    plt.xlabel('Division length (um)')
    plt.ylabel('Density')
    plt.xlim(xmin, xmax)
    
def plot_ln_division_length_pdf(data, fig, distribution='norm', xmin=0, xmax=np.log(20)):
    """
    Plots a histogram of division length as probability density function.
    Currently, the function assumes that 'length' is the total length of the cell;
    CellModeller takes the length as L_total - 2*r
    
    @param data     Analyzed data from CelProfiler as dataframe
    @param fig      pyplot figure
    @param xmin     Start plotting the distribution at xmin
    @param xmax     Stop plotting the distribution at xmax
    """
    # Fetch division length of each cell
    div_length_list = [] # The list of dataframes to fill
    for data_set in data:
        div_length_list.append(np.log(data_set.groupby(['id'])['targetVol'].mean() - 2*data_set.groupby(['id'])['radius'].mean()))
    division_length = pd.concat(div_length_list)
    
    # plot historgram
    plt.hist(division_length, bins = 'auto', density = True, facecolor = 'b', alpha = 1)
    fit_dist(division_length, xmin, xmax, distribution)
    
    # labels
    plt.xlabel('ln(Division Length (um))')
    plt.ylabel('Density')
    plt.xlim(xmin, xmax)
    

def plot_growth_rate_pdf(data, fig, filter_min, filter_max, distribution='norm'):
    """
    Plots a histogram of growth rate as probability density function
    
    @param  data        analyzed data from CelProfiler as dataframe
    @param  fig         pyplot figure
    @param  filter_min  minimum growth rate required to be included in analysis
    @param  filter_max  maxiimum growth rate required to be included in analysis
    """
    # Fetch growth rate of each cell
    growth_list = []
    for data_set in data:
        # Filter between min and max values
        growth_rate_filter = (
            data_set.
            groupby("id", dropna=True).
            filter(lambda x: x['growthRate'].mean() >= filter_min and x['growthRate'].mean() <= filter_max)
        )
        growth_list.append(growth_rate_filter.groupby(['id'])['growthRate'].mean())
        
    growth_rate = pd.concat(growth_list)
    growth_rate = growth_rate*60 # Convert to 1/h

    # plot historgram
    xmin = filter_min*60; xmax = filter_max*60
    plt.hist(growth_rate, bins = 25, density = True, facecolor = 'b', alpha = 1)
    fit_dist(growth_rate, xmin, xmax, distribution)

    # labels
    plt.xlabel('Growth Rate (1/h)')
    plt.ylabel('Density')
    plt.xlim(xmin, xmax)
      
def fit_dist(data, xmin, xmax, distribution, pts=100):
    """
    Fits a normal distribution to data and plots
    
    @param data Experiental data
    @param xmin Start plotting the distribution at xmin
    @param xmax Stop plotting the distribution at xmax
    @param distribution The type of distribution to fit
    """
    samples = len(data)
    if distribution == 'norm':
        mu, std = norm.fit(data)
        x = np.linspace(xmin, xmax, pts)
        p = norm.pdf(x, mu, std)
        title = "Fit results: mean = %.3f,  std = %.3f,  n = %d" % (mu, std, samples)
    elif distribution == 'lognorm':
        shape, loc, scale = lognorm.fit(data) # shape = std, scale = exp(mean)
        x = np.linspace(xmin, xmax, pts)
        p = lognorm.pdf(x, shape, loc, scale)
        title = "Fit results: mean = %.2f,  std = %.2f, loc = %.2f,  n = %d" % (scale, shape, loc, samples)
    elif distribution == 'weibull_min':
        shape, loc, scale = weibull_min.fit(data) # shape = std, scale = exp(mean)
        x = np.linspace(xmin, xmax, pts)
        p = weibull_min.pdf(x, shape, loc, scale)
        title = "Fit results: mean = %.2f,  std = %.2f, loc = %.2f,  n = %d" % (scale, shape, loc, samples)
        
    plt.plot(x, p, 'k', linewidth=2)
    
    plt.title(title)
    
if __name__ == "__main__":
    # Read in data
    print('Select desired analysis csv file(s)')
    csv_files = askopenfilenames()
    print(csv_files)
    data = []
    for i, csv_file in enumerate(csv_files):
        df = pd.read_csv(csv_file, usecols=['growthRate','radius','length','targetVol'])
        data.append(pd.read_csv(csv_file))
        data[i].dropna()
    
    # Plot growth rates
    growth_filter_min = 0 #0.25/60
    growth_filter_max = 3/60 #2.1/60 #doubling time of 0.33 h (20 min) [growth_rate = ln(2)/doubling_time]
    growth_fig = plt.figure(1)
    plot_growth_rate_pdf(data, growth_fig, growth_filter_min, growth_filter_max, distribution='norm')
        
    # Plot radii
    radius_fig = plt.figure(2)
    plot_radius_pdf(data, radius_fig, distribution='norm')
    
    # Plot division length
    division_length_fig = plt.figure(3)
    plot_division_length_pdf(data, division_length_fig, distribution='lognorm')
    
    # Plot ln(division length)
    division_length_fig = plt.figure(4)
    plot_ln_division_length_pdf(data, division_length_fig, distribution='norm')
    
    plt.show()