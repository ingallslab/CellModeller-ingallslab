import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

def plot_line_default_params(df):
    """
    Show plot of effect of gamma on summary stats when reg_param = 1/gamma
    """    
    #summary_stats = ['aspect_ratio', 'density_parameter', 'order_parameter', 'gr_vs_centroid']
    summary_stats = ['aspect_ratio', 'density_parameter', 'order_parameter', 'growth_parameter']

    # For each value of gamma, plot summary_stat vs. gamma where reg_param = 1/gamma
    n = 1
    for stat in summary_stats:
        x = []
        y = []
        for gamma in df['gamma'].unique():
            fig = plt.figure(n)
            x.append(gamma)
            y.append(df.loc[(df['reg_param'] == 1.0/gamma) & (df['gamma'] == gamma), stat].values[0])          
        plt.plot(x, y, marker='o', markersize=6, color='k')
        plt.ylabel(stat)
        n += 1
            
        # Set axis to log scale
        plt.xscale('log')
    
        # Add labels
        plt.xlabel('gamma')
        
        # Save
        #plt.savefig("gamma%s.png" % (key,t), bbox_inches='tight')
        
        plt.show()       
        
pickle_full_path = 'sensitivity_analysis_results.pickle'
df = pickle.load(open(pickle_full_path, 'rb'))  
plot_line_default_params(df)
