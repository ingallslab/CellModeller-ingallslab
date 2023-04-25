import matplotlib.pyplot as plt
import pandas as pd
import re

"""
plot the summary stats of exp and sim data
"""

def get_ensemble_stats(df, params):
    """
    Calculate mean and standard deviation from an ensemble of repeated simulations
    @param  df: dataframe of summary stats from simulation
    @param params a list of parameters
    @return two list of dataframes storing mean, std of simulations grouped by parameters
    """

    # drop unnamed:0 column
    df = df.drop('Unnamed: 0', axis=1)

    # Group by each parameter, then calculate mean and std
    # Return as a DataFrame to enable another groupby downstream
    df_mean_list = []
    df_std_list = []
    for param_i in params:
        df_ensemble_mean = df.groupby([param_i]).mean(numeric_only = True).reset_index()
        df_ensemble_std = df.groupby([param_i]).std(numeric_only = True).reset_index()
        df_mean_list.append(df_ensemble_mean)
        df_std_list.append(df_ensemble_std)
    #print(df_std_list)
    
    return df_mean_list, df_std_list
    
def plot_ensemble_data(df_mean_list, df_std_list, params, exp_df, summary_statistic_method_list):
    """
    plot simulation sum stats based on varying parameters, also plot experiment stats distribution
    @param df_mean_list: list of dataframes that stores sum stats mean from simulation grouped by parameters
    @param df_std_list: list of dataframes that stores sum stats std from simulation grouped by parameters
    @param list of varying parameters
    @param exp_df: a dataframe of the experimental data
    @param summary_statistic_method_list: a list of sum stats
    """    

    # For each value of each parameter, plot summary_stat vs itself
    n = 1
    for i in range(len(params)):
        param_i = params[i] # the parameter we are looking at
        df_mean_i = df_mean_list[i] # the dataframe grouped by param_i
        df_std_i = df_std_list[i]
        for stat in summary_statistic_method_list:
            x = []
            y = []
            yerr = []
            exp_y = [] # experimantal data
            #exp_yerr = []
            for values in df_mean_i[param_i].unique():
                x.append(values)
                y.append(df_mean_i.loc[(df_mean_i[param_i] == values), stat].values[0])    
                yerr.append(df_std_i.loc[(df_std_i[param_i] == values), stat].values[0])

            # exp data
            exp_y = exp_df[stat].values[0:-2]
            print(stat)
            print(exp_y)

            # plot sim data
            plt.errorbar(x, y, yerr=yerr, capsize=6, marker='o', markersize=6, color='k', label="sim")
            plt.ylabel(stat)
            n += 1

            # plot exp data
            plt.scatter([1 for i in range(len(exp_y))], exp_y, color='b', label="exp")
                
            # Set axis to log scale
            plt.xscale('log')
        
            # Add labels and legend
            plt.xlabel(param_i)
            plt.legend()
            
            # Save
            plt.savefig("plots/" + param_i + "_" + stat + ".png") # stats_plots/
              
            plt.close()    
    
def main():
    # from get_summary_stats.py
    summary_statistic_method_list = ["Aspect Ratio", 
                                     "Anisotropy", 
                                     "Density", 
                                     "growth_rate_exp_deviation",
                                     "cell_orientaion_on_boundary", 
                                     "AgeDistanceDistribution", 
                                     "dyadStructure",
                                     "fourier_descriptor",
                                     "convexity"
                                     ]
    
    # 
    sim_sum_stats_path = "sim_summary_stats.csv"
    df = pd.read_csv(sim_sum_stats_path)

    # remove mean and std
    df = df.drop(df.index[-1])
    df = df.drop(df.index[-1])

    # extract the parameter values for each simulation
    params = [  'adhesion',
                'gamma',
                'reg_param',
                'targetVol']
    gamma_default = 100
    reg_param_default = 1/100
    adh_default = 0
    targetVol_default = 5
    adh = []
    gamma = []
    reg = []
    vol = []
    for ind in df.index:
        name_str = df['name'][ind]
        name_list = re.split(r'_', name_str)
        #print(name_list)
        if name_list[0] == 'adh':
            adh.append(float(name_list[1]))
            gamma.append(gamma_default)
            reg.append(reg_param_default)
            vol.append(targetVol_default)
        elif name_list[0] == 'gamma':
            adh.append(adh_default)
            gamma.append(float(name_list[1]))
            reg.append(reg_param_default)
            vol.append(targetVol_default)
        elif name_list[0] == 'reg':
            adh.append(adh_default)
            gamma.append(gamma_default)
            reg.append(float(name_list[2]))
            vol.append(targetVol_default)
        elif name_list[0] == 'targetVol':
            adh.append(adh_default)
            gamma.append(gamma_default)
            reg.append(reg_param_default)
            vol.append(float(name_list[1]))
    
    # add the parameter value columns to df
    df['adhesion'] = adh
    df['gamma'] = gamma
    df['reg_param'] = reg
    df['targetVol'] = vol
    df.to_csv('processed_summary_stats.csv', index=False)

    # gourp summary stats by parameter values    
    df_ensemble_mean_list, df_ensemble_std_list = get_ensemble_stats(df, params)
    #print("mean:\n", df_ensemble_mean_list, "\nstd:\n",  df_ensemble_std_list)
    #df_ensemble_mean_list[0].to_csv('1.csv', index=False)
    #df_ensemble_std_list[0].to_csv('2.csv', index=False)

    # now process experimental data
    exp_sum_stats_path = "exp_summary_stats.csv"
    exp_df = pd.read_csv(exp_sum_stats_path, index_col=0)
    #print(mean_std_row)
    
    # plot
    plot_ensemble_data(df_ensemble_mean_list, df_ensemble_std_list, params, exp_df, summary_statistic_method_list)

if __name__ == '__main__':
    main()
