#import sys
#sys.path.append('../../../')
import os
import pandas as pd
from os import listdir
from os.path import isfile, join
import glob

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates

# scripts that we are testing
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AspectRatio import calc_aspect_ratio
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.Anisotropy import get_global_order_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.density_calculation import get_density_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.growth_rate_exp_deviation import get_exp_deviation
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.CellOrientationOnBoundary import calc_cell_orientation_on_boundary
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AgeDistanceDistribution import calcAgeDistanceDistribution
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.dyadStructure import calcDyadStructure
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.fourier_descriptor import calc_fourier_descriptor
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.convexity_smart import cal_convexity


"""
produce a csv containing summary stats of scenes
Note: this is a single code that is capable of calculating summary statistics of both exp. and sim.
"""
# helper functions
def get_max_cell_num(pickle_files_directory):
    path = pickle_files_directory + "/*.pickle"
    pickle_files_list = [filename for filename in sorted(glob.glob(path))]
    length = len(pickle_files_list)
    picklefile = pickle_files_list[length-1]
    cs = load_cellStates("", picklefile)
    return len(cs), length

def create_stat_dict(summary_statistic_method_list):
    # store summary statistics over time-lapse (for each individual time-step)
    stats_value = {}
    for summary_stat in summary_statistic_method_list:
        stats_value[summary_stat] = []

    return stats_value


if __name__ == '__main__':

     # summary statistics: "Aspect Ratio", "Anisotropy", "Density" , "dist_vs_growth_rate", "growth_rate_exp_deviation, 
     #                      AgeDistanceDistribution, cell_orientaion_on_boundary, dyadStructure"
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

    # optional arguments
    # 1. mode: local (calculate for each micro colony) or global(calculate over time lapse without considering
    # micro-colonies) (default: 'global')
    mode = 'global'

    # 2. minimum size of micro colony (unit: um - default: 2)
    min_size_of_micro_colony = 2

    # 3. maximum distance between bacteria that can be neighbours. (unit: pixel - default: 3.4)
    max_distance_between_cells = 3.4

    # 4. convert micro meter to pixel (default: 0.144)
    um_pixel_ratio = 0.144

    # interval time: 3 min, 0.05=3/60
    dt = 0.05  #input("interval time (0.05 is experiment default): ")

    # alpha shape for fourier descriptor
    shape = 0.048

    # margin int pixel representing the margin to add to the generated image.
    margin=10

    # plot export path
    fig_path = ''

    # create dataframe to store summary stats    
    dict = {}
    dict['name'] = []
    dict['sim/exp'] = []
    dict['cell_num'] = []
    dict['time'] = []
    dict.update(create_stat_dict(summary_statistic_method_list))
    df = pd.DataFrame.from_dict(dict)
    #print(df)

    # pickle directories path
    dirpath = "../unit test"  #input("input experimental data directory: ")
    onlydirs = [f for f in listdir(dirpath) if not isfile(join(dirpath, f))]
    
    # recursively extract  pickles from each directory            
    print('start running experimental data')
    # recursively extract data from each directory
    for dir in onlydirs:
        print(dir)
        pickle_files_directory = os.path.join(dirpath, dir)
        path = pickle_files_directory + "/*.pickle"
        pickle_files_list = [filename for filename in sorted(glob.glob(path))]
        if len(pickle_files_list) == 0: # which means this is not a pickle folder
            continue

        # get_max_cell_num and last timestep
        cell_num, time = get_max_cell_num(pickle_files_directory)

        # create new row to store 
        new_row = pd.DataFrame({'name': dir,
                                'sim/exp': 1 ,
                                'cell_num': cell_num,
                                'time': time},
                                index=[0])

        # calculate summary stats
        ## summary stats that work on the last phase
        cs_last = load_cellStates("", pickle_files_list[-1])

        # calculation of aspect ratio
        if "Aspect Ratio" in summary_statistic_method_list:
            aspect_ratio = calc_aspect_ratio(cs_last)
            # store
            new_row['Aspect Ratio'] = aspect_ratio
            
        # calculation of Anisotropy
        if "Anisotropy" in summary_statistic_method_list:
            mean_anisotropy = get_global_order_parameter(cs_last)
            # store
            new_row['Anisotropy'] = mean_anisotropy

        # calculation of Density
        if "Density" in summary_statistic_method_list:
            density_parameter = get_density_parameter(cs_last, um_pixel_ratio)
            # store
            new_row['Density'] = density_parameter

        """# calculation of Correlate growth penalty
        if "dist_vs_growth_rate" in summary_statistic_method_list:
            growth_rate_parameter = get_growth_parameter(cs_last, dt, show_plots=False)
            # store anisotropy
            stats_value["dist_vs_growth_rate"].append(growth_rate_parameter)"""

        # calculation of summary statistic for deviation from exponential growth
        if "growth_rate_exp_deviation" in summary_statistic_method_list:
            exp_deviation = get_exp_deviation(pickle_files_directory, dt)
            # store
            new_row['growth_rate_exp_deviation'] = exp_deviation

        # calculation of summary statistic for cell_orientaion_on_boundary
        if "cell_orientaion_on_boundary" in summary_statistic_method_list :
            mean_cell_orientaion = calc_cell_orientation_on_boundary(cs_last, fig_path, dir, False) # to show the plot, use True
            # store
            new_row['cell_orientaion_on_boundary'] = mean_cell_orientaion
            
        # calculation of summary statistic for AgeDistanceDistribution
        if "AgeDistanceDistribution" in summary_statistic_method_list :
            age_distance = calcAgeDistanceDistribution(cs_last, fig_path, dir, False)
            # store
            new_row['AgeDistanceDistribution'] = age_distance

        
        # calculation of summary statistic for dyadStructure
        if "dyadStructure" in summary_statistic_method_list:
            # find the last phase with 2 cells for dyad structure
            for cnt, filename in enumerate(pickle_files_list):
                cs = load_cellStates("", filename)
                if len(cs)==2:
                    next_len = len(load_cellStates("", pickle_files_list[cnt+1]))
                    if next_len > 2 :
                        cs_dyad = cs
                        dyad_stats = calcDyadStructure(cs_dyad)
                        # store
                        new_row['dyadStructure'] = dyad_stats
                        break
    
            
        # calculation of fourier_descriptor
        if "fourier_descriptor" in summary_statistic_method_list:
            fourierDescriptor = calc_fourier_descriptor(cs_last, fig_name = dir, shape = shape, margin=margin)
            # store fourier descriptor
            new_row['fourier_descriptor'] = fourierDescriptor

        # calculation of convexity
        if "convexity" in summary_statistic_method_list:
            convexity, alpha = cal_convexity(cs_last, fig_path, dir, shape, margin)
            # store 
            new_row['convexity'] = convexity
            new_row['alpha'] = alpha

        # add summary stat to dataframe
        df = pd.concat([df,new_row])
    
    # calculating averages
    df_mean = df.groupby(['sim/exp'], as_index=False).mean(numeric_only=True)
    df_mean.at[len(df_mean.index)-1,'name']="mean"
    df_std = df.groupby(['sim/exp'], as_index=False).std(numeric_only=True)
    df_std.at[len(df_std.index)-1,'name']="std"
    #print("mean:\n", df_mean)
    df = pd.merge(df,df_mean, how='outer')
    df = pd.merge(df,df_std, how='outer')
    # rename mean sim and mean exp
    """df.at[len(df.index)-4,'name']="sim_mean"
    df.at[len(df.index)-3,'name']="exp_mean"
    df.at[len(df.index)-2,'name']="sim_std"
    df.at[len(df.index)-1,'name']="exp_std" """

    print(df)
    df.to_csv('exp_summary_stats.csv')






