import os
import sys
import pandas as pd
import CellModeller
sys.path.append(r'C:\Users\aaron\CellModeller-ingallslab\Scripts')
import CellProfilerAnalysis  # required for processing experimental data
from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates_full_path, read_time_step, create_pickle_list_full_path
# Colony shape metrics
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AspectRatio import calc_aspect_ratio
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.convexity import cal_convexity
# Interior structure metrics
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.Anisotropy import get_global_order_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.density_calculation import get_density_parameter
# Growth metric
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.growth_rate_exp_deviation import get_exp_deviation, get_norm_growth_rate

'''
USER INPUTS
'''
start_from_step = 60
cells_to_analyze = 200
processing_dir = r'C:\Users\aaron\OneDrive - University of Waterloo\Research\Papers\Biophysics Calibration\experimental_data\e_coli\MG1655\FOSBE_data\processed_data\Output-ImageProcessing+PostProcessing' # path to directory where all experiments are stored
dt = 3/60 #h

def run_summary_stats(experiment, stepNum, n_cells, cells, data_dir, dt):
    summary_stats = {}
    summary_stats["experiment"] = experiment
    summary_stats["stepNum"] = stepNum
    summary_stats["n_cells"] = n_cells
    summary_stats["aspect_ratio"] = calc_aspect_ratio(cells)
    summary_stats['convexity'] = cal_convexity(cells)
    summary_stats["order_parameter"] = get_global_order_parameter(cells)
    summary_stats["density_parameter"] = get_density_parameter(cells)
    summary_stats["R2_exp_fit"] = get_exp_deviation(data_dir, dt)
    summary_stats["growth_rate"] = get_norm_growth_rate(data_dir, dt)
    return summary_stats
    

if __name__ == '__main__':
    # Get list of experiments to process
    experiments = os.listdir(processing_dir)
    print(f'processing {len(experiments)} experiments: {experiments}')

    main_storage_list = []
    for experiment in experiments:       
        # Set input image path
        pickles_dir = 'pickles'
        data_dir = os.path.abspath(os.path.join(processing_dir, experiment, pickles_dir))
        
        # Read pickle file
        pickle_list = create_pickle_list_full_path(data_dir)
        
        for i, pickle in enumerate(pickle_list):
            stepNum = read_time_step(pickle) - 1
            if stepNum < 60:
                continue
                
            cells = load_cellStates_full_path(pickle) # load last pickle file
            n_cells = len(cells)
            
            if n_cells < cells_to_analyze and i < len(pickle_list)-1:
                continue
            elif n_cells < cells_to_analyze and i == len(pickle_list)-1: # Case if there are not enough cells by end of experiment
                print(f'Processing {experiment} {pickle}')
                summary_stats = run_summary_stats(experiment, stepNum, n_cells, cells, data_dir, dt)
                main_storage_list.append(summary_stats)
                break
            else:
                print(f'Processing {experiment} {pickle}')
                summary_stats = run_summary_stats(experiment, stepNum, n_cells, cells, data_dir, dt)
                main_storage_list.append(summary_stats)
                break
            '''
            # Run summary stats
            summary_stats = {}
            summary_stats["experiment"] = experiment
            summary_stats["stepNum"] = stepNum
            summary_stats["n_cells"] = n_cells
            summary_stats["aspect_ratio"] = calc_aspect_ratio(cells)
            summary_stats['convexity'] = cal_convexity(cells)
            summary_stats["order_parameter"] = get_global_order_parameter(cells)
            summary_stats["density_parameter"] = get_density_parameter(cells)
            summary_stats["R2_exp_fit"] = get_exp_deviation(data_dir, dt)
            
            # Append to storage list
            main_storage_list.append(summary_stats)
            
            if n_cells >= cells_to_analyze:
                break
            '''
        
    df = pd.DataFrame(main_storage_list)
    print(f'Exporting {os.path.join(processing_dir, "summary_stats.csv")}')
    df.to_csv(os.path.join(processing_dir, 'summary_stats.csv'), index=False)