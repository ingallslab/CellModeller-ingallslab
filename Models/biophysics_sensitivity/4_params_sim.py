import os
import pickle
import sys

import numpy as np
import pandas as pd
import re
from CellModeller.Simulator import Simulator
import pyopencl as cl

# Local/custom modules
sys.path.append(os.path.join(os.path.expanduser('~'), 'CellModeller-ingallslab/Scripts/summaryStatistics/biophysics/'))
sys.path.append(os.path.join(os.path.expanduser('~'), 'CellModeller-ingallslab/Scripts/summaryStatistics/'))

cellmodeller_module = "adhesion_reg_param_sim.py"
max_cs = 250

'''
This scripts produce simulation results varying parameters while keeping others default

Instructions:
0. Setup the simulation module accordingly; ensure there is a time step variable dt = x.y; 
                                            ensure corresponding parameters are add in setparams() 
1. Define parametric space in main()
2. Run. This produces pickles of the simulation
3. use get_sum_stat.py to get summary stats from the pickles (just like simulation)
4. run process_summary_stats.py to plot sum stats

Note:

'''

def run_model(parameters, file_name):
    # Defining input/output locations
    export_path = "data/" + file_name
    sys.stdout = open(os.devnull, 'w') # Disable printing from simulations
    
    # Run CellModeller simulation
    simulate(cellmodeller_module, parameters, export_path)
    sys.stdout = sys.__stdout__ # Re-enable printing
    print("Simulation completed")
    
    # Load cstates
    #pickle_list = helperFunctions.create_pickle_list(export_path)
    #cs = helperFunctions.load_cellStates(export_path, pickle_list[-1]) # load last pickle file, and get its cell status
    
    # Extract time parameters from simulation (could hard-code to avoid unnecessary repetition)
    #sim_file = os.path.abspath(cellmodeller_module)
    #dt = sim_time_parameters(sim_file)
    
    """# Calculate summary statistics
    summary_stats = {}
    summary_stats['aspect_ratio'] = calc_aspect_ratio(cs)
    summary_stats['anisotropy'] = get_global_order_parameter(cs)
    summary_stats['density'] = get_density_parameter(cs)
    summary_stats['growth_rate_exp_deviation'] = get_exp_deviation(cs)
    #summary_stats['convexity'] = calc_aspect_ratio(cs)
    summary_stats['jaggedness'] = calc_fourier_descriptor(cs)
    summary_stats['cell_orientaion_on_boundary'] = calc_cell_orientation_on_boundary(cs)
    summary_stats['AgeDistanceDistribution'] = calcAgeDistanceDistribution(cs)
    #summary_stats['dyadStructure'] = calcDyadStructure(cs)"""
    
    
    return export_path 

def simulate(modfilename, params, export_path):
    """
    Same as batch.py with modifications to allow external parameters and termination based on stepNum
    
    @param modfilename  module name (e.g., "Tutorial1.py")
    @param export_path  output directory name
    @param params       dictionary of parameter values passed to the simulation
    """
    # Get module file name and append module to system path
    (path,name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    
    # Extract time parameters from simulation
    sim_file = os.path.abspath(modfilename)
    dt = sim_time_parameters(sim_file)
    
    # Run simulation
    sim = Simulator(modname, dt, clPlatformNum=0, clDeviceNum=0, saveOutput=True, outputDirName=export_path, psweep=True, params=params)
    while len(sim.cellStates) <= max_cs:
        sim.step()
        
    # Export pickle at final timestep
    sim.writePickle()
    
def sim_time_parameters(file_path):
    """
    Extract the simulation time step from the module file. 
    Must be written as dt = X.Y in the module file
    
    @param  file_path   str of file name containing parameters to be modified
    @return dt          Simulation time step specified in CellModeller module file
    """
    search_dt = f'dt = (\d+\.\d+)'
    with open(file_path, "r+") as file:
        file_contents = file.read()
        dt = re.findall(search_dt, file_contents)
    return float(dt[0])
    
def write_to_pickle(data, file_name):
    # Open a file, where you want to store the data
    file = open(file_name, 'wb')

    # Dump information to that file
    pickle.dump(data, file)

    # Close the file
    file.close()

def main():
    # Default parameter values
    gamma_default = 100
    reg_param_default = 1/100
    adh_default = 0
    targetVol_default = 5

    # Define parametric space
    gamma = [10, 20, 40, 100]
    reg_param = [1, 1/10, 1/100]
    adh = [0.01, 0.1, 1, 10, 100] 
    #targetVol = [2, 4, 5, 8]

    """# Sensitivity analysis
    df = pd.DataFrame(columns=[ 'gamma', 
                                'reg_param', 
                                'adhesion',
                                'targetVol',
                                'aspect_ratio',
                                'anisotropy',
                                'density',
                                'growth_rate_exp_deviation',
                                #'convexity',
                                'jaggedness',
                                'cell_orientaion_on_boundary',
                                'AgeDistanceDistribution',
                                #'dyadStructure'
                                ])"""
 
    """# Platform
    platforms = cl.get_platforms()
    print("Select OpenCL platform:")
    for i in range(len(platforms)):
        print(('press '+str(i)+' for '+str(platforms[i])))
    platnum = int(eval(input('Platform Number: ')))

    # Device
    devices = platforms[platnum].get_devices()
    print("Select OpenCL device:")
    for i in range(len(devices)):
        print(('press '+str(i)+' for '+str(devices[i])))
    devnum = int(eval(input('Device Number: ')))"""
    

    
    samples = 4 # how many simulation rounds per variation
    count = 0 # for indexing the file

    for n in range(samples):
        for gamma_i in gamma:
            # Run model
            parameters = {'gamma': gamma_i, 'reg_param': reg_param_default, 'adhesion': adh_default, 'targetVol': targetVol_default}
            export_path = run_model(parameters, "gamma_"+str(gamma_i)+"_"+str(n+count))
            print(export_path)

        for reg_param_i in reg_param:
            # Run model
            parameters = {'gamma': gamma_default, 'reg_param': reg_param_i, 'adhesion': adh_default, 'targetVol': targetVol_default}
            export_path = run_model(parameters, "reg_param_"+str(reg_param_i)+"_"+str(n+count))
            print(export_path)

        for adh_i in adh:
            # Run model
            parameters = {'gamma': gamma_default, 'reg_param': reg_param_default, 'adhesion': adh_i, 'targetVol': targetVol_default}
            export_path = run_model(parameters, "adh_"+str(adh_i)+"_"+str(n+count))
            print(export_path)
        
        """for targetVol_i in targetVol:
            # Run model
            parameters = {'gamma': gamma_default, 'reg_param': reg_param_default, 'adhesion': adh_default, 'targetVol': targetVol_i}
            export_path = run_model(parameters, "targetVol_"+str(targetVol_i)+"_"+str(n+count))
            print(export_path)"""


# Make sure we are running as a script
if __name__ == "__main__": 
    main()
