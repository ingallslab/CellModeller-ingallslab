import os
import sys
import re
import uuid
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import pyabc
import GPUtil
import random
import time
from CellModeller.Simulator import Simulator

# Summary statistics imports
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AspectRatio import calc_aspect_ratio
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.Anisotropy import get_global_order_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.density_calculation import get_density_parameter
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.growth_rate_exp_deviation import get_exp_deviation
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.CellOrientationOnBoundary import calc_cell_orientation_on_boundary
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.AgeDistanceDistribution import calcAgeDistanceDistribution
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.dyadStructure import calcDyadStructure
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.fourier_descriptor import calc_fourier_descriptor
from Scripts.run_ss_on_exp_sim.scripts.summary_statistics.convexity_smart import cal_convexity
from Scripts.run_ss_on_exp_sim.scripts.helper_functions import helperFunctions 

# Module
cellmodeller_module = "simulation_module.py"

# GPU parallelization settings
first_population = True
n_gpu = len(GPUtil.getGPUs())
sim_per_gpu = 3
max_load = 0.15 #conjugation_sim takes up 15% max 
max_mem = 0.0225 #conjugation_sim takes up 6% of 12 GB; 2.25% of 32 GB

# Misc
np.seterr(all="ignore") #suppress warnings for cleaner output

def model(parameters):
    summary_stats = {}
    
    # Wait on first population to ensure even loading of GPUs
    global first_population
    if first_population:
        wait = random.uniform(1, 20)
        time.sleep(wait)
        first_population = False
    
    # Get device number
    global n_gpu                     
    devnum = get_devnum(n_gpu)
    
    global max_cells
    global max_time
    # Defining input/output locations
    sim_id = str(uuid.uuid4())
    export_path = "data/" + sim_id
    
    # Extract time parameters from simulation (could hard-code to avoid unnecessary repetition)
    sim_file = os.path.abspath(cellmodeller_module)
    dt = sim_time_parameters(sim_file)
    
    # Run CellModeller simulation
    print(f"Simulation {sim_id} started on device {devnum}")
    sys.stdout = open(os.devnull, 'w') # Disable printing from simulations
    start = time.time()
    try:
        simulate(cellmodeller_module, parameters, export_path, devnum, dt, max_cells=max_cells, max_time=max_time)
        end = time.time()
        total_time = (end - start)/3600
        print(f"Simulation {sim_id} completed on device {devnum} in {total_time} h")
        
        # Load cellStates
        pickle_list = helperFunctions.create_pickle_list_full_path(export_path)
        cells = helperFunctions.load_cellStates_full_path(pickle_list[-1]) # load last pickle file
        
        # Calculate summary statistics
        summary_stats["aspect_ratio"] = calc_aspect_ratio(cells)
        summary_stats["anisotropy"] = get_global_order_parameter(cells)
        summary_stats["density"] = get_density_parameter(cells)
        summary_stats["growth_rate_exp_deviation"] = get_exp_deviation(export_path, dt)
        summary_stats['convexity'] = cal_convexity(cells)
        #summary_stats["fourier_descriptor"] = calc_fourier_descriptor(cells)
        #summary_stats['cell_orientaion_on_boundary'] = calc_cell_orientation_on_boundary(cells)
        #summary_stats['AgeDistanceDistribution'] = calcAgeDistanceDistribution(cells)
        #summary_stats['dyadStructure'] = calcDyadStructure(cells)
    except:
        sys.stdout = sys.__stdout__ # Re-enable printing
        print(f"Simulation {sim_id} crashed on device {devnum}")
        summary_stats["aspect_ratio"] = 0
        summary_stats["anisotropy"] = 0
        summary_stats["density"] = 0
        summary_stats["growth_rate_exp_deviation"] = 0
        summary_stats['convexity'] = 0

    # Write summary stats to file for convenience (optional) 
    file = open(export_path + "/summary_stats.txt","w")
    for key, value in summary_stats.items(): 
        file.write('%s: %.3f\n' % (key, value)) 
    file.close()

    return summary_stats
   
def distance(sim_stats, exp_stats):
    """
    Calculate Euclidean distance between simulations and experimental observations
    
    @param  sim_stats   Dict of simulated summary statistics
    @param  exp_stats   Dict of experimental summary statistics
    @return distance    Euclidean distance between simuluations and experiment
    """
    distance_list = []
    for x in sim_stats.keys():
        diff = np.abs(sim_stats[x] - exp_stats[x])/exp_stats[x]
        distance_list.append(diff)

    distance = np.linalg.norm(distance_list)

    return distance   

def simulate(modfilename, params, export_path, devnum, dt, max_cells=None, max_time=None):
    """
    Same as batch.py with modificationsto allow external parameters and termination based on stepNum
    
    @param modfilename  module name (e.g., "Tutorial1.py")
    @param export_path  output directory name
    @param params       dictionary of parameter values passed to the simulation
    """
    (path, name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    
    # Extract time parameters from simulation
    sim_file = os.path.abspath(modfilename)
    (dt, sim_time) = sim_time_parameters(sim_file)
    
    sim = Simulator(modname, dt, saveOutput=True, outputDirName=export_path, psweep=True, params=params, clDeviceNum=devnum)
        
    # Run simulation and terminate when simulation reaches a certain number of cells
    if max_cells:
        while len(sim.cellStates) <= max_cells:
            sim.step()
            
    if max_time:
        while sim.stepNum*dt <= max_time:
            sim.step()
        
    # Export pickle at final timestep
    sim.writePickle()
    
    # Delete sim object to prevent OOMing?? I guess it works, but I got cell index error!?
    del sim    
        
def sim_time_parameters(file_path):
    """
    Extract the simulation time step from the module file. 
    Must be written as dt = X.Y and sim_time = A.B in the module file
    
    @param  file_path   str of file name containing parameters to be modified
    """
    search_dt = f'dt = (\d+\.\d+)'
    search_sim_time = f'sim_time = (\d+\.\d+)'
    with open(file_path, "r+") as file:
        file_contents = file.read()
        dt = re.findall(search_dt, file_contents)
        sim_time = re.findall(search_sim_time, file_contents)
    return float(dt[0]), float(sim_time[0])
    
def get_devnum(n_gpu):
    """
    Chooses GPU to use. 
    Option 1: choose randomly when all GPUs are available
    Option 2: when not all GPUs available (i.e., when running), chooses based on minimum load
    Option 3: if no GPUs "available" due to miscalculation of max_load or max_mem, chooses randomly
    """
    global max_load, max_mem, sim_per_gpu
    
    # Get list of available GPUs
    gpu_list = GPUtil.getAvailable(order='random', 
                                   limit=n_gpu, 
                                   maxLoad=max_load*sim_per_gpu, 
                                   maxMemory=max_mem*sim_per_gpu
                                  )
    # Choose by lowest load if sims already running
    if len(gpu_list) < n_gpu:
       gpu_list = GPUtil.getAvailable(order='load', 
                                   limit=n_gpu, 
                                   maxLoad=max_load*sim_per_gpu, 
                                   maxMemory=max_mem*sim_per_gpu
                                  )                                       
    
    # Select GPU
    if gpu_list:
        devnum = gpu_list[0]
    else:
        devnum = random.randrange(n_gpu)    
    
    return devnum
    
if __name__ == '__main__':
    # Define simulation termination condition; keep undesired option as None
    max_cells = 200
    max_time = None

    # Define ABC-SMC settings
    n_cores = 8
    population_size = 40 # Number of simulations we need to accept to complete one calibration round
    min_epsilon = 0.05  # Stop if epsilon becomes smaller than this  
    max_populations = 5 #5 #  Number of calibration rounds(stages)
    n_sim_parallel = sim_per_gpu*n_gpu  

    # Define experimental data
    exp_summary_stats = {}
    exp_summary_stats["aspect_ratio"] = 0.706879319359634
    exp_summary_stats["anisotropy"] = 0.75535426433106
    exp_summary_stats["density"] = 0.953001564059211
    exp_summary_stats["growth_rate_exp_deviation"] = 0.985259710208264
    exp_summary_stats["fourier_descriptor"] = 0.0206358837797498
    exp_summary_stats['convexity'] = 0.6 # 0.6 is made up for now.

    # Define prior distribution for each parameter [lb, ub]
    param_config = {'gamma': [10, 1000], 'reg_param': [1/1000, 1/10]}
    prior_distributions = {}
    for parameter_name in param_config.keys():
        param_low = param_config[parameter_name][0]
        param_hi = param_config[parameter_name][1]
        width = abs(param_hi - param_low)
        prior_distributions[parameter_name] = {"type": "uniform", "args": (param_low, width), "kwargs": {}}

    # create instance of `Distribution` class
    prior = pyabc.Distribution()
    prior = prior.from_dictionary_of_dictionaries(prior_distributions)
   
    abc = pyabc.ABCSMC(model, 
                       prior, 
                       distance, 
                       population_size=population_size, 
                       #eps=pyabc.epsilon.MedianEpsilon(initial_epsilon, median_multiplier),
                       sampler=pyabc.sampler.MulticoreEvalParallelSampler(n_sim_parallel))
    db_path = "results.db"
    abc.new("sqlite:///" + db_path, {"data": observation})
    
    # Run
    history = abc.run(minimum_epsilon=min_epsilon, 
                      max_nr_populations=max_populations)
