import os
import sys
import re
import uuid
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import pyabc
from CellModeller.Simulator import Simulator

cellmodeller_module = "simulation_module.py"    
max_cells = 100

def model(parameters):
    # Defining input/output locations
    export_path = "data/" + str(uuid.uuid4())
    sys.stdout = open(os.devnull, 'w') # Disable printing from simulations
    
    # Run CellModeller simulation
    simulate(cellmodeller_module, parameters, export_path)
    sys.stdout = sys.__stdout__ # Re-enable printing
    print("Simulation completed")
    
    # Obtain summary statistics from simulation
    dummy_result = np.array([parameters["gamma"], parameters["reg_param"]]) #Just return parameter values until summary statistics are working
    return {"data": dummy_result}

def simulate(modfilename, params, export_path):
    """
    Same as batch.py with modificationsto allow external parameters and termination based on stepNum
    
    @param modfilename  module name (e.g., "Tutorial1.py")
    @param export_path  output directory name
    @param params       dictionary of parameter values passed to the simulation
    """
    (path,name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    
    # Extract time parameters from simulation
    sim_file = os.path.abspath(modfilename)
    dt = sim_time_parameters(sim_file)
    
    # Run simulation
    sim = Simulator(modname, dt, saveOutput=True, outputDirName=export_path, psweep=True, params=params)
    while len(sim.cellStates) <= max_cells:
        sim.step()
        
    # Export pickle at final timestep
    sim.writePickle()
    
def sim_time_parameters(file_path):
    """
    Extract the simulation time step from the module file. 
    Must be written as dt = X.Y in the module file
    
    @param  file_path   str of file name containing parameters to be modified
    """
    search_dt = f'dt = (\d+\.\d+)'
    with open(file_path, "r+") as file:
        file_contents = file.read()
        dt = re.findall(search_dt, file_contents)
    return float(dt[0])
    
if __name__ == '__main__':
    # Define experimental data
    sumstat1 = 10
    sumstat2 = 0.1
    observation = {"data": np.array([sumstat1, sumstat2])} # currently using dummy values
    
    # Define prior distribution for each parameter
    lower_bound = 0
    scale = 1
    prior = pyabc.Distribution(gamma=pyabc.RV("uniform", 1, 99), #lower_bound, scale
                               reg_param=pyabc.RV("uniform", 1e-4, 1)
                               )
    
    # Define distance metric
    distance = pyabc.PNormDistance(p=2) # p=2: Euclidean distance
    
    # Define ABC-SMC settings
    n_cores = 8
    abc = pyabc.ABCSMC(model, 
                       prior, 
                       distance, 
                       population_size=10, 
                       sampler=pyabc.sampler.MulticoreEvalParallelSampler(n_cores))
    
    # Create database                   
    db_path = "results.db"
    abc.new("sqlite:///" + db_path, observation)
    
    # Run
    history = abc.run(minimum_epsilon=0.05, 
                      max_nr_populations=10) 
                      
