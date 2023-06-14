import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium_reg_param import CLBacterium_reg_param
from CellModeller.GUI import Renderers
import numpy as np

# Calibration parameters
gamma = 100
reg_param = 0.01 

# Physiological parameters
growth_mu = 1.397
growth_sigma = 0.426
division_mu = 5.785
division_sigma = 2.526

def setup(sim):
    # Set biophysics module
    biophys = CLBacterium_reg_param(sim, jitter_z=False, max_cells=5000, reg_param=reg_param, gamma=gamma)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0), dir=(1,0,0))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    
    # Specify how often data is saved
    sim.pickleSteps = 100
    sim.dt = 0.025 #h

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = np.random.normal(division_mu, division_sigma)
    # Specify growth rate of cells
    
    cell.growthRate = np.random.normal(growth_mu, growth_sigma)
    cell.color = (0.0,1.0,0.0)

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.strainRate < 0.01:
            cell.color = (1,0,0)
            #print(cell.cellAge)
        else:
            cell.color = (0,1,0)
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = np.random.normal(division_mu, division_sigma)
    d2.targetVol = np.random.normal(division_mu, division_sigma)
    d1.growthRate = np.random.normal(growth_mu, growth_sigma)
    d2.growthRate = np.random.normal(growth_mu, growth_sigma)
    
def setparams(param_dict):
    global gamma, reg_param, adhesion
    gamma = param_dict['gamma']
    reg_param = param_dict['reg_param']
    #adhesion = param_dict['adhesion']

