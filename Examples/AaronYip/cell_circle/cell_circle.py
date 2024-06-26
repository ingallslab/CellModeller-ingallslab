import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

R = 10
n_cells = 15
theta_array = numpy.linspace(0, 2*math.pi, n_cells)

def setup(sim):
    # Set biophysics module
    biophys = CLBacterium(sim, jitter_z=False)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    for theta in theta_array:
        cell_pos = (R*math.cos(theta), R*math.sin(theta), 0)
        cell_dir = (-math.sin(theta), math.cos(theta), 0)
        sim.addCell(cellType=0, pos=cell_pos, dir=cell_dir)

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    
    # Specify how often data is saved
    sim.pickleSteps = 100

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 4
    # Specify growth rate of cells
    cell.growthRate = 0
    cell.color = (0.0,1.0,0.0)

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

