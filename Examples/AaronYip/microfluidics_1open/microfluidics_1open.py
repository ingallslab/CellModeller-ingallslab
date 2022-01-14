'''
This is a demo simulation for growth in a microfluidic device. We will solve a reaction-diffusion scenario, where cells consume a nutrient field and grow proportionally based on Monod kinetics. The geometry is as follows:

|--------------open-----------------|
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|___________________________________|

Domain size is 20 X 20 um.
'''

import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math
import os #needed for outputting files in PDEsolver

cell_cols = {0:[0,1.0,0], 1:[1.0,0,0], 2:[0,0,1.0]} #RGB cell colours
cell_lens = {0:3.5, 1:3.5, 2:3.5} #target cell lengths
cell_growr = {0:1, 1:1, 2:1.05} #growth rates

#Geometry

len_x = 20.0
len_y = 20.0
num_x = 10.0
num_y = 10.0

#Picklesteps
save_every = 10

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, max_planes=3)
    solverParams = dict(N_x = num_x, N_y = num_y, L_x = len_x, L_y = len_y, mesh_type = 'crossed', rel_tol = 1e-6, pickleSteps = save_every)

    # use this file for reg too
    #regul = ModuleRegulator(sim, sim.moduleName)
    regul = ModuleRegulator(sim, __file__)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None, solverParams)
 
    planeWeight = 1.0
    biophys.addPlane((0,0,0), (1,0,0), planeWeight) 		    # left side
    biophys.addPlane((len_x,0,0), (-1,0,0), planeWeight) 		# right side
    biophys.addPlane((0,0,0), (0,1,0), planeWeight) 			# base
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,1,0), dir=(1,0,0))
    sim.addCell(cellType=1, pos=(6,1,0), dir=(1,0,0))
    sim.addCell(cellType=2, pos=(12,1,0), dir=(1,0,0))


    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = save_every

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell_lens[cell.cellType] + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = cell_growr[cell.cellType]
    cell.color = cell_cols[cell.cellType]

def update(cells): 
    #Iterate through each cell; remove cells exiting open boundaries; flag cells for division;
   
    for (id, cell) in cells.items():
        if cell.pos[1] > (len_y-5.0): #cells must be removed from sim while stil within FEniCS domain. TODO: only include non-removed cells, or remove cells slightly before opening
            cell.deathFlag = True
        if cell.volume > cell.targetVol and cell.deathFlag == False:
            cell.divideFlag = True
            
def solvePDEandGrowth(sim,cells):
    #cells = sim.cellStates 
    #info available in sim.cellStates is in biophys.updateCellState
    '''
    This function requires the simulation object as an input.
    
    get cell centers and volumes
    call PDE solver; return nutrient field u_local 
    update growthrates for each cell based on their u_local; okay to loop through like in Tutorial 3
    TODO: 1) write a function in ModuleRegulator.py similar to kill()
          2) call the regulator.solvePDEandGrowth() somewhere in simulator.py
    '''
    # Reproducing work of WPJS
    phys = sim.phys
    sim.phys.get_cells()
    centers = phys.cell_centers[0:phys.n_cells]
    vols = phys.cell_vols_dev.get()[0:phys.n_cells] #TODO: vols in phys are actually wrong; update that at some point.
    
    filename = os.path.join(sim.outputDirPath, 'step-%05i.pvd' % sim.stepNum)
    directory = ""
    
    u_local = sim.solver.SolvePDE(centers, vols, filename, directory, sim.stepNum)

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = cell_lens[parent.cellType] + random.uniform(0.0,0.5)
    d2.targetVol = cell_lens[parent.cellType] + random.uniform(0.0,0.5)
    
def kill(cell):
    """
    This function specifies how to change a cell's state when it dies.
    """ 
    # define what happens to a cell's state when it dies
    cell.growthRate = 0.0       # dead cells can't grow any more
    cell.divideFlag = False     # dead cells can't divide
