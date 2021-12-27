'''
ConjugationCrossFeeding.py is a combination of Conjugation.py and Tutorial_3.py

This file simulates spread of a conjugative plasmid through a two-strain auxotrophic community growing in a microcolony.
The expectation is that conjugation in a mutualism will lead to a higher conjugation efficiency than in a neutral community.
For now, the simulation starts with one population as donors, and one as recipients. In the future, allow for having a fraction of one initial population as donor.

Cell type = 0 -> recipient (red)
Cell type = 1 -> donor (green)
Cell type = 2 -> transconjugant (blue), conjugating plasmid

TODO: the transconjugants seem to be growing much faster! check their growth rate... Dec 3, 2021 AY

Next steps: create this simulation in a microfluidic chamber with flow (i.e. constant concentration) boundary conditions.
'''

import random
import numpy
import math

#Import basic biophysics
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium

#Import diffusion-related libraries
from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator

#Import GUI functionality
from CellModeller.GUI import Renderers

max_cells = 10000

#Cell length parameters
delta = 1.9
delta_sig = 0.45

#Specify parameter for solving diffusion dynamics
grid_dim = (80, 80, 8) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-160, -160, -16) # where to place the diffusion space onto simulation space

n_signals = 2
n_species = 2

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, max_cells=10000,gamma=10.0, cgs_tol=1E-5,compNeighbours=True)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [10.0, 10.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(-3,0,0), dir=(1,0,0), length=delta,rad=0.4) #acceptor
    sim.addCell(cellType=1, pos=(3,0,0), dir=(1,0,0), length=delta,rad=0.4) #donor

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    sim.addRenderer(sigrend) #Add
    
    # Specify how often data is saved
    sim.pickleSteps = 50

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell.length + random.gauss(delta,delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 0.0
    if cell.cellType == 0: cell.color = (1.0,0.0,0.0)
    elif cell.cellType == 1: cell.color = (0.0,1.0,0.0)
    elif cell.cellType == 2: cell.color = (0.0,0.0,1.0)
    
    # Specify initial concentration of chemical species
    cell.species[:] = [0.0]*n_species
    # Specify initial concentration of signaling molecules
    cell.signals[:] = [0.0]*n_signals

cl_prefix = \
    '''
        const float Da = 1.0f;
        const float Db = 1.0f;
        const float ka = 1.f;
        const float kb = 1.f;
        
        float  alpha_in = species[0];
        float  alpha = signals[0];
        
        float beta_in = species[1];
        float beta = signals[1];
        '''


# Da = diffusion rate of alpha through the cell membrane
# Db = diffusion rate of beta through the cell membrane

def specRateCL(): # Add if/else, new species
    global cl_prefix
    return cl_prefix + '''
        if (cellType!=1){
        rates[0] = ka + Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = Db*(beta-beta_in)*area/gridVolume;
        
        } else {
        rates[0] = Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = kb + Db*(beta-beta_in)*area/gridVolume;
        }
        '''


def sigRateCL(): #Add
    global cl_prefix
    return cl_prefix + '''
        rates[0] = -Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = -Db*(beta-beta_in)*area/gridVolume;
        
        '''

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    v_max = 0.9
    Km = 0.1 #This can be though of as the degree of mutualism - turns mutualism off when set to 0.
    
    for (id, cell) in cells.items():
        #Growth and division
        if cell.cellType==0:
            cell.growthRate = 0.1 + v_max * cell.species[1] / (Km + cell.species[1])
        else:
            cell.growthRate = 0.0 + v_max * cell.species[0] / (Km + cell.species[0])
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
        
        #Conjugation
        infect_chances = 0
        if cell.cellType == 0: #if a cell is an acceptor,
            for index in cell.neighbours:#loop through all contacts
                if cells[index].cellType != 0: #if donor or transconjugant is in contact
                    #if random.random() < 0.01: # constant probability per unit time
                    if random.random() < cells[index].effGrowth/10.0: # probability of infection is proportional to donor growth rate
                        cell.cellType = 2 #become transconjugant

        if cell.cellType == 0: cell.color = (1.0,0.0,0.0)
        elif cell.cellType == 1: cell.color = (0.0,1.0,0.0)
        elif cell.cellType == 2: cell.color = (0.0,0.0,1.0)

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)

