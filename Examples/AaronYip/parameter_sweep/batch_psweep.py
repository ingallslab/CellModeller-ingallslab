import os
import sys
import subprocess
import string
import shutil

from CellModeller.Simulator import Simulator

dt = 0.025
max_cells = 100
pickleSteps = 10
cell_buffer = 0

def simulate(modfilename, export_path):
    (path,name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, dt, pickleSteps=pickleSteps, saveOutput=True, outputDirName=export_path, psweep=True)
    while len(sim.cellStates) < max_cells-cell_buffer:
        sim.step()

def main():
    if len(sys.argv)<2:
        print("Please specify a model (.py) file")
        exit(0)
    else:
        moduleName = sys.argv[1]
    
    export_path = sys.argv[2]
    simulate(moduleName, export_path)

if __name__ == "__main__": 
    main()
