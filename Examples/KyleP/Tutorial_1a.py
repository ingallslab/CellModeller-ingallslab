import random

import numpy as np

from CellModeller.GUI.Renderers.BacterialRenderers.BacillusRenderer import (
    BacillusRenderer,
)
from CellModeller.Modules.Biophysics.BacterialModels.BacillusModel import BacillusModel
from CellModeller.Modules.CMModule import UserModule


class MySimulation(UserModule):
    # Set biophysics module
    loaded_modules = [BacillusModel()]

    # Specify the initial cell and its location in the simulation
    initial_cells = [{"pos": (0, 0, 0), "dir": (1, 0, 0)}]

    # Set scaling factor between real time and simulation time
    time_factor = 400

    # Add some objects to draw the models
    renderers = [BacillusRenderer()]

    # Specify amount of console logging
    verbosity = 5

    # Specify additional cell attributes used by this simulation
    cell_attrs = np.dtype([("target_len", "f")])

    def on_sim_ready(self) -> None:
        return

    def sim_step(self, dt, cell_states) -> None:
        # Iterate through cells and flag those that reach target size for division
        for cell in cell_states:
            if cell.length > cell.target_len:
                cell.divide_flag = True

    def new_cell(self, cell) -> None:
        # Specify mean and distribution of target cell size
        cell.target_len = 3.5 + random.uniform(-0.5, 0.5)
        # Specify growth rate of cells in %/min
        cell.growth_rate = 0.05
        cell.color = (0.0, 1.0, 0.0)
