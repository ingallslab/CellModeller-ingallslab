from time import perf_counter

import numpy

from CellModeller.GUI.Renderers.BacterialRenderers.BacillusRenderer import (
    BacillusRenderer,
)
from CellModeller.Modules.Biophysics.BacterialModels.FilamentingModel import (
    FilamentingModel,
)
from CellModeller.Modules.CMModule import UserModule, cell_event

"""
This example demonstrates the behaviour of the newly implemented filamenting
cells in a simulation with non-filamenting cells and their interactions. The
restoring_strength attribute can be tuned between 0 and 1 to control how quickly
the filamenting cell straightens out. Currently filaments, joints, and
neighbouring_filaments have to be set manually at initialization, but eventually
add some sort of "generate_filament_flag" cell event to FilamentingModel that
will automate this process and allow for exploring different filament
segmentation patterns.

Change max_gen to control when the simulation pauses after reaching a cell count
of 2**max_gen.
"""


class NewBenchmark(UserModule):
    loaded_modules = [FilamentingModel()]
    renderers = [BacillusRenderer(wireframe=True)]
    initial_cells = [
        {
            "pos": (-0.25, 0, 0),
            "dir": (1, 0, 0),
            "filamenting": True,
            "joints": [1, 3],
            "neighbouring_filaments": [2, 4],
            "growth_rate": 0,
            "color": (1, 0, 0),
        },
        {
            "pos": (2, 0, 0),
            "dir": (-1, 0, 0),
            "length": 0,
            "filamenting": True,
            "joint": True,
            "restoring_strength": 0.1,
            "filaments": [0, 2],
            "growth_rate": 0,
            "color": (1, 0, 0),
        },
        {
            "pos": (2, -2.25, 0),
            "dir": (0, -1, 0),
            "filamenting": True,
            "joints": [0, 1],
            "neighbouring_filaments": [-1, -1],
            "growth_rate": 0,
            "color": (1, 0, 0),
        },
        {
            "pos": (-2.5, 0, 0),
            "dir": (1, 0, 0),
            "length": 0,
            "filamenting": True,
            "joint": True,
            "restoring_strength": 0.1,
            "filaments": [0, 4],
            "growth_rate": 0,
            "color": (1, 0, 0),
        },
        {
            "pos": (-2.5, 2.25, 0),
            "dir": (0, -1, 0),
            "filamenting": True,
            "joints": [0, 3],
            "neighbouring_filaments": [-1, -1],
            "growth_rate": 0,
            "color": (1, 0, 0),
        },
        {
            "pos": (5, -1, 0),
            "dir": (1, 1.5, 0),
            "growth_rate": 0.15,
            "color": (0, 1, 0),
        },
        {
            "pos": (1, 3, 0),
            "dir": (1, 0.5, 0),
            "growth_rate": 0.15,
            "color": (0, 1, 0),
        },
        {
            "pos": (-1, -3, 0),
            "dir": (1, 1, 0),
            "growth_rate": 0.15,
            "color": (0, 1, 0),
        },
    ]
    cell_attrs = numpy.dtype([("target_vol", "f")])

    def __init__(self) -> None:
        self.max_gen = 6
        self.max_gen = 2**self.max_gen
        self.paused_once = False

    def on_sim_ready(self) -> None:
        self.start = perf_counter()

    def new_cell(self, cell) -> None:
        cell.target_vol = 3.5

    def sim_step(self, dt, cell_states) -> None:
        del dt
        for cell in cell_states:
            if cell.length > cell.target_vol:
                cell.divide_flag = True

    @cell_event("divide_flag")
    def stop_sim(self, cells):
        parent, d1, d2 = cells
        d1.growth_rate = parent.growth_rate
        d2.growth_rate = parent.growth_rate
        d1.color = parent.color
        d2.color = parent.color
        if parent.cell_count > self.max_gen and not self.paused_once:
            self.pause()
            print(f"Elapsed time: {perf_counter() - self.start:.2f} seconds")
            self.paused_once = True
