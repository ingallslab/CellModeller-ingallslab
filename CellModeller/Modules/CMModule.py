"""
@brief Namespace for classes related to module system.
@details Modules form the core of the new CellModeller framework which allows
        for easily defining and sharing cell attributes between different parts
        of the application, synchronizing simulator update steps between various
        modules, and defining and responding to custom cell events.
"""

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import Any, Callable, TypeVar

import numpy as np
import pyopencl as cl

from CellModeller.CellState import CellState
from CellModeller.GUI.Renderers.Renderer import Renderer


class SimParams:
    """
    @brief A custom datatype for sharing sim info across modules.
    @note Properties are intended to be read-only! Modifying these from within
            a module may lead to undefined behaviour.
    """

    ## @brief Whether or not to print to console.
    verbosity_level: Callable[[int], bool]
    ## @brief Number of sim steps = cell.time / (dt * time_factor)
    sim_steps: int = 0
    ## @brief May differ from user-defined value.
    max_cells: int
    ## @brief Primarily for biophysics modules, does not affect mem alloc.
    max_contacts: int
    ## @brief Merged dtype of all modules, for use in reading OpenCL array.
    dtype: np.dtype
    ## @brief Indicates whether OpenCL is being used.
    using_svm: bool
    ## @brief Only available if using_svm is True.
    cl_context: cl.Context | None = None
    ## @brief see above.
    cl_queue: cl.CommandQueue | None = None


class CMModule(ABC):
    """
    @brief Abstract base class for all CellModeller modules.
    @details Exposes methods to Simulator that initialize SimParams, step
            forward the simulation, and initialize new cells. Also sets the
            cell_attrs property that Simulator reads from to assign cell
            attribute fields.
    """

    ## @brief Optional cell attributes that this module defines.
    cell_attrs: np.dtype | None = None

    ## @brief Read-only global simulation values.
    sim_params: SimParams
    ## @brief Function to pause simulation from within a module.
    pause: Callable

    def set_sim_attrs(self, sim_params: SimParams, pause_func: Callable) -> None:
        """
        @brief Called by Simulator.
        @details Calls abstract method on_sim_ready() to initialize module with
                values from SimParams.
        """
        self.sim_params = sim_params
        self.pause = pause_func
        self.on_sim_ready()

    @abstractmethod
    def on_sim_ready(self) -> None:
        """
        @brief Abstract initialization function for modules.
        @details Defer init to here instead of __init__() when appropriate.
        """
        ...

    def sim_step_cl(self, dt: float, cell_array: cl.MemoryObject) -> None:
        """
        @brief Step function for OpenCL kernels.
        @details Called before regular step function. Not required to implement.
        @param dt Timestep in simulation seconds.
        @param array pyopencl buffer to entire cell array as defined in
                SimParams.dtype.
        """
        del dt, cell_array

    @abstractmethod
    def sim_step(self, dt: float, cell_states: Sequence[CellState]) -> None:
        """
        @brief Abstract step function for all modules in simulation.
        @details Called each time step before cell events are handled.
        @param dt Timestep in simulation seconds.
        @param states Limited access to CellArrays via Sequence interface.
        """
        ...

    def new_cell(self, cell: CellState) -> None:
        """
        @brief Initializes a cell in this module.
        @details Can be used to set default cell attributes required by this
                module or update internal states. Not required to implement.
        """
        del cell


## @brief Generic TypeVar for classes implementing cell events.
SelfT = TypeVar("SelfT", bound=CMModule)
## @brief Type signature for functions that respond to cell events.
# @param cells A tuple of cells involved in the event.
EventHandler = Callable[[SelfT, tuple[CellState, ...]], None]
## @brief Type signature for bound functions that do not require a self param.
# @param cells A tuple of cells involved in the event.
BoundEventHandler = Callable[[tuple[CellState, ...]], None]


def cell_event(event: str, priority: int = 0) -> Callable[[EventHandler], EventHandler]:
    """
    @brief Function decorator for custom cell events.
    @details Mark module methods as responding to a specific cell event. Events
            are triggered when the corresponding cell attribute is truthy.
    @param event The cell attribute that acts as a flag for the handled event.
    @param priority Higher values indicate that an event should be processed
            before lower priority events. If modules differ in prioriy for a
            given event, the highest priority is used and all modules handle
            that event at the same priority level.
    """

    def decorator(handler: EventHandler) -> EventHandler:
        setattr(handler, "_cell_event", (event, priority))
        return handler

    return decorator


class UserModule(CMModule):
    """
    @brief Abstract base class for modules that configure simulations.
    @details These are the files users select to run simulations. Creates
            instances of other modules and renderers to attach to Simulator.
    """

    ## @brief Distribution of cells specified as a list of dict.
    # @note Do not try to set cell indices here as that is undefined behaviour.
    #       Cells are added in order, starting from index 0, if that is
    #       relevant to setting up the simulation.
    initial_cells: list[dict[str, Any]] = []
    ## @brief Renderer objects used by the GUI, if running in a GUI.
    renderers: list[Renderer] = []
    ## @brief Modules.
    # @details Can technically include another UserModule too...
    loaded_modules: list[CMModule] = []
    ## @brief Not yet implemented but frequency of saving sim state.
    pickle_steps: int = 50
    ## @brief Sets simulation level of console printing.
    verbosity: int = 0
    ## @brief Sets simulation frequency of console printing.
    output_freq: int = 20
    ## @brief Number of cells at which to pause simulation.
    # @details Also determines how much memory is allocated for cell attributes
    #       so setting it too high might result in an error!
    max_cells: int = 1000
    ## @brief Number of contacts used by biophysics modules.
    max_contacts: int = 10000
    ## @brief Conversion of real time to simulation time.
    # @details Ex. The default value of 400 means that every 0.05 seconds of
    #       real time increments 20 seconds of simulation time. In practice,
    #       3 second real time => 20 minutes simulation => e. coli generation
    #       time at @ growth rate of 5%/min
    time_factor: float = 400
