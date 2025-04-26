"""
@brief Main simulator class for CellModeller.
@details Integrates all components of the application:
        - OpenGL GUI and renderers
        - modules and cell events
        - OpenCL arrays and cell attributes
"""

import importlib
import importlib.util
import inspect
import os
import sys
from collections.abc import Callable
from time import perf_counter

import numpy as np

from CellModeller.CellState import CellArrays, CellManager, CellState, StateSequence
from CellModeller.GUI.Renderers.Renderer import Renderer, WriteFunc
from CellModeller.Modules.CMModule import (
    BoundEventHandler,
    CMModule,
    UserModule,
)

## @brief Type signature for functions that are called before cell events.
# @param sim The CellManager instance for creating/deleting and accessing
#       cell states.
# @param cell The cell that triggered the event.
# @return A tuple of cells involved in the event, which may simply be the
#       triggering cell by itself, other cells it interacts with, or even new
#       cells created as a result of the event.
PreEvent = Callable[[CellManager[StateSequence], CellState], tuple[CellState, ...]]
## @brief Type signature for functions that are called after cell events.
# @param sim The CellManager instance for creating/deleting and accessing
#       cell states.
# @param cells A tuple of cells involved in the event.
PostEvent = Callable[[CellManager[StateSequence], tuple[CellState, ...]], None]


def event_decorator(event_name: str, pre_or_post: str) -> Callable:
    def decorator(handler: PreEvent | PostEvent):
        setattr(handler, pre_or_post, event_name)
        return handler

    return decorator


def pre_event(event_name: str) -> Callable[[PreEvent], PreEvent]:
    """
    @brief Function decorator for simulator events.
    @details Once an event is triggered, it is called for all modules that it
            occurs in, but the flag should be turned off in at least on of them,
            or in a simulator event handler, to avoid unintentionally triggering
            again in the next time step.
    @remark pre and post events are not called if a cell event does not occur
            in any loaded_modules!
    @note Unlike cell event handlers, simulator events are declared at the
            top level of files and do not accept a self parameter. There cannot
            be multiple PreEvent or PostEvent handlers for a cell event. Instead,
            the last found instance of each is the only one that will be called.
            The search order is as follows:
                - Simulator this class
                - CMModule in order of loaded_modules
                - UserModule for this simulation
    """
    return event_decorator(event_name, "_pre_event")


def post_event(event_name: str) -> Callable[[PostEvent], PostEvent]:
    """
    @brief Function decorator for simulator events.
    @details see above.
    """
    return event_decorator(event_name, "_post_event")


class ModuleTypeError(Exception):
    """
    @brief Empty error class for passing a specific error to PyGLCMViewer.
    """

    pass


class Simulator(CellManager[CellArrays]):
    """
    @brief Runs the simulation.
    @details Does all of the wonderful stuff it's supposed to do!
    @todo Write some better documentation for this class maybe? (incl. params)
    """

    def load_user_module(self, module_name: str, module_str: str) -> None:
        """
        @brief Loads UserModule.
        @details Retains ability to load from string for when pickles get
                implemented again.
        @param module_name Name of UserModule to load.
        @param module_str Contents of file as string.
        """
        if module_str:
            print("Importing model %s from string" % (module_name))
            spec = importlib.util.spec_from_loader(module_name, loader=None)
            if spec is None:
                raise RuntimeError()
            module = importlib.util.module_from_spec(spec)
            # is this uhhh... safe?
            exec(module_str, module.__dict__)
        else:
            # setup model module
            # In case the module_name is a path to a python file:
            # Get path and file name
            (path, name) = os.path.split(module_name)
            # path = blank; name = module (without the .py)
            # Append path to PYTHONPATH, if no path do nothing
            if path:
                if path not in sys.path:
                    sys.path.append(path)
            # Remove .py extension if present
            module_name = str(name).split(".")[0]
            print("Importing model %s" % (module_name))

            if module_name in sys.modules:
                module = sys.modules[module_name]
                importlib.reload(module)
            else:
                module = __import__(module_name, globals(), locals(), [], 0)

        # search for UserModule subclass in module file
        module_class = None
        for name in dir(module):
            attr = getattr(module, name)
            if (
                inspect.isclass(attr)
                and issubclass(attr, UserModule)
                and not inspect.isabstract(attr)
            ):
                module_class = attr
        if not module_class:
            raise ModuleTypeError()
        self.user_module = module_class()

    user_module: UserModule
    cell_states: CellArrays
    cell_events: dict[str, tuple[int, list[BoundEventHandler]]]

    def __init__(
        self,
        moduleName: str,
        moduleStr: str = "",
        clPlatformNum: int = 0,
        clDeviceNum: int = 0,
        verbosity: int = 0,
        is_gui: bool = False,
        **_,
    ) -> None:
        self.clPlatformNum = clPlatformNum
        self.clDeviceNum = clDeviceNum

        self.stepNum: int = 0
        self.paused: bool = False
        self.write_since_last_buffer_copy: bool = True
        self.start_time: float = perf_counter()

        self.load_user_module(moduleName, moduleStr)
        if not self.user_module.loaded_modules:
            print("No simulation modules loaded!")
        if not self.user_module.renderers and is_gui:
            print("No renderers added to GUI!")
        if not self.user_module.initial_cells:
            print("No cells added to initial distribution!")
        self.pickleSteps = self.user_module.pickle_steps

        self.verbosity = max(verbosity, self.user_module.verbosity)

        # handle edge case of every cell dividing in a single timestep
        self.orig_max_cells = self.user_module.max_cells
        self.user_module.max_cells = 2 ** int(np.ceil(np.log2(self.orig_max_cells)) + 1)

        # merge defined cell attributes and custom cell events
        # also do runtime checks that classes implement ABC valid
        cell_attrs: dict[str, np.dtype] = {}
        self.cell_events = {}
        # renderers first, only cell attributes
        for renderer in self.user_module.renderers:
            if not isinstance(renderer, Renderer):
                raise RuntimeError(f"{renderer.__class__.__name__} not a Renderer")
            if renderer.cell_attrs is not None:
                cell_attrs |= {renderer.__class__.__name__: renderer.cell_attrs}
        # include usermodule last in list (for update loop as well!)
        self.user_module.loaded_modules += [self.user_module]
        for module in self.user_module.loaded_modules:
            if not isinstance(module, CMModule):
                raise RuntimeError(f"{module.__class__.__name__} is not CMModule")
            if module.cell_attrs is not None:
                cell_attrs |= {module.__class__.__name__: module.cell_attrs}
            # check for cell events in module methods
            for name in dir(module):
                attr = getattr(module, name)
                if hasattr(attr, "_cell_event"):
                    event_name, priority = getattr(attr, "_cell_event")
                    if event_name not in self.cell_events:
                        self.cell_events[event_name] = (-priority, [])
                    prev_priority, event_handlers = self.cell_events[event_name]
                    priority = min(-priority, prev_priority)
                    event_handlers += [attr]
                    self.cell_events[event_name] = (prev_priority, event_handlers)
        # sort cell events by priority
        self.cell_events = dict(sorted(self.cell_events.items(), key=lambda x: x[1][0]))

        # find pre and post events, search order this module < loaded < user
        self.pre_handlers: dict[str, PreEvent] = {}
        self.post_handlers: dict[str, PostEvent] = {}
        search_dirs = [self] + self.user_module.loaded_modules
        namespaces = set(module.__class__.__module__ for module in search_dirs)
        for namespace in namespaces:
            for name in dir(sys.modules[namespace]):
                attr = getattr(sys.modules[namespace], name)
                if hasattr(attr, "_pre_event"):
                    self.pre_handlers[getattr(attr, "_pre_event")] = attr
                if hasattr(attr, "_post_event"):
                    self.post_handlers[getattr(attr, "_post_event")] = attr

        # instantiate cell array to construct final dtype
        self.cell_states = CellArrays(
            self.verbosity,
            self.user_module.max_cells,
            self.clPlatformNum,
            self.clDeviceNum,
            cell_attrs=cell_attrs,
        )

        # set global params and pass to all modules
        params = self.user_module.get_sim_attrs()
        params.verbosity = self.verbosity
        params.dtype = self.cell_states.dtype
        params.using_svm = self.cell_states.using_svm
        if self.cell_states.using_svm:
            params.cl_context = self.cell_states.context
            params.cl_queue = self.cell_states.queue
        # this also calls on_sim_ready for each module
        for module in self.user_module.loaded_modules:
            module.set_sim_attrs(params, self.pause)

        # add starting cells to simulation
        with self.cell_states:
            for cell in self.user_module.initial_cells:
                # initializes cell in each module as well
                view = self.new_cell()
                # dict values override defaults however
                for key, val in cell.items():
                    setattr(view, key, val)

    def pause(self) -> None:
        """
        @brief Handle to pause_func for loaded_modules.
        """
        self.paused = True

    def step(self, dt: float) -> None:
        """
        @brief Advance the simulation by dt seconds real time.
        @details This is the main event loop and also most time consuming part
                of the simulation. Converts dt to simulation time before calling
                sim_step functions of each module in the order they are defined
                in the active UserModule, and finally in the UserModule itself.
        """
        dt = self.user_module.time_factor * dt
        self.start_time = perf_counter()

        # OpenCL kernels if available
        if self.cell_states.using_svm:
            for module in self.user_module.loaded_modules:
                module.sim_step_cl(dt, self.cell_states.get_buffer())

        # context manager for mapping cell arrays from OpenCL memory
        with self.cell_states as sim_cells:
            # regular step function for each module
            for module in self.user_module.loaded_modules:
                module.sim_step(dt, sim_cells)

            # check each cell for configured events
            for cell in sim_cells:
                for event in self.cell_events:
                    if getattr(cell, event):
                        # if no pre event found only apply to single cell
                        if event in self.pre_handlers:
                            cells = self.pre_handlers[event](self, cell)
                        else:
                            cells = (cell,)
                        for event_handler in self.cell_events[event][1]:
                            event_handler(cells)
                        if event in self.post_handlers:
                            self.post_handlers[event](self, cells)

            # if exceed user-defined max cell count pause
            if sim_cells.cell_count >= self.orig_max_cells:
                self.paused = True
                print(f"Reached max cell count of {self.orig_max_cells}!")

        self.cell_states.step(dt)
        self.write_since_last_buffer_copy = True
        self.stepNum += 1

        if self.verbosity > 2:
            print(f"sim step took {perf_counter() - self.start_time:.2f}s")

    def write_to_buffer(self, buffers_to_write: dict[str, list[WriteFunc]]) -> None:
        """
        @brief Method to pass cell data to OpenGL renderers.
        @param buffers_to_write Dict of field names and corresponding buffer
                copy functions to pass memoryview of blocks to.
        """
        if self.write_since_last_buffer_copy:
            with self.cell_states.read_blocks(list(buffers_to_write.keys())) as blocks:
                for block, write_funcs in zip(blocks, buffers_to_write.values()):
                    for write_func in write_funcs:
                        write_func(block)

            self.write_since_last_buffer_copy = False

    def new_cell(self) -> CellState:
        """
        @brief Required by CellManager abstract class.
        @details Adds a new cell to cell arrays and calls new_cell() on each
                attached module to initialize it.
        @return The newly created CellState.
        """
        cell = self.cell_states.new_cell()
        for module in self.user_module.loaded_modules:
            module.new_cell(cell)
        return cell

    def del_cell(self, cell: CellState) -> None:
        """
        @brief Required by CellManager abstract class.
        @details Does not interact with modules as modules are not required to
                define a method for removing cells. However, they may handle
                the remove_flag event if they need to.
        """
        self.cell_states.del_cell(cell)


@pre_event("divide_flag")
def divide(sim: CellManager[StateSequence], cell: CellState) -> tuple[CellState, ...]:
    """
    @brief PreEvent for generating daughter cells for divide_flag.
    @note Simply marking a flag will not trigger pre or post events unless at
            least 1 module handles that event, even if it is an empty function.
    """
    daughter1 = sim.new_cell()
    daughter2 = sim.new_cell()
    cell.remove_flag = True
    return (cell, daughter1, daughter2)


@post_event("remove_flag")
def remove(sim: CellManager[StateSequence], cells: tuple[CellState, ...]):
    """
    @brief PostEvent for removing cells due to remove_flag.
    @note see above.
    """
    for cell in cells:
        sim.del_cell(cell)
