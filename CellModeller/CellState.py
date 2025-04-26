"""
@brief Namespace for classes that manage cell data.
@details Cells define attributes that different components of the application
        use to share details of cells. Can be extended to include additional
        and extend the behaviour of certain modules with cell events that are
        triggered upon setting associated flags.
"""

import warnings
from abc import ABC, abstractmethod
from collections.abc import Sequence
from contextlib import contextmanager
from typing import Any, Callable, ContextManager, Generic, Iterator, TypeVar, overload

import numpy as np
import pyopencl as cl
from numpy.typing import NDArray


class CellState:
    """
    @brief Creates a view to a specific index in the underlying numpy array.
    @details Known cell attributes are accessed as fields of the array, attempt
            to access anything not defined in modules will store the attribute
            in a dictionary and raise a warning. Attributes accessed in this
            way cannot be used in OpenCL.
    @param array Numpy void dtype used to access main array.
    @param extra_array Array of python objects holding a dict for each cell.
    @param index Index of cell to view in array.
    @param time Global simulation time, can be accessed per cell. Not the same
            as cell age!
    @param cell_count Global cell count, can be accessed per cell.
    @note Index is technically a mutable property to allow accessing different
            cells without creating a new CellState instance. Avoid using this
            feature in practice for clarity.
    """

    ## @brief Attributes part of this CellState instance, not numpy array.
    __slots__ = ("_array", "_extra_array", "index", "time", "cell_count")
    ## @brief Private dict of defined attribute fields.
    _fields: dict[str, Callable[["CellState", Any], None]] = {}

    def __init__(
        self,
        array: np.void,
        extra_array: np.ndarray,
        index: int,
        time: float,
        cell_count: int,
    ) -> None:
        self._array = array
        self._extra_array = extra_array
        self.index = index
        self.time = time
        self.cell_count = cell_count

    def __getattr__(self, name: str) -> Any:
        """
        @brief Override default python attribute access method.
        @details Only called if __getattribute__() raises an AttributeError,
                which happens if no @property getter is found for a key.
        """
        return self._extra_array[self.index][name]

    def __setattr__(self, name: str, value: Any) -> None:
        """
        @brief Override default python attribute setter.
        @details Unlike __getattr__(), this is called on all keys, so must
                explicitly check for attribute fields and __slots__, fallback
                to dict for keys that are not found in either.
        """
        if name in self._fields:
            self._fields[name](self, value)
        elif name in self.__slots__:
            object.__setattr__(self, name, value)
        else:
            if not isinstance(self._extra_array[self.index], dict):
                self._extra_array[self.index] = dict()
            self._extra_array[self.index][name] = value
            warnings.warn(
                f"{name} is not a known key for CellState, writing anyway. See CellState docs for more info.",
                stacklevel=2,
            )

    def __str__(self) -> str:
        out = {"index": self.index, "time": self.time}
        for field in self._fields:
            out |= {field: getattr(self, field)}

        return str(out)

    def __repr__(self) -> str:
        return f"Cell {self.index}"

    @staticmethod
    def array_get(field: str) -> Callable[["CellState"], Any]:
        return lambda self: self._array[field][self.index]

    @staticmethod
    def array_set(field: str) -> Callable[["CellState", Any], None]:
        return lambda self, val: self._array[field].__setitem__(self.index, val)

    @classmethod
    def make_with_fields(cls, *field_names: str) -> type["CellState"]:
        """
        @brief Class factory method to bake dtype fields into access structure.
        @details Create and attach a property getter for each name, and put the
                setter in a class-level attribute dictionary for fast lookup.
        @param field_names Unpacked as a list of strings.
        @return A subclass of CellState that can be instantiated with included
                access to _array with defined fields.
        """
        for field in field_names:
            setattr(cls, field, property(cls.array_get(field)))
            cls._fields |= {field: cls.array_set(field)}

        return cls


## @brief Generic TypeVar for classes implementing collections.abc.Sequence.
StateSequence = TypeVar("StateSequence", bound=Sequence[CellState])


class CellManager(ABC, Generic[StateSequence]):
    """
    @brief Abstract type for classes that modify simulation cells.
    @details Extends Sequence[CellState] to include methods that can add and
            remove cells from the simulation. For use in sim event handlers.
            See PreEvent and PostEvent docs for more details.
    """

    cell_states: StateSequence

    @abstractmethod
    def new_cell(self) -> CellState:
        """
        @brief Call necessary methods to initialize a new cell and return it.
        """
        ...

    @abstractmethod
    def del_cell(self, cell: CellState) -> None:
        """
        @brief Call necessary methods to remove a cell from the simulation.
        """
        ...


class CellArrays(Sequence[CellState]):
    """
    @brief Manage allocation and access to numpy array for cell data.
    @details Initializes OpenCL context and queue, implements a context manager
            for mapping memory between host and compute layer. While inside of
            the context manager, behaves as a python sequence of cells, to
            allow indexed access of cells and iterating over all active cells.
    @param max_cells To allocate memory for.
    @param platform_id To initialize the OpenCL context.
    @param device_id see above.
    @param cell_attrs Cell attributes by name of module defining them, in order
            to allow contiguous block access to renderer-defined attributes.
    @param use_svm Whether or not to attempt to create OpenCL arrays using
            shared virtual memory (SVM), which may not be available on all
            platforms in which case a message will be printed and numpy will
            be used to manage cell data instead.
    @note Any index can be accessed with __getitem__() but ONLY ACTIVE CELLS
            are returned by __iter__(). Therefore, there is a discrepancy
            between the value returned by len() and the length of the Iterator.
            Use cell_states[0].cell_count to get the latter.
    @bug __getitem__() not correctly implemented for slice indices yet.
    """

    _attrs_ndarray: np.ndarray
    _attrs_svmptr: cl.SVM
    _attrs_ctxmgr: ContextManager[NDArray[np.void]]
    _attrs_array: np.void | None
    _extra_attrs: np.ndarray
    _block_locs: dict[str, tuple[int, int]]
    _zero_vals: dict[str, np.void]
    _makecellstate: Callable[[int, np.void], CellState]

    ## @brief The constructed CellState subclass with fields.
    cellstatetype: type[CellState]
    ## @brief Default attributes to use for all simulations.
    # @details Referenced in Simulator or used internally.
    default_dtype = {
        "CellArrays": np.dtype(
            [
                ("_active_flag", "?"),
                ("divide_flag", "?"),
                ("remove_flag", "?"),
            ]
        )
    }

    def __init__(
        self,
        max_cells: int,
        platform_id: int,
        device_id: int,
        cell_attrs: dict[str, np.dtype] = {},
        use_svm: bool = True,
    ) -> None:
        # global sim values passed to cellstates
        self.time = 0.0
        self.cell_count = 0

        # private array management
        self._free_indices: list[int] = []
        self._last_index: int = 0

        # merge module attributes into a single np.dtype
        # tuple with start/end of contiguous attributes
        self._block_locs = {}
        # default value used to reset new cells
        self._zero_vals = {}
        # dtype is constructed with list[tuple[name, format, shape]]
        attrs_dtype = []
        offset = 0
        for module, dtype in (self.default_dtype | cell_attrs).items():
            self._block_locs |= {module: (offset, dtype.itemsize * max_cells)}
            for name, fmt, *rest in dtype.descr:
                shape = rest[0] if rest else ()
                attrs_dtype += [(name, fmt, (max_cells,) + shape)]
                self._zero_vals[name] = np.zeros(1, dtype=fmt)[0]
                self._block_locs |= {name: (offset, np.dtype(fmt).itemsize * max_cells)}
            offset += dtype.itemsize * max_cells
        # this will raise an error if attempt to define duplicate attribute names
        self.dtype = np.dtype(attrs_dtype)

        # create cellstate access class using dtype fields
        if self.dtype.names is not None:
            self.cellstatetype = CellState.make_with_fields(*self.dtype.names)
        else:
            raise RuntimeError("No fields in CellState.")
        self._makecellstate = lambda i, arr: self.cellstatetype(
            arr, self._extra_attrs, i, self.time, self.cell_count
        )

        # create OpenCL context
        self.platform: cl.Platform = cl.get_platforms()[platform_id]
        self.device: cl.Device = self.platform.get_devices()[device_id]
        self.context: cl.Context = cl.Context(
            devices=[self.device],
            properties=[(cl.context_properties.PLATFORM, self.platform)],
        )
        self.queue: cl.CommandQueue = cl.CommandQueue(self.context)
        self.using_svm = use_svm
        if use_svm and not (
            self.device.svm_capabilities
            & cl.device_svm_capabilities.COARSE_GRAIN_BUFFER
        ):
            self.using_svm = False
            print(
                "OpenCL device does not support SVM (shared virtual memory).",
                # "Didn't bother to implement another way of doing it you can if you want to."
                # + "\\_(0.0)_/",
                "Defaulting to numpy.ndarray backed cell data, OpenCL kernels may not work.",
            )

        # create arrays
        self._attrs_array = None
        self._extra_attrs = np.zeros(max_cells, dtype=np.object_)
        if self.using_svm:
            # cannot create elements with size != a power of 2 upto some value
            # so initialize as array of void type and read with np.frombuffer
            alignment: int = self.device.preferred_global_atomic_alignment
            self._attrs_svmptr = cl.SVM(
                cl.csvm_empty(
                    self.context,
                    int(np.ceil(self.dtype.itemsize / alignment)),
                    dtype=np.dtype(f"{alignment}V"),
                )
            )
        else:
            self._attrs_ndarray = np.zeros(1, self.dtype)

    def __enter__(self) -> "CellArrays":
        if self.using_svm:
            self._attrs_ctxmgr = self._attrs_svmptr.map_rw(self.queue)
            array = self._attrs_ctxmgr.__enter__()
            self._attrs_array = np.frombuffer(array.data, dtype=self.dtype, count=1)[0]
        else:
            self._attrs_array = self._attrs_ndarray[0]
        return self

    def __exit__(self, type, value, traceback) -> None:
        self._attrs_array = None
        if self.using_svm:
            self._attrs_ctxmgr.__exit__(type, value, traceback)

    @overload
    def __getitem__(self, index: int) -> CellState: ...

    @overload
    def __getitem__(self, index: slice) -> Sequence[CellState]: ...

    def __getitem__(self, index: int | slice) -> CellState | Sequence[CellState]:
        if self._attrs_array is None:
            raise ValueError("SVM not mapped, use with statement to access cells!")
        if isinstance(index, int):
            return self._makecellstate(index, self._attrs_array)
        elif isinstance(index, slice):
            # FIXME: unsure how to store slices without duplicating entire obj
            return self

    def __setitem__(self, index: int, new_cell: CellState) -> None:
        old_cell = self[index]
        for field in old_cell._fields:
            setattr(old_cell, field, getattr(new_cell, field))

    def __len__(self) -> int:
        return self._last_index

    def __iter__(self) -> Iterator[CellState]:
        if self._attrs_array is None:
            raise ValueError("SVM not mapped, use with statement to access cells!")
        view = self._makecellstate(0, self._attrs_array)
        for i in range(self._last_index):
            view.index = i
            if view._active_flag:
                yield view

    def __str__(self) -> str:
        out = "\n".join(str(cell) for cell in self)
        return out

    def __repr__(self) -> str:
        return str([repr(cell) for cell in self])

    def get_buffer(self) -> cl.MemoryObject:
        """
        @brief Get pyopencl MemoryObject for use in kernel computations.
        """
        return self._attrs_svmptr.as_buffer(self.context)

    @contextmanager
    def read_blocks(self, fields: list[str]) -> Iterator[list[memoryview]]:
        """
        @brief Return python buffers corresponding to blocks of cell data.
        @details Can be used to access individual cell attributes or entire
                sets of attributes defined together in a single module.
        @param fields List of fields or module names to look up.
        @return A list of memoryviews as a context manager so that OpenCL
                memory can be synced afterwards.
        """

        def gen_blocks(arr: np.ndarray):
            view = arr.data.cast("B")
            return [view[slice(*self._block_locs[field])] for field in fields]

        if self.using_svm:
            with self._attrs_svmptr.map_ro(self.queue) as ro_array:
                yield gen_blocks(ro_array)
        else:
            yield gen_blocks(self._attrs_ndarray)

    def step(self, dt: float) -> None:
        """
        @brief Step global time.
        """
        self.time += dt

    def new_cell(self) -> CellState:
        """
        @brief Add a cell to the array.
        @details First attempt to access any indices left behind by removed
                cells, otherwise extend used portion of allocated array. Assume
                all values are zeroed out, either by initial allocation or when
                removing cells with del_cell().
        """
        if len(self._free_indices) > 0:
            index = self._free_indices.pop()
        else:
            index = self._last_index
            self._last_index += 1

        cell = self[index]
        cell._active_flag = True
        self.cell_count += 1

        return cell

    def del_cell(self, cell: CellState) -> None:
        """
        @brief Remove a cell from the array.
        @details Mark a cell as inactive and add its index to be reused by new
                cells. Also clean up all values so that new cells are zeroed.
        """
        cell._active_flag = False
        self.cell_count -= 1

        for field in cell._fields:
            setattr(cell, field, self._zero_vals[field])
        self._extra_attrs[cell.index] = {}

        if cell.index == self._last_index - 1:
            self._last_index -= 1
        else:
            self._free_indices += [cell.index]
