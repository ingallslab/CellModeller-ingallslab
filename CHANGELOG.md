# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [5.0.0a1] - 2025-04-25

Major rewrite and breaking changes to API. Restructure application to make it
more flexible and amenable to future development. Still in unstable alpha,
as more features are added it might eventually make it to next major version.

**NOTE:** Important refactoring changes that affect *users*
 - Changed cell attributes to use Python convention snake case i.e. `divide_flag`
 instead of `divideFlag`
 - `cell.growth_rate` is now measured in %/min not %/second. See `time_factor`
 documentation in `CellModeller.Modules.CMModule` for more details. Essentially,
 a growth rate of 0.05 corresponds to a growth rate of 1 in the previous version.
 - `cell.index` is not a unique identifier. It is guaranteed to be constant
 through the duration of any individual cell's life, but new indices are not
 assigned in a particular order. If a simulation requires unique ids assign in
 new_cell() method.

### Added

 - New guidelines to `contributing.md`
   - Doxygen doctrings
   - Black code formatter
   - Conventional commits format
   - Semantic versioning usage
   - Up-to-date TODOs
   - Maintaining the changelog
 - Support for keyboard controls in GUI
 - OpenGL retained mode classes and shaders
   - Allow renderers to define cell attributes for intended cell types
   - Base wrapper class that handles shader compiling and memory management
   - Default class for GUI background grid
   - Default class for rod-shaped cells
 - New classes to manage cell data efficiently
   - Define custom cell attributes using `numpy.dtype` declaration
   - Preallocate memory for known cell attributes
   - Manage centralized OpenCL shared memory between modules
 - Default biophysics implementation for rod-shaped cells
   - Performs the same role as in the previous version of CellModeller
   - Extensible interface makes it easier to implement different cell types

### Changed

 - From `setup.py` to `pyproject.toml`
   - pre-commit hook for automatic code formatting with black
   - Raise Python requirement to 3.10+ for improved `typing` support
 - DoxygenConfig to version 1.10.0
   - Include documentation from python docstrings
   - Hide undocumented members from output
 - Main GUI implementation
   - Use retained mode instead of immediate mode for OpenGL rendering
   - Improve performance by separating UI and sim in different threads
     - Simulator syncs with real time, processing speed permitting
     - Simulation speed independently controlled by `time_factor` variable
 - Underlying `CellState` implementation
   - Warn users about using undefined cell attributes
 - Main simulator implementation
   - Support for custom cell events and flexible sim loop updates
   - Use any number of default and user-defined modules simultaneously
   - Allows OpenGL renderers to read cell attributes directly

### Deprecated

 - (Temporarily) save to pickles features while rewriting application
 - Old biophysics classes
 - Old integration classes until reimplemented
 - Old signalling classes until reimplemented
 - Regulation/`ModuleRegulator.py`, replaced by `CMModule`

### Removed

 - `PyGLWidget` class
 - Immediate mode renderers in `CellModeller/GUI/Renderers.py`

### Fixed

 - Address underlying array fields with object.attribute syntax
   - Better syncing between OpenCL arrays and cell data in modules
 - Use SciPy sparse library to improve clarity of mathematical implementation

### BREAKING CHANGES

 - All old scripts and modules that rely on:
   - Specific attributes of `Simulator.py`
   - The old fixed-type 4-module system
   - Cells being stored as a list of Python objects (incl. pickles)
 - Everything else not tested or included as part of v5.0.0
