# CellModeller TODOs

Multicellular agent-based modelling framework

Please move finished items to unreleased section of `CHANGELOG.md`

### Planned

 - [ ] Improve `README.md` to reflect recent updates to application
   - [ ] Include instructions for setup and installation
   - [ ] Link to tutorial examples for writing a basic simulation
   - [ ] How to generate documentation with Doxygen
   - [ ] Brief overview of module structure and function
     - [ ] Biophysics implementation
 - [ ] Add `BoxContact` entity to mimic original plane contacts
 - [ ] Add OpenGL renderer for filamenting cells
 - [ ] Add tutorial examples that demonstrate use cases of new features
   - [ ] Moving planes?
   - [ ] Come up with more ideas, ideally 2-3?
 - [ ] Investigate sources of **numerical instability** in biophysics solver...
 - [ ] Fix `CellArrays` slice implementation
 - [ ] Add/port useful features from original application
   - [ ] Saving to pickles/periodically writing simulation state
   - [ ] Right-click on a cell in GUI to display cell attributes
   - [ ] Integration with PDESolver and visualize in ParaView
 - [ ] Fix OpenCL kernels in contact finding
   - This is **not a trivial problem!!!**
   - Heavily platform-dependent optimizations required
   - See linked resources for examples
     - [OpenCL matrix-multiplication SGEMM tutorial (Cedric Nugteren, 2014)](https://cnugteren.github.io/tutorial/pages/page1.html)
     - [NVIDIA GPU Computing Webinars Best Practices for OpenCL Computing (2009)](https://developer.download.nvidia.com/CUDA/training/NVIDIA_GPU_Computing_Webinars_Best_Practises_For_OpenCL_Programming.pdf)
     - [Intel® FPGA SDK for OpenCL™ Standard Edition: Best Practices Guide (2018)](https://www.intel.com/content/www/us/en/docs/programmable/683176/18-1/introduction-to-standard-edition-best.html)
   - Even with newer technology, parallelization is still highly nuanced
   - Preferably prior experience in high-performance computing

### Working On

 - [ ] Add tasks to `TODO.md`
   - [x] Housekeeping guidelines
   - [x] Render pipeline overhaul
   - [x] New flexible module framework
   - [x] Filamenting support + module examples
   - [x] Starting points for full app migration
 - [ ] Add new biophysics classes using `ModuleProtocol`
   - [x] Default bacterial model
     - [ ] Document physics translated from original application
   - [ ] Filamenting model demonstrating module flexibility

### Done

 - [x] Add new guidelines to `contributing.md`
   - [x] Doxygen doctrings
   - [x] Black code formatter
   - [x] Conventional commits format
   - [x] Semantic versioning usage
   - [x] Up-to-date TODOs
   - [x] Maintaining the changelog
 - [x] Migrate from `setup.py` to `pyproject.toml`
   - [x] pre-commit hook for automatic code formatting with black
   - [x] Raise Python requirement to 3.10+ for improved `typing` support
 - [x] Update DoxygenConfig to version 1.10.0
   - [x] Include documentation from python docstrings
   - [x] Hide undocumented members from output
 - [x] Refactor `Scripts/CellModellerGUI.py` to use `if __name__ == "__main__":`
 - [x] Remove `PyGLWidget` class
 - [x] Change `CellModeller.GUI.PyGLCMViewer` implementation
   - [x] Use retained mode instead of immediate mode for OpenGL rendering
   - [x] Improve performance by separating UI and sim in different threads
     - [x] Simulator syncs with real time, processing speed permitting
     - [x] Simulation speed independently controlled by `time_factor` variable
   - [x] Support for keyboard controls
 - [x] Remove save to pickles features while rewriting application
 - [x] Change `CellModeller.CellState` implementation
   - [x] Address underlying `numpy.ndarray` fields with object.attribute syntax
   - [x] Warn users about using undefined cell attributes
 - [x] Remove immediate mode renderers in `CellModeller/GUI/Renderers.py`
 - [x] Add OpenGL retained mode classes and shaders
   - [x] Base wrapper class that handles shader compiling and memory management
     - [x] Default class for GUI background grid
     - [x] Default class for rod-shaped cells
   - [x] Allow renderers to define cell attributes for intended cell types
 - [x] Add `CellArrays` class to manage cell data efficiently
   - [x] Define custom cell attributes using `numpy.dtype` declaration
   - [x] Preallocate memory for known cell attributes
   - [x] Manage centralized OpenCL shared memory between modules
     - [x] Require a context manager to synchronize
   - [x] Conform to `collection.abc.Sequence` specification to access cells
   - [x] Provide `memoryview`s to contiguous blocks of cell data for renderers
 - [x] Add new `CellModeller.Modules.ModuleProtocols`
   - [x] Support for custom cell events and flexible sim loop updates
 - [x] Change `CellModeller.Simulator` implementation
   - [x] Integrate with `ModuleProtocol` and `UserModule` classes
   - [x] Check for custom cell events in update loop
   - [x] Forward access layer for renderers to cell data

---

### Guideline for changelog verbs: [^1]

Verb    | Use for...
--------|-----------
Add     | New features, files, or functionality
Remove  | Deprecated or obsolete features/code/files
Change  | Behavior changes or refactors that alter how something works
Update  | Version bumps, metadata changes, or refreshed content (e.g. docs, deps)
Fix     | Bug fixes
Migrate | Moving from one system to another (e.g. `setup.py` -> `pyproject.toml`)
Improve | Enhancements to performance, UX, maintainability, etc.
Refactor | Internal restructuring without changing external behavior

[^1]: This table was made with GenAI: ChatGPT using the GPT-4-turbo model, April
2025 release.
