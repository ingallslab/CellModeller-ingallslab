# CellModeller TODOs

Multicellular agent-based modelling framework

Please move finished items to unreleased section of `CHANGELOG.md`

### Planned

 - [ ] Add `BoxContact` entity to mimic original plane contacts
 - [ ] Add tutorial examples that demonstrate use cases of new features
   - [ ] Moving planes?
   - [ ] Conjugation events?
   - [ ] Diffusion-based toxins?
 - [ ] Investigate sources of **numerical instability** in biophysics solver...
   - [ ] Explore possible preconditioners for contact matrix
 - [ ] Fix `CellArrays` slice implementation
 - [ ] Add/port useful features from original application
   - [ ] Saving to pickles/periodically writing simulation state
   - [ ] Right-click on a cell in GUI to display cell attributes
   - [ ] Integration with PDESolver and visualize in ParaView
   - [ ] and... lots of other features that are definitely possible to migrate but will take time
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

 - [x] Filamenting model demonstrating module flexibility
   - [ ] Cell events for creation of new segments in filamenting cell
   - [ ] Add OpenGL renderer for filamenting cells
 - [ ] Improve `README.md` to reflect recent updates to application
   - [x] Include instructions for setup and installation
   - [x] Link to tutorial examples for writing a basic simulation
   - [x] How to generate documentation with Doxygen
   - [ ] Brief overview of module structure and function
   - [ ] Actually write something that isn't just AI boilerplate...
   - [ ] Throw in some pictures/screenshots!

### Done

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
