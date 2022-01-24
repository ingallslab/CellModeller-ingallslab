CellModeller4
=============
Multicellular modelling framework, created by Tim Rudge, PJ Steiner, and Jim Haseloff, University of Cambridge

## Updates:
This fork of CellModeller can simulate growth in microfluidic traps. This is permitted by the following functionalities:

- Cell removal from the simulation
- PDE solver with user-defined boundary conditions

Cell removal is handled by deleting the cell id from sim.CellStates and excluding these cells from the biophysics handling. See [this tutorial](Examples/AaronYip/cell_removal).

Scalar fields are handled by [FEniCS](https://fenicsproject.org/) 2019.1.0, an open-source PDE solver. The CellModeller engine and FEniCS solver are fully coupled. Please download FEniCS from their website.

Cells and fields can be visualized simultaneously in [Paraview](https://www.paraview.org/) by converting .pickle files to .vtp format. Some parts of this are still in progress.

All the above updates are based on the work of @WilliamPJSmith.

/!\ Simulations interfaced with FEniCS cannot be run in the GUI.
