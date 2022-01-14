"""
Use with: microfluidics_1open.py; eventually each simulation will have its own PDE solver file, but first need to dynamically import in Simulator.py

FEniCS tutorial demo program: Poisson equation with Dirichlet and Neumann conditions.
  laplace(u) = f  
            f = 1 
            grad(u) = at walls
            u = 0 at open boundary
"""

'''
This is a demo simulation for growth in a microfluidic device. We will solve a reaction-diffusion scenario, where a antibiotic field degrades at a constant rate. The cells will grow independently of the antibiotic field.

|--------------open-----------------|
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|                                   |
|___________________________________|

Domain size is 20 X 20 um.
'''

from dolfin import *
import matplotlib.pyplot as plt
import numpy
import math

class DolfinSolver:
    
    def __init__(self, solverParams):
        """ 
        Initialise the dolfin solver using a dictionary of params.
        """
        # extract fixed params from params dictionary
        self.N_x = int(solverParams['N_x'])
        self.N_y = int(solverParams['N_y'])
        self.L_x = solverParams['L_x']
        self.L_y = solverParams['L_y']
        self.mesh_type = solverParams['mesh_type']
        self.rel_tol = solverParams['rel_tol']
        self.pickleSteps = solverParams['pickleSteps']

        self.zeroGradient = Constant(0.0) #use this for zeroGradient BCs
        
        # some attributes that we'll update on the fly: set them to None for now
        self.boundaryCondition = None
        self.mesh = None
        self.V = None
        self.solution = None
        
    def SolvePDE(self, centers, areas, filename, dir, stepNum=0):
        """
        High-level function to be called during the function.
        """
        global L_y
        global L_x
        
        L_x = self.L_x
        L_y = self.L_y
        
        # we'll need to extend this to something more general later
        # also: clear up the file I/O stuff, it's messy!
        self.SetMesh()
        self.V = FunctionSpace(self.mesh, "CG", 1)
        
        # choose boundary conditions
        self.set_bcs()
        
        # Use cell centres to evaluate volume occupancy of mesh
        # self.GetVolFracsCrossMesh(centers, areas)
        
        # call a solver and save the solution
        self.NewtonIterator()
        if stepNum % self.pickleSteps == 0:
            self.WriteFieldToFile(dir+filename+'.pvd', self.solution)
        
        # interpolate solution to cell centres
        u_local = self.InterpolateToCenters2D(centers)
        
        return u_local
        
    def SetMesh(self):
        '''
        SetMesh defines a rectangular mesh starting from (0,0). Length and number of computational cells are defined in the CellModeller simulation. 
        '''    
        mesh_type = self.mesh_type
        self.mesh = RectangleMesh(Point(0.0,0.0), Point(self.L_x,self.L_y), self.N_x, self.N_y)
    
    def set_bcs(self):
        '''
        Initialise boundary conditions on the mesh.
        /!\ Assumes that the function V is up-to-date
        '''
        u_outlet = Constant(0.0) #set outlet BC to a constant
        dbc = TopDirichletBoundary()
        self.boundaryCondition = DirichletBC(self.V, u_outlet, dbc)
        
    def SourceTerm(self, u):
        ''' 
        Source term (e.g. Monod for nutrients, Michaelis-Menten for antibiotics)
        '''
        source = Constant(1.0) #Currently set to a constant for testing
        
        return source
        
    def NewtonIterator(self):
        """
        A Newton iterator for solving non-linear problems.
        /!\ Assumes that function space (V), boundaryCondition, vol_fracs are up-to-date.
        """
        # Define variational problem       
        u = Function(self.V, name = "Species A")
        v = TestFunction(self.V)
        
        #This formulation depends on the governing equation
        a = dot(grad(u), grad(v))*dx
        L = self.SourceTerm(u)*v*dx - self.zeroGradient*v*ds
        F = a - L

        # Set some sort of logger if desired. Format below seems out of date
        #set_log_level(PROGRESS) # very detailed info for testing
        #set_log_level(WARNING) # near-silent info for simulations
        #set_log_active(False) # suppress solver text
        
        # Call built-in Newton solver
        solve(F == 0, u, self.boundaryCondition, solver_parameters={"newton_solver":{"relative_tolerance":self.rel_tol}})
        
        self.solution = u    
        
    def InterpolateToCenters2D(self, centers):
        """
        Interpolate a solution object u onto a list of cell coordinates
        """
        u = self.solution
        data_t = tuple(map(tuple, centers)) # Convert to tuple format
        u_local = numpy.zeros(len(data_t),dtype=numpy.float64)  # preallocate solution array
        for i in range(0,len(data_t)):  # loop over all cells
            u_local[i] = numpy.maximum(u(data_t[i][0:2]), 0.0) # extrapolate solution value at cell centre; 
        return u_local    

    def WriteFieldToFile(self, filename, u):
        '''
        Export the PDE solution as a pvd mesh.
        '''
        print('Writing fields...') 
        File(filename) << u
        print('Done.') 
        
class TopDirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        global L_y
        return bool(near(x[1], L_y) and on_boundary)
