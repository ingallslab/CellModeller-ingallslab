"""
Use with: microfluidics_1open.py

FEniCS tutorial demo program: Poisson equation with Dirichlet and Neumann conditions.
  laplace(u) = f  
            f = 0 (no generation or consumption) 
            grad(u) = 0 at walls
            u = 1 at open boundary
            
The solution to the diffusion equation with these BCs will be a linear concentration gradient.
"""

'''
This is a demo simulation for growth in a microfluidic device. The goal is to demonstrate that the fields from FEniCS can be used to adjust CellModeller growth rates. 

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
        extend_elements = 2 #Extend domain by two cells in y
        extend_length_y = solverParams['L_y']/solverParams['N_y']*extend_elements #Extend domain 
        
        # extract fixed params from params dictionary
        self.N_x = int(solverParams['N_x'])
        self.N_y = int(solverParams['N_y']) + extend_elements
        self.L_x = solverParams['L_x'] 
        self.L_y = solverParams['L_y'] + extend_length_y
        self.D = solverParams['D']
        self.u0 = solverParams['u0']
        self.mu_max = solverParams['mu_max']
        self.K = solverParams['K']
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
        global N_x
        global N_y
        
        L_x = self.L_x
        L_y = self.L_y
        N_x = self.N_x
        N_y = self.N_y
        
        # we'll need to extend this to something more general later
        # also: clear up the file I/O stuff, it's messy!
        self.SetMesh()
        self.V = FunctionSpace(self.mesh, "CG", 1)
        
        # choose boundary conditions
        self.set_bcs()
        
        # Use cell centres to evaluate volume occupancy of mesh
        self.GetVolFracs2D(centers, areas)
        #self.GetVolFracsCrossMesh(centers, areas)
               
        # call a solver and save the solution
        set_log_active(False)
        #set_log_level(50) #INFO=20, ERROR=40, #CRITICAL=50
        self.NewtonIterator()
        if stepNum % self.pickleSteps == 0:
            self.WriteFieldToFile(dir+filename+'_field.pvd', self.solution)
            # Export a meshfile showing element occupancies (for testing purposes)
            G = self.VolumeFractionOnElements()
            g = Function(self.V, name = "Volume fraction")
            g.interpolate(G)
            self.WriteFieldToFile(dir+filename+'_VolFracsCheck'+'.pvd', g)
            
        # interpolate solution to cell centres
        u_local = self.InterpolateToCenters2D(centers)
        
        return u_local
        
        
    def SetMesh(self):
        '''
        SetMesh defines a rectangular mesh starting from (0,0). Length and number of computational cells are defined in the CellModeller simulation. 
        '''    
        mesh_type = self.mesh_type
        self.mesh = RectangleMesh(Point(0.0,0.0), Point(self.L_x,self.L_y), self.N_x, self.N_y, mesh_type)
        
    
    def set_bcs(self):
        '''
        Initialise boundary conditions on the mesh.
        /!\ Assumes that the function V is up-to-date
        '''
        tbc = TopDirichletBoundary()
        
        bc0 = DirichletBC(self.V, Constant(self.u0), tbc)
        
        self.boundaryCondition = bc0
                     
        
    def NewtonIterator(self):
        """
        A Newton iterator for solving non-linear problems.
        /!\ Assumes that function space (V), boundaryCondition, vol_fracs are up-to-date.
        """
        # Define variational problem       
        u = Function(self.V, name = "Species A")
        v = TestFunction(self.V)
        
        #This formulation depends on the governing equation
        a = Constant(self.D)*dot(grad(u), grad(v))*dx
        L = self.SourceTerm(u)*v*dx - self.zeroGradient*v*ds
        F = a - L

        # Set some sort of logger if desired. Format below seems out of date
        #set_log_level(PROGRESS) # very detailed info for testing
        #set_log_level(WARNING) # near-silent info for simulations
        #set_log_active(False) # suppress solver text
        
        # Call built-in Newton solver
        solve(F == 0, u, self.boundaryCondition, solver_parameters={"newton_solver":{"relative_tolerance":self.rel_tol,
                                    "convergence_criterion":"incremental",
                                    "linear_solver":"lu",
                                    "maximum_iterations":500}})        
        self.solution = u    
        
        
    def SourceTerm(self, u):
        ''' 
        Source term (e.g. Monod for nutrients, Michaelis-Menten for antibiotics)
        '''
        mu_max = Constant(self.mu_max)
        K = Constant(self.K)
        
        return -1*mu_max*u*VolumeFraction()/(u + K)
        
    
    def GetVolFracsCrossMesh(self, centers, vols):
        """
        Create a global list of the cell volume fractions in mesh elements.
        Assumes that self.mesh is up to date.
        'Volumes' are equivalent to areas in 2D.
        /!\ Exports the array vol_fracs as a global array, for use by VolumeFraction.
        /!\ Takes 
        """
        global vol_fracs, L_y, N_y

        # assign elements of cells
        elements = self.AssignElementsToDataCrossMesh(centers)
        
        # need to define volume fraction for every element in the mesh
        # (not just the canonical elements for counting)
        num_elements = self.mesh.num_cells()

        # sum cell volumes over each element
        a = (self.L_x*L_y) / (4.0*float(self.N_x)*float(N_y))
        vol_fracs = numpy.bincount(elements,vols,num_elements) / a
        
        
    def AssignElementsToDataCrossMesh(self, centers):
        """
        Sort cell centres into their respective mesh elements.
        """
        N = centers.shape[0]
        elements = numpy.zeros((N), numpy.int32)
        for i in range(0,N):
            point = centers[i]
            elements[i] = self.GetElementIndexCrossMesh(point)
            
        return elements
        
        
    def GetElementIndexCrossMesh(self, point):
        """
        Get tetrahedron and cube indices and calculate global element index.
        """
        [s, sqOrigin] = self.GetSquareIndex(point)
        t = self.GetTriangleIndexCrossMesh(point,sqOrigin)
        
        return t + 4*s
        
        
    def GetSquareIndex(self, Point):
        """
        Given mesh dimensions, assign which square a point is in.
        /!\ Assumes L_y and N_y are supplied as global variables.
        """
        global L_y, N_y
        
        p_x = Point[0]; p_y = Point[1]
        p = int(numpy.floor(p_x*self.N_x / float(self.L_x)))          # index along x
        q = int(numpy.floor(p_y*N_y / float(L_y)))           # index along y

        s = p + q*self.N_x                                             # global index of this square
        sqOrigin = [p*self.L_x / float(self.N_x),\
                    q*L_y / float(N_y)]                        # coordinates of this square's origin                
        
        return int(s), sqOrigin
        
        
    def GetTriangleIndexCrossMesh(self, Point, Origin):
        """
        Given mesh square, assign which triangle a point is in.
        /!\ Assumes that we're using a crossed rectilinear mesh.
        /!\ Assumes L_y and N_y are supplied as global variables.
        """
        global L_y, N_y
        
        p_x = Point[0]; p_y = Point[1]
        a_x = Origin[0]; a_y = Origin[1]
        dx = p_x - a_x
        dy = p_y - a_y
        hx = self.L_x / self.N_x
        hy = L_y / N_y
        gr = hy / hx; 
        
        return 1*(dy-gr*dx>0) + 2*(dy+gr*dx-hy>0);
        
        
    def InterpolateToCenters2D(self, centers):
        """
        Interpolate a solution object u onto a list of cell coordinates
        """
        u = self.solution
        data_t = tuple(map(tuple, centers)) # Convert to tuple format
        u_local = numpy.zeros(len(data_t),dtype=numpy.float64)  # preallocate solution array
        for i in range(len(data_t)):  # loop over all cells
            u_local[i] = numpy.maximum(u(data_t[i][0:2]), 0.0) # extrapolate solution value at cell centre; 
        return u_local   
        
        
    def GetVolFracs2D(self, centers, vols):
        """
        Create a global list of the cell volume fractions in mesh elements.
        Assumes that self.mesh and self.h are up-to-date.
        'Volumes' are equivalent to areas in 2D.
        NOTE: this is pretty slow for >1000 cells, so I've written GetVolFracsCrossMechs as a faster alternative.
        /!\ Exports the array vol_fracs as a global array, for use by VolumeFraction.
        """
        global vol_fracs

        # assign elements of cells
        elements = self.AssignElementsToData2D(centers)
    
        # compute element volumes
        num_elements = self.mesh.num_cells()
        element_vols = numpy.zeros((num_elements,)) 
        for i in range(0,num_elements):
            element_vols[i] = Cell(self.mesh,i).volume()
        
        # Here, we need to:
        #  >  define volume fraction for every element in the mesh
        #  >  normalise volume totals by mesh element volumes
        binned_vol_fracs = numpy.bincount(elements,vols,num_elements)
        vol_fracs = numpy.divide(binned_vol_fracs, element_vols)


    def AssignElementsToData2D(self, centers):
        """
        Modified cell assignment function using the mesh cell 'collides' method.
        Unlike ASSIGNELEMENTSTODATA2D, this will work on unstructured [2D] meshes.
        """
        N = centers.shape[0]
        elements = numpy.zeros((N), numpy.int32)
    
        for i in range(0,N):
            for j in range(0,self.mesh.num_cells()):
                if Cell(self.mesh, j).collides(Point(float(centers[i][0]),\
                                                     float(centers[i][1]),\
                                                     0)):
                    elements[i] = j
                    break
        return elements
         

    def WriteFieldToFile(self, filename, u):
        '''
        Export the PDE solution as a pvd mesh.
        '''
        print('Writing fields...') 
        File(filename) << u
        print('Done.') 
        
    def VolumeFractionOnElements(self):
        """ 
        Returns the bacteria volume fraction for each element
        """
        
        return VolumeFraction()        
        
'''
Supporting classes
'''

class VolumeFraction(UserExpression):
    """
    Supporting class defining element-wise volume fraction function, for nutrient PDEs.
    
    TODO: define shape_value as scalar to get rid of warning message.
    """

    def eval_cell(self, value, x, ufc_cell):
        """
        Evaluate the cell volume fraction for this mesh element.
        /!\ Assumes vol_fracs is being supplied as a global variable.
        """
        global vol_fracs
        value[0] = min(vol_fracs[ufc_cell.index], 1.0) #Vol frac cannot exceed 1
        
        
    def value_shape(self):
        return ()
        
class TopDirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        global L_y
        return bool(near(x[1], L_y) and on_boundary)
