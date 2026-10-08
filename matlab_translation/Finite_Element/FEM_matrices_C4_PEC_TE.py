import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def FEM_matrices_C4_PEC_TE(p,t,xy_box,xy_inc):
    
    # Inputs:
    #   p:      double, (2,n_nodes) array of mesh nodes
    #   t:      int, (n_elements, n_dof) triangular element matrix
    #   xy_box: the whole rectangular computational domain
    #   xy_inc: grid points for the inclusion
    
    # Outputs:
    #   M: double, (sparse), finite element mass matrix int_D u' * v' dx
    #   Mx: double, (sparse), finite element mass matrix int_D (u')_x * v dx
    #   My: double, (sparse), finite element mass matrix int_D (u')_y * v dx
    #   Me: double, (sparse), finite element stiffness matrix int_D u' * v' dx
    
    n_nodes = np.size(p,1)      # number of nodes
    n_elements = np.size(t,1)   # number of elements
    n_dof = np.size(t,2)        # number of degree of freedom per element
    
    x_min = xy_box[0]; x_max = xy_box[1]
    y_min = xy_box[2]; y_max = xy_box[3]
    
    # ==============
    # Label Nodes   
    # ==============
    tol = 1e-5                  # Tolerance, "should be smaller than mesh width"
    idx = np.zeros((n_nodes,1)) # Node IDs?
    n_equations = 0             # Number of equations
    
    # label the nodes initially, not including nodes on the right and top boundaries
    
    for i in range(1,n_nodes):
        IN, ON = inpolygon(p[i,1],p[i,2],xy_inc[:,1],xy_inc[:,2]) # Check if node is inside or on the inclusion boundary
        
