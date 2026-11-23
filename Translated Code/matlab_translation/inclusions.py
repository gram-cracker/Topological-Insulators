import sys
from pathlib import Path

import matplotlib
#matplotlib.use("Agg")
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh as scipy_eigsh

#Cupy
# import cupy as cp 
# from cupyx.scipy.sparse import csr_matrix as cupy_csr_matrix
# from cupyx.scipy.sparse.linalg import eigsh as cupy_eigsh


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))  

from Finite_Element.FEM_matrices_C4_PEC_TE import FEM_matrices_C4_PEC_TE
from Mesh.mesh_generation import mesh_generation_C4_PEC

import time





def geo_PEC_inclusions(lat_c=1.0, PEC_radius=0.1):
    """Return the polygonal circular PEC boundary used by the MATLAB model."""
    # radius = 0.1 * lat_c
    center = np.array([0.5 * lat_c, 0.5 * lat_c])
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    return center + PEC_radius * np.column_stack((np.cos(theta), np.sin(theta)))

def geo_PEC_inclusions_squircle(lat_c=1.0, PEC_radius=0.1, Square_constant=1):
    """Return the polygonal square PEC boundary used by the MATLAB model."""
    #radius = 0.1 * lat_c
    center = np.array([0.5 * lat_c, 0.5 * lat_c])
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    s = Square_constant
    ro = (PEC_radius * np.sqrt(2))/(s * np.abs(np.sin(2 * (theta + 0.0001))) ) * np.sqrt(1 - np.sqrt(1 - s**2 * (np.sin(2*(theta+ 0.0001))**2)))
    return center +   np.column_stack((ro * np.cos(theta), ro * np.sin(theta)))

def _assemble_cell(lat_c=1.0, mesh_size=0.02, PEC_radius=0.1):
    xy_box = np.array([0.0, lat_c, 0.0, lat_c])
    xy_cell = np.array([
        [0.0, 0.0],
        [lat_c, 0.0],
        [lat_c, lat_c],
        [0.0, lat_c],
    ])
    #xy_inc = geo_PEC_inclusions(lat_c=lat_c, PEC_radius=PEC_radius)
    xy_inc = geo_PEC_inclusions_squircle(lat_c=lat_c, PEC_radius=PEC_radius)
    p, t, _ = mesh_generation_C4_PEC(xy_cell, xy_inc, mesh_size=mesh_size)
    M, Mx, My, Me, _ = FEM_matrices_C4_PEC_TE(p, t, xy_box, xy_inc)
    return M, Mx, My, Me