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


def _lowest_characteristic_values_numpy(matrices, kappa, n_omega):
    M, Mx, My, Me = matrices

    operator = M - 2.0j * (kappa[0] * Mx + kappa[1] * My)
    operator = operator + np.dot(kappa, kappa) * Me

    
    if n_omega >= operator.shape[0]:
        raise ValueError("n_omega must be smaller than the number of FEM equations.")

    # This is the generalized sparse eigensolve used by MATLAB's eigs(...,'smallestabs').
    if n_omega >= operator.shape[0] - 1:
        eigenvalues = eigh(
            operator.toarray(),
            Me.toarray(),
            eigvals_only=True,
            subset_by_index=(0, n_omega - 1),
        )
    else:
        
        eigenvalues = scipy_eigsh(
            operator,
            k=n_omega,
            M=Me,
            sigma=0.0,
            which="LM",
            return_eigenvectors=False,
        )
    eigenvalues = np.real_if_close(eigenvalues, tol=1000)
    if np.iscomplexobj(eigenvalues):
        raise ArithmeticError("The FEM eigenvalues are not real within numerical tolerance.")
    eigenvalues = np.real(eigenvalues)
    if np.any(eigenvalues < -1e-10):
        raise ArithmeticError("The FEM system produced negative eigenvalues.")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    return np.sort(eigenvalues)


def _lowest_characteristic_values_cupy(matrices, kappa, n_omega):
    """Compute characteristic values using CuPy and the GPU eigensolver."""
    M, Mx, My, Me = matrices
    M_gpu = cupy_csr_matrix(M)
    Mx_gpu = cupy_csr_matrix(Mx)
    My_gpu = cupy_csr_matrix(My)
    Me_gpu = cupy_csr_matrix(Me)
    kappa_gpu = cp.asarray(kappa)

    operator_gpu = (
        M_gpu
        - 2.0j * (kappa_gpu[0] * Mx_gpu + kappa_gpu[1] * My_gpu)
        + cp.dot(kappa_gpu, kappa_gpu) * Me_gpu
    )
    if n_omega >= operator_gpu.shape[0] - 1:
        raise ValueError(
            "The CuPy sparse eigensolver requires n_omega < matrix size - 1."
        )

    eigenvalues_gpu = cupy_eigsh(
        operator_gpu,
        k=n_omega,
        M=Me_gpu,
        sigma=0.0,
        which="LM",
        return_eigenvectors=False,
    )
    eigenvalues = cp.asnumpy(eigenvalues_gpu)
    eigenvalues = np.real_if_close(eigenvalues, tol=1000)
    if np.iscomplexobj(eigenvalues):
        raise ArithmeticError("The FEM eigenvalues are not real within numerical tolerance.")
    eigenvalues = np.real(eigenvalues)
    if np.any(eigenvalues < -1e-10):
        raise ArithmeticError("The FEM system produced negative eigenvalues.")
    return np.sort(np.maximum(eigenvalues, 0.0))


def _lowest_characteristic_values(matrices, kappa, n_omega):
    # Comment out one line and uncomment the other to switch backends.
    return _lowest_characteristic_values_numpy(matrices, kappa, n_omega)
    #return _lowest_characteristic_values_cupy(matrices, kappa, n_omega)

def _lowest_eigenfrequencies(matrices, kappa, n_omega):
    return np.sqrt(_lowest_characteristic_values(matrices, kappa, n_omega))

def compare_to_analytic_estimate(lat_c=1.0, radii=[0.25, 0.1, 0.01]):
    ''' Compute touching characteristic values at M for several PEC inclusion radii and compare to analyic estimate. '''
    from bandpath import calc_m_point_lowest_touching

    C_0 = 8* (np.pi**3)/(lat_c**4)
    # print("Normalized characteristic values omega * a / (2 pi)")
    for eps in radii:
        estimate = np.sqrt(2*np.pi**2+C_0*(eps**2))
        # *lat_c/(2*np.pi)
        touching = calc_m_point_lowest_touching(lat_c=lat_c, PEC_radius=eps)
        # *lat_c/(2*np.pi)
        avg = np.mean(touching)
        print(f"Radius = {eps:.2f}   Touch @ M: {touching[0]:.6f}, {touching[1]:.6f}   Averaged: {avg:.6f}   Error Gap: {touching[1] - touching[0]:.6f}")
        # print(f"                Averaged:   {avg:.6f}")
        print(f"                Analytic:  {estimate:.6f}")
        print(f"                Diffs:     {touching[0]-estimate:.6f}, {touching[1]-estimate:.6f}   Averaged Diff: {avg-estimate:.6f}\n")
