#import cupy as np
import numpy as np
from scipy import sparse

from Mesh.mesh_generation import points_in_polygon


def _triangle_gradients(xy):
    """Return the constant gradients of the 3 linear basis functions for one triangle."""
    x1, y1 = xy[0]
    x2, y2 = xy[1]
    x3, y3 = xy[2]
    area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if abs(area2) < 1e-14:
        raise ValueError('Degenerate triangle in FEM assembly.')

    b = np.array([y2 - y3, y3 - y1, y1 - y2], dtype=float)
    c = np.array([x3 - x2, x1 - x3, x2 - x1], dtype=float)
    grad = np.vstack((b, c)) / area2
    return grad, 0.5 * abs(area2)


def FEM_matrices_C4_PEC_TE(p, t, xy_box, xy_inc):
    """Assemble the FEM matrices for the square-lattice PEC periodic cell.

    The implementation follows the MATLAB logic but uses vectorized local-element
    assembly with NumPy/SciPy sparse matrices.
    """
    p = np.asarray(p, dtype=float)
    t = np.asarray(t, dtype=int)
    n_nodes = p.shape[0]
    n_elements = t.shape[0]
    n_dof = t.shape[1]

    x_min, x_max, y_min, y_max = xy_box
    tol = 1e-8

    idx = np.zeros(n_nodes, dtype=int)
    n_equations = 0

    is_inclusion = points_in_polygon(p, xy_inc)
    for i in range(n_nodes):
        x, y = p[i]
        if abs(x - x_max) < tol and abs(y - y_max) < tol:
            idx[i] = -1
        elif abs(x - x_max) < tol:
            idx[i] = -2
        elif abs(y - y_max) < tol:
            idx[i] = -3
        elif is_inclusion[i]:
            idx[i] = 0
        else:
            n_equations += 1
            idx[i] = n_equations

    corner = np.argmin(np.abs(p[:, 0] - x_min) + np.abs(p[:, 1] - y_min))
    left = np.where(np.abs(p[:, 0] - x_min) < tol)[0]
    bottom = np.where(np.abs(p[:, 1] - y_min) < tol)[0]

    if left.size == 0 or bottom.size == 0:
        raise RuntimeError('Periodic boundary nodes were not found for the mesh.')

    idx_left = idx[left]
    idx_bottom = idx[bottom]

    for i in range(n_nodes):
        if idx[i] == -1:
            idx[i] = idx[corner]
        elif idx[i] == -2:
            matches = np.where(np.abs(p[i, 1] - p[left, 1]) < tol)[0]
            if matches.size:
                idx[i] = idx_left[matches[0]]
        elif idx[i] == -3:
            matches = np.where(np.abs(p[i, 0] - p[bottom, 0]) < tol)[0]
            if matches.size:
                idx[i] = idx_bottom[matches[0]]

    if np.any(idx < 0):
        raise RuntimeError('Periodic boundary mapping was not successful; adjust the mesh size.')

    rows = []
    cols = []
    mats = []
    matx = []
    maty = []
    mate = []

    for ne in range(n_elements):
        nodes_loc = t[ne]
        xy_loc = p[nodes_loc]
        grad, area = _triangle_gradients(xy_loc)

        mass_loc = np.array([
            [2.0, 1.0, 1.0],
            [1.0, 2.0, 1.0],
            [1.0, 1.0, 2.0],
        ], dtype=float) * area / 12.0

        mx_loc = np.outer(np.ones(3), grad[0]) * area / 3.0
        my_loc = np.outer(np.ones(3), grad[1]) * area / 3.0
        stiffness_loc = area * (grad.T @ grad)

        for ni in range(n_dof):
            g_i = idx[nodes_loc[ni]]
            if g_i <= 0:
                continue
            for nj in range(n_dof):
                g_j = idx[nodes_loc[nj]]
                if g_j <= 0:
                    continue
                rows.append(g_i - 1)
                cols.append(g_j - 1)
                mats.append(stiffness_loc[ni, nj])
                matx.append(mx_loc[ni, nj])
                maty.append(my_loc[ni, nj])
                mate.append(mass_loc[ni, nj])

    M = sparse.coo_matrix((mats, (rows, cols)), shape=(n_equations, n_equations)).tocsr()
    Mx = sparse.coo_matrix((matx, (rows, cols)), shape=(n_equations, n_equations)).tocsr()
    My = sparse.coo_matrix((maty, (rows, cols)), shape=(n_equations, n_equations)).tocsr()
    Me = sparse.coo_matrix((mate, (rows, cols)), shape=(n_equations, n_equations)).tocsr()
    return M, Mx, My, Me, idx
