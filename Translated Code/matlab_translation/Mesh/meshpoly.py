import numpy as np
from scipy.spatial import Delaunay

from Mesh.mesh_generation import _boundary_vertices_inside_triangles, points_in_polygon


def meshpoly(node, edge, qtree, p, options=None):
    """Simple triangular mesh generator for polygonal regions.

    This lightweight replacement keeps the MATLAB/NumPy terminology while using the
    robust triangulation from SciPy. It is enough for the square-lattice FEM band
    structure translation used in the topological-insulators project.
    """
    if options is None:
        options = {'output': False}

    node = np.asarray(node, dtype=float)
    if node.ndim == 1:
        node = node.reshape(-1, 2)

    if node.shape[0] < 3:
        raise ValueError('At least three vertices are required to build a polygon mesh.')

    # Use a uniform grid of candidate points inside the polygon and triangulate.
    xs = np.linspace(node[:, 0].min(), node[:, 0].max(), 25)
    ys = np.linspace(node[:, 1].min(), node[:, 1].max(), 25)
    X, Y = np.meshgrid(xs, ys, indexing='xy')
    pts = np.column_stack((X.ravel(), Y.ravel()))

    points = np.unique(
        np.vstack((node, pts[points_in_polygon(pts, node, include_boundary=False)])),
        axis=0,
    )
    tri = Delaunay(points).simplices
    centroids = points[tri].mean(axis=1)
    keep = points_in_polygon(centroids, node, include_boundary=False)
    edge_midpoints = np.vstack((
        (points[tri[:, 0]] + points[tri[:, 1]]) / 2.0,
        (points[tri[:, 1]] + points[tri[:, 2]]) / 2.0,
        (points[tri[:, 2]] + points[tri[:, 0]]) / 2.0,
    ))
    n_triangles = tri.shape[0]
    midpoint_outside = ~points_in_polygon(edge_midpoints, node)
    keep &= ~(
        midpoint_outside[:n_triangles]
        | midpoint_outside[n_triangles:2 * n_triangles]
        | midpoint_outside[2 * n_triangles:]
    )
    keep &= ~_boundary_vertices_inside_triangles(tri, points, node)
    tri = tri[keep]
    if tri.shape[0] == 0:
        raise ValueError("Triangulation produced no elements inside the polygon.")

    p = points
    t = tri.astype(int, copy=False)
    return p, t
