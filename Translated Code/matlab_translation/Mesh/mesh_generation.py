#import cupy as np
import numpy as np
from scipy.spatial import Delaunay


def points_in_polygon(points, polygon, include_boundary=True):
    """Vectorized point-in-polygon test for a closed polygon."""
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    polygon = np.asarray(polygon, dtype=float).reshape(-1, 2)
    if polygon.shape[0] < 3:
        raise ValueError("A polygon must contain at least three vertices.")

    x = points[:, 0]
    y = points[:, 1]
    inside = np.zeros(points.shape[0], dtype=bool)
    on_boundary = np.zeros(points.shape[0], dtype=bool)
    for i in range(polygon.shape[0]):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % polygon.shape[0]]
        dx = x2 - x1
        dy = y2 - y1
        cross = (x - x1) * dy - (y - y1) * dx
        scale = max(1.0, abs(dx), abs(dy))
        on_segment = (
            (np.abs(cross) <= 1e-12 * scale)
            & (x >= min(x1, x2) - 1e-12)
            & (x <= max(x1, x2) + 1e-12)
            & (y >= min(y1, y2) - 1e-12)
            & (y <= max(y1, y2) + 1e-12)
        )
        on_boundary |= on_segment
        crosses_ray = (y1 > y) != (y2 > y)
        intersection_x = x1 + dx * (y - y1) / (dy + np.finfo(float).tiny)
        inside ^= crosses_ray & (x < intersection_x)

    if include_boundary:
        return inside | on_boundary
    return inside & ~on_boundary


def _boundary_vertices_inside_triangles(triangles, points, vertices):
    """Identify Delaunay triangles that enclose a vertex of the hole."""
    a = points[triangles[:, 0]]
    b = points[triangles[:, 1]]
    c = points[triangles[:, 2]]
    ab = b - a
    ac = c - a
    denominator = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    enclosed = np.zeros(triangles.shape[0], dtype=bool)

    for point in vertices:
        candidates = (
            (point[0] > np.minimum(np.minimum(a[:, 0], b[:, 0]), c[:, 0]))
            & (point[0] < np.maximum(np.maximum(a[:, 0], b[:, 0]), c[:, 0]))
            & (point[1] > np.minimum(np.minimum(a[:, 1], b[:, 1]), c[:, 1]))
            & (point[1] < np.maximum(np.maximum(a[:, 1], b[:, 1]), c[:, 1]))
        )
        if not np.any(candidates):
            continue
        relative = point - a[candidates]
        u = (relative[:, 0] * ac[candidates, 1] - relative[:, 1] * ac[candidates, 0]) / denominator[candidates]
        v = (ab[candidates, 0] * relative[:, 1] - ab[candidates, 1] * relative[:, 0]) / denominator[candidates]
        enclosed[np.flatnonzero(candidates)] |= (
            (u > 1e-10) & (v > 1e-10) & (u + v < 1.0 - 1e-10)
        )
    return enclosed


def mesh_generation_C4_PEC(xy_cell=None, xy_inc=None, mesh_size=0.02):
    """Create a triangular mesh for a square periodic cell with a PEC hole.

    The boundary sampling and maximum spacing follow the MATLAB
    ``mesh_generation_C4_PEC`` routine.  SciPy's Delaunay triangulation is
    filtered so elements inside or crossing the polygonal inclusion are removed.
    
    
    """
    if mesh_size <= 0:
        raise ValueError("mesh_size must be positive.")
    if xy_cell is None:
        xy_cell = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    xy_cell = np.asarray(xy_cell, dtype=float).reshape(-1, 2)
    if xy_cell.shape[0] < 4:
        raise ValueError("The periodic cell must have at least four boundary points.")
    if xy_inc is None:
        # Default to a circular inclusion with radius 10% of the cell size, centered in the cell.
        center = (xy_cell.min(axis=0) + xy_cell.max(axis=0)) / 2.0
        radius = 0.1 * min(np.ptp(xy_cell[:, 0]), np.ptp(xy_cell[:, 1]))
        theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
        xy_inc = center + radius * np.column_stack((np.cos(theta), np.sin(theta)))
    xy_inc = np.asarray(xy_inc, dtype=float).reshape(-1, 2)
    if xy_inc.shape[0] < 3:
        raise ValueError("The inclusion must have at least three boundary points.")

    edge_lengths = np.linalg.norm(np.roll(xy_inc, -1, axis=0) - xy_inc, axis=1)
    subdivisions = np.maximum(1, np.ceil(edge_lengths / mesh_size).astype(int))
    inclusion_boundary = np.vstack([
        xy_inc[i] + np.arange(subdivisions[i])[:, None] / subdivisions[i]
        * (xy_inc[(i + 1) % xy_inc.shape[0]] - xy_inc[i])
        for i in range(xy_inc.shape[0])
    ])

    x_min, y_min = xy_cell.min(axis=0)
    x_max, y_max = xy_cell.max(axis=0)
    grid_spacing = mesh_size / np.sqrt(2.0)
    nx = max(2, int(np.ceil((x_max - x_min) / grid_spacing)) + 1)
    ny = max(2, int(np.ceil((y_max - y_min) / grid_spacing)) + 1)
    xs = np.linspace(x_min, x_max, nx)
    ys = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    grid_points = np.column_stack((X.ravel(), Y.ravel()))

    outer_boundary = np.vstack((
        np.column_stack((xs, np.full(xs.shape, y_min))),
        np.column_stack((xs, np.full(xs.shape, y_max))),
        np.column_stack((np.full(ys.shape, x_min), ys)),
        np.column_stack((np.full(ys.shape, x_max), ys)),
    ))
    grid_points = grid_points[~points_in_polygon(grid_points, xy_inc)]
    points = np.unique(
        np.vstack((grid_points, outer_boundary, xy_cell, inclusion_boundary)),
        axis=0,
    )
    triangles = Delaunay(points).simplices

    centroids = points[triangles].mean(axis=1)
    keep = ~points_in_polygon(centroids, xy_inc, include_boundary=False)

    edge_midpoints = np.vstack((
        (points[triangles[:, 0]] + points[triangles[:, 1]]) / 2.0,
        (points[triangles[:, 1]] + points[triangles[:, 2]]) / 2.0,
        (points[triangles[:, 2]] + points[triangles[:, 0]]) / 2.0,
    ))
    edge_midpoints_inside = points_in_polygon(edge_midpoints, xy_inc, include_boundary=False)
    keep &= ~(
        edge_midpoints_inside[:triangles.shape[0]]
        | edge_midpoints_inside[triangles.shape[0]:2 * triangles.shape[0]]
        | edge_midpoints_inside[2 * triangles.shape[0]:]
    )
    keep &= ~_boundary_vertices_inside_triangles(triangles, points, inclusion_boundary)
    triangles = triangles[keep]
    if triangles.size == 0:
        raise ValueError("Mesh generation produced no elements; check the cell and inclusion.")

    region_numbers = np.zeros(triangles.shape[0], dtype=int)
    return points, triangles.astype(int, copy=False), region_numbers
