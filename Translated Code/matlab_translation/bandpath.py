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
from solver import _lowest_characteristic_values, _lowest_eigenfrequencies
from inclusions import _assemble_cell





def _band_path(lat_c, n_kappa, path_type="TXMT"):
    if path_type == "TXMT":
        """Return the T to X to M to T Bloch path in MATLAB coordinates."""
        points = [
            (np.linspace(0.0, np.pi, n_kappa),np.zeros(n_kappa)),
            (np.pi * np.ones(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
            (np.linspace(np.pi, 0.0, n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
        ]
    elif path_type == "MXTX'M":
        """Return the M to X to T to X' to M Bloch path in MATLAB coordinates."""
        points = [
            (np.pi * np.ones(n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
            (np.linspace(np.pi, 0.0, n_kappa), np.zeros(n_kappa)),
            (np.zeros(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
            (np.linspace(0.0, np.pi, n_kappa), np.pi * np.ones(n_kappa)),
        ]
    elif path_type =="TX'M'T":
        """Return the T to X' to M' to T Bloch path in MATLAB coordinates."""
        points = [
            (np.zeros(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
            (np.linspace(0.0, -np.pi, n_kappa), np.pi * np.ones(n_kappa)),
            (np.linspace(-np.pi, 0.0, n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
        ]
    else:
        raise ValueError("path_type unsupported.")
    return np.vstack([
        np.column_stack((kx, ky)) / lat_c
        for kx, ky in points
    ])


def solve_band_structure(lat_c=1.0, n_kappa=100, n_omega=7, mesh_size=0.02, PEC_radius=0.1, path_type="TXMT"):
    """Compute bands along the MATLAB path_type path."""
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=PEC_radius)
    return _solve_band_structure(matrices, lat_c, n_kappa, n_omega, path_type=path_type)


def _solve_band_structure(matrices, lat_c, n_kappa, n_omega, path_type="TXMT"):
    kappa_vec = _band_path(lat_c, n_kappa, path_type)
    omega = np.empty((n_omega, kappa_vec.shape[0]), dtype=float)
    for n, kappa in enumerate(kappa_vec):
        omega[:, n] = (
            # _lowest_characteristic_values(matrices, kappa, n_omega)
            _lowest_eigenfrequencies(matrices, kappa, n_omega) * lat_c / (2.0 * np.pi)
        )
    return omega, kappa_vec

def solve_frequency_surfaces(
    lat_c=1.0,
    grid_size=21,
    n_omega=5,
    mesh_size=0.02,
    radius_fraction=0.1,
    kappa_center=(0.0, 0.0),
    PEC_radius = 0.1,
):
    """Compute the lowest bands on the rectangular Bloch-vector grid in MATLAB."""
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=PEC_radius)
    return _solve_frequency_surfaces(matrices, lat_c, grid_size, n_omega, radius_fraction, kappa_center)


def _solve_frequency_surfaces(matrices, lat_c, grid_size, n_omega, radius_fraction, kappa_center):
    kx = kappa_center[0] + (radius_fraction * 2 * np.pi / lat_c) * np.linspace(-1.0, 1.0, grid_size)
    ky = kappa_center[1] + (radius_fraction * 2 * np.pi / lat_c) * np.linspace(-1.0, 1.0, grid_size)
    kappa_x, kappa_y = np.meshgrid(kx, ky)
    omega = np.empty((grid_size, grid_size, n_omega), dtype=float)
    for row in range(grid_size):
        for column in range(grid_size):
            omega[row, column, :] = _lowest_eigenfrequencies(
                matrices,
                np.array([kappa_x[row, column], kappa_y[row, column]]),
                n_omega,
            ) * lat_c / (2.0 * np.pi)
    return omega, kappa_x, kappa_y


# Depricated
# def solve_m_point_bands(
#     lat_c=1.0,
#     grid_size=21,
#     n_bands=2,
#     mesh_size=0.02,
#     radius_fraction=0.1,
#     estimate_samples=11,
#     PEC_radius=0.1,
# ):
#     """Solve the two lowest characteristic-value bands in a neighborhood of M."""
#     matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=PEC_radius)
#     return _solve_m_point_bands(
#         matrices, lat_c, grid_size, n_bands, radius_fraction, estimate_samples
#     )


# maybe depricated?
def _solve_m_point_bands(
    matrices, lat_c, grid_size, n_bands, radius_fraction, estimate_samples
):
    values, kappa_x, kappa_y, estimate = _solve_lattice_point_bands(
        matrices,
        np.array([np.pi / lat_c, np.pi / lat_c]),
        lat_c,
        grid_size,
        n_bands,
        radius_fraction,
        estimate_samples,
    )
    estimate["m_point"] = estimate["kappa_center"]
    return values, kappa_x, kappa_y, estimate


def solve_lattice_point_bands(
    kappa_center,
    lat_c=1.0,
    grid_size=21,
    n_bands=2,
    mesh_size=0.02,
    radius_fraction=0.1,
    estimate_samples=11,
    PEC_radius=0.1,
):
    """Compute nearby bands around a Bloch point, e.g. Gamma (0,0), X (pi/a,0), M (pi/a,pi/a)."""
    if lat_c <= 0:
        raise ValueError("lat_c must be positive.")
    center = np.asarray(kappa_center, dtype=float)
    if center.shape != (2,) or not np.isfinite(center).all():
        raise ValueError("kappa_center must contain two finite Bloch-vector coordinates.")
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=PEC_radius)
    return _solve_lattice_point_bands(
        matrices,
        center,
        lat_c,
        grid_size,
        n_bands,
        radius_fraction,
        estimate_samples,
    )


def _solve_lattice_point_bands(
    matrices,
    kappa_center,
    lat_c,
    grid_size,
    n_bands,
    radius_fraction,
    estimate_samples,
):
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer of at least 3.")
    if n_bands < 2:
        raise ValueError("n_bands must be at least 2 to estimate the band gap.")
    if lat_c <= 0:
        raise ValueError("lat_c must be positive.")
    if not 0.0 < radius_fraction:
        raise ValueError("radius_fraction must be positive.")
    if estimate_samples < 3:
        raise ValueError("estimate_samples must be at least 3.")

    kappa_center = np.asarray(kappa_center, dtype=float)
    if kappa_center.shape != (2,) or not np.isfinite(kappa_center).all():
        raise ValueError("kappa_center must contain two finite Bloch-vector coordinates.")
    radius = radius_fraction * np.pi / lat_c
    offsets = np.linspace(-radius, radius, grid_size)
    kappa_x, kappa_y = np.meshgrid(
        kappa_center[0] + offsets, kappa_center[1] + offsets
    )
    characteristic_values = np.empty((grid_size, grid_size, n_bands), dtype=float)
    for row in range(grid_size):
        for column in range(grid_size):
            characteristic_values[row, column, :] = _lowest_characteristic_values(
                matrices,
                np.array([kappa_x[row, column], kappa_y[row, column]]),
                n_bands,
            ) * lat_c / (2.0 * np.pi)

    #values_at_center = _lowest_characteristic_values(
       # matrices, kappa_center, n_bands
    #)
    values_at_center = (
    _lowest_characteristic_values(matrices, kappa_center, n_bands)
    * lat_c
    / (2.0 * np.pi)
)
    directions = {
        "kx": np.array([1.0, 0.0]),
        "ky": np.array([0.0, 1.0]),
        "diagonal": np.array([1.0, 1.0]) / np.sqrt(2.0),
        "anti-diagonal": np.array([1.0, -1.0]) / np.sqrt(2.0),
    }
    radii = np.linspace(0.0, radius, estimate_samples)
    fit_count = max(2, (estimate_samples - 1) // 2)
    directional_slopes = {}
    for name, direction in directions.items():
        gaps = np.empty(estimate_samples, dtype=float)
        for index, distance in enumerate(radii):
            values = _lowest_characteristic_values(
                matrices, kappa_center + distance * direction, n_bands
            )
            gaps[index] = values[1] - values[0]
        radial_distance = radii[1:fit_count + 1]
        gap_increase = (
            gaps[1:fit_count + 1]
            - values_at_center[1]
            + values_at_center[0]
        )
        directional_slopes[name] = float(
            np.dot(radial_distance, gap_increase)
            / np.dot(radial_distance, radial_distance)
        )

    estimate = {
        "kappa_center": kappa_center,
        "characteristic_values": values_at_center,
        "direct_gap": [float(values_at_center[i+1] - values_at_center[i]) for i in range(n_bands-1)],
        "mean_gap_slope": float(np.mean(list(directional_slopes.values()))),
        "directional_gap_slopes": directional_slopes,
        "fit_radius": float(radii[fit_count]),
    }
    return characteristic_values, kappa_x, kappa_y, estimate


def plot_band_structure(omega, path_type=""):
    fig, ax = plt.subplots()
    for band in range(omega.shape[0]):
        ax.plot(omega[band, :], linewidth=1.8, label=f"band {band + 1}")
    ax.set_xlabel("Bloch point index")
    ax.set_ylabel(r"Normalized Frequency ($\omega\frac{a}{2\pi}$)")
    
    if path_type == "TXMT":
        ax.set_xticks([0, omega.shape[1] // 3, 2 * omega.shape[1] // 3, omega.shape[1] - 1])
        ax.set_xticklabels([r"$\Gamma$", "X", "M", r"$\Gamma$"])
        title = r"$\Gamma$-X-M-$\Gamma$ band structure"
    elif path_type == "MXTxM":
        ax.set_xticks([0, omega.shape[1] // 4, 2 * omega.shape[1] // 4, 3 * omega.shape[1] // 4, omega.shape[1] - 1])
        ax.set_xticklabels(["M", "X", r"$\Gamma$", "X'", "M"])
        title = r"M-X-$\Gamma$-X'-M band structure"
    elif path_type == "TX'M'T":
        ax.set_xticks([0, omega.shape[1] // 3, 2 * omega.shape[1] // 3, omega.shape[1] - 1])
        ax.set_xticklabels([r"$\Gamma$", "X'", "M'", r"$\Gamma$"])
        title = r"$\Gamma$-X'-M'-$\Gamma$ band structure"
    else:
        title = "Band structure"
        
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig, ax


def plot_band_comparison(band_structures, bands=None):
    """Overlay band results from multiple runs on shared axes.

    ``band_structures`` maps case labels to arrays shaped
    ``(number_of_bands, number_of_path_points)``. All cases must use the same
    Bloch-path sampling so a given x-coordinate represents the same path point.
    ``bands`` optionally selects 1-based band numbers, e.g. ``(1, 3)``.
    """
    if not band_structures:
        raise ValueError("Provide at least one labeled band structure to compare.")

    normalized = {}
    n_path_points = None
    for label, values in band_structures.items():
        values = np.asarray(values, dtype=float)
        if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
            raise ValueError(
                f"Band structure {label!r} must be a non-empty 2-D array."
            )
        if not np.isfinite(values).all():
            raise ValueError(f"Band structure {label!r} contains non-finite values.")
        if n_path_points is None:
            n_path_points = values.shape[1]
        elif values.shape[1] != n_path_points:
            raise ValueError(
                "All band structures must have the same number of path points."
            )
        normalized[str(label)] = values

    if bands is None:
        selected_bands = sorted({
            band_index
            for values in normalized.values()
            for band_index in range(values.shape[0])
        })
    else:
        selected_bands = sorted(set(bands))
        if not selected_bands or any(
            not isinstance(band, (int, np.integer)) or band < 1
            for band in selected_bands
        ):
            raise ValueError("bands must contain positive, 1-based integer band numbers.")
        missing = [
            band
            for band in selected_bands
            if any(band > values.shape[0] for values in normalized.values())
        ]
        if missing:
            raise ValueError(
                f"Selected band(s) {missing} are not present in every band structure."
            )

    fig, ax = plt.subplots()
    colors = plt.get_cmap("tab10")
    line_styles = ("-", "--", ":", "-.")
    for case_index, (label, values) in enumerate(normalized.items()):
        for band_number in selected_bands:
            band = band_number - 1
            ax.plot(
                values[band],
                color=colors(band % colors.N),
                linestyle=line_styles[case_index % len(line_styles)],
                linewidth=1.8,
                label=f"{label} - band {band_number}",
            )

    ax.set_xlabel("Bloch path point index")
    ax.set_ylabel("Frequency")
    ax.set_title("Band-structure comparison")
    ax.legend()
    fig.tight_layout()
    return fig, ax


def plot_frequency_surfaces(omega, kappa_x, kappa_y, lat_c=1.0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for band in range(omega.shape[2]):
        ax.plot_surface(
            kappa_x / np.pi,
            kappa_y / np.pi,
            omega[:, :, band],
            alpha=0.5,
            linewidth=0.3,
            label=f"Band {band + 1}"
        )
    
    MIN = np.min(omega)
    MAX = np.max(omega)
    T = np.array([0.0, 0.0])
    X = np.array([np.pi/lat_c, 0.0])
    M = np.array([np.pi/lat_c, np.pi/lat_c])
    ax.stem([T[0] / np.pi],[T[1] / np.pi],[MAX])
    ax.scatter(
        2*[T[0] / np.pi],
        2*[T[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$\Gamma$",
        s=55,
        # label=r"$\Gamma$",
    )
    ax.stem([X[0] / np.pi],[X[1] / np.pi],[MAX])
    ax.scatter(
        2*[X[0] / np.pi],
        2*[X[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$X$",
        s=55,
        # label=r"X",
    )
    ax.stem([M[0] / np.pi],[M[1] / np.pi],[MAX])
    ax.scatter(
        2*[M[0] / np.pi],
        2*[M[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$M$",
        s=55,
        # label=r"M",
    )
    
    ax.set_xlim([kappa_x[0].min() / np.pi, kappa_x[0].max() / np.pi])
    ax.set_ylim([kappa_y[:,0].min() / np.pi, kappa_y[:,0].max() / np.pi])
    ax.set_xlabel(r"$\kappa_x/\pi$")
    ax.set_ylabel(r"$\kappa_y/\pi$")
    ax.set_zlabel(r"Frequency ($\omega\frac{a}{2\pi}$)")
    ax.set_title("Square-lattice bands")
    ax.legend()
    ax.view_init(elev=20, azim=60)
    
    fig.tight_layout()
    return fig, ax


def _plot_lattice_point_bands(
    characteristic_values, kappa_x, kappa_y, estimate, point_label, bands=[n for n in range(2)]
):
    """Plot the lowest characteristic-value surfaces around a Bloch point."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cmaps = ["viridis", "plasma", "cividis"]
    for band in bands:
        ax.plot_surface(
            kappa_x / np.pi,
            kappa_y / np.pi,
            characteristic_values[:, :, band],
            # color=colors[band],
            alpha=0.78,
            linewidth=0,
            label=f"Band {band + 1}",
            cmap = cmaps[band % len(cmaps)]
        )

    point = estimate["kappa_center"]
    values_at_m = np.array(estimate["characteristic_values"])[bands]
    # values_at_m = estimate["characteristic_values"]
    ax.scatter(
        [point[0] / np.pi],
        [point[1] / np.pi],
        [np.mean(values_at_m)],
        color="black",
        marker="x",
        s=55,
        label="M point",
    )
    ax.set_xlabel(r"$\kappa_x/\pi$")
    ax.set_ylabel(r"$\kappa_y/\pi$")
    ax.set_zlabel(r"Normalized characteristic value ($\omega\frac{a}{2\pi}$)")
    ax.set_title(f"Characteristic-value bands near {point_label}")
    ax.legend(loc="lower right")
    ax.view_init(elev=20, azim=60)
    
    fig.tight_layout()
    return fig, ax


def plot_lattice_point_bands(
    kappa_center,
    lat_c=1.0,
    grid_size=21,
    n_bands=2,
    mesh_size=0.02,
    radius_fraction=0.1,
    estimate_samples=11,
    point_label=None,
    PEC_radius=0.1,
):
    """Compute and plot bands near a Bloch point; coordinates are in inverse length.

    For lattice constant ``lat_c``, common points are Gamma=(0,0),
    X=(pi/lat_c,0), and M=(pi/lat_c,pi/lat_c). Returns ``(fig, ax, estimate)``.
    """
    center = np.asarray(kappa_center, dtype=float)
    values, kappa_x, kappa_y, estimate = solve_lattice_point_bands(
        center,
        lat_c=lat_c,
        grid_size=grid_size,
        n_bands=n_bands,
        mesh_size=mesh_size,
        radius_fraction=radius_fraction,
        estimate_samples=estimate_samples,
        PEC_radius=PEC_radius,
    )
    if point_label is None:
        point_label = f"({center[0] / np.pi:.3g} pi, {center[1] / np.pi:.3g} pi)"
    fig, ax = _plot_lattice_point_bands(
        values, kappa_x, kappa_y, estimate, point_label
    )
    return fig, ax, estimate


def plot_m_point_bands(characteristic_values, kappa_x, kappa_y, estimate):
    """Plot precomputed characteristic-value bands close to M."""
    return _plot_lattice_point_bands(
        characteristic_values, kappa_x, kappa_y, estimate, "M"
    )


def plot_m_point_touching_bands(characteristic_values, kappa_x, kappa_y, estimate):
    """Plot precomputed characteristic-value bands close to M."""
    estimate["direct_gap"] = estimate["direct_gap"][1:3]
    return _plot_lattice_point_bands(
        characteristic_values, kappa_x, kappa_y, estimate, "M", bands=[1,2]
    )
    
    
def calc_m_point_lowest_touching(lat_c=1.0, PEC_radius=0.1, mesh_size=0.02):
    M = np.array([np.pi / lat_c, np.pi / lat_c])
    mesh_size = min(mesh_size, 0.5 * PEC_radius)
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=PEC_radius)
    characteristic_vals = _lowest_eigenfrequencies(matrices, M, 3)
    touching = characteristic_vals[1:3]
    return touching


def compare_bands_2_and_3(lat_c=1.0, grid_size=21, mesh_size=0.02, PEC_radius=0.1):
    omega1, kappa_x, kappa_y = solve_frequency_surfaces(lat_c=lat_c, grid_size=2*grid_size-1, n_omega=3, mesh_size=mesh_size, radius_fraction=1.0, PEC_radius=PEC_radius)
    
    omega = np.empty((grid_size, grid_size, 2), dtype=float)
    SHIFT = (2*grid_size-2)//4
    # MEDIAN = np.median(omega1[0:grid_size,0:grid_size,2])
    MEDIAN = np.median((omega1[0:grid_size,0:grid_size,2]+omega1[0:grid_size,0:grid_size,1])/2)
    omega[:,:,0] = omega1[SHIFT:grid_size+SHIFT, SHIFT:grid_size+SHIFT, 1]
    omega[:,:,1] = 2*MEDIAN - omega1[0:grid_size,0:grid_size,2]
    # omega[:,:,2] = 2*MEDIAN-omega[:,:,2]
    # shift = [np.pi/lat_c, np.pi/lat_c]
    # for row in range(grid_size):
    #     for column in range(grid_size):
    #         shift_omega[row,column,0] = omega[row+SHIFT,column+SHIFT,1]
    #         shift_omega[row,column,1] = 2*MEDIAN - omega[row,column,2]
    # omega = shift_omega
    kappa_x = kappa_x[SHIFT:grid_size+SHIFT, SHIFT:grid_size+SHIFT]
    kappa_y = kappa_y[SHIFT:grid_size+SHIFT, SHIFT:grid_size+SHIFT]
    
    diff = omega[:,:,1]-omega[:,:,0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for band in range(omega.shape[2]):
        ax.plot_surface(
            kappa_x / np.pi,
            kappa_y / np.pi,
            omega[:, :, band],
            alpha=0.5,
            linewidth=0.3,
            label=f"Band {band + 1}"
        )
    
    ax.plot_surface(
        kappa_x / np.pi,
        kappa_y / np.pi,
        diff,
        alpha=0.5,
        linewidth=0.3,
        label=f"Difference"
    )
    MIN = np.min(omega)
    MAX = np.max(omega)
    T = np.array([0.0, 0.0])
    X = np.array([np.pi/lat_c, 0.0])
    M = np.array([np.pi/lat_c, np.pi/lat_c])
    ax.stem([T[0] / np.pi],[T[1] / np.pi],[MAX])
    ax.scatter(
        2*[T[0] / np.pi],
        2*[T[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$\Gamma$",
        s=55,
        # label=r"$\Gamma$",
    )
    ax.stem([X[0] / np.pi],[X[1] / np.pi],[MAX])
    ax.scatter(
        2*[X[0] / np.pi],
        2*[X[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$X$",
        s=55,
        # label=r"X",
    )
    ax.stem([M[0] / np.pi],[M[1] / np.pi],[MAX])
    ax.scatter(
        2*[M[0] / np.pi],
        2*[M[1] / np.pi],
        [0,1.1*MAX],
        color="black",
        marker=r"$M$",
        s=55,
        # label=r"M",
    )
    
    ax.set_xlim([kappa_x[0].min() / np.pi, kappa_x[0].max() / np.pi])
    ax.set_ylim([kappa_y[:,0].min() / np.pi, kappa_y[:,0].max() / np.pi])
    ax.set_xlabel(r"$\kappa_x/\pi$")
    ax.set_ylabel(r"$\kappa_y/\pi$")
    ax.set_zlabel("Frequency")
    ax.set_title("Square-lattice bands near the Gamma point")
    ax.legend()
    ax.view_init(elev=20, azim=60)
    
    fig.tight_layout()
    
    print(f"After shifting by (pi/a,/pi/a) and flipping band 3 about median midpoint between bands 2 and 3 ---")
    print(f"Min Difference: {np.min(diff):.6f}, Max Difference: {np.max(diff):.6f}, Median Difference: {np.median(diff):.6f}")
    return fig, ax