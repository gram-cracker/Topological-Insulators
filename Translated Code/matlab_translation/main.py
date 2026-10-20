import sys
from pathlib import Path

import matplotlib
# matplotlib.use("Agg")
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Finite_Element.FEM_matrices_C4_PEC_TE import FEM_matrices_C4_PEC_TE
from Mesh.mesh_generation import mesh_generation_C4_PEC


def geo_PEC_inclusions(lat_c=1.0):
    """Return the polygonal circular PEC boundary used by the MATLAB model."""
    radius = 0.1 * lat_c
    center = np.array([0.5 * lat_c, 0.5 * lat_c])
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    return center + radius * np.column_stack((np.cos(theta), np.sin(theta)))


def _assemble_cell(lat_c=1.0, mesh_size=0.02):
    xy_box = np.array([0.0, lat_c, 0.0, lat_c])
    xy_cell = np.array([
        [0.0, 0.0],
        [lat_c, 0.0],
        [lat_c, lat_c],
        [0.0, lat_c],
    ])
    xy_inc = geo_PEC_inclusions(lat_c)
    p, t, _ = mesh_generation_C4_PEC(xy_cell, xy_inc, mesh_size=mesh_size)
    M, Mx, My, Me, _ = FEM_matrices_C4_PEC_TE(p, t, xy_box, xy_inc)
    return M, Mx, My, Me


def _lowest_characteristic_values(matrices, kappa, n_omega):
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
        eigenvalues = eigsh(
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


def _lowest_eigenfrequencies(matrices, kappa, n_omega):
    return np.sqrt(_lowest_characteristic_values(matrices, kappa, n_omega))


# def _band_path_TXMT(lat_c, n_kappa):
#     """Return the T-to-X-to-M-to-T Bloch path in MATLAB coordinates."""
#     points = [
#         (np.linspace(0.0, np.pi, n_kappa),np.zeros(n_kappa)),
#         (np.pi * np.ones(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
#         (np.linspace(np.pi, 0.0, n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
#     ]
#     return np.vstack([
#         np.column_stack((kx, ky)) / lat_c
#         for kx, ky in points
#     ])
    
    
def _band_path(lat_c, n_kappa, path_type="MXTxM"):
    if path_type == "TXMT":
        """Return the T-to-X-to-M-to-T Bloch path in MATLAB coordinates."""
        points = [
            (np.linspace(0.0, np.pi, n_kappa),np.zeros(n_kappa)),
            (np.pi * np.ones(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
            (np.linspace(np.pi, 0.0, n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
        ]
    elif path_type == "MXTxM":
        """Return the M-to-X-to-T-to-X'-to-M Bloch path in MATLAB coordinates."""
        points = [
            (np.pi * np.ones(n_kappa), np.linspace(np.pi, 0.0, n_kappa)),
            (np.linspace(np.pi, 0.0, n_kappa), np.zeros(n_kappa)),
            (np.zeros(n_kappa), np.linspace(0.0, np.pi, n_kappa)),
            (np.linspace(0.0, np.pi, n_kappa), np.pi * np.ones(n_kappa)),
        ]
    else:
        raise ValueError("path_type must be either 'TXMT' or 'MXTxM'.")
    return np.vstack([
        np.column_stack((kx, ky)) / lat_c
        for kx, ky in points
    ])


def solve_band_structure(lat_c=1.0, n_kappa=100, n_omega=7, mesh_size=0.02):
    # """Compute bands along the MATLAB M-to-X-to-T-to-M Bloch path."""
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size)
    return _solve_band_structure(matrices, lat_c, n_kappa, n_omega)


def _solve_band_structure(matrices, lat_c, n_kappa, n_omega, path_type="MXTxM"):
    # if path_type == "TXMT":
    #     kappa_vec = _band_path_TXMT(lat_c, n_kappa)
    #     # title = "$\Gamma$-X-M-$\Gamma$"
    # else:
    #     kappa_vec = _band_path_MXTxM(lat_c, n_kappa)
    #     # title = "M-X-$\Gamma$-X'-M"
    kappa_vec = _band_path(lat_c, n_kappa, path_type)
    omega = np.empty((n_omega, kappa_vec.shape[0]), dtype=float)
    for n, kappa in enumerate(kappa_vec):
        omega[:, n] = (
            _lowest_eigenfrequencies(matrices, kappa, n_omega)
            * lat_c / (2.0 * np.pi)
        )
    return omega, kappa_vec


def solve_frequency_surfaces(
    lat_c=1.0,
    grid_size=21,
    n_omega=5,
    mesh_size=0.02,
    kappa_center=(0.0, 0.0),
):
    """Compute the lowest bands on the rectangular Bloch-vector grid in MATLAB."""
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size)
    return _solve_frequency_surfaces(matrices, lat_c, grid_size, n_omega, kappa_center)


def _solve_frequency_surfaces(matrices, lat_c, grid_size, n_omega, kappa_center):
    kx = kappa_center[0] + (0.1 * np.pi / lat_c) * np.linspace(-1.0, 1.0, grid_size)
    ky = kappa_center[1] + (0.1 * np.pi / lat_c) * np.linspace(-1.0, 1.0, grid_size)
    kappa_x, kappa_y = np.meshgrid(kx, ky)
    omega = np.empty((grid_size, grid_size, n_omega), dtype=float)
    for row in range(grid_size):
        for column in range(grid_size):
            omega[row, column, :] = _lowest_eigenfrequencies(
                matrices,
                np.array([kappa_x[row, column], kappa_y[row, column]]),
                n_omega,
            )
    return omega, kappa_x, kappa_y


def solve_m_point_bands(
    lat_c=1.0,
    grid_size=21,
    n_bands=2,
    mesh_size=0.02,
    radius_fraction=0.1,
    estimate_samples=11,
):
    """Solve the two lowest characteristic-value bands in a neighborhood of M."""
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size)
    return _solve_m_point_bands(
        matrices, lat_c, grid_size, n_bands, radius_fraction, estimate_samples
    )


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
):
    """Compute nearby bands around a Bloch point, e.g. Gamma (0,0), X (pi/a,0), M (pi/a,pi/a)."""
    if lat_c <= 0:
        raise ValueError("lat_c must be positive.")
    center = np.asarray(kappa_center, dtype=float)
    if center.shape != (2,) or not np.isfinite(center).all():
        raise ValueError("kappa_center must contain two finite Bloch-vector coordinates.")
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size)
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
            )

    values_at_center = _lowest_characteristic_values(
        matrices, kappa_center, n_bands
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
    ax.set_ylabel("Frequency")
    
    if path_type == "TXMT":
        ax.set_xticks([0, omega.shape[1] // 3, 2 * omega.shape[1] // 3, omega.shape[1] - 1])
        ax.set_xticklabels(["$\Gamma$", "X", "M", "$\Gamma$"])
        title = "$\Gamma$-X-M-$\Gamma$ band structure"
    elif path_type == "MXTxM":
        ax.set_xticks([0, omega.shape[1] // 4, 2 * omega.shape[1] // 4, 3 * omega.shape[1] // 4, omega.shape[1] - 1])
        ax.set_xticklabels(["M", "X", "$\Gamma$", "X'", "M"])
        title = "M-X-$\Gamma$-X'-M band structure"
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


def plot_frequency_surfaces(omega, kappa_x, kappa_y):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for band in range(omega.shape[2]):
        ax.plot_surface(
            kappa_x / np.pi,
            kappa_y / np.pi,
            omega[:, :, band],
            alpha=0.75,
            linewidth=0,
        )
    ax.set_xlabel(r"$\kappa_x/\pi$")
    ax.set_ylabel(r"$\kappa_y/\pi$")
    ax.set_zlabel("Frequency")
    ax.set_title("Square-lattice bands near the Gamma point")
    fig.tight_layout()
    return fig, ax


def _plot_lattice_point_bands(
    characteristic_values, kappa_x, kappa_y, estimate, point_label
):
    # """Plot the two lowest characteristic-value surfaces around a Bloch point."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    colors = ("tab:blue", "tab:orange")
    for band in range(characteristic_values.shape[2]):
        ax.plot_surface(
            kappa_x / np.pi,
            kappa_y / np.pi,
            characteristic_values[:, :, band],
            # color=colors[band],
            alpha=0.78,
            linewidth=0,
            label=f"Band {band + 1}",
        )

    point = estimate["kappa_center"]
    values_at_m = estimate["characteristic_values"]
    ax.scatter(
        [point[0] / np.pi],
        [point[1] / np.pi],
        [np.mean(values_at_m)],
        color="black",
        marker="x",
        s=55,
        label="M point",
    )
    slopes = estimate["directional_gap_slopes"]
    slope_text = "\n".join(
        f"{direction}: {slope:.4g}"
        for direction, slope in slopes.items()
    )
    ax.text2D(
        0.02,
        0.98,
        (
            f"At {point_label}: lambda1={values_at_m[0]:.6g}, "
            f"lambda2={values_at_m[1]:.6g}\n"
            f"Direct gap={[f'{gap:.4g}' for gap in estimate['direct_gap']]}\n"
            f"Mean near-M gap slope={estimate['mean_gap_slope']:.4g}\n"
            f"Directional slopes (d gap / d|k|):\n{slope_text}"
        ),
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    ax.set_xlabel(r"$\kappa_x/\pi$")
    ax.set_ylabel(r"$\kappa_y/\pi$")
    ax.set_zlabel("Characteristic value (lambda)")
    ax.set_title(f"Lowest characteristic-value bands near {point_label}")
    ax.legend(loc="lower right")
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
    estimate["characteristic_values"] = estimate["characteristic_values"][1:3]
    estimate["direct_gap"] = estimate["direct_gap"][1:3]
    return _plot_lattice_point_bands(
        characteristic_values[:,:,1:3], kappa_x, kappa_y, estimate, "M"
    )

def main():
    print(matplotlib.get_backend())
    # Graph images are written here; change this path to choose another location.
    output_dir = ROOT / "outputs"
    output_dir.mkdir(exist_ok=True)

    # Shared lattice/FEM setup: lattice constant a=1.0 and mesh spacing h=0.02.
    # Smaller h refines the mesh (and can shift results) at higher runtime cost.
    lat_c = 1.0
    mesh_size = 0.02
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size)

    # # Band path: 100 samples on each of the path segments,
    # # with the lowest 7 bands computed at every Bloch vector.
    # path_type = "TXMT"  # or "MXTxM"
    # path_samples = 100
    # path_band_count = 7
    # omega, _ = _solve_band_structure(
    #     matrices, lat_c, path_samples, path_band_count, path_type
    # )
    # band_fig, _ = plot_band_structure(omega, path_type)
    # band_fig.savefig(output_dir / "band_structure.png", dpi=180)
    # plt.close(band_fig)

    # # Frequency surface centered at Gamma (0,0), spanning +/-0.1*pi/a
    # # on each reciprocal-lattice axis, sampled on a 21-by-21 grid.
    # # Increase grid_size for finer sampling; change center to inspect elsewhere.
    # gamma_grid_size = 21
    # gamma_band_count = 5
    # gamma_center = (0.0, 0.0)
    # surface_omega, kappa_x, kappa_y = _solve_frequency_surfaces(
    #     matrices,
    #     lat_c,
    #     gamma_grid_size,
    #     gamma_band_count,
    #     gamma_center,
    # )
    # surface_fig, _ = plot_frequency_surfaces(surface_omega, kappa_x, kappa_y)
    # surface_fig.savefig(output_dir / "frequency_surfaces.png", dpi=180)
    # plt.close(surface_fig)

    # M-point plot: center is (pi/a, pi/a); radius_fraction=0.1 samples
    # +/-0.1*pi/a along each axis. grid_size controls the 2-D plot resolution.
    # Two lowest characteristic-value bands are shown; estimate_samples controls
    # the number of radial samples used for the near-M gap-slope estimate.
    m_grid_size = 31
    m_band_count = 3
    m_radius_fraction = 0.1
    m_estimate_samples = 11
    m_values, m_kx, m_ky, m_estimate = _solve_m_point_bands(
        matrices,
        lat_c,
        m_grid_size,
        m_band_count,
        m_radius_fraction,
        m_estimate_samples,
    )
    # m_fig, _ = plot_m_point_bands(m_values, m_kx, m_ky, m_estimate)
    m_fig, _ = plot_m_point_touching_bands(m_values, m_kx, m_ky, m_estimate)
    plt.show()
    m_fig.savefig(output_dir / "m_point_lower_bands.png", dpi=180)
    plt.close(m_fig)

    values_at_m = m_estimate["characteristic_values"]
    normalized_frequencies = np.sqrt(values_at_m) / (2.0 * np.pi)
    print(f"Saved band-structure graphs to {output_dir}")
    print(
        "M-point estimate: "
        f"lambda1={values_at_m[0]:.8g}, lambda2={values_at_m[1]:.8g}, "
        f"direct gap={[f'{gap:.8g}' for gap in m_estimate['direct_gap']]}"
    )
    print(
        "M-point normalized frequencies: "
        f"f1={normalized_frequencies[0]:.8g}, f2={normalized_frequencies[1]:.8g}"
    )
    print(
        "Mean near-M characteristic-value gap slope "
        f"d(lambda2-lambda1)/d|k|={m_estimate['mean_gap_slope']:.8g}"
    )


if __name__ == "__main__":
    main()
