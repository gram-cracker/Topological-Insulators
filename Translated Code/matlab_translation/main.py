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
from solver import compare_to_analytic_estimate
from bandpath import (
    _solve_band_structure,
    _solve_lattice_point_bands,
    plot_band_structure,
    plot_m_point_touching_bands,
)
from inclusions import _assemble_cell, geo_PEC_inclusions_squircle, geo_PEC_inclusions, geo_PEC_inclusions_triangle




# import warnings
# warnings.filterwarnings("ignore")



def main():
    print("MatPlotLib backend: " + matplotlib.get_backend())
    # Graph images are written here; change this path to choose another location.
    output_dir = ROOT / "outputs"
    output_dir.mkdir(exist_ok=True)
    
    time_stamp = time.perf_counter()
    
    # Shared lattice/FEM setup: lattice constant a=1.0 and mesh spacing h=0.02.
    # Smaller h refines the mesh (and can shift results) at higher runtime cost.
    lat_c = 1.0
    mesh_size = 0.02
    eps = 0.1  # inclusion radius parameter
    #inc1 = geo_PEC_inclusions_squircle(lat_c, PEC_radius=eps)
    inc1 = geo_PEC_inclusions_triangle(lat_c, PEC_radius=eps)
    #inc1 = geo_PEC_inclusions(lat_c, PEC_radius=eps)
    plt.plot(inc1[:,0],inc1[:,1])
    plt.show()

    # Compute touching characteristic values at M for several PEC inclusion radii and compare to analyic estimate.
    compare_to_analytic_estimate(lat_c=lat_c, radii=[0.2,0.1,0.05,0.01])
    #exit() for when you want to just gain eigen estimates
    matrices = _assemble_cell(lat_c=lat_c, mesh_size=mesh_size, PEC_radius=eps)
    
    # Surface plots!    
    # omega, kappa_x, kappa_y = solve_frequency_surfaces(grid_size=41, n_omega=3, radius_fraction=1)
    # surface_fig, _ = plot_frequency_surfaces(omega, kappa_x, kappa_y)
    # surface_fig.savefig(output_dir / "surface_plot.png", dpi=180)
    # plt.show()
    # surface_fig.savefig(output_dir / "surface_plot_snapshot.png", dpi=180)
    # plt.close(surface_fig)

    # Band path: 100 samples on each of the path segments,
    # with the lowest 7 bands computed at every Bloch vector.
    path_type = "TXMT"  # or "MXTxM" or "TXMT" or "TX'M'T"
    path_samples = 200
    path_band_count = 7
    omega, _ = _solve_band_structure(
        matrices, lat_c, path_samples, path_band_count, path_type
    )
    band_fig, _ = plot_band_structure(omega, path_type)
    band_fig.savefig(output_dir / f"Triangle_band_structure_1_eps=.1{path_type}.png", dpi=180)
    plt.close(band_fig)
    
    exit()
    # M-point plot: center is (pi/a, pi/a); radius_fraction=0.1 samples
    # +/-0.1*pi/a along each axis. grid_size controls the 2-D plot resolution.
    # Two lowest characteristic-value bands are shown; estimate_samples controls
    # the number of radial samples used for the near-M gap-slope estimate.
    M = np.array([np.pi / lat_c, np.pi / lat_c])
    m_grid_size = 41
    m_band_count = 3
    m_radius_fraction = 0.1
    m_estimate_samples = 11
    # m_values, m_kx, m_ky, m_estimate = _solve_m_point_bands(
    #     matrices,
    #     lat_c,
    #     m_grid_size,
    #     m_band_count,
    #     m_radius_fraction,
    #     m_estimate_samples,
    # )
    m_values, m_kx, m_ky, m_estimate = _solve_lattice_point_bands(
        matrices=matrices,
        kappa_center=M,
        lat_c=lat_c,
        grid_size=m_grid_size,
        n_bands=m_band_count,
        radius_fraction=m_radius_fraction,
        estimate_samples=m_estimate_samples
    )
    # m_fig, _ = plot_m_point_bands(m_values, m_kx, m_ky, m_estimate)
    m_fig, _ = plot_m_point_touching_bands(m_values, m_kx, m_ky, m_estimate)
    m_fig.savefig(output_dir / "m_point_touch.png", dpi=180)
    plt.show()
    m_fig.savefig(output_dir / "m_point_touch_snapshot.png", dpi=180)
    plt.close(m_fig)

    values_at_m = m_estimate["characteristic_values"]
    normalized_frequencies = np.sqrt(values_at_m) / (2.0 * np.pi)
    print(f"Saved band-structure graphs to {output_dir}")
    print(
        "M-point estimate: "
        f"lambda_1={values_at_m[0]:.8g}, lambda_2={values_at_m[1]:.8g}, "
        f"direct gap={[f'{gap:.8g}' for gap in m_estimate['direct_gap']]}"
    )
    print(
        "M-point normalized frequencies: "
        f"f1={normalized_frequencies[0]:.8g}, f2={normalized_frequencies[1]:.8g}"
    )
    print(
        "Mean near-M characteristic-value gap slope "
        f"d(lambda_2-lambda_1)/d|k|={m_estimate['mean_gap_slope']:.8g}"
    )
    print(f"Program took {time.perf_counter() - time_stamp} seconds")


if __name__ == "__main__":
    main()
