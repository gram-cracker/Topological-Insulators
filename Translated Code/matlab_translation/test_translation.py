import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from Finite_Element.FEM_matrices_C4_PEC_TE import FEM_matrices_C4_PEC_TE
from Mesh.mesh_generation import mesh_generation_C4_PEC, points_in_polygon
from Mesh.meshpoly import meshpoly
from main import (
    _lowest_eigenfrequencies,
    _solve_lattice_point_bands,
    _solve_m_point_bands,
    main,
    plot_band_structure,
    plot_band_comparison,
    plot_frequency_surfaces,
    plot_lattice_point_bands,
    plot_m_point_bands,
    solve_band_structure,
    solve_frequency_surfaces,
    solve_lattice_point_bands,
)


class MeshTranslationTests(unittest.TestCase):
    def test_pec_mesh_respects_cell_and_inclusion(self):
        xy_inc = np.array([
            [0.4, 0.5],
            [0.5, 0.4],
            [0.6, 0.5],
            [0.5, 0.6],
        ])
        points, triangles, regions = mesh_generation_C4_PEC(
            mesh_size=0.1, xy_inc=xy_inc
        )

        self.assertEqual(triangles.shape[1], 3)
        self.assertEqual(regions.shape, (triangles.shape[0],))
        self.assertGreater(triangles.shape[0], 0)
        centroids = points[triangles].mean(axis=1)
        self.assertFalse(np.any(points_in_polygon(centroids, xy_inc)))
        element_points = points[triangles]
        edge_lengths = np.linalg.norm(
            element_points - np.roll(element_points, -1, axis=1), axis=2
        )
        self.assertLessEqual(float(edge_lengths.max()), 0.1 + 1e-12)
        self.assertTrue(np.all(points >= 0.0))
        self.assertTrue(np.all(points <= 1.0))

    def test_polygon_mesh_helper_returns_valid_triangles(self):
        polygon = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        points, triangles = meshpoly(polygon, None, None, None)
        self.assertEqual(triangles.shape[1], 3)
        self.assertGreater(points.shape[0], 3)
        self.assertTrue(np.all(triangles >= 0))
        self.assertTrue(np.all(triangles < points.shape[0]))
        self.assertTrue(all(np.any(np.all(points == vertex, axis=1)) for vertex in polygon))
        centroids = points[triangles].mean(axis=1)
        self.assertTrue(np.all(points_in_polygon(centroids, polygon, include_boundary=False)))


class FEMTranslationTests(unittest.TestCase):
    def test_assembled_mass_and_stiffness_match_linear_triangle(self):
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, -1.0],
            [-1.0, 2.0],
            [2.0, -1.0],
            [2.0, 2.0],
        ])
        triangles = np.array([[0, 1, 2]])
        inclusion = np.array([[1.4, 1.4], [1.6, 1.4], [1.6, 1.6], [1.4, 1.6]])
        stiffness, mx, my, mass, _ = FEM_matrices_C4_PEC_TE(
            points, triangles, np.array([-1.0, 2.0, -1.0, 2.0]), inclusion
        )
        expected_stiffness = np.array([
            [1.0, -0.5, -0.5],
            [-0.5, 0.5, 0.0],
            [-0.5, 0.0, 0.5],
        ])
        expected_mass = np.array([
            [1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0],
            [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
            [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0],
        ])
        np.testing.assert_allclose(stiffness.toarray()[:3, :3], expected_stiffness)
        np.testing.assert_allclose(mass.toarray()[:3, :3], expected_mass)
        np.testing.assert_allclose(mx.toarray()[:3, :3], np.tile([-1.0 / 6, 1.0 / 6, 0.0], (3, 1)))
        np.testing.assert_allclose(my.toarray()[:3, :3], np.tile([-1.0 / 6, 0.0, 1.0 / 6], (3, 1)))

    def test_generalized_eigenfrequencies_are_sorted_and_real(self):
        matrices = (
            sparse.diags([1.0, 4.0], format="csr"),
            sparse.csr_matrix((2, 2)),
            sparse.csr_matrix((2, 2)),
            sparse.eye(2, format="csr"),
        )
        frequencies = _lowest_eigenfrequencies(matrices, np.zeros(2), 1)
        np.testing.assert_allclose(frequencies, [1.0])


class PlotAndSolverTests(unittest.TestCase):
    def test_reduced_band_and_grid_solvers_and_plots(self):
        bands, kappa = solve_band_structure(
            n_kappa=2, n_omega=2, mesh_size=0.25
        )
        self.assertEqual(bands.shape, (2, 8))
        self.assertEqual(kappa.shape, (8, 2))
        self.assertTrue(np.isfinite(bands).all())

        surfaces, kappa_x, kappa_y = solve_frequency_surfaces(
            grid_size=2, n_omega=2, mesh_size=0.25
        )
        self.assertEqual(surfaces.shape, (2, 2, 2))
        self.assertEqual(kappa_x.shape, (2, 2))
        self.assertEqual(kappa_y.shape, (2, 2))
        self.assertTrue(np.isfinite(surfaces).all())

        band_fig, _ = plot_band_structure(bands)
        comparison_fig, comparison_ax = plot_band_comparison(
            {"baseline": bands, "adjusted": bands + 0.1}
        )
        selected_fig, selected_ax = plot_band_comparison(
            {
                "baseline": np.arange(15.0).reshape(3, 5),
                "adjusted": np.arange(15.0, 30.0).reshape(3, 5),
            },
            bands=(1, 3),
        )
        surface_fig, _ = plot_frequency_surfaces(surfaces, kappa_x, kappa_y)
        self.assertEqual(len(band_fig.axes), 1)
        self.assertEqual(len(comparison_fig.axes), 1)
        self.assertEqual(len(comparison_ax.lines), 4)
        self.assertEqual(len(selected_ax.lines), 4)
        self.assertEqual(
            [line.get_label() for line in selected_ax.lines],
            [
                "baseline - band 1",
                "baseline - band 3",
                "adjusted - band 1",
                "adjusted - band 3",
            ],
        )
        self.assertEqual(len(surface_fig.axes), 1)
        import matplotlib.pyplot as plt

        plt.close(band_fig)
        plt.close(comparison_fig)
        plt.close(selected_fig)
        plt.close(surface_fig)

    def test_band_comparison_rejects_different_path_sampling(self):
        with self.assertRaisesRegex(ValueError, "same number of path points"):
            plot_band_comparison({
                "coarse": np.ones((2, 5)),
                "fine": np.ones((2, 7)),
            })

    def test_band_comparison_rejects_missing_selected_band(self):
        with self.assertRaisesRegex(ValueError, "not present in every band structure"):
            plot_band_comparison(
                {"three bands": np.ones((3, 5)), "two bands": np.ones((2, 5))},
                bands=(1, 3),
            )

    def test_m_point_characteristic_values_and_dirac_gap_estimate(self):
        def conical_values(_matrices, kappa, _n_bands):
            center = np.array([np.pi, np.pi])
            distance = float(np.linalg.norm(kappa - center))
            return np.array([4.0 - 0.5 * distance, 4.0 + 0.5 * distance])

        with patch("main._lowest_characteristic_values", side_effect=conical_values):
            values, kappa_x, kappa_y, estimate = _solve_lattice_point_bands(
                (), np.array([np.pi, np.pi]), 1.0, 5, 2, 0.1, 7
            )

        self.assertEqual(values.shape, (5, 5, 2))
        self.assertEqual(kappa_x.shape, (5, 5))
        self.assertEqual(kappa_y.shape, (5, 5))
        np.testing.assert_allclose(estimate["characteristic_values"], [4.0, 4.0])
        self.assertAlmostEqual(estimate["direct_gap"], 0.0)
        for slope in estimate["directional_gap_slopes"].values():
            self.assertAlmostEqual(slope, 1.0)

        fig, _ = plot_m_point_bands(
            values,
            kappa_x,
            kappa_y,
            {**estimate, "m_point": estimate["kappa_center"]},
        )
        self.assertEqual(len(fig.axes), 1)
        plt.close(fig)

    def test_custom_lattice_point_plot_uses_requested_center(self):
        center = np.array([np.pi, 0.0])

        def conical_values(_matrices, kappa, _n_bands):
            distance = float(np.linalg.norm(kappa - center))
            return np.array([3.0 - 0.25 * distance, 3.0 + 0.25 * distance])

        with patch("main._assemble_cell", return_value=()):
            with patch("main._lowest_characteristic_values", side_effect=conical_values):
                fig, ax, estimate = plot_lattice_point_bands(
                    center,
                    grid_size=5,
                    radius_fraction=0.1,
                    estimate_samples=5,
                    point_label="X",
                )

        np.testing.assert_allclose(estimate["kappa_center"], center)
        np.testing.assert_allclose(estimate["characteristic_values"], [3.0, 3.0])
        self.assertAlmostEqual(estimate["direct_gap"], 0.0)
        self.assertEqual(ax.get_title(), "Two lowest characteristic-value bands near X")
        plt.close(fig)

    def test_main_saves_all_graphs(self):
        bands = np.array([[0.2, 0.3], [0.4, 0.5]])
        path = np.zeros((2, 2))
        surfaces = np.ones((2, 2, 2))
        grid = np.array([[0.0, 1.0], [0.0, 1.0]])
        characteristic_values = np.ones((3, 3, 2))
        m_grid = np.arange(3.0).reshape(1, 3) + np.arange(3.0).reshape(3, 1)
        estimate = {
            "m_point": np.array([np.pi, np.pi]),
            "kappa_center": np.array([np.pi, np.pi]),
            "characteristic_values": np.array([1.0, 1.1]),
            "direct_gap": 0.1,
            "mean_gap_slope": 0.5,
            "directional_gap_slopes": {"kx": 0.5},
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch("main.ROOT", Path(temp_dir)):
                with patch("main._assemble_cell", return_value=()):
                    with patch("main._solve_band_structure", return_value=(bands, path)):
                        with patch(
                            "main._solve_frequency_surfaces",
                            return_value=(surfaces, grid, grid),
                        ):
                            with patch(
                                "main._solve_m_point_bands",
                                return_value=(
                                    characteristic_values,
                                    m_grid,
                                    m_grid,
                                    estimate,
                                ),
                            ):
                                main()
            self.assertTrue((Path(temp_dir) / "outputs" / "band_structure.png").is_file())
            self.assertTrue((Path(temp_dir) / "outputs" / "frequency_surfaces.png").is_file())
            self.assertTrue((Path(temp_dir) / "outputs" / "m_point_lower_bands.png").is_file())


if __name__ == "__main__":
    unittest.main()
