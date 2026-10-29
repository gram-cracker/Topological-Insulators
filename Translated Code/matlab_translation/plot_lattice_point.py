from pathlib import Path

import matplotlib.pyplot as plt
#import cupy as np
import numpy as np

from main import plot_lattice_point_bands

fig, ax, estimate = plot_lattice_point_bands(
    (0, 0),  # X point; try (0, 0) for Gamma or (np.pi, np.pi) for M
    point_label="Gamma",
)

output_dir = Path(__file__).resolve().parent / "outputs"
output_dir.mkdir(exist_ok=True)
fig.savefig(output_dir / "Gamma_point_bands.png", dpi=180)

print(estimate)
plt.show()