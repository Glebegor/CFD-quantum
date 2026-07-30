import sys
import unittest
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from lbm_solver import D2Q9Solver, equilibrium, physical_fields  # noqa: E402
from macrocomputations import reconstruct_fields  # noqa: E402


class LbmSolverTests(unittest.TestCase):
    def test_equilibrium_reconstructs_requested_fields(self):
        density = np.full((3, 4), 1.02)
        velocity = np.zeros((3, 4, 2))
        velocity[..., 0] = 0.05
        reconstructed_density, reconstructed_velocity, _ = reconstruct_fields(
            equilibrium(density, velocity)
        )
        self.assertTrue(np.allclose(reconstructed_density, density))
        self.assertTrue(np.allclose(reconstructed_velocity, velocity))

    def test_inlet_is_one_metre_per_second_after_conversion(self):
        solver = D2Q9Solver(8, 4, inlet_velocity_lattice=0.05)
        solver.step()
        fields = physical_fields(solver.f, 1.225, 101325.0, 20.0)
        self.assertTrue(np.allclose(fields["velocity"][:, 0, 0], 1.0))
        self.assertTrue(np.allclose(fields["velocity"][:, 0, 1], 0.0))
        self.assertTrue(np.all(np.isfinite(fields["pressure"])))


if __name__ == "__main__":
    unittest.main()
