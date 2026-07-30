import sys
import unittest
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from lbm_solver import apply_channel_boundaries  # noqa: E402
from macrocomputations import reconstruct_fields  # noqa: E402


class ChannelBoundaryTests(unittest.TestCase):
    def test_left_inlet_and_right_outlet_are_both_applied(self):
        rng = np.random.default_rng(42)
        f = rng.random((9, 4, 8))
        expected_outlet = f[:, :, -2].copy()

        bounded = apply_channel_boundaries(f, 0.03)
        density, velocity, _ = reconstruct_fields(bounded)

        np.testing.assert_allclose(density[:, 0], 1.0, atol=1e-14)
        np.testing.assert_allclose(velocity[:, 0, 0], 0.03, atol=1e-14)
        np.testing.assert_allclose(velocity[:, 0, 1], 0.0, atol=1e-14)
        np.testing.assert_allclose(bounded[:, :, -1], expected_outlet)


if __name__ == "__main__":
    unittest.main()
