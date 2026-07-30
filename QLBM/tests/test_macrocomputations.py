import sys
import unittest
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from macrocomputations import reconstruct_fields  # noqa: E402


class MacrocomputationTests(unittest.TestCase):
    def test_reconstruct_fields(self):
        f = np.zeros((9, 4, 4))
        f[1, 2, 1] = 0.25
        f[3, 2, 1] = 0.25
        density, velocity, pressure = reconstruct_fields(f)
        self.assertEqual(density.shape, (4, 4))
        self.assertEqual(velocity.shape, (4, 4, 2))
        self.assertEqual(pressure.shape, (4, 4))
        self.assertAlmostEqual(density[2, 1], 0.5)
        self.assertTrue(np.allclose(velocity[2, 1], (0.0, 0.0)))
        self.assertAlmostEqual(pressure[2, 1], 1.0 / 6.0)


if __name__ == "__main__":
    unittest.main()
