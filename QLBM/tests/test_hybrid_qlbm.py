import sys
import unittest
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Statevector

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

import gates  # noqa: E402
from hybrid_qlbm import (  # noqa: E402
    classical_stream,
    populations_to_amplitudes,
    probabilities_to_populations,
)


class HybridQlbmTests(unittest.TestCase):
    def test_amplitude_round_trip(self):
        generator = np.random.default_rng(7)
        f = generator.random((9, 4, 4))
        amplitudes, total = populations_to_amplitudes(f, 2, 2, 4)
        reconstructed, valid = probabilities_to_populations(
            np.abs(amplitudes) ** 2,
            total,
            2,
            2,
            4,
        )
        self.assertAlmostEqual(valid, 1.0)
        self.assertTrue(np.allclose(reconstructed, f))

    def test_quantum_stream_matches_classical_stream(self):
        generator = np.random.default_rng(11)
        f = generator.random((9, 4, 4))
        amplitudes, total = populations_to_amplitudes(f, 2, 2, 4)
        circuit = gates.qk.QuantumCircuit(8)
        from qiskit.circuit.library import StatePreparation

        circuit.append(StatePreparation(amplitudes), circuit.qubits)
        circuit.append(gates.StreamG(2, 2, 4), circuit.qubits)
        probabilities = Statevector.from_instruction(circuit).probabilities()
        actual, valid = probabilities_to_populations(
            probabilities,
            total,
            2,
            2,
            4,
        )
        self.assertAlmostEqual(valid, 1.0)
        self.assertTrue(np.allclose(actual, classical_stream(f), atol=1e-10))


if __name__ == "__main__":
    unittest.main()
