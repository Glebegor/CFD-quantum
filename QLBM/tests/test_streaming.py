"""Exhaustive basis-state tests for the 4x4 D2Q9 streaming circuit."""

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from qiskit.quantum_info import Statevector

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

import gates  # noqa: E402
from functions import DecodeIndex, EncodeIndex  # noqa: E402
from qlbm import QLBM  # noqa: E402

D2Q9 = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
]


class StreamingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Statevector tests intentionally use 4x4. A 512x512 statevector
        # circuit has 22 qubits and is not suitable for exhaustive unit tests.
        cls.test_params = SimpleNamespace(
            CELLS_SIZE=4,
            REG_X_SIZE=2,
            REG_Y_SIZE=2,
            REG_D_SIZE=4,
            D=D2Q9,
            D_NUMBER=9,
        )
        cls.model = QLBM(cls.test_params)
        cls.reg_x_size = cls.test_params.REG_X_SIZE
        cls.reg_y_size = cls.test_params.REG_Y_SIZE
        cls.reg_d_size = cls.test_params.REG_D_SIZE
        cls.nx = 1 << cls.reg_x_size
        cls.ny = 1 << cls.reg_y_size
        cls.stream_gate = gates.StreamG(
            cls.reg_x_size,
            cls.reg_y_size,
            cls.reg_d_size,
        )

    def simulate(self, x, y, direction, stream=True):
        circuit = self.model.PrepareBasisState(x, y, direction)
        if stream:
            circuit.append(self.stream_gate, circuit.qubits)

        statevector = Statevector.from_instruction(circuit)
        probabilities = statevector.probabilities()
        index = int(np.argmax(probabilities))
        probability = float(probabilities[index])
        decoded = DecodeIndex(
            index,
            self.reg_x_size,
            self.reg_y_size,
            self.reg_d_size,
        )
        return decoded, probability, statevector

    def expected_stream(self, x, y, direction):
        if direction < self.test_params.D_NUMBER:
            cx, cy = self.test_params.D[direction]
            return (
                (x + cx) % self.nx,
                (y + cy) % self.ny,
                direction,
            )
        return x, y, direction

    def assert_basis_output(self, actual, probability, statevector, expected):
        self.assertEqual(actual, expected)
        self.assertTrue(np.isclose(probability, 1.0, atol=1e-10))
        self.assertTrue(
            np.isclose(np.linalg.norm(statevector.data), 1.0, atol=1e-10)
        )

    def test_encode_decode_round_trip_all_256_states(self):
        for direction in range(16):
            for y in range(self.ny):
                for x in range(self.nx):
                    with self.subTest(x=x, y=y, direction=direction):
                        index = EncodeIndex(
                            x,
                            y,
                            direction,
                            self.reg_x_size,
                            self.reg_y_size,
                            self.reg_d_size,
                        )
                        self.assertEqual(
                            DecodeIndex(
                                index,
                                self.reg_x_size,
                                self.reg_y_size,
                                self.reg_d_size,
                            ),
                            (x, y, direction),
                        )

    def test_full_512_grid_register_structure_without_statevector(self):
        gate = gates.StreamG(9, 9, 4)
        self.assertEqual(gate.num_qubits, 22)

        representative_states = [
            (0, 0, 0),
            (1, 2, 1),
            (511, 256, 1),
            (0, 256, 3),
            (511, 511, 5),
            (0, 0, 7),
            (123, 456, 15),
        ]
        for state in representative_states:
            with self.subTest(state=state):
                index = EncodeIndex(*state, 9, 9, 4)
                self.assertEqual(DecodeIndex(index, 9, 9, 4), state)

    def test_prepare_basis_state_before_streaming(self):
        actual, probability, statevector = self.simulate(1, 2, 1, stream=False)
        self.assert_basis_output(
            actual,
            probability,
            statevector,
            (1, 2, 1),
        )

    def test_specialized_basis_stream_matches_general_stream(self):
        for direction in range(9):
            with self.subTest(direction=direction):
                circuit = self.model.PrepareBasisState(1, 2, direction)
                specialized = gates.BasisStreamG(
                    self.reg_x_size,
                    self.reg_y_size,
                    direction,
                    self.test_params.D,
                )
                circuit.append(
                    specialized,
                    [
                        *self.model.qreg_x,
                        *self.model.qreg_y,
                    ],
                )
                statevector = Statevector.from_instruction(circuit)
                index = int(np.argmax(statevector.probabilities()))
                actual = DecodeIndex(
                    index,
                    self.reg_x_size,
                    self.reg_y_size,
                    self.reg_d_size,
                )
                self.assertEqual(
                    actual,
                    self.expected_stream(1, 2, direction),
                )

    def test_all_d2q9_directions_from_middle_cell(self):
        for direction in range(9):
            with self.subTest(direction=direction):
                actual, probability, statevector = self.simulate(
                    1, 2, direction
                )
                self.assert_basis_output(
                    actual,
                    probability,
                    statevector,
                    self.expected_stream(1, 2, direction),
                )

    def test_periodic_boundaries(self):
        cases = [
            ((3, 2, 1), (0, 2, 1)),
            ((0, 2, 3), (3, 2, 3)),
            ((3, 3, 5), (0, 0, 5)),
            ((0, 0, 7), (3, 3, 7)),
        ]
        for inputs, expected in cases:
            with self.subTest(inputs=inputs):
                actual, probability, statevector = self.simulate(*inputs)
                self.assert_basis_output(
                    actual,
                    probability,
                    statevector,
                    expected,
                )

    def test_unused_directions_are_identity(self):
        for direction in range(9, 16):
            with self.subTest(direction=direction):
                expected = (1, 2, direction)
                actual, probability, statevector = self.simulate(*expected)
                self.assert_basis_output(
                    actual,
                    probability,
                    statevector,
                    expected,
                )

    def test_exhaustive_streaming_all_256_basis_states(self):
        for direction in range(16):
            for y in range(self.ny):
                for x in range(self.nx):
                    with self.subTest(x=x, y=y, direction=direction):
                        actual, probability, statevector = self.simulate(
                            x, y, direction
                        )
                        self.assert_basis_output(
                            actual,
                            probability,
                            statevector,
                            self.expected_stream(x, y, direction),
                        )


if __name__ == "__main__":
    unittest.main()
