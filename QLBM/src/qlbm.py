"""Core QLBM register preparation and streaming smoke test."""

import csv
import json
import numpy as np
import qiskit as qk
from qiskit.quantum_info import Statevector
import qiskit_ibm_runtime as qk_runtime
from qiskit_aer import AerSimulator
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_ibm_runtime import EstimatorV2 as Estimator, SamplerV2 as Sampler
import gates
from animation import create_time_series_animation
from functions import DecodeIndex, EncodeIndex
from macrocomputations import reconstruct_fields
from pathlib import Path
from time import perf_counter


class QLBM:
    """QLBM circuit with register order [x][y][direction]."""

    def __init__(self, params):
        self.params = params
        expected_cells = 1 << params.REG_X_SIZE
        if params.REG_X_SIZE != params.REG_Y_SIZE:
            raise ValueError("This configuration requires equal x and y registers")
        if params.CELLS_SIZE != expected_cells:
            raise ValueError(
                "CELLS_SIZE must equal 2**REG_X_SIZE "
                f"({expected_cells}), got {params.CELLS_SIZE}"
            )
        if (1 << params.REG_D_SIZE) < params.D_NUMBER:
            raise ValueError("REG_D_SIZE cannot encode every D2Q9 direction")

        self.qreg_x = qk.QuantumRegister(params.REG_X_SIZE, "x")
        self.qreg_y = qk.QuantumRegister(params.REG_Y_SIZE, "y")
        self.qreg_d = qk.QuantumRegister(params.REG_D_SIZE, "direction")


        # Pipeline parameters
        self.circuit = None
        self.isa_circuit = None
        self.backend = None
        self.job = None
        self.result = None
        self.counts = None
        self.input_state = None
        self.expected_state = None

    def PrepareBasisState(self, x, y, d):
        """Prepare ``|x,y,d>`` with qubit zero as each register's LSB."""
        max_x = 1 << self.params.REG_X_SIZE
        max_y = 1 << self.params.REG_Y_SIZE
        max_d = 1 << self.params.REG_D_SIZE

        if not 0 <= x < max_x:
            raise ValueError(f"x must be in [0, {max_x - 1}]")
        if not 0 <= y < max_y:
            raise ValueError(f"y must be in [0, {max_y - 1}]")
        if not 0 <= d < max_d:
            raise ValueError(f"d must be in [0, {max_d - 1}]")

        circuit = qk.QuantumCircuit(
            self.qreg_x,
            self.qreg_y,
            self.qreg_d,
        )

        for bit in range(self.params.REG_X_SIZE):
            if (x >> bit) & 1:
                circuit.x(self.qreg_x[bit])

        for bit in range(self.params.REG_Y_SIZE):
            if (y >> bit) & 1:
                circuit.x(self.qreg_y[bit])

        for bit in range(self.params.REG_D_SIZE):
            if (d >> bit) & 1:
                circuit.x(self.qreg_d[bit])

        return circuit

    # Python-style alias while preserving the existing public method.
    prepare_basis_state = PrepareBasisState

    def prepareCircuit(self, input_state=None):
        print("START OF QLBM")
        if input_state is None:
            input_state = (
                getattr(self.params, "INPUT_X", 1),
                getattr(self.params, "INPUT_Y", 1),
                getattr(self.params, "INPUT_D", 8),
            )
        input_x, input_y, input_d = input_state
        if not 0 <= input_d < self.params.D_NUMBER:
            raise ValueError(
                f"INPUT_D must be a valid D2Q9 direction (0-{self.params.D_NUMBER - 1})"
            )

        self.input_state = (input_x, input_y, input_d)
        dx, dy = self.params.D[input_d]
        self.expected_state = (
            (input_x + dx) % (1 << self.params.REG_X_SIZE),
            (input_y + dy) % (1 << self.params.REG_Y_SIZE),
            input_d,
        )
        self.circuit = self.PrepareBasisState(input_x,input_y,input_d)

        if getattr(self.params, "SPECIALIZE_BASIS_STREAM", False):
            streaming_gate = gates.BasisStreamG(
                self.params.REG_X_SIZE,
                self.params.REG_Y_SIZE,
                input_d,
                self.params.D,
            )
            self.circuit.append(
                streaming_gate,
                [*self.qreg_x, *self.qreg_y],
            )
            stream_kind = "specialized basis-state"
        else:
            streaming_gate = gates.StreamG(
                self.params.REG_X_SIZE,
                self.params.REG_Y_SIZE,
                self.params.REG_D_SIZE,
            )
            self.circuit.append(streaming_gate, self.circuit.qubits)
            stream_kind = "general superposition"

        print(f"1. Circuit prepared ({stream_kind} streaming).")
        print(f"   Input: {self.input_state}; expected: {self.expected_state}")

    def circuitOptimization(self):
        execution_circuit = self.circuit.copy()

        # Sampler returns bit-string counts and therefore requires explicit
        # measurements and an output classical register. Keep self.circuit
        # measurement-free so it remains usable by Statevector tests.
        if self.params.MODE is False:
            execution_circuit.measure_all()

        pass_manager = generate_preset_pass_manager(
            backend=self.backend,
            optimization_level=getattr(self.params, "OPTIMIZATION_LEVEL", 3),
        )
        self.isa_circuit = pass_manager.run(execution_circuit)
        print("3. Circuit optimized.")
        print(f"   ISA qubits: {self.isa_circuit.num_qubits}")
        print(f"   ISA depth: {self.isa_circuit.depth()}")
        print(f"   ISA operations: {dict(self.isa_circuit.count_ops())}")

    def setupBackend(self, service: qk_runtime.QiskitRuntimeService):
        if self.params.TYPE == "prod":
            if service is None:
                raise ValueError("IBM Runtime service is required in prod mode")
            self.backend = service.least_busy(
                simulator=False,
                operational=True,
                min_num_qubits=(
                    self.params.REG_X_SIZE
                    + self.params.REG_Y_SIZE
                    + self.params.REG_D_SIZE
                ),
            )
        else:
            self.backend = AerSimulator()

        print(f"2. Backend prepared: {self.backend.name}")

    def execute(self):
        # Sampler
        if self.params.MODE is False:
            mode = Sampler(self.backend)
        # Estimator
        else:
            raise NotImplementedError(
                "Estimator mode requires an observable; use MODE=False "
                "for measured streaming results."
            )

        self.job = mode.run(
            [self.isa_circuit],
            shots=self.params.SHOTS
        )
        self.result = self.job.result()
        self.counts = self.result[0].data.meas.get_counts()
        print(f"4. Execution complete. Counts: {self.counts}")
        self.printResultSummary()

    def printResultSummary(self):
        """Print QPU correctness and invalid-direction diagnostics."""
        if not self.counts:
            raise RuntimeError("No Sampler counts available. Run execute() first.")

        total = sum(self.counts.values())
        expected_count = 0
        valid_direction_count = 0
        for bitstring, count in self.counts.items():
            state = DecodeIndex(
                int(bitstring.replace(" ", ""), 2),
                self.params.REG_X_SIZE,
                self.params.REG_Y_SIZE,
                self.params.REG_D_SIZE,
            )
            if state[2] < self.params.D_NUMBER:
                valid_direction_count += count
            if state == self.expected_state:
                expected_count += count

        print(
            "   Expected-state probability: "
            f"{expected_count / total:.3%} ({expected_count}/{total})"
        )
        print(
            "   Valid D2Q9-direction probability: "
            f"{valid_direction_count / total:.3%} "
            f"({valid_direction_count}/{total})"
        )


    def showResults(
        self,
        output_path="../../Images/qlbm_results.png",
        show=True,
        step=None,
        time_seconds=None,
        wall_time_seconds=None,
    ):
        """Decode Sampler counts and plot the measured QLBM state.

        The figure contains spatial probability, direction probabilities,
        the mean D2Q9 velocity at measured cells, and the most frequent
        computational-basis states.
        """
        if not self.counts:
            raise RuntimeError("No Sampler counts available. Run execute() first.")

        import matplotlib.pyplot as plt

        nx = 1 << self.params.REG_X_SIZE
        ny = 1 << self.params.REG_Y_SIZE
        total_shots = sum(self.counts.values())

        spatial_probability = np.zeros((ny, nx), dtype=float)
        direction_probability = np.zeros(1 << self.params.REG_D_SIZE, dtype=float)
        populations = np.zeros((self.params.D_NUMBER, ny, nx), dtype=float)

        decoded_rows = []
        for bitstring, count in self.counts.items():
            clean_bits = bitstring.replace(" ", "")
            basis_index = int(clean_bits, 2)
            x, y, direction = DecodeIndex(
                basis_index,
                self.params.REG_X_SIZE,
                self.params.REG_Y_SIZE,
                self.params.REG_D_SIZE,
            )
            probability = count / total_shots
            spatial_probability[y, x] += probability
            direction_probability[direction] += probability

            if direction < self.params.D_NUMBER:
                populations[direction, y, x] += probability

            decoded_rows.append(
                (bitstring, count, probability, x, y, direction)
            )

        density, velocity, pressure = reconstruct_fields(populations)
        occupied = density > 0.0

        figure, axes = plt.subplots(2, 2, figsize=(15, 11))

        heatmap = axes[0, 0].imshow(
            spatial_probability,
            origin="lower",
            cmap="viridis",
            interpolation="nearest",
        )
        axes[0, 0].set_title("Spatial measurement probability")
        axes[0, 0].set_xlabel("x")
        axes[0, 0].set_ylabel("y")
        figure.colorbar(heatmap, ax=axes[0, 0], label="probability")

        direction_indices = np.arange(len(direction_probability))
        colors = [
            "tab:blue" if direction < self.params.D_NUMBER else "tab:gray"
            for direction in direction_indices
        ]
        axes[0, 1].bar(
            direction_indices,
            direction_probability,
            color=colors,
        )
        axes[0, 1].set_title("Direction probability")
        axes[0, 1].set_xlabel("D2Q9 direction")
        axes[0, 1].set_ylabel("probability")
        axes[0, 1].set_xticks(direction_indices)

        measured_y, measured_x = np.nonzero(occupied)
        axes[1, 0].quiver(
            measured_x,
            measured_y,
            velocity[..., 0][occupied],
            velocity[..., 1][occupied],
            density[occupied],
            cmap="plasma",
            angles="xy",
            scale_units="xy",
            scale=1,
        )
        axes[1, 0].set_title("Mean measured D2Q9 velocity")
        axes[1, 0].set_xlabel("x")
        axes[1, 0].set_ylabel("y")
        axes[1, 0].set_xlim(-0.5, nx - 0.5)
        axes[1, 0].set_ylim(-0.5, ny - 0.5)
        axes[1, 0].set_aspect("equal")
        axes[1, 0].grid(alpha=0.25)

        top_rows = sorted(
            decoded_rows,
            key=lambda row: row[1],
            reverse=True,
        )[:16]
        labels = [
            f"({x},{y},d={direction})"
            for _, _, _, x, y, direction in top_rows
        ]
        probabilities = [row[2] for row in top_rows]
        axes[1, 1].barh(labels[::-1], probabilities[::-1])
        axes[1, 1].set_title("Most frequent decoded states")
        axes[1, 1].set_xlabel("probability")

        time_label = (
            "" if time_seconds is None else f" — t={time_seconds:.3f} s"
        )
        figure.suptitle(
            f"QLBM Sampler results — {total_shots} shots{time_label}",
            fontsize=15,
        )
        figure.tight_layout()

        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=160, bbox_inches="tight")

        expected_count = sum(
            count
            for _, count, _, x, y, direction in decoded_rows
            if (x, y, direction) == self.expected_state
        )
        correct_direction_count = sum(
            count
            for _, count, _, _, _, direction in decoded_rows
            if direction == self.expected_state[2]
        )
        valid_direction_count = sum(
            count
            for _, count, _, _, _, direction in decoded_rows
            if direction < self.params.D_NUMBER
        )

        data_path = output.with_suffix(".npz")
        np.savez_compressed(
            data_path,
            f=populations,
            density=density,
            velocity=velocity,
            pressure=pressure,
            spatial_probability=spatial_probability,
            direction_probability=direction_probability,
            bitstrings=np.asarray(list(self.counts), dtype=str),
            counts=np.asarray(list(self.counts.values()), dtype=np.int64),
            input_state=np.asarray(self.input_state, dtype=np.int64),
            expected_state=np.asarray(self.expected_state, dtype=np.int64),
            step=-1 if step is None else step,
            time_seconds=np.nan if time_seconds is None else time_seconds,
            wall_time_seconds=(
                np.nan if wall_time_seconds is None else wall_time_seconds
            ),
        )

        metadata = {
            "schema": "qlbm-d2q9-snapshot-v1",
            "step": step,
            "time_seconds": time_seconds,
            "wall_time_seconds": wall_time_seconds,
            "grid": [ny, nx],
            "measurement_shots": total_shots,
            "input_state": list(self.input_state),
            "expected_state": list(self.expected_state),
            "expected_state_probability": expected_count / total_shots,
            "correct_direction_probability": (
                correct_direction_count / total_shots
            ),
            "valid_direction_probability": (
                valid_direction_count / total_shots
            ),
            "invalid_direction_probability": (
                1.0 - valid_direction_count / total_shots
            ),
            "population_normalization": "measurement_probability",
            "pressure_model": "isothermal_lattice_p=rho/3",
            "physical_lbm_warning": (
                "This snapshot represents measurements of one encoded "
                "basis-state stream, not a complete domain-wide LBM ensemble."
            ),
            "files": {
                "png": output.name,
                "npz": data_path.name,
                "csv": output.with_suffix(".csv").name,
            },
        }
        metadata_path = output.with_suffix(".json")
        metadata_path.write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )

        csv_path = output.with_suffix(".csv")
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["bitstring", "count", "probability", "x", "y", "direction"]
            )
            writer.writerows(decoded_rows)

        print(f"5. Snapshot plot saved to: {output}")
        print(f"   Numerical fields saved to: {data_path}")
        print(f"   Metadata saved to: {metadata_path}")
        print(f"   Decoded counts saved to: {csv_path}")

        if show:
            plt.show()
        else:
            plt.close(figure)

        return {
            "spatial_probability": spatial_probability,
            "direction_probability": direction_probability,
            "f": populations,
            "density": density,
            "velocity": velocity,
            "pressure": pressure,
            "decoded_states": decoded_rows,
            "figure": figure,
            "metadata": metadata,
        }

    show_results = showResults

    def Run(self, service: qk_runtime.QiskitRuntimeService = None):
        self.prepareCircuit()
        self.setupBackend(service)
        self.circuitOptimization()
        self.execute()
        if getattr(self.params, "SHOW_RESULTS", True):
            self.showResults(
                output_path=getattr(
                    self.params,
                    "RESULT_PLOT",
                    "output/qlbm_results.png",
                ),
                show=True,
            )

    def RunTimeSeries(self, service: qk_runtime.QiskitRuntimeService = None):
        """Run and persist one measured streaming snapshot per physical step."""
        steps = int(
            getattr(
                self.params,
                "SNAPSHOT_COUNT",
                getattr(self.params, "SHOTS_NUMBER", 40),
            )
        )
        dt = float(getattr(self.params, "TIME_STEP", 0.1))
        if steps <= 0 or dt <= 0.0:
            raise ValueError("SHOTS_NUMBER and TIME_STEP must be positive")

        output_dir = Path(
            getattr(self.params, "TIME_SERIES_DIR", "output/time_series")
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        self.setupBackend(service)
        current_state = (
            getattr(self.params, "INPUT_X", 1),
            getattr(self.params, "INPUT_Y", 1),
            getattr(self.params, "INPUT_D", 8),
        )
        snapshots = []
        series_started = perf_counter()

        for step in range(1, steps + 1):
            snapshot_started = perf_counter()
            time_seconds = step * dt
            print(f"\n=== Snapshot {step}/{steps}, t={time_seconds:.3f} s ===")
            self.prepareCircuit(current_state)
            self.circuitOptimization()
            self.execute()
            wall_time_seconds = perf_counter() - snapshot_started
            stem = output_dir / f"step_{step:04d}"
            result = self.showResults(
                output_path=stem.with_suffix(".png"),
                show=False,
                step=step,
                time_seconds=time_seconds,
                wall_time_seconds=wall_time_seconds,
            )
            snapshots.append(result["metadata"])
            current_state = self.expected_state

        total_wall_time_seconds = perf_counter() - series_started
        manifest = {
            "schema": "qlbm-d2q9-time-series-v1",
            "time_step_seconds": dt,
            "snapshot_count": steps,
            "duration_seconds": steps * dt,
            "total_wall_time_seconds": total_wall_time_seconds,
            "mean_wall_time_per_snapshot_seconds": (
                total_wall_time_seconds / steps
            ),
            "grid": [
                1 << self.params.REG_Y_SIZE,
                1 << self.params.REG_X_SIZE,
            ],
            "snapshots": snapshots,
        }
        manifest_path = output_dir / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"\nTime-series manifest saved to: {manifest_path}")
        gif_path, mp4_path = create_time_series_animation(output_dir, dt)
        print(f"GIF animation saved to: {gif_path}")
        if mp4_path:
            print(f"MP4 animation saved to: {mp4_path}")
        return manifest





def TestQLBM(input_x, input_y, input_d):
    """Run and assert one streaming basis-state transition."""
    import params_test

    model = QLBM(params_test)
    nx = 1 << params_test.REG_X_SIZE
    ny = 1 << params_test.REG_Y_SIZE

    circuit = model.PrepareBasisState(input_x, input_y, input_d)
    streaming_gate = gates.StreamG(
        params_test.REG_X_SIZE,
        params_test.REG_Y_SIZE,
        params_test.REG_D_SIZE,
    )
    circuit.append(streaming_gate, circuit.qubits)

    statevector = Statevector.from_instruction(circuit)
    probabilities = statevector.probabilities()
    actual_index = int(np.argmax(probabilities))
    actual_probability = float(probabilities[actual_index])
    actual = DecodeIndex(
        actual_index,
        params_test.REG_X_SIZE,
        params_test.REG_Y_SIZE,
        params_test.REG_D_SIZE,
    )

    if input_d < params_test.D_NUMBER:
        cx, cy = params_test.D[input_d]
        expected = (
            (input_x + cx) % nx,
            (input_y + cy) % ny,
            input_d,
        )
    else:
        expected = (input_x, input_y, input_d)

    print("Input:      ", (input_x, input_y, input_d))
    print("Expected:   ", expected)
    print("Actual:     ", actual)
    print("Probability:", actual_probability)

    assert actual == expected
    assert np.isclose(actual_probability, 1.0, atol=1e-10)
    assert np.isclose(np.linalg.norm(statevector.data), 1.0, atol=1e-10)

    # The resulting index must use the same [x][y][direction] layout.
    assert actual_index == EncodeIndex(
        *expected,
        params_test.REG_X_SIZE,
        params_test.REG_Y_SIZE,
        params_test.REG_D_SIZE,
    )

    return actual, actual_probability


if __name__ == "__main__":
    TestQLBM(1, 2, 1)
