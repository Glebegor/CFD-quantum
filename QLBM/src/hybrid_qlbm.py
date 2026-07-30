"""Hybrid D2Q9 solver: classical BGK collision and quantum streaming."""

import argparse
import json
import shutil
from pathlib import Path
from time import perf_counter

import numpy as np
import qiskit as qk
from qiskit.circuit.library import StatePreparation
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

import gates
from animation import create_time_series_animation
from functions import DecodeIndex, EncodeIndex
from lbm_solver import (
    apply_channel_boundaries,
    equilibrium,
    physical_fields,
    save_snapshot,
)
from macrocomputations import reconstruct_fields


def populations_to_amplitudes(f, reg_x_size, reg_y_size, reg_d_size=4):
    """Amplitude-encode non-negative D2Q9 populations."""
    populations = np.asarray(f, dtype=float)
    ny = 1 << reg_y_size
    nx = 1 << reg_x_size
    if populations.shape != (9, ny, nx):
        raise ValueError(f"f must have shape (9, {ny}, {nx})")
    if np.any(populations < -1e-14):
        raise ValueError("amplitude encoding requires non-negative populations")

    populations = np.clip(populations, 0.0, None)
    total_population = float(populations.sum())
    if total_population <= 0.0:
        raise ValueError("total population must be positive")

    amplitudes = np.zeros(
        1 << (reg_x_size + reg_y_size + reg_d_size),
        dtype=complex,
    )
    for direction in range(9):
        for y in range(ny):
            for x in range(nx):
                index = EncodeIndex(
                    x,
                    y,
                    direction,
                    reg_x_size,
                    reg_y_size,
                    reg_d_size,
                )
                amplitudes[index] = np.sqrt(
                    populations[direction, y, x] / total_population
                )
    return amplitudes, total_population


def probabilities_to_populations(
    probabilities,
    total_population,
    reg_x_size,
    reg_y_size,
    reg_d_size=4,
):
    """Decode probabilities, postselect valid D2Q9 directions and restore mass."""
    probabilities = np.asarray(probabilities, dtype=float)
    nx = 1 << reg_x_size
    ny = 1 << reg_y_size
    f = np.zeros((9, ny, nx), dtype=float)
    valid_probability = 0.0
    for index, probability in enumerate(probabilities):
        if probability <= 0.0:
            continue
        x, y, direction = DecodeIndex(
            index,
            reg_x_size,
            reg_y_size,
            reg_d_size,
        )
        if direction < 9:
            f[direction, y, x] += probability
            valid_probability += probability
    if valid_probability <= 0.0:
        raise RuntimeError("no valid D2Q9 measurements were observed")
    f *= total_population / valid_probability
    return f, valid_probability


def classical_stream(f):
    """Reference periodic D2Q9 streaming used for verification."""
    from macrocomputations import CX, CY

    streamed = np.empty_like(f)
    for direction, (dx, dy) in enumerate(zip(CX.astype(int), CY.astype(int))):
        streamed[direction] = np.roll(
            f[direction],
            shift=(dy, dx),
            axis=(0, 1),
        )
    return streamed


class HybridQLBM:
    """Classical collision/boundaries with amplitude-encoded quantum streaming."""

    def __init__(self, params, qpu=False, service=None):
        self.params = params
        self.qpu = qpu
        self.service = service
        self.cells = int(getattr(params, "HYBRID_CELLS_SIZE", 4))
        qubits = int(np.log2(self.cells))
        if 1 << qubits != self.cells:
            raise ValueError("HYBRID_CELLS_SIZE must be a power of two")
        self.reg_x_size = qubits
        self.reg_y_size = qubits
        self.reg_d_size = 4
        self.backend = None
        self.last_job_id = None
        self.last_circuit_metrics = None

    def setup_backend(self):
        if not self.qpu:
            return
        if self.service is None:
            raise ValueError("IBM Runtime service is required for QPU mode")
        self.backend = self.service.least_busy(
            simulator=False,
            operational=True,
            min_num_qubits=(
                self.reg_x_size + self.reg_y_size + self.reg_d_size
            ),
        )
        print(f"Hybrid QPU backend: {self.backend.name}")

    def build_streaming_circuit(self, f, measure):
        amplitudes, total_population = populations_to_amplitudes(
            f,
            self.reg_x_size,
            self.reg_y_size,
            self.reg_d_size,
        )
        circuit = qk.QuantumCircuit(
            self.reg_x_size + self.reg_y_size + self.reg_d_size
        )
        circuit.append(StatePreparation(amplitudes), circuit.qubits)
        circuit.append(
            gates.StreamG(
                self.reg_x_size,
                self.reg_y_size,
                self.reg_d_size,
            ),
            circuit.qubits,
        )
        if measure:
            circuit.measure_all()
        return circuit, total_population

    def quantum_stream(self, f, shots):
        if self.qpu:
            circuit, total_population = self.build_streaming_circuit(
                f,
                measure=True,
            )
            pass_manager = generate_preset_pass_manager(
                backend=self.backend,
                optimization_level=3,
            )
            isa_circuit = pass_manager.run(circuit)
            self.last_circuit_metrics = {
                "depth": isa_circuit.depth(),
                "operations": dict(isa_circuit.count_ops()),
                "physical_qubits": isa_circuit.num_qubits,
            }
            sampler = Sampler(self.backend)
            job = sampler.run([isa_circuit], shots=shots)
            self.last_job_id = job.job_id()
            print(
                f"Hybrid job {self.last_job_id}: "
                f"depth={self.last_circuit_metrics['depth']}"
            )
            counts = job.result()[0].data.meas.get_counts()
            probabilities = np.zeros(1 << circuit.num_qubits, dtype=float)
            for bitstring, count in counts.items():
                probabilities[int(bitstring.replace(" ", ""), 2)] += (
                    count / shots
                )
        else:
            circuit, total_population = self.build_streaming_circuit(
                f,
                measure=False,
            )
            simulator = AerSimulator(method="statevector")
            circuit.save_statevector()
            compiled = qk.transpile(
                circuit,
                simulator,
                optimization_level=3,
            )
            result = simulator.run(compiled).result()
            statevector = result.get_statevector(compiled)
            probabilities = statevector.probabilities()
            self.last_circuit_metrics = {
                "depth": circuit.decompose(reps=5).depth(),
                "logical_qubits": circuit.num_qubits,
            }

        streamed, valid_probability = probabilities_to_populations(
            probabilities,
            total_population,
            self.reg_x_size,
            self.reg_y_size,
            self.reg_d_size,
        )
        return streamed, valid_probability

    def collide(self, f):
        density, velocity, _ = reconstruct_fields(f)
        f_eq = equilibrium(density, velocity)
        omega = 1.0 / float(
            getattr(self.params, "HYBRID_RELAXATION_TIME", 0.8)
        )
        collided = f - omega * (f - f_eq)
        negative_mass = float(-collided[collided < 0.0].sum())
        collided = np.clip(collided, 0.0, None)
        return collided, negative_mass

    def apply_boundaries(self, f):
        lattice_velocity = float(
            getattr(self.params, "HYBRID_LATTICE_INLET_VELOCITY", 0.03)
        )
        return apply_channel_boundaries(f, lattice_velocity)

    def run(self):
        self.setup_backend()
        snapshots = int(
            getattr(
                self.params,
                "HYBRID_QPU_SNAPSHOT_COUNT" if self.qpu else "HYBRID_SNAPSHOT_COUNT",
                1 if self.qpu else 40,
            )
        )
        shots = int(getattr(self.params, "HYBRID_SHOTS", 32768))
        output_dt = float(getattr(self.params, "TIME_STEP", 0.1))
        output_dir = Path(
            getattr(
                self.params,
                "HYBRID_QPU_OUTPUT_DIR" if self.qpu else "HYBRID_AER_OUTPUT_DIR",
                "output/hybrid_qpu" if self.qpu else "output/hybrid_aer",
            )
        )
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        density = np.ones((self.cells, self.cells), dtype=float)
        velocity = np.zeros((self.cells, self.cells, 2), dtype=float)
        f = equilibrium(density, velocity)
        f = self.apply_boundaries(f)

        velocity_scale = (
            float(getattr(self.params, "PHYSICAL_INLET_VELOCITY_M_S", 1.0))
            / float(
                getattr(self.params, "HYBRID_LATTICE_INLET_VELOCITY", 0.03)
            )
        )
        reference_density = float(
            getattr(self.params, "REFERENCE_DENSITY_KG_M3", 1.225)
        )
        reference_pressure = float(
            getattr(self.params, "REFERENCE_PRESSURE_PA", 101325.0)
        )
        history = []
        started = perf_counter()

        for step in range(1, snapshots + 1):
            collided, negative_mass = self.collide(f)
            streamed, valid_probability = self.quantum_stream(
                collided,
                shots,
            )
            f = self.apply_boundaries(streamed)
            fields = physical_fields(
                f,
                reference_density,
                reference_pressure,
                velocity_scale,
            )
            plotted_density = fields["density"][:, 1:]
            plotted_pressure = fields["pressure"][:, 1:]
            plotted_speed = np.linalg.norm(fields["velocity"][:, 1:], axis=-1)
            time_seconds = step * output_dt
            metadata = {
                "schema": "hybrid-qlbm-snapshot-v1",
                "engine": "IBM QPU" if self.qpu else "Aer statevector",
                "snapshot": step,
                "time_seconds": time_seconds,
                "grid": [self.cells, self.cells],
                "boundaries": {
                    "left": "fixed velocity inlet",
                    "right": "zero-gradient outlet",
                    "top_bottom": "periodic",
                },
                "plot_excludes_x_cells": [0],
                "measurement_shots": shots if self.qpu else None,
                "valid_direction_probability": valid_probability,
                "postselection_probability": valid_probability,
                "negative_population_mass_clipped": negative_mass,
                "job_id": self.last_job_id,
                "circuit": self.last_circuit_metrics,
                "mean_density_kg_m3": float(plotted_density.mean()),
                "max_velocity_m_s": float(plotted_speed.max()),
                "pressure_range_pa": float(
                    plotted_pressure.max() - plotted_pressure.min()
                ),
            }
            history.append(metadata)
            save_snapshot(
                output_dir,
                step,
                time_seconds,
                f,
                fields,
                metadata,
                history,
            )
            print(
                f"Hybrid {'QPU' if self.qpu else 'Aer'} "
                f"{step}/{snapshots}: valid={valid_probability:.3%}, "
                f"max|u|={plotted_speed.max():.4f} m/s"
            )

        manifest = {
            "schema": "hybrid-qlbm-series-v1",
            "engine": "IBM QPU" if self.qpu else "Aer statevector",
            "grid": [self.cells, self.cells],
            "snapshot_count": snapshots,
            "time_step_seconds": output_dt,
            "duration_seconds": snapshots * output_dt,
            "wall_time_seconds": perf_counter() - started,
            "snapshots": history,
        }
        manifest_path = output_dir / "manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
        gif, mp4 = create_time_series_animation(output_dir, output_dt)
        print(f"Hybrid manifest: {manifest_path}")
        print(f"Hybrid GIF: {gif}")
        if mp4:
            print(f"Hybrid MP4: {mp4}")
        return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qpu", action="store_true")
    arguments = parser.parse_args()

    import params

    service = None
    if arguments.qpu:
        import credentials

        service = QiskitRuntimeService(
            token=credentials.TOKEN,
            instance=credentials.CRN,
            channel="ibm_quantum_platform",
        )
    HybridQLBM(params, qpu=arguments.qpu, service=service).run()


if __name__ == "__main__":
    main()
