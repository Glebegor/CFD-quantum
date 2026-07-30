"""Bit-accurate encoding helpers for the |x, y, direction> registers."""


def EncodeIndex(
    x,
    y,
    direction,
    reg_x_size,
    reg_y_size,
    direction_qubits=4,
):
    """Encode ``|x,y,d>`` into Qiskit's little-endian basis index."""
    max_x = 1 << reg_x_size
    max_y = 1 << reg_y_size
    max_direction = 1 << direction_qubits

    if not 0 <= x < max_x:
        raise ValueError(f"x must be in [0, {max_x - 1}]")
    if not 0 <= y < max_y:
        raise ValueError(f"y must be in [0, {max_y - 1}]")
    if not 0 <= direction < max_direction:
        raise ValueError(
            f"direction must be in [0, {max_direction - 1}]"
        )

    return (
        x
        | (y << reg_x_size)
        | (direction << (reg_x_size + reg_y_size))
    )


def DecodeIndex(
    index,
    reg_x_size,
    reg_y_size,
    direction_qubits=4,
):
    """Decode Qiskit's basis index into ``(x, y, direction)``."""
    total_qubits = reg_x_size + reg_y_size + direction_qubits
    if not 0 <= index < (1 << total_qubits):
        raise ValueError(
            f"index must be in [0, {(1 << total_qubits) - 1}]"
        )

    x_mask = (1 << reg_x_size) - 1
    y_mask = (1 << reg_y_size) - 1
    direction_mask = (1 << direction_qubits) - 1

    x = index & x_mask
    y = (index >> reg_x_size) & y_mask
    direction = (
        index >> (reg_x_size + reg_y_size)
    ) & direction_mask

    return x, y, direction


# Python-style aliases while preserving the project's public names.
encode_index = EncodeIndex
decode_index = DecodeIndex


def getCircuitDeep(circuit, decompose_reps=0):
    """Return the depth of a Qiskit circuit or custom gate.

    Parameters
    ----------
    circuit : qiskit.QuantumCircuit or qiskit.circuit.Gate
        Circuit whose depth should be calculated. For a custom gate, its
        definition is used automatically.
    decompose_reps : int, optional
        Number of decomposition passes applied before calculating depth.
        Zero returns the logical circuit depth.
    """
    if not isinstance(decompose_reps, int) or decompose_reps < 0:
        raise ValueError("decompose_reps must be a non-negative integer")

    target = circuit
    if not hasattr(target, "depth"):
        target = getattr(circuit, "definition", None)

    if target is None or not hasattr(target, "depth"):
        raise TypeError("circuit must be a Qiskit QuantumCircuit or Gate")

    if decompose_reps:
        target = target.decompose(reps=decompose_reps)

    return target.depth()


# Correctly-spelled Python-style alias; keep the requested public name.
getCircuitDepth = getCircuitDeep
get_circuit_depth = getCircuitDeep
