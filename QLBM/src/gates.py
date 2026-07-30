import qiskit as qk
from qiskit.circuit.library import MCXGate, XGate


# 0 → 1
# 1 → 2
# 2 → 3
# 3 → 0
def IncrementG(num_qubits):
    '''
    |x> -> x+1 mod 2^n
    '''
    
    circuit = qk.QuantumCircuit(num_qubits, name=f"INT_{num_qubits}")
    
    for target in range(num_qubits-1,0,-1):
        controls=list(range(target))
        
        circuit.append(
            MCXGate(num_ctrl_qubits=len(controls)),
            controls + [target],
        )
    
    circuit.x(0)
    
    return circuit.to_gate()

# 3 → 2
# 2 → 1
# 1 → 0
# 0 → 3
def DecrementG(num_qubits):
    '''
    |x> -> x-1 mod 2^n
    '''
    
    circuit = qk.QuantumCircuit(num_qubits, name=f"DEC_{num_qubits}")
    
    for target in range(num_qubits-1,0,-1):
        controls=list(range(target))
        
        for control in controls:
            circuit.x(control)
            
        circuit.append(
            MCXGate(num_ctrl_qubits=len(controls)),
            controls + [target],
        )
        for control in controls:
            circuit.x(control)
                
    circuit.x(0)
    
    return circuit.to_gate()

# Change direction
def DirectionG(circuit, op, reg_d, reg_target, d_value):
    '''
    Apply operation to target_register only when:
        reg_d == d_value
    
    reg_d[0] is the least sign bit.
    '''
    
    direction_bits = len(reg_d)
    
    if not 0 <= d_value < 2**direction_bits:
        raise ValueError(
            f"d_value must be between 0 and {2**direction_bits - 1}"
        )
        
    for bit, qubit in enumerate(reg_d):
        expected_bit = (d_value >> bit) & 1
        
        if expected_bit == 0:
            circuit.x(qubit)
            
    controlled_operation = op.control(
        num_ctrl_qubits=direction_bits,
    )
    
    circuit.append(
        controlled_operation,
        [
            *reg_d,
            *reg_target
        ],
    )
    
    for bit, qubit in enumerate(reg_d):
        expected_bit = (d_value >>bit) & 1
        
        if expected_bit == 0:
            circuit.x(qubit)

# Stream gate
def StreamG(qx_n, qy_n, qd_n=4):
    """
    Create periodic D2Q9 quantum streaming
    |x,y,d> -> |x+cx[d], y+cy[d], d>
    """
    
    if qx_n <= 0 or qy_n <= 0 or qd_n <= 0:
        raise ValueError("register sizes must be positive")
    
    x_reg = qk.QuantumRegister(qx_n, "x")
    y_reg = qk.QuantumRegister(qy_n, "y")
    d_reg = qk.QuantumRegister(qd_n, "d")
    
    circuit = qk.QuantumCircuit(
        x_reg,
        y_reg,
        d_reg,
        name="D2Q9_STREAM"
    )
    
    increment_x = IncrementG(qx_n)
    increment_y = IncrementG(qy_n)
    
    decrement_x = DecrementG(qx_n)
    decrement_y = DecrementG(qy_n)
    
    # Movement
    
    
    # XX
    for direction in (1, 5, 8):
        DirectionG(
            circuit=circuit, 
            op=increment_x, 
            reg_d=d_reg,
            reg_target=x_reg,
            d_value=direction
            )
    
    for direction in (3, 6, 7):
            DirectionG(
                circuit=circuit, 
                op=decrement_x, 
                reg_d=d_reg,
                reg_target=x_reg,
                d_value=direction
                )
            
            
     # YY
    for direction in (2, 5, 6):
            DirectionG(
                circuit=circuit, 
                op=increment_y, 
                reg_d=d_reg,
                reg_target=y_reg,
                d_value=direction
                )
        
    for direction in (4, 7, 8):
            DirectionG(
                circuit=circuit, 
                op=decrement_y, 
                reg_d=d_reg,
                reg_target=y_reg,
                d_value=direction
                )
    return circuit.to_gate()


def BasisStreamG(qx_n, qy_n, direction, velocities):
    """Streaming gate specialized for one known computational-basis direction.

    This is equivalent to ``StreamG`` for a basis-state direction, but omits
    all direction comparisons and inactive arithmetic branches.  It is meant
    for small QPU smoke tests, not for a superposition of directions.
    """
    if qx_n <= 0 or qy_n <= 0:
        raise ValueError("register sizes must be positive")
    if not 0 <= direction < len(velocities):
        raise ValueError(
            f"direction must be between 0 and {len(velocities) - 1}"
        )

    x_reg = qk.QuantumRegister(qx_n, "x")
    y_reg = qk.QuantumRegister(qy_n, "y")
    circuit = qk.QuantumCircuit(x_reg, y_reg, name=f"BASIS_STREAM_D{direction}")
    dx, dy = velocities[direction]

    if dx == 1:
        circuit.append(IncrementG(qx_n), x_reg)
    elif dx == -1:
        circuit.append(DecrementG(qx_n), x_reg)
    elif dx != 0:
        raise ValueError("only D2Q9 velocities -1, 0 and 1 are supported")

    if dy == 1:
        circuit.append(IncrementG(qy_n), y_reg)
    elif dy == -1:
        circuit.append(DecrementG(qy_n), y_reg)
    elif dy != 0:
        raise ValueError("only D2Q9 velocities -1, 0 and 1 are supported")

    return circuit.to_gate()
    
    
def DirectionCircuit(d_value, direction_qubits=4):
    """Build a visual example of a gate conditioned on direction ``d``."""
    reg_d = qk.QuantumRegister(direction_qubits, "d")
    reg_target = qk.QuantumRegister(1, "target")
    circuit = qk.QuantumCircuit(
        reg_d,
        reg_target,
        name=f"DIR_{d_value}",
    )

    DirectionG(
        circuit=circuit,
        op=XGate(),
        reg_d=reg_d,
        reg_target=reg_target,
        d_value=d_value,
    )

    return circuit

def ShowAllCircuits():
    """Show increment, decrement, direction and streaming circuits."""
    import matplotlib.pyplot as plt
    import params

    register_sizes = {
        params.REG_X_SIZE,
        params.REG_Y_SIZE,
        params.REG_D_SIZE,
    }

    for num_qubits in sorted(register_sizes, reverse=True):
        increment_circuit = IncrementG(num_qubits).definition
        increment_figure = increment_circuit.draw(output="mpl", fold=-1)
        increment_figure.suptitle(f"Increment gate — {num_qubits} qubits")

        decrement_circuit = DecrementG(num_qubits).definition
        decrement_figure = decrement_circuit.draw(output="mpl", fold=-1)
        decrement_figure.suptitle(f"Decrement gate — {num_qubits} qubits")

    for d_value, (dx, dy) in enumerate(params.D):
        direction_circuit = DirectionCircuit(
            d_value=d_value,
            direction_qubits=params.REG_D_SIZE,
        )
        direction_figure = direction_circuit.draw(output="mpl", fold=-1)
        direction_figure.suptitle(
            f"Direction gate d={d_value}, velocity=({dx}, {dy})"
        )

    stream_circuit = StreamG(
        qx_n=params.REG_X_SIZE,
        qy_n=params.REG_Y_SIZE,
        qd_n=params.REG_D_SIZE,
    ).definition
    stream_figure = stream_circuit.draw(output="mpl", fold=-1)
    stream_figure.suptitle(
        "D2Q9 stream gate: |x,y,d> → |x+cx[d],y+cy[d],d>"
    )

    plt.show()


if __name__ == "__main__":
    ShowAllCircuits()
    
