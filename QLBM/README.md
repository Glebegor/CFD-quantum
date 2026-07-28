# Shared specification


grid: 512x512
lattice: d2q9
flow: left to right
Geometry: wing
Boundary conditions, periodic first, then inlet/ outlet and bounce -back
test grid: 4 x 4
test qubits: 2qx, 2qy, 4qd
primary data: f[direction, y, x]
Outputs: density, velocity, pressure

6  2  5
 \ | /
3  0  1
 / | \
7  4  8

velocity vectors:
C = [
    ( 0,  0),
    ( 1,  0),
    ( 0,  1),
    (-1,  0),
    ( 0, -1),
    ( 1,  1),
    (-1,  1),
    (-1, -1),
    ( 1, -1),
]


interfaces

f.shape == (9, ny, nx)
density.shape == (ny, nx)
velocity.shape == (ny, nx, 2)
pressure.shape == (ny, nx)
solid_mask.shape == (ny, nx)



registers:
q[0:2] = x
q[2:4] = y
q[4:8] = direction


functions conversion
encode_index(x, y, direction)
decode_index(index)
prepare_basis_state(x, y, direction)
