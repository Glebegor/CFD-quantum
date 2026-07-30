TYPE="prod"
# If true use estimator, if false use sampler
MODE=False
SHOTS=1000
SHOW_RESULTS = True
RESULT_PLOT = "output/qlbm_results.png"

# Cells
CELLS_SIZE = 16

# Directions
# 6  2  5
#  \ | /
# 3  0  1
#  / | \
# 7  4  8
D = [
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
D_NUMBER = 9

# Registers
REG_X_SIZE = 4
REG_Y_SIZE = 4
REG_D_SIZE = 4

# QPU basis-state streaming smoke test.  Since this input is classical and
# known before execution, specialize the streaming circuit to its direction
# instead of executing all twelve expensive direction-controlled branches.
INPUT_X = 3
INPUT_Y = 3
INPUT_D = 1
SPECIALIZE_BASIS_STREAM = True
OPTIMIZATION_LEVEL = 3

# Shape
WING_X = 1.0
WING_Y = 1.0
WING_CHORD = 1.0
WING_ATTACK_ANGLE = 0.0

# Time
TIME_STEP = 0.1
SNAPSHOT_COUNT = 40
# Backward-compatible name used by older scripts.
SHOTS_NUMBER = SNAPSHOT_COUNT
TIME_SERIES_DIR = "output/time_series"

# Classical D2Q9 reference CFD model.  A stable internal LBM step is used and
# eight internal steps are grouped into every requested 0.1-second snapshot.
CFD_CELLS_SIZE = 128
PHYSICAL_CELL_SIZE_M = 0.25
PHYSICAL_INLET_VELOCITY_M_S = 1.0
LBM_INLET_VELOCITY = 0.05
LBM_RELAXATION_TIME = 0.6
REFERENCE_DENSITY_KG_M3 = 1.225
REFERENCE_PRESSURE_PA = 101325.0
CFD_TIME_SERIES_DIR = "output/cfd_series"

# Full hybrid QLBM proof of concept. Collision and boundaries are classical;
# the complete f_i field is amplitude-encoded and streamed by StreamG.
HYBRID_CELLS_SIZE = 32
HYBRID_LATTICE_INLET_VELOCITY = 0.03
HYBRID_RELAXATION_TIME = 0.8
HYBRID_SHOTS = 32768
HYBRID_SNAPSHOT_COUNT = 40
HYBRID_QPU_SNAPSHOT_COUNT = 20
HYBRID_AER_OUTPUT_DIR = "output/hybrid_aer"
HYBRID_QPU_OUTPUT_DIR = "output/hybrid_qpu"

# Constants
