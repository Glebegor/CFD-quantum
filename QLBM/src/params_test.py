"""Small exhaustive-test configuration for D2Q9 streaming."""

TYPE="test"
# If true use estimator, if false use sampler
MODE=False
SHOTS=1000
SHOW_RESULTS = True
RESULT_PLOT = "output/qlbm_results_local.png"

# Cells
CELLS_SIZE = 4

# Directions
# 6  2  5
#  \ | /
# 3  0  1
#  / | \
# 7  4  8
D = [
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
D_NUMBER = 9

REG_X_SIZE = 2
REG_Y_SIZE = 2
REG_D_SIZE = 4

INPUT_X = 1
INPUT_Y = 2
INPUT_D = 1
SPECIALIZE_BASIS_STREAM = True
OPTIMIZATION_LEVEL = 3

WING_X = 1.0
WING_Y = 1.0
WING_CHORD = 1.0
WING_ATTACK_ANGLE = 0.0

TIME_STEP = 0.1
SNAPSHOT_COUNT = 40
# Backward-compatible name used by older scripts.
SHOTS_NUMBER = SNAPSHOT_COUNT
TIME_SERIES_DIR = "output/time_series"
