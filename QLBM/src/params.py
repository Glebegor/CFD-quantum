import math



# Cells
CELLS_SIZE = 512

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
REG_X_SIZE = 9
REG_Y_SIZE = 9
REG_D_SIZE = 4

# Shape
WING_X = 1.0
WING_Y = 1.0
WING_CHORD = 1.0
WING_ATTACK_ANGLE = 0.0

# Time
TIME_STEP = 0.1
SHOTS_NUMBER = 40

# Constants

