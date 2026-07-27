# ==========================================
# ISLBM Cylinder CFD Solver
# C++20 + OpenMP
# ==========================================
#
# macOS note:
#   Apple Clang does not ship libgomp, so a real GCC is used.
#   Override the compiler if your GCC has a different suffix, e.g.
#       make CXX=g++-15
#   or point at any OpenMP-capable compiler.

CXX = /opt/homebrew/bin/g++-16

TARGET = lbm

# Compiler flags.
#   -ffast-math is deliberately NOT used: it lets the compiler assume
#   no NaN/Inf, which would silently defeat the isfinite() stability
#   guards in the solver.
CXXFLAGS = \
	-std=c++20 \
	-O3 \
	-march=native \
	-fopenmp \
	-Wall \
	-Wextra \
	-Iinclude

# Directories
SRC_DIR = src
APP_DIR = app
BUILD_DIR = build

# Source files
SOURCES = \
	$(APP_DIR)/main.cpp \
	$(SRC_DIR)/Mesh.cpp \
	$(SRC_DIR)/Boundary.cpp \
	$(SRC_DIR)/Collision.cpp \
	$(SRC_DIR)/Streaming.cpp \
	$(SRC_DIR)/LBM.cpp \
	$(SRC_DIR)/VTKWriter.cpp

# Object files
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# ==========================================
# Default target
# ==========================================
all: $(TARGET)

# Linking
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==========================================
# Run simulation
# ==========================================
run: $(TARGET)
	./$(TARGET)

# ==========================================
# Clean / rebuild
# ==========================================
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)

rebuild: clean all

# ==========================================
# Debug build (no optimisation, assertions live)
# ==========================================
debug:
	$(MAKE) CXXFLAGS="-std=c++20 -g -O0 -Wall -Wextra -fopenmp -Iinclude"

.PHONY: all run clean rebuild debug
