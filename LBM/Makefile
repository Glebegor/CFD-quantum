# ==========================================
# ISLBM Cylinder CFD Solver
# C++20 + OpenMP
# ==========================================


CXX = /opt/homebrew/bin/g++-16

TARGET = lbm

# Compiler flags
CXXFLAGS = \
	-std=c++20 \
	-O3 \
	-march=native \
	-fopenmp \
	-ffast-math \
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

	$(CXX) \
	$(CXXFLAGS) \
	$^ \
	-o $@




# Compilation

$(BUILD_DIR)/%.o: %.cpp

	mkdir -p $(dir $@)

	$(CXX) \
	$(CXXFLAGS) \
	-c $< \
	-o $@



# ==========================================
# Run simulation
# ==========================================

run: $(TARGET)

	./$(TARGET)



# ==========================================
# Clean
# ==========================================

clean:

	rm -rf $(BUILD_DIR)

	rm -f $(TARGET)



# ==========================================
# Debug build
# ==========================================

debug:

	$(MAKE) \
	CXXFLAGS="-std=c++20 -g -O0 -Wall -Wextra -Iinclude -fopenmp"
