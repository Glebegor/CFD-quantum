# Makefile for the OpenFOAM settings utility
# Built for macOS + Apple Clang.

CXX ?= clang++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -I$(SRC_DIR)
LDFLAGS ?=

TARGET ?= settings_app
BUILD_DIR ?= build
SRC_DIR ?= src

SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/mesh_generator.cpp

OBJECTS := \
	$(BUILD_DIR)/main.o \
	$(BUILD_DIR)/mesh_generator.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/mesh_generator.h
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/mesh_generator.o: $(SRC_DIR)/mesh_generator.cpp $(SRC_DIR)/mesh_generator.h
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all run clean
