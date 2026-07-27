
#include "mesh_generator.h"

#include <exception>
#include <iostream>

int main()
{
    try {
        meshgen::Config config;
        meshgen::generate_mesh(config);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Mesh generation failed: " << error.what() << std::endl;
        return 1;
    }
}