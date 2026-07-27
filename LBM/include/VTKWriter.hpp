#pragma once

#include "Mesh.hpp"
#include "LBM.hpp"

#include <string>

class VTKWriter
{

public:
    static void write(

        const std::string &filename,

        const Mesh &mesh,

        const LBMSolver &solver

    );
};