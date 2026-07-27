#include "VTKWriter.hpp"

#include <fstream>
#include <stdexcept>
#include <cmath>

void VTKWriter::write(
    const std::string &filename,
    const Mesh &mesh,
    const LBMSolver &solver
)
{
    std::ofstream file(filename);
    if (!file)
    {
        throw std::runtime_error(
            "Cannot open VTK file");
    }

    auto &f =
        solver.distributions();

    const size_t N =
        mesh.nodes.size();

    /*
        VTK legacy format.

        POINT_DATA obsahuje:

        velocity
        density
        speed
        solid

    */

    file
        << "# vtk DataFile Version 3.0\n"
        << "ISLBM cylinder flow\n"
        << "ASCII\n\n";

    file
        << "DATASET POLYDATA\n";

    file
        << "POINTS "
        << N
        << " float\n";

    for (
        const auto &node :
        mesh.nodes)
    {

        file
            << node.x
            << " "
            << node.y
            << " 0\n";
    }

    file
        << "\nPOINT_DATA "
        << N
        << "\n";

    /*
        Velocity vector
    */

    file
        << "VECTORS velocity float\n";

    for (size_t id = 0;
         id < N;
         id++)
    {

        double rho = 0.0;

        double ux = 0.0;

        double uy = 0.0;

        for (int q = 0; q < 9; q++)
        {

            double fq =
                f[q][id];

            rho += fq;

            ux +=
                fq *
                LBMConstants::cx[q];

            uy +=
                fq *
                LBMConstants::cy[q];
        }

        ux /= rho;

        uy /= rho;

        file
            << ux
            << " "
            << uy
            << " 0\n";
    }

    /*
        Density
    */

    file
        << "\nSCALARS density float 1\n";

    file
        << "LOOKUP_TABLE default\n";

    for (size_t id = 0;
         id < N;
         id++)
    {

        double rho = 0;

        for (int q = 0; q < 9; q++)
            rho += f[q][id];

        file
            << rho
            << "\n";
    }

    /*
        Solid mask
    */

    file
        << "\nSCALARS solid int 1\n";

    file
        << "LOOKUP_TABLE default\n";

    for (
        const auto &node :
        mesh.nodes)
    {

        file
            << (node.solid ? 1 : 0)
            << "\n";
    }
}