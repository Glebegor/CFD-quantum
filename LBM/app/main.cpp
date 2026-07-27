#include "Mesh.hpp"
#include "LBM.hpp"
#include "VTKWriter.hpp"
#include <filesystem>

#include <iostream>

int main()
{

  std::filesystem::create_directories("output");

  Mesh mesh;

  mesh.generate(
      64,
      0.08);

  mesh.buildInterpolation();

  LBMSolver solver(
      mesh,
      100.0);

  solver.initialize();

  for (int step = 0;
       step < 5000;
       step++)
  {

    solver.step();

    if (step % 500 == 0)
    {

      std::cout
          << "step "
          << step
          << std::endl;

      VTKWriter::write(
          "output/output_" + std::to_string(step) + ".vtk",
          mesh,
          solver);
    }
  }
}