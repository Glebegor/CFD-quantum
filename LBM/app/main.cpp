#include "Parameters.hpp"
#include "Mesh.hpp"
#include "LBM.hpp"
#include "VTKWriter.hpp"

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

namespace
{
  // Zero-padded step index so the VTK snapshots sort correctly as a
  // ParaView time series (output_000000.vtk, output_001000.vtk, ...).
  std::string snapshotName(const std::string &dir, int step)
  {
    std::ostringstream name;
    name << dir << "/output_"
         << std::setw(6) << std::setfill('0') << step << ".vtk";
    return name.str();
  }

  void printStartupSummary(
      const Parameters &params,
      const Mesh &mesh,
      const LBMSolver &solver)
  {
    std::cout << std::fixed << std::setprecision(6)
              << "==================================================\n"
              << " ISLBM D2Q9 Cylinder Flow Solver\n"
              << "==================================================\n"
              << " Resolution        : " << params.resolution << " x "
              << params.resolution << "\n"
              << " Mesh nodes        : " << mesh.nodeCount() << "\n"
              << " Mesh cells (quad) : " << mesh.cellCount() << "\n"
              << " Reynolds number   : " << params.reynoldsNumber << "\n"
              << " Inlet velocity    : " << params.inletVelocity
              << " (Mach " << params.machNumber() << ")\n"
              << " Char. length L    : " << params.characteristicLength()
              << " cells\n"
              << " Viscosity nu      : " << params.viscosity() << "\n"
              << " Relaxation tau    : " << params.relaxationTime() << "\n"
              << " Relaxation omega  : " << solver.relaxationFrequency() << "\n"
              << " Max steps         : " << params.maxSteps << "\n"
              << " VTK interval      : " << params.vtkInterval << "\n"
              << " Output directory  : " << params.outputDirectory << "\n"
              << "==================================================\n"
              << std::flush;
  }
}

int main()
{
  Parameters params;

  std::filesystem::create_directories(params.outputDirectory);

  // --- build mesh + interpolation topology ---
  Mesh mesh;
  mesh.generate(
      params.resolution,
      params.cylinderRadius,
      params.cylinderX,
      params.cylinderY);
  mesh.buildInterpolation();

  // --- solver ---
  LBMSolver solver(mesh, params);
  solver.initialize();

  printStartupSummary(params, mesh, solver);

  // --- time loop ---
  for (int step = 0; step <= params.maxSteps; step++)
  {
    solver.step();

    if (step % params.logInterval == 0)
      std::cout << "  step " << std::setw(7) << step
                << " / " << params.maxSteps << "\n"
                << std::flush;

    if (step % params.vtkInterval == 0)
      VTKWriter::write(snapshotName(params.outputDirectory, step), mesh, solver);
  }

  std::cout << "Simulation complete. VTK output in '"
            << params.outputDirectory << "'.\n";
  return 0;
}
