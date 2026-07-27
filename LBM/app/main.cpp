#include "Parameters.hpp"
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "LBM.hpp"
#include "VTKWriter.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  // Zero-padded step index so the VTK snapshots sort correctly as a
  // ParaView time series (output_000000.vtk, output_001000.vtk, ...).
  std::string snapshotName(int step)
  {
    std::ostringstream name;
    name << "output_"
         << std::setw(6) << std::setfill('0') << step << ".vtk";
    return name.str();
  }

  // ParaView collection file: loads every snapshot as one animated
  // time series with an explicit PHYSICAL time value (seconds) per
  // frame. This is what makes the velocity field play back as motion
  // (and enables the Particle Tracer / Temporal Pathlines filters).
  void writeCollection(const std::string &dir,
                       const std::vector<int> &steps,
                       double dtPhysical)
  {
    std::ofstream pvd(dir + "/simulation.pvd");
    pvd.precision(9);
    pvd << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n"
        << "  <Collection>\n";
    for (int step : steps)
      pvd << "    <DataSet timestep=\"" << step * dtPhysical
          << "\" file=\"" << snapshotName(step) << "\"/>\n";
    pvd << "  </Collection>\n"
        << "</VTKFile>\n";
  }

  void printStartupSummary(
      const Parameters &params,
      const Mesh &mesh,
      const LBMSolver &solver)
  {
    const bool airfoil = params.obstacle == Obstacle::Airfoil;
    std::cout << std::fixed << std::setprecision(6)
              << "==================================================\n"
              << " ISLBM D2Q9 Flow Solver\n"
              << "==================================================\n"
              << " Resolution        : " << params.resolution << " x "
              << params.resolution << "\n"
              << " Mesh nodes        : " << mesh.nodeCount() << "\n"
              << " Mesh cells (quad) : " << mesh.cellCount() << "\n"
              << " Obstacle          : "
              << (airfoil ? "airfoil" : "cylinder") << "\n";
    if (airfoil)
      std::cout << " Chord / AoA       : " << params.chord << " / "
                << params.angleOfAttack << " deg\n";
    std::cout
              << " Reynolds number   : " << params.reynoldsNumber << "\n"
              << " Inlet velocity    : " << params.inletVelocity
              << " (Mach " << params.machNumber() << ")\n"
              << " Char. length L    : " << params.characteristicLength()
              << " cells\n"
              << " Viscosity nu      : " << params.viscosity() << "\n"
              << " Relaxation tau    : " << params.relaxationTime() << "\n"
              << " Relaxation omega  : " << solver.relaxationFrequency() << "\n"
              << "--------------------------------------------------\n"
              << " Physical chord    : " << params.physicalChord << " m\n"
              << " Physical velocity : " << params.physicalVelocity << " m/s\n"
              << " Cell size dx      : " << params.dxPhysical() << " m\n"
              << " Timestep dt       : " << params.dtPhysical() << " s\n"
              << " Domain size       : " << params.lengthScale() << " m\n"
              << " Frame rate        : " << params.framesPerSecond << " fps\n"
              << " Duration          : " << params.durationSeconds << " s\n"
              << " Steps / frame     : " << params.snapshotStride() << "\n"
              << " Total steps       : " << params.totalSteps() << "\n"
              << " Snapshots         : "
              << (params.totalSteps() / params.snapshotStride() + 1) << "\n"
              << " Output directory  : " << params.outputDirectory << "\n"
              << "==================================================\n"
              << std::flush;
  }
}

int main()
{
  Parameters params;

  std::filesystem::create_directories(params.outputDirectory);

  // --- build obstacle geometry (cylinder, airfoil, or any imported
  //     outline; see makeGeometry / Geometry.hpp to add shapes) ---
  const Geometry geometry = makeGeometry(params);

  // --- build mesh + interpolation topology ---
  Mesh mesh;
  mesh.generate(
      params.resolution,
      [&geometry](double x, double y)
      { return geometry.isSolid(x, y); });
  mesh.buildInterpolation();

  // --- solver ---
  LBMSolver solver(mesh, params);
  solver.initialize();

  printStartupSummary(params, mesh, solver);

  // --- time loop (driven by the physical duration / frame rate) ---
  const int totalSteps = params.totalSteps();
  const int stride = params.snapshotStride();

  std::vector<int> snapshots;
  for (int step = 0; step <= totalSteps; step++)
  {
    solver.step();

    if (step % params.logInterval == 0)
      std::cout << "  step " << std::setw(7) << step
                << " / " << totalSteps
                << "   (t = " << step * params.dtPhysical() << " s)\n"
                << std::flush;

    if (step % stride == 0)
    {
      VTKWriter::write(
          params.outputDirectory + "/" + snapshotName(step),
          mesh, solver, params);
      snapshots.push_back(step);
    }
  }

  // Emit the ParaView collection describing the animated series.
  writeCollection(params.outputDirectory, snapshots, params.dtPhysical());

  std::cout << "Simulation complete. VTK output in '"
            << params.outputDirectory << "'.\n"
            << "Open '" << params.outputDirectory
            << "/simulation.pvd' in ParaView for the animation.\n";
  return 0;
}
