#include "VTKWriter.hpp"

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace
{
  [[noreturn]] void fail(const std::string &what, std::size_t id)
  {
    std::cerr << "\nVTK EXPORT ABORTED: " << what
              << " at node " << id << "\n"
              << std::flush;
    std::abort();
  }
}

void VTKWriter::write(
    const std::string &filename,
    const Mesh &mesh,
    const LBMSolver &solver)
{
  std::ofstream file(filename);
  if (!file)
    throw std::runtime_error("Cannot open VTK file: " + filename);

  const auto &f = solver.distributions();
  const std::size_t N = mesh.nodeCount();
  const std::size_t C = mesh.cellCount();

  // Full precision so ParaView reads exact values.
  file.precision(9);

  // ---------------------------------------------------------------
  // Header + geometry
  // ---------------------------------------------------------------
  file << "# vtk DataFile Version 3.0\n"
       << "ISLBM cylinder flow\n"
       << "ASCII\n"
       << "DATASET UNSTRUCTURED_GRID\n";

  file << "POINTS " << N << " float\n";
  for (const auto &node : mesh.nodes)
    file << node.x << " " << node.y << " 0\n";

  // ---------------------------------------------------------------
  // Topology: one quad (4 nodes) per cell.
  // CELLS line count = C, total integers = C * (1 + 4).
  // ---------------------------------------------------------------
  file << "\nCELLS " << C << " " << C * 5 << "\n";
  for (const auto &cell : mesh.cells)
    file << "4 " << cell[0] << " " << cell[1] << " "
         << cell[2] << " " << cell[3] << "\n";

  file << "\nCELL_TYPES " << C << "\n";
  for (std::size_t c = 0; c < C; c++)
    file << "9\n"; // VTK_QUAD

  // ---------------------------------------------------------------
  // Point data: reconstruct macroscopic fields once, validate, cache.
  // ---------------------------------------------------------------
  std::vector<double> ux(N), uy(N), rho(N);

  for (std::size_t id = 0; id < N; id++)
  {
    if (mesh.nodes[id].solid)
    {
      // Masked obstacle interior: report a quiescent, unit-density
      // state so contours/glyphs stay well defined.
      ux[id] = 0.0;
      uy[id] = 0.0;
      rho[id] = 1.0;
      continue;
    }

    double r, u, v;
    LBMConstants::computeMoments(f, id, r, u, v);

    if (!std::isfinite(r) || r <= 0.0)
      fail("invalid density", id);
    if (!std::isfinite(u) || !std::isfinite(v))
      fail("invalid velocity", id);

    rho[id] = r;
    ux[id] = u;
    uy[id] = v;
  }

  file << "\nPOINT_DATA " << N << "\n";

  // Velocity vector field.
  file << "VECTORS velocity float\n";
  for (std::size_t id = 0; id < N; id++)
    file << ux[id] << " " << uy[id] << " 0\n";

  // Density scalar field.
  file << "\nSCALARS density float 1\n"
       << "LOOKUP_TABLE default\n";
  for (std::size_t id = 0; id < N; id++)
    file << rho[id] << "\n";

  // Velocity-magnitude scalar (convenient for contours in ParaView).
  file << "\nSCALARS velocity_magnitude float 1\n"
       << "LOOKUP_TABLE default\n";
  for (std::size_t id = 0; id < N; id++)
    file << std::sqrt(ux[id] * ux[id] + uy[id] * uy[id]) << "\n";

  // Solid mask.
  file << "\nSCALARS solid int 1\n"
       << "LOOKUP_TABLE default\n";
  for (const auto &node : mesh.nodes)
    file << (node.solid ? 1 : 0) << "\n";
}
