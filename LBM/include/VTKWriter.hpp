#pragma once

#include "Mesh.hpp"
#include "LBM.hpp"
#include "Parameters.hpp"

#include <string>

/*
    VTK Legacy UNSTRUCTURED_GRID exporter.

    Writes points, quad-cell connectivity, cell types and the
    macroscopic fields (velocity vector, density scalar, solid mask)
    so ParaView can render the surface, contours, glyphs and the
    stream tracer. Output is guaranteed to contain no NaN/Inf values;
    numerical corruption aborts with a diagnostic.
*/
class VTKWriter
{

public:
  // Coordinates are written in metres and velocity in m/s, using the
  // physical scaling in `params`, so ParaView animates real units and
  // the Particle Tracer integrates consistently.
  static void write(
      const std::string &filename,
      const Mesh &mesh,
      const LBMSolver &solver,
      const Parameters &params);
};
