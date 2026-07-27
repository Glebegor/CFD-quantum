#pragma once

#include "Mesh.hpp"
#include "Parameters.hpp"
#include "Streaming.hpp"
#include "Collision.hpp"
#include "D2Q9.hpp"

/*
    D2Q9 lattice-Boltzmann solver (BGK collision, ISLBM streaming).

    Owns the two distribution buffers (structure-of-arrays) and the
    collision / streaming operators. The mesh and configuration are
    referenced, not owned.
*/
class LBMSolver
{

public:
  LBMSolver(Mesh &mesh, const Parameters &params);

  // Seed the domain with the equilibrium of the inlet flow.
  void initialize();

  // Advance one lattice-Boltzmann time step:
  //   collide -> stream -> obstacle bounce-back -> domain conditions.
  void step();

  // Read-only access to the current distributions (for output).
  [[nodiscard]] const LBMConstants::Distributions &distributions() const
  {
    return f_;
  }

  [[nodiscard]] double relaxationFrequency() const
  {
    return collision_.relaxationFrequency();
  }

private:
  Mesh &mesh_;
  const Parameters &params_;

  // Double buffer for streaming.
  LBMConstants::Distributions f_;
  LBMConstants::Distributions fNext_;

  BGKCollision collision_;
  ISLBMStreaming streaming_;
};
