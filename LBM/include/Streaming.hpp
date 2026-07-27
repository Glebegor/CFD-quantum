#pragma once

#include "Mesh.hpp"
#include "D2Q9.hpp"

/*
    Interpolation-supplemented (off-lattice) streaming.

    Instead of the classical index shift f_q(x) <- f_q(x - e_q), each
    population is reconstructed at its departure point from the
    per-direction interpolation stencil built by Mesh. Uses an explicit
    double buffer: the streamed state is written into `fNext` and then
    swapped into `f`.
*/
class ISLBMStreaming
{

public:
  void apply(
      const Mesh &mesh,
      LBMConstants::Distributions &f,
      LBMConstants::Distributions &fNext) const;
};
