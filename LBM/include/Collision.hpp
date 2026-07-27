#pragma once

#include "D2Q9.hpp"

#include <cstddef>

/*
    BGK single-relaxation-time collision operator.

    The relaxation frequency omega is supplied by the configuration
    layer (Parameters) so that no viscosity / Reynolds physics is
    hardcoded here.
*/
class BGKCollision
{

public:
  explicit BGKCollision(double omega);

  // Relax the distributions at node `id` towards local equilibrium.
  // Performs numerical sanity checks and aborts loudly on corruption.
  void apply(
      LBMConstants::Distributions &f,
      std::size_t id) const;

  [[nodiscard]] double relaxationFrequency() const { return omega_; }

private:
  double omega_;
};
