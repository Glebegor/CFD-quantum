#pragma once

#include "Mesh.hpp"
#include "D2Q9.hpp"

/*
    Boundary conditions.

    Two independent conditions are applied after streaming:

      * bounceBack  - no-slip on the cylinder (obstacle) nodes.
      * applyDomain - external flow conditions on the four domain
                      edges: velocity inlet (left), zero-gradient
                      outflow (right), free-stream far field (top,
                      bottom).
*/
class Boundary
{

public:
  // Halfway bounce-back: reflect populations on solid nodes.
  static void bounceBack(
      const Mesh &mesh,
      LBMConstants::Distributions &f);

  // Domain-edge conditions driven by the inlet velocity.
  static void applyDomain(
      const Mesh &mesh,
      LBMConstants::Distributions &f,
      double inletVelocity);
};
