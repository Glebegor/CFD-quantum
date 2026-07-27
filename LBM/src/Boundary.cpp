#include "Boundary.hpp"

void Boundary::bounceBack(
    const Mesh &mesh,
    LBMConstants::Distributions &f)
{
  const long n = static_cast<long>(mesh.nodes.size());

#pragma omp parallel for schedule(static)
  for (long id = 0; id < n; id++)
  {
    if (!mesh.nodes[id].solid)
      continue;

    // No-slip: every population is reflected into its opposite
    // direction, so momentum incident on the obstacle is reversed.
    for (int q = 1; q < LBMConstants::Q; q++)
    {
      const int qOpp = LBMConstants::opposite[q];
      if (q < qOpp) // swap each pair once
        std::swap(f[q][id], f[qOpp][id]);
    }
  }
}

void Boundary::applyDomain(
    const Mesh &mesh,
    LBMConstants::Distributions &f,
    double inletVelocity)
{
  const int res = mesh.resolution;

#pragma omp parallel for schedule(static)
  for (long id = 0; id < static_cast<long>(mesh.nodes.size()); id++)
  {
    const Node &node = mesh.nodes[id];

    // Obstacle nodes are handled by bounceBack().
    if (node.solid)
      continue;

    const int i = node.i;
    const int j = node.j;

    // Velocity inlet on the left edge: impose equilibrium at the
    // prescribed free-stream velocity (Dirichlet).
    if (i == 0)
    {
      for (int q = 0; q < LBMConstants::Q; q++)
        f[q][id] = LBMConstants::equilibrium(q, 1.0, inletVelocity, 0.0);
      continue;
    }

    // Zero-gradient (convective) outflow on the right edge: copy the
    // upstream neighbour so disturbances leave the domain cleanly.
    if (i == res - 1)
    {
      const int upstream = mesh.index(i - 1, j);
      for (int q = 0; q < LBMConstants::Q; q++)
        f[q][id] = f[q][upstream];
      continue;
    }

    // Free-stream far field on the top and bottom edges.
    if (j == 0 || j == res - 1)
    {
      for (int q = 0; q < LBMConstants::Q; q++)
        f[q][id] = LBMConstants::equilibrium(q, 1.0, inletVelocity, 0.0);
    }
  }
}
