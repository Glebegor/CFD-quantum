#include "Mesh.hpp"

#include <cmath>
#include <random>
#include <algorithm>
#include <array>

/*
    Node generation.

    A regular resolution x resolution grid over the unit square is
    perturbed by a small jitter so that the mesh is genuinely
    off-lattice (the premise of ISLBM). Logical (i, j) indices are
    retained so topology queries stay O(1).
*/
void Mesh::generate(
    int res,
    const SolidPredicate &isSolid)
{
  resolution = res;

  nodes.clear();
  nodes.reserve(static_cast<std::size_t>(res) * res);

  std::mt19937 rng(42);

  // Jitter amplitude: a few percent of the cell size, small enough
  // that the interpolation stencil stays well conditioned.
  const double h = 1.0 / static_cast<double>(res - 1);
  std::uniform_real_distribution<double> jitter(-0.05 * h, 0.05 * h);

  for (int j = 0; j < res; j++)
  {
    for (int i = 0; i < res; i++)
    {
      Node node;
      node.i = i;
      node.j = j;

      node.x = static_cast<double>(i) * h;
      node.y = static_cast<double>(j) * h;

      // Off-lattice perturbation (interior nodes only; the domain
      // border is kept aligned so boundary conditions stay clean).
      if (i > 0 && i < res - 1 && j > 0 && j < res - 1)
      {
        node.x += jitter(rng);
        node.y += jitter(rng);
      }

      // Solid mask: delegate the obstacle test to the geometry.
      node.solid = isSolid(node.x, node.y);

      nodes.push_back(node);
    }
  }

  /*
      Quad-cell connectivity (VTK_QUAD).

      Each logical cell (i, j) connects the four surrounding nodes in
      counter-clockwise order. Cells whose four corners are all solid
      lie fully inside the cylinder and are dropped, leaving a clean
      hole in the exported surface.
  */
  cells.clear();
  cells.reserve(static_cast<std::size_t>(res - 1) * (res - 1));

  for (int j = 0; j < res - 1; j++)
  {
    for (int i = 0; i < res - 1; i++)
    {
      const int a = index(i, j);
      const int b = index(i + 1, j);
      const int c = index(i + 1, j + 1);
      const int d = index(i, j + 1);

      const bool fullySolid =
          nodes[a].solid && nodes[b].solid &&
          nodes[c].solid && nodes[d].solid;

      if (fullySolid)
        continue;

      cells.push_back({a, b, c, d});
    }
  }
}

/*
    Build the ISLBM interpolation stencils.

    For every node and every lattice direction q we reconstruct the
    population arriving from the departure point

        x_dep = x_node - e_q * h.

    Because the mesh is a lightly perturbed grid, the departure point
    lies very close to the logical node (i - cx, j - cy). We therefore
    gather the clamped 3x3 logical block around that node as candidates
    and build inverse-distance-squared weights to the nearest four.
    This keeps the search O(1) per direction (no global neighbour
    query) while remaining a genuine off-lattice interpolation.

    The comment on the original brute-force O(N^2) version noted that a
    KD-tree would be needed for production; exploiting the known logical
    topology removes that cost entirely for structured node clouds.
*/
void Mesh::buildInterpolation()
{
  const int res = resolution;
  const double h = 1.0 / static_cast<double>(res - 1);

  constexpr int maxStencil = 4;

#pragma omp parallel for schedule(static)
  for (long id = 0; id < static_cast<long>(nodes.size()); id++)
  {
    Node &node = nodes[id];

    for (int q = 0; q < LBMConstants::Q; q++)
    {
      node.stencil[q].clear();
      node.stencilWeights[q].clear();

      // Rest velocity: population stays put.
      if (LBMConstants::cx[q] == 0 && LBMConstants::cy[q] == 0)
      {
        node.stencil[q].push_back(static_cast<int>(id));
        node.stencilWeights[q].push_back(1.0);
        continue;
      }

      // Departure point for this direction.
      const double xDep = node.x - LBMConstants::cx[q] * h;
      const double yDep = node.y - LBMConstants::cy[q] * h;

      // Logical node closest to the departure point.
      const int baseI = node.i - LBMConstants::cx[q];
      const int baseJ = node.j - LBMConstants::cy[q];

      // Gather the clamped 3x3 candidate block and keep the nearest.
      std::array<std::pair<double, int>, 9> candidates;
      int count = 0;

      for (int dj = -1; dj <= 1; dj++)
      {
        for (int di = -1; di <= 1; di++)
        {
          int ci = baseI + di;
          int cj = baseJ + dj;

          if (ci < 0)
            ci = 0;
          if (cj < 0)
            cj = 0;
          if (ci > res - 1)
            ci = res - 1;
          if (cj > res - 1)
            cj = res - 1;

          const int cid = index(ci, cj);
          const double ddx = nodes[cid].x - xDep;
          const double ddy = nodes[cid].y - yDep;

          candidates[count++] = {ddx * ddx + ddy * ddy, cid};
        }
      }

      // Nearest candidates first.
      std::sort(
          candidates.begin(),
          candidates.begin() + count,
          [](const auto &lhs, const auto &rhs)
          { return lhs.first < rhs.first; });

      const int k = std::min(maxStencil, count);

      // Inverse-distance-squared weighting. The near-coincident node
      // dominates, so the scheme reduces to standard lattice streaming
      // on a uniform grid and adds only sub-cell interpolation when the
      // mesh is perturbed.
      double sum = 0.0;
      for (int n = 0; n < k; n++)
      {
        const int cid = candidates[n].second;

        // Skip duplicates from clamping at the border.
        if (std::find(node.stencil[q].begin(), node.stencil[q].end(), cid) !=
            node.stencil[q].end())
          continue;

        const double w = 1.0 / (candidates[n].first + 1e-12);
        node.stencil[q].push_back(cid);
        node.stencilWeights[q].push_back(w);
        sum += w;
      }

      for (double &w : node.stencilWeights[q])
        w /= sum;
    }
  }
}
