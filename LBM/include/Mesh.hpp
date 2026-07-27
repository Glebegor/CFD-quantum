#pragma once

#include "D2Q9.hpp"

#include <array>
#include <vector>
#include <cstddef>

/*
    Off-lattice node used by the interpolation-supplemented LBM.

    Node positions are a regular grid perturbed by a small jitter, so
    the departure point x - e_q of a population never lands exactly on
    another node. Streaming is therefore reconstructed by interpolation
    from a small neighbour stencil, built once per direction.
*/
struct Node
{
  double x = 0.0;
  double y = 0.0;

  // Logical (i, j) grid coordinates. Retained so the solver can reason
  // about topology (edges, obstacle neighbourhood) without a search.
  int i = 0;
  int j = 0;

  bool solid = false;

  /*
     ISLBM interpolation stencils, one per lattice direction q.

     The post-collision population that must arrive at this node along
     direction q originates at the departure point

         x_dep = x_node - e_q * h

     which is reconstructed as

         f_q(x_dep) = sum_k stencilWeights[q][k] * f_q(stencil[q][k])

     For the rest velocity (q = 0) the stencil is simply the node itself.
  */
  std::array<std::vector<int>, LBMConstants::Q> stencil;
  std::array<std::vector<double>, LBMConstants::Q> stencilWeights;
};

/*
    Unstructured quadrilateral mesh over the unit square.

    Although the nodes are generated from a logical grid, the topology
    is stored explicitly as VTK-compatible quad cells so the result can
    be exported as an UNSTRUCTURED_GRID (surface, contours, streamlines
    in ParaView) rather than a bare point cloud.
*/
class Mesh
{

public:
  std::vector<Node> nodes;

  // Quad cell connectivity: four node indices per cell, CCW ordering.
  std::vector<std::array<int, 4>> cells;

  // Grid resolution (nodes per axis). Zero until generate() is called.
  int resolution = 0;

  // Build the node cloud and the quad-cell connectivity.
  void generate(
      int resolution,
      double cylinderRadius,
      double cylinderX,
      double cylinderY);

  // Build the per-direction interpolation stencils for ISLBM streaming.
  void buildInterpolation();

  // Flat index of the logical grid node (i, j).
  [[nodiscard]] int index(int i, int j) const
  {
    return j * resolution + i;
  }

  [[nodiscard]] std::size_t nodeCount() const { return nodes.size(); }
  [[nodiscard]] std::size_t cellCount() const { return cells.size(); }
};
