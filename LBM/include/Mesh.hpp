#pragma once

#include <vector>

struct Node
{

  double x;
  double y;

  bool solid = false;

  /*
     ISLBM interpolation stencil.

     For arbitrary node:

     f(x-e,t)

     becomes:

     sum(
        weight_i *
        f(neighbour_i)
     )

  */

  std::vector<int> neighbours;

  std::vector<double> weights;
};

class Mesh
{

public:
  std::vector<Node> nodes;

  void generate(
      int resolution,
      double cylinderRadius);

  void buildInterpolation();
};