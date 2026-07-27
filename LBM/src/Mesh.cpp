#include "Mesh.hpp"

#include <cmath>
#include <random>
#include <algorithm>

void Mesh::generate(
    int resolution,
    double cylinderRadius)
{

  nodes.reserve(
      resolution * resolution);

  /*
      Off-lattice node generation.

      Nejde o klasický LBM grid.

      Každý bod má pouze:
          x
          y
          interpolation neighbours

      Streaming bude později řešen
      pomocí interpolace.
  */

  std::mt19937 rng(42);

  std::uniform_real_distribution<double>
      jitter(-0.0008, 0.0008);

  for (int j = 0; j < resolution; j++)
  {

    for (int i = 0; i < resolution; i++)
    {

      double x =
          static_cast<double>(i) /
          (resolution - 1);

      double y =
          static_cast<double>(j) /
          (resolution - 1);

      /*
          Small perturbation:

          vytváří off-lattice distribuci.
      */

      x += jitter(rng);
      y += jitter(rng);

      Node node;

      node.x = x;
      node.y = y;

      double dx =
          x - 0.25;

      double dy =
          y - 0.5;

      if (
          std::sqrt(
              dx * dx +
              dy * dy) <
          cylinderRadius)
      {
        node.solid = true;
      }

      nodes.push_back(node);
    }
  }
}

void Mesh::buildInterpolation()
{

  /*
      ISLBM:

      potřebujeme najít body:

          x - ei

      v našem prostoru.


      Produkční verze:

          KD-tree
          k-nearest neighbours
          MLS weights


      Tady používáme jednoduchý
      lokální search jako baseline.
  */

  constexpr int stencilSize = 8;

  for (size_t id = 0;
       id < nodes.size();
       id++)
  {

    auto &node =
        nodes[id];

    std::vector<
        std::pair<double, int>>
        distance;

    distance.reserve(
        nodes.size());

    for (size_t j = 0;
         j < nodes.size();
         j++)
    {

      if (id == j)
        continue;

      double dx =
          nodes[j].x - node.x;

      double dy =
          nodes[j].y - node.y;

      double d =
          dx * dx + dy * dy;

      distance.push_back(
          {d,
           static_cast<int>(j)});
    }

    std::partial_sort(

        distance.begin(),

        distance.begin() +
            stencilSize,

        distance.end()

    );

    for (int k = 0; k < stencilSize; k++)
    {

      node.neighbours.push_back(
          distance[k].second);
    }

    /*
        Simple inverse distance weighting.

        Produkce:

        MLS / RBF

    */

    double sum = 0;

    for (int n :
         node.neighbours)
    {

      double dx =
          nodes[n].x - node.x;

      double dy =
          nodes[n].y - node.y;

      double w =
          1.0 /
          (std::sqrt(
               dx * dx +
               dy * dy) +
           1e-12);

      node.weights.push_back(w);

      sum += w;
    }

    for (auto &w :
         node.weights)
    {
      w /= sum;
    }
  }
}