#include "Streaming.hpp"

#include <omp.h>

void ISLBMStreaming::apply(

    Mesh &mesh,

    std::array<
        std::vector<double>,
        9> &f,

    std::array<
        std::vector<double>,
        9> &fNext

)
{

#pragma omp parallel for

  for (long id = 0;
       id < (long)mesh.nodes.size();
       id++)
  {

    Node &node =
        mesh.nodes[id];

    if (node.solid)
      continue;

    /*
        ISLBM streaming


        místo:

        f[q][x-e]


        děláme:


        sum(
            interpolationWeight *
            neighbourDistribution
        )


    */

    for (int q = 0; q < 9; q++)
    {

      double value = 0.0;

      for (size_t k = 0;
           k < node.neighbours.size();
           k++)
      {

        int neighbour =
            node.neighbours[k];

        value +=

            node.weights[k]

            *

            f[q][neighbour];
      }

      fNext[q][id] = value;
    }
  }

  std::swap(
      f,
      fNext);
}