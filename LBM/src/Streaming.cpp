#include "Streaming.hpp"

#include <utility>

void ISLBMStreaming::apply(
    const Mesh &mesh,
    LBMConstants::Distributions &f,
    LBMConstants::Distributions &fNext) const
{
  const long n = static_cast<long>(mesh.nodes.size());

#pragma omp parallel for schedule(static)
  for (long id = 0; id < n; id++)
  {
    const Node &node = mesh.nodes[id];

    // Solid nodes do not stream; carry their state across the buffer
    // swap so the bounce-back step operates on a consistent state.
    if (node.solid)
    {
      for (int q = 0; q < LBMConstants::Q; q++)
        fNext[q][id] = f[q][id];
      continue;
    }

    /*
        Off-lattice streaming:

            fNext_q(x) = sum_k w_k * f_q(neighbour_k)

        where the stencil for direction q interpolates the departure
        point x - e_q. The reconstruction reads only from `f`
        (post-collision) and writes only to `fNext`, so there is no
        read/write aliasing across nodes.
    */
    for (int q = 0; q < LBMConstants::Q; q++)
    {
      const auto &nb = node.stencil[q];
      const auto &w = node.stencilWeights[q];

      double value = 0.0;
      for (std::size_t k = 0; k < nb.size(); k++)
        value += w[k] * f[q][nb[k]];

      fNext[q][id] = value;
    }
  }

  // Double-buffer swap: `f` now holds the streamed populations.
  std::swap(f, fNext);
}
