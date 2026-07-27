#include "Boundary.hpp"

void Boundary::bounceBack(

    Mesh &mesh,

    std::array<
        std::vector<double>,
        9> &f

)
{

  for (size_t id = 0;
       id < mesh.nodes.size();
       id++)
  {

    if (!mesh.nodes[id].solid)
      continue;

    for (int q = 0;
         q < 9;
         q++)
    {

      int opposite =
          LBMConstants::opposite[q];

      std::swap(

          f[q][id],

          f[opposite][id]

      );
    }
  }
}