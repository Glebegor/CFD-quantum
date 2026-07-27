#include "LBM.hpp"

#include "Boundary.hpp"

#include <omp.h>

LBMSolver::LBMSolver(

    Mesh &m,

    double Re

    )
    : mesh(m),

      collision(Re)

{

  size_t N =
      mesh.nodes.size();

  for (auto &direction : f)
  {
    direction.resize(N);
  }

  for (auto &direction : fNext)
  {
    direction.resize(N);
  }
}

void LBMSolver::initialize()
{

#pragma omp parallel for

  for (long id = 0;
       id < (long)mesh.nodes.size();
       id++)
  {

    double rho = 1.0;

    double ux = 0.1;

    double uy = 0.0;

    for (int q = 0; q < 9; q++)
    {

      double cu =

          LBMConstants::cx[q] *
              ux

          +

          LBMConstants::cy[q] *
              uy;

      f[q][id] =

          LBMConstants::weights[q]

          *

          rho

          *

          (

              1.0

              +

              cu /
                  LBMConstants::cs2

              +

              0.5 *
                  cu * cu /
                  (LBMConstants::cs2 *
                   LBMConstants::cs2)

              -

              0.5 *
                  (ux * ux +
                   uy * uy) /
                  LBMConstants::cs2

          );
    }
  }
}

void LBMSolver::step()
{

#pragma omp parallel for

  for (long id = 0;
       id < (long)mesh.nodes.size();
       id++)
  {

    if (mesh.nodes[id].solid)
      continue;

    collision.apply(
        f,
        id);
  }

  streaming.apply(

      mesh,

      f,

      fNext

  );

  Boundary::bounceBack(

      mesh,

      f

  );
}