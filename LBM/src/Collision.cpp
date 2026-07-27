#include "Collision.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <cstddef>
#include <omp.h>

BGKCollision::BGKCollision(
    double reynolds)
{

  /*
      Navier-Stokes:

          Re = U L / nu


      LBM viscosity:

          nu = cs² (tau - 0.5)


      omega = 1/tau


      Pro jednoduchost:

          U = 0.1
          L = 1

  */

  double U = 0.1;
  double L = 1.0;

  double nu =
      U * L / reynolds;

  double tau =
      0.5 +
      nu / LBMConstants::cs2;

  omega =
      1.0 / tau;
}

void BGKCollision::apply(

    std::array<
        std::vector<double>,
        9> &f,

    size_t id

)
{

  double rho = 0.0;

  double ux = 0.0;

  double uy = 0.0;

  /*
      Macroscopic quantities

      rho = sum(fi)

      u = sum(fi*ci)/rho

  */

  for (int q = 0; q < 9; q++)
  {

    double fq =
        f[q][id];

    rho += fq;

    ux +=
        fq *
        LBMConstants::cx[q];

    uy +=
        fq *
        LBMConstants::cy[q];
  }

  ux /= rho;

  uy /= rho;

  /*
      BGK collision:


      f = f - omega(f-feq)

  */

  for (int q = 0; q < 9; q++)
  {

    double cu =
        LBMConstants::cx[q] * ux +
        LBMConstants::cy[q] * uy;

    double u2 =
        ux * ux +
        uy * uy;

    double feq =

        LBMConstants::weights[q] *
        rho *
        (1.0 +
         cu /
             LBMConstants::cs2

         +

         0.5 *
             cu * cu /
             (LBMConstants::cs2 *
              LBMConstants::cs2)

         -

         0.5 *
             u2 /
             LBMConstants::cs2);

    f[q][id] +=
        omega *
        (feq -
         f[q][id]);
  }
}