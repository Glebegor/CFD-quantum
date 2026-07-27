#pragma once

#include <array>
#include <vector>
#include <cstddef>
#include <cmath>

/*
    D2Q9 lattice definition and shared lattice-Boltzmann helpers.

    Lattice velocity numbering:

        6   2   5
          \ | /
        3 - 0 - 1
          / | \
        7   4   8
*/

namespace LBMConstants
{

  constexpr int Q = 9;

  // Discrete lattice velocities.
  constexpr std::array<int, Q> cx =
      {
          0, 1, 0, -1, 0,
          1, -1, -1, 1};

  constexpr std::array<int, Q> cy =
      {
          0, 0, 1, 0, -1,
          1, 1, -1, -1};

  // Index of the velocity pointing in the opposite direction
  // (used by the bounce-back boundary condition).
  constexpr std::array<int, Q> opposite =
      {
          0, 3, 4, 1, 2,
          7, 8, 5, 6};

  // Lattice weights for the equilibrium distribution.
  constexpr std::array<double, Q> weights =
      {
          4.0 / 9.0,

          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,

          1.0 / 36.0,
          1.0 / 36.0,
          1.0 / 36.0,
          1.0 / 36.0};

  // Speed of sound squared, and speed of sound.
  constexpr double cs2 = 1.0 / 3.0;
  const double cs = std::sqrt(cs2);

  // Numerical stability ceiling for the macroscopic velocity
  // magnitude (lattice units). Well above any physical inflow but
  // low enough to catch a diverging solution before it produces
  // NaN/Inf. Corresponds to Mach ~= 0.7.
  constexpr double maxVelocity = 0.4;

  // Structure-of-arrays distribution storage: one contiguous
  // vector per lattice direction.
  using Distributions = std::array<std::vector<double>, Q>;

  // -----------------------------------------------------------
  // Discrete equilibrium distribution (second order in velocity):
  //
  //   feq_q = w_q * rho * [ 1 + (c.u)/cs^2
  //                           + (c.u)^2 / (2 cs^4)
  //                           - u^2 / (2 cs^2) ]
  // -----------------------------------------------------------
  inline double equilibrium(
      int q,
      double rho,
      double ux,
      double uy)
  {
    const double cu = cx[q] * ux + cy[q] * uy;
    const double u2 = ux * ux + uy * uy;

    return weights[q] * rho *
           (1.0 + cu / cs2 + 0.5 * (cu * cu) / (cs2 * cs2) - 0.5 * u2 / cs2);
  }

  // -----------------------------------------------------------
  // Macroscopic moments at a single node:
  //
  //   rho    = sum_q f_q
  //   rho*u  = sum_q f_q * c_q
  //
  // Kept in one place so the collision, boundary and output code
  // never re-derive the reconstruction differently.
  // -----------------------------------------------------------
  inline void computeMoments(
      const Distributions &f,
      std::size_t id,
      double &rho,
      double &ux,
      double &uy)
  {
    rho = 0.0;
    double momentumX = 0.0;
    double momentumY = 0.0;

    for (int q = 0; q < Q; q++)
    {
      const double fq = f[q][id];
      rho += fq;
      momentumX += fq * cx[q];
      momentumY += fq * cy[q];
    }

    ux = momentumX / rho;
    uy = momentumY / rho;
  }

}
