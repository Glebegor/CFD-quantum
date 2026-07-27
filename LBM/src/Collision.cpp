#include "Collision.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace
{
  // Emit a detailed diagnostic and terminate. Silent numerical
  // corruption is never allowed to propagate.
  [[noreturn]] void fail(
      const char *what,
      std::size_t id,
      const LBMConstants::Distributions &f)
  {
    std::cerr << "\nLBM INSTABILITY: " << what << "\n"
              << "  node: " << id << "\n";
    for (int q = 0; q < LBMConstants::Q; q++)
      std::cerr << "  f[" << q << "] = " << f[q][id] << "\n";
    std::cerr << std::flush;
    std::abort();
  }
}

BGKCollision::BGKCollision(double omega)
    : omega_(omega)
{
  // A valid BGK relaxation frequency lies in (0, 2). Values outside
  // this range mean the derived viscosity / tau is unphysical.
  if (!std::isfinite(omega_) || omega_ <= 0.0 || omega_ >= 2.0)
  {
    std::cerr << "LBM CONFIG ERROR: relaxation frequency omega = "
              << omega_ << " is outside the stable range (0, 2).\n"
              << "Check Reynolds number, inlet velocity and resolution.\n";
    std::abort();
  }
}

void BGKCollision::apply(
    LBMConstants::Distributions &f,
    std::size_t id) const
{
  // --- reconstruct macroscopic moments ---
  for (int q = 0; q < LBMConstants::Q; q++)
    if (!std::isfinite(f[q][id]))
      fail("non-finite distribution before collision", id, f);

  double rho, ux, uy;
  LBMConstants::computeMoments(f, id, rho, ux, uy);

  // --- sanity checks: density must be positive & finite ---
  if (!std::isfinite(rho) || rho <= 0.0)
    fail("non-positive or non-finite density", id, f);

  if (!std::isfinite(ux) || !std::isfinite(uy))
    fail("non-finite velocity", id, f);

  // --- velocity ceiling: catch divergence before it becomes NaN ---
  if (ux * ux + uy * uy >
      LBMConstants::maxVelocity * LBMConstants::maxVelocity)
    fail("velocity exceeded stability limit (Mach too high)", id, f);

  // --- BGK relaxation towards equilibrium: f += omega (feq - f) ---
  for (int q = 0; q < LBMConstants::Q; q++)
  {
    const double feq = LBMConstants::equilibrium(q, rho, ux, uy);
    f[q][id] += omega_ * (feq - f[q][id]);

    if (!std::isfinite(f[q][id]))
      fail("collision produced a non-finite value", id, f);
  }
}
