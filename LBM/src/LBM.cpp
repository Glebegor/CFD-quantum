#include "LBM.hpp"
#include "Boundary.hpp"

LBMSolver::LBMSolver(Mesh &mesh, const Parameters &params)
    : mesh_(mesh),
      params_(params),
      collision_(params.omega())
{
  const std::size_t n = mesh_.nodes.size();

  for (auto &direction : f_)
    direction.assign(n, 0.0);

  for (auto &direction : fNext_)
    direction.assign(n, 0.0);
}

void LBMSolver::initialize()
{
  // Uniform flow at the inlet velocity: every population is set to the
  // corresponding equilibrium value (rho = 1, u = U_inlet, v = 0).
  const double u = params_.inletVelocity;

#pragma omp parallel for schedule(static)
  for (long id = 0; id < static_cast<long>(mesh_.nodes.size()); id++)
  {
    for (int q = 0; q < LBMConstants::Q; q++)
      f_[q][id] = LBMConstants::equilibrium(q, 1.0, u, 0.0);
  }
}

void LBMSolver::step()
{
  // 1. Collision (BGK) — in place on the fluid nodes.
#pragma omp parallel for schedule(static)
  for (long id = 0; id < static_cast<long>(mesh_.nodes.size()); id++)
  {
    if (mesh_.nodes[id].solid)
      continue;
    collision_.apply(f_, id);
  }

  // 2. Streaming (ISLBM interpolation) with a double-buffer swap.
  streaming_.apply(mesh_, f_, fNext_);

  // 3. No-slip on the cylinder.
  Boundary::bounceBack(mesh_, f_);

  // 4. External flow conditions on the domain edges.
  Boundary::applyDomain(mesh_, f_, params_.inletVelocity);
}
