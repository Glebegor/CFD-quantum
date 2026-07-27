#pragma once

#include "D2Q9.hpp"

#include <string>

/*
    Central configuration / physics layer.

    All physical and numerical constants live here so that no
    magic numbers are hardcoded inside the solver source files.

    Convention
    ----------
    The domain is the unit square [0,1] x [0,1] discretised on a
    `resolution x resolution` node grid. The lattice spacing is

        h = 1 / (resolution - 1)      [physical length per cell]

    LBM works in lattice units (cells / timesteps). The inlet
    velocity is therefore already expressed in lattice units
    (cells per timestep) and must stay well below the lattice
    speed of sound cs = sqrt(1/3) for the low-Mach / incompressible
    assumption to hold.
*/

struct Parameters
{
    // =========================
    // Mesh
    // =========================

    // Number of nodes per axis (square domain).
    int resolution = 512;

    // =========================
    // Geometry (obstacle)
    // =========================

    // Cylinder radius in domain units (fraction of unit square).
    double cylinderRadius = 0.08;

    // Cylinder centre in domain units.
    double cylinderX = 0.25;
    double cylinderY = 0.50;

    // =========================
    // Physics
    // =========================

    // Reynolds number, Re = U * L / nu.
    double reynoldsNumber = 100.0;

    // Inlet velocity in lattice units (cells / timestep).
    double inletVelocity = 0.05;

    // =========================
    // Simulation control
    // =========================

    int maxSteps = 20000;

    // VTK export interval (steps between snapshots).
    int vtkInterval = 1000;

    // Console logging interval (steps between log lines).
    int logInterval = 1000;

    // =========================
    // Output
    // =========================

    std::string outputDirectory = "output";

    // ---------------------------------------------------------
    // Derived lattice quantities.
    //
    // Kept here (not in the solver) so the physics definition
    // has a single source of truth.
    // ---------------------------------------------------------

    // Lattice spacing (physical length represented by one cell).
    [[nodiscard]] double gridSpacing() const
    {
        return 1.0 / static_cast<double>(resolution - 1);
    }

    // Characteristic length of the obstacle in lattice cells
    // (the cylinder diameter). This is the L used in Re.
    [[nodiscard]] double characteristicLength() const
    {
        return (2.0 * cylinderRadius) / gridSpacing();
    }

    // Kinematic viscosity in lattice units: nu = U * L / Re.
    [[nodiscard]] double viscosity() const
    {
        return (inletVelocity * characteristicLength()) / reynoldsNumber;
    }

    // BGK relaxation time: tau = 0.5 + nu / cs^2.
    [[nodiscard]] double relaxationTime() const
    {
        return 0.5 + viscosity() / LBMConstants::cs2;
    }

    // BGK relaxation frequency: omega = 1 / tau.
    [[nodiscard]] double omega() const
    {
        return 1.0 / relaxationTime();
    }

    // Mach number of the inflow (should stay well below ~0.3).
    [[nodiscard]] double machNumber() const
    {
        return inletVelocity / LBMConstants::cs;
    }
};
