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

// Obstacle immersed in the flow.
enum class Obstacle
{
    Cylinder,
    Airfoil
};

struct Parameters
{
    // =========================
    // Mesh
    // =========================

    // Number of nodes per axis (square domain).
    int resolution = 128;

    // =========================
    // Geometry (obstacle)
    // =========================

    // Which obstacle to immerse in the flow.
    Obstacle obstacle = Obstacle::Airfoil;

    // --- Cylinder ---
    // Cylinder radius in domain units (fraction of unit square).
    double cylinderRadius = 0.08;

    // Cylinder centre in domain units.
    double cylinderX = 0.25;
    double cylinderY = 0.50;

    // --- Airfoil (NACA) ---
    // Coordinate file (relative to the run directory). If it cannot be
    // opened, the geometry falls back to the analytic NACA 0012 shape.
    std::string airfoilFile = "../wing_naca_0012.txt";

    // Chord length in domain units.
    double chord = 0.35;

    // Leading-edge reference point in domain units.
    double airfoilX = 0.25;
    double airfoilY = 0.50;

    // Angle of attack in degrees (positive = nose up into the flow).
    double angleOfAttack = 10.0;

    // =========================
    // Physics
    // =========================

    // Reynolds number, Re = U * L / nu.
    double reynoldsNumber = 100.0;

    // Inlet velocity in lattice units (cells / timestep). Kept small
    // and fixed for numerical stability; the physical velocity is set
    // separately below and mapped onto this value.
    double inletVelocity = 0.05;

    // =========================
    // Physical scaling (real-world units)
    // =========================
    //
    // LBM is dimensionless. To obtain physical time / length we pin one
    // length scale (the chord) and one velocity scale, from which the
    // duration of a timestep follows. Everything else (snapshot stride,
    // total steps, output units) is then derived, not hand-tuned.

    // Chord length of the real wing. Default: Airbus A320 mean
    // aerodynamic chord (~4.2 m).
    double physicalChord = 4.2;

    // Free-stream speed in metres per second.
    double physicalVelocity = 1.0;

    // Animation cadence and length.
    double framesPerSecond = 10.0;
    double durationSeconds = 4.0;

    // Console logging interval (steps between log lines).
    int logInterval = 500;

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

    // Characteristic length of the obstacle in lattice cells: the
    // cylinder diameter, or the airfoil chord. This is the L used in Re.
    [[nodiscard]] double characteristicLength() const
    {
        const double lengthInDomain =
            (obstacle == Obstacle::Airfoil) ? chord : (2.0 * cylinderRadius);
        return lengthInDomain / gridSpacing();
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

    // ---------------------------------------------------------
    // Physical-unit conversion.
    //
    //   u_lattice = u_phys * dt / dx   =>   dt = u_lattice * dx / u_phys
    //
    // with dx fixed by the chord (metres per cell).
    // ---------------------------------------------------------

    // Physical size of one lattice cell [m].
    [[nodiscard]] double dxPhysical() const
    {
        return physicalChord / characteristicLength();
    }

    // Physical duration of one timestep [s].
    [[nodiscard]] double dtPhysical() const
    {
        return inletVelocity * dxPhysical() / physicalVelocity;
    }

    // Factor converting lattice velocity -> m/s (for output).
    [[nodiscard]] double velocityScale() const
    {
        return physicalVelocity / inletVelocity;
    }

    // Factor converting a domain-unit coordinate -> metres (for output).
    // One domain unit spans the whole square, i.e. chord/chordFraction.
    [[nodiscard]] double lengthScale() const
    {
        return physicalChord / chord;
    }

    // Steps between snapshots to hit the requested frame rate.
    [[nodiscard]] int snapshotStride() const
    {
        const double s = (1.0 / framesPerSecond) / dtPhysical();
        return s < 1.0 ? 1 : static_cast<int>(s + 0.5);
    }

    // Total number of steps covering the requested duration.
    [[nodiscard]] int totalSteps() const
    {
        return static_cast<int>(durationSeconds / dtPhysical() + 0.5);
    }
};
