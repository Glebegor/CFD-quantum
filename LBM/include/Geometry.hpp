#pragma once

#include "Parameters.hpp"

#include <array>
#include <memory>
#include <string>
#include <vector>

/*
    Extensible obstacle geometry.

    Every immersed body implements a single query: "is this point inside
    the solid?". The mesh only ever sees that predicate, so new shapes
    can be added without touching the solver.

    Provided shapes
    ---------------
      Circle   - analytic cylinder cross-section.
      Polygon  - arbitrary closed outline; an airfoil is just a polygon
                 loaded from a coordinate file (NACA / Lednicer format)
                 or generated analytically.

    Several shapes can be combined in a Geometry (logical union), so a
    scene may contain e.g. a wing plus a second body.
*/

using Point = std::array<double, 2>;

// Base class for an immersed solid body.
struct Shape
{
  virtual ~Shape() = default;

  // World-coordinate solid test.
  [[nodiscard]] virtual bool contains(double x, double y) const = 0;
};

// Analytic circle (cylinder cross-section).
class Circle : public Shape
{
public:
  Circle(double centreX, double centreY, double radius);
  [[nodiscard]] bool contains(double x, double y) const override;

private:
  double cx_;
  double cy_;
  double r2_;
};

// Closed polygon, tested with even-odd ray casting.
class Polygon : public Shape
{
public:
  explicit Polygon(std::vector<Point> worldPoints);
  [[nodiscard]] bool contains(double x, double y) const override;

  [[nodiscard]] const std::vector<Point> &points() const { return pts_; }

  // ---- section factories (return a unit-chord outline, x in ~[0,1]) ----

  // Read a 2-column coordinate file. Handles both a single closed loop
  // and the two-surface NACA/Lednicer layout. Returns an empty vector
  // if the file cannot be read or parsed.
  static std::vector<Point> loadSection(const std::string &path);

  // Analytic NACA 4-digit symmetric section (thickness as a fraction,
  // e.g. 0.12 for NACA 0012).
  static std::vector<Point> nacaSection(double thickness, int samplesPerSide);

  // Place a unit-chord section into the domain: scale to `chord`, rotate
  // by the angle of attack (degrees, positive = nose up) about the
  // leading edge, then move the leading edge to (refX, refY).
  static Polygon place(const std::vector<Point> &section,
                       double chord, double angleOfAttackDeg,
                       double refX, double refY);

private:
  std::vector<Point> pts_;
};

// Logical union of shapes: solid if inside any member.
class Geometry
{
public:
  void add(std::unique_ptr<Shape> shape) { shapes_.push_back(std::move(shape)); }

  [[nodiscard]] bool isSolid(double x, double y) const
  {
    for (const auto &shape : shapes_)
      if (shape->contains(x, y))
        return true;
    return false;
  }

  [[nodiscard]] bool empty() const { return shapes_.empty(); }

private:
  std::vector<std::unique_ptr<Shape>> shapes_;
};

// Build the geometry described by the configuration (single entry point
// to change what is immersed in the flow).
Geometry makeGeometry(const Parameters &params);
