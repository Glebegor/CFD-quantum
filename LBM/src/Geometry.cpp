#include "Geometry.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

// =====================================================================
// Circle
// =====================================================================
Circle::Circle(double centreX, double centreY, double radius)
    : cx_(centreX), cy_(centreY), r2_(radius * radius) {}

bool Circle::contains(double x, double y) const
{
  const double dx = x - cx_;
  const double dy = y - cy_;
  return dx * dx + dy * dy < r2_;
}

// =====================================================================
// Polygon
// =====================================================================
Polygon::Polygon(std::vector<Point> worldPoints)
    : pts_(std::move(worldPoints)) {}

// Even-odd ray-casting point-in-polygon test.
bool Polygon::contains(double x, double y) const
{
  const std::size_t n = pts_.size();
  if (n < 3)
    return false;

  bool inside = false;
  for (std::size_t i = 0, j = n - 1; i < n; j = i++)
  {
    const double xi = pts_[i][0], yi = pts_[i][1];
    const double xj = pts_[j][0], yj = pts_[j][1];

    const bool straddles = (yi > y) != (yj > y);
    if (straddles &&
        x < (xj - xi) * (y - yi) / (yj - yi) + xi)
      inside = !inside;
  }
  return inside;
}

std::vector<Point> Polygon::loadSection(const std::string &path)
{
  std::ifstream file(path);
  if (!file)
    return {};

  // Read every parseable coordinate pair, rejecting header lines and
  // the point-count line (values well outside the normalised range).
  std::vector<Point> raw;
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream ss(line);
    double x, y;
    if (!(ss >> x >> y))
      continue; // title / non-numeric line
    if (x < -1.5 || x > 1.5 || y < -1.5 || y > 1.5)
      continue; // e.g. the "66.  66." count line
    raw.push_back({x, y});
  }

  if (raw.size() < 3)
    return {};

  // Detect a mid-file reset of x (two-surface NACA/Lednicer layout,
  // where the upper and lower surfaces both run leading -> trailing).
  std::size_t split = 0;
  for (std::size_t i = 1; i < raw.size(); i++)
    if (raw[i][0] < raw[i - 1][0] - 0.5)
    {
      split = i;
      break;
    }

  if (split == 0)
    return raw; // already a single ordered loop

  // Stitch into one loop: upper (LE -> TE) then lower reversed (TE -> LE).
  std::vector<Point> loop;
  loop.reserve(raw.size());
  for (std::size_t i = 0; i < split; i++)
    loop.push_back(raw[i]);
  for (std::size_t i = raw.size(); i-- > split;)
    loop.push_back(raw[i]);
  return loop;
}

std::vector<Point> Polygon::nacaSection(double thickness, int samplesPerSide)
{
  // Symmetric NACA 4-digit half-thickness distribution.
  const auto yt = [thickness](double x)
  {
    return 5.0 * thickness *
           (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
            0.2843 * x * x * x - 0.1015 * x * x * x * x);
  };

  std::vector<Point> upper, lower;
  upper.reserve(samplesPerSide);
  lower.reserve(samplesPerSide);

  // Cosine spacing clusters points at the leading and trailing edges.
  for (int k = 0; k < samplesPerSide; k++)
  {
    const double beta = M_PI * k / (samplesPerSide - 1);
    const double x = 0.5 * (1.0 - std::cos(beta));
    upper.push_back({x, yt(x)});
    lower.push_back({x, -yt(x)});
  }

  std::vector<Point> loop;
  loop.reserve(2 * samplesPerSide);
  for (const auto &p : upper)
    loop.push_back(p);
  for (auto it = lower.rbegin(); it != lower.rend(); ++it)
    loop.push_back(*it);
  return loop;
}

Polygon Polygon::place(const std::vector<Point> &section,
                       double chord, double angleOfAttackDeg,
                       double refX, double refY)
{
  // Normalise the section to its x-extent so any outline (not just a
  // unit-chord airfoil) scales correctly.
  double xmin = section.front()[0], xmax = section.front()[0];
  for (const auto &p : section)
  {
    xmin = std::min(xmin, p[0]);
    xmax = std::max(xmax, p[0]);
  }
  const double extent = (xmax - xmin) > 1e-12 ? (xmax - xmin) : 1.0;
  const double scale = chord / extent;

  // Positive angle of attack = nose up against a +x free stream, i.e. a
  // clockwise rotation of the section about its leading edge.
  const double theta = -angleOfAttackDeg * M_PI / 180.0;
  const double c = std::cos(theta);
  const double s = std::sin(theta);

  std::vector<Point> world;
  world.reserve(section.size());
  for (const auto &p : section)
  {
    const double sx = (p[0] - xmin) * scale;
    const double sy = p[1] * scale;
    world.push_back({refX + (sx * c - sy * s),
                     refY + (sx * s + sy * c)});
  }
  return Polygon(std::move(world));
}

// =====================================================================
// Factory
// =====================================================================
Geometry makeGeometry(const Parameters &params)
{
  Geometry geometry;

  if (params.obstacle == Obstacle::Cylinder)
  {
    geometry.add(std::make_unique<Circle>(
        params.cylinderX, params.cylinderY, params.cylinderRadius));
    return geometry;
  }

  // Airfoil: prefer the coordinate file, fall back to analytic NACA 0012.
  std::vector<Point> section = Polygon::loadSection(params.airfoilFile);
  if (section.empty())
  {
    std::cout << "[geometry] could not read '" << params.airfoilFile
              << "'; using analytic NACA 0012.\n";
    section = Polygon::nacaSection(0.12, 80);
  }
  else
  {
    std::cout << "[geometry] loaded airfoil '" << params.airfoilFile
              << "' (" << section.size() << " points).\n";
  }

  geometry.add(std::make_unique<Polygon>(Polygon::place(
      section, params.chord, params.angleOfAttack,
      params.airfoilX, params.airfoilY)));
  return geometry;
}
