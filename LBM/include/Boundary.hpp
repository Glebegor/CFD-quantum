#pragma once

#include "Mesh.hpp"
#include "D2Q9.hpp"

#include <array>
#include <vector>

class Boundary
{

public:
  static void bounceBack(
      Mesh &mesh,
      std::array<std::vector<double>, 9> &f);
};