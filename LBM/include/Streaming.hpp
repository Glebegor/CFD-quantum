#pragma once

#include "Mesh.hpp"
#include "D2Q9.hpp"

#include <array>
#include <vector>

class ISLBMStreaming
{

public:
  void apply(

      Mesh &mesh,

      std::array<
          std::vector<double>, 9> &f,

      std::array<
          std::vector<double>, 9> &fNext

  );
};