#pragma once

#include "D2Q9.hpp"

#include <array>
#include <vector>
#include <cstddef>

class BGKCollision
{

private:
  double omega;

public:
  BGKCollision(
      double reynolds);

  void apply(

      std::array<
          std::vector<double>,
          9> &f,

      size_t id

  );
};