#pragma once

#include <array>

namespace LBMConstants
{

  constexpr int Q = 9;

  constexpr std::array<int, Q> cx =
      {
          0, 1, 0, -1, 0,
          1, -1, -1, 1};

  constexpr std::array<int, Q> cy =
      {
          0, 0, 1, 0, -1,
          1, 1, -1, -1};

  constexpr std::array<int, Q> opposite =
      {
          0, 3, 4, 1, 2,
          7, 8, 5, 6};

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

  constexpr double cs2 = 1.0 / 3.0;

}