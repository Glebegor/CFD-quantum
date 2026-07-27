#pragma once

#include "Mesh.hpp"
#include "Streaming.hpp"
#include "Collision.hpp"

class LBMSolver
{

  Mesh &mesh;

  std::array<
      std::vector<double>,
      9>
      f;

  std::array<
      std::vector<double>,
      9>
      fNext;

  BGKCollision collision;

  ISLBMStreaming streaming;

public:
  LBMSolver(
      Mesh &mesh,
      double Re);

  void initialize();

  void step();

  auto &distributions() const
  {
    return f;
  }
};