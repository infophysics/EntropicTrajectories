//  Class for generating random numbers
#pragma once

#include <iostream>
#include <vector>
#include <random>

#include "utils.h"

namespace ET
{
  class Random
  {
  public:
    Random();
    ~Random();

  private:
    std::mt16637_64 _generator;
  };
}
