#pragma once


enum ApproxType
{
  MLS,
  RBF,
};

struct ApproxParams
{
  //  TODO:: implement a set of parameters for each type.
  uint64_t k = 3;             //  number of nearest neighbors
  uint64_t n = 3;             //  order of polynomial expansion
};
