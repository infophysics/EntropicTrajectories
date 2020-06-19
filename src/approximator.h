//  Approximator for fields
#pragma once


#include <vector>
#include <iostream>
#include <map>
#include "grid.h"
#include "params.h"
#include "utils.h"
#include "matrix.h"


namespace ET
{
  template<typename T> class ScalarField;
}
#include "scalar.h"
namespace ET
{
  template<typename T>
  class Approximator
  {
  public:
    Approximator();
    ~Approximator();

    //  Setters
    void setDerivative(std::string type);

    //  Gradient functions
    std::vector<T> gradient(Grid<T>* grid, ScalarField<T>* field, uint64_t index);
    std::vector<T> gradientMLS(Grid<T>* grid, ScalarField<T>* field, uint64_t index);

    Matrix<T> constructBMatrix(Grid<T>* grid, std::vector<uint64_t>* neighbors,
      uint64_t index, uint64_t order);
  private:
    enum ApproxType _type;
    struct ApproxParams _params;
  };


  template class Approximator<double>;
}
