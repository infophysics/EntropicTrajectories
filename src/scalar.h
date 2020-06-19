//  Scalar field class
#pragma once

#include <vector>
#include <iostream>

#include "grid.h"

namespace ET
{
  template<typename T> class Approximator;
}
#include "approximator.h"
namespace ET
{
  template<typename T>
  class ScalarField
  {
  public:
    ScalarField();
    ~ScalarField();
    ScalarField(Grid<T>* micro, std::vector<T> field);

    //  Getters
    Approximator<T>* getApproximator();

    //  Setters
    void setDerivative(std::string type);
    //  Methods for calculating derivatives
    Matrix<T> constructBMatrix();


  private:
    Grid<T>* _micro;
    std::vector<T> _field;
    Approximator<T>* _approx;

  };

  template class ScalarField<double>;
}
