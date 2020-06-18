//  Approximator c++
#include "approximator.h"


namespace ET
{
  //  map for a string to enum of derivative type
  std::map<std::string, ApproxType> ApproxTypeMap =
  {
    { "MLS", ApproxType::MLS },
    { "RBF", ApproxType::RBF }
  };

  template<typename T>
  Approximator<T>::Approximator()
  {

  }

  template<typename T>
  Approximator<T>::~Approximator()
  {

  }

  //  Setters
  template<typename T>
  void Approximator<T>::setDerivative(std::string type)
  {
    _type = ApproxTypeMap[type];
  }

  template<typename T>
  std::vector<T> Approximator<T>::gradient(Grid<T>* grid, ScalarField<T>* field, uint64_t index)
  {
    if (_type == 0)
      return gradientMLS(grid, field, index);
  }
  template<typename T>
  std::vector<T> Approximator<T>::gradientMLS(Grid<T>* grid, ScalarField<T>* field, uint64_t index)
  {
    std::vector<T> result;
    //  TODO:: implement MLS here.
    return result;
  }
}
