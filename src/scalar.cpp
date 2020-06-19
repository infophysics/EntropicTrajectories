//  Scalar field cpp file



#include "scalar.h"


namespace ET
{
  template<typename T>
  ScalarField<T>::ScalarField()
  {

  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  {

  }

  template<typename T>
  Approximator<T>* ScalarField<T>::getApproximator()
  {
    return _approx;
  }

  template<typename T>
  void ScalarField<T>::setDerivative(std::string type)
  {
    _approx->setDerivative(type);
  }

  template<typename T>
  ScalarField<T>::ScalarField(Grid<T>* micro, std::vector<T> field)
  : _micro(micro), _field(field)
  {
    _approx = new Approximator<T>();
  }
}
