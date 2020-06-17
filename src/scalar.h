//  Scalar field class
#pragma once


namespace ET
{
  template<typename T>
  class ScalarField
  {
  public:
    ScalarField();
    ~ScalarField();

  private:
    Grid<T>* _micro;
    std::vector<T> _field;

  };
}
