//  Vector cpp for Entropic Trajectories and EDMC
#include <vector>
#include <iostream>
#include <stdio.h>
#include <complex>
#include "vector.h"

namespace ET
{

  template<typename T>
  Vector<T>::Vector()
  {

  }
  template<typename T>
  Vector<T>::~Vector()
  {

  }
  template<typename T>
  Vector<T>::Vector(unsigned int n) : _n(n)
  {

  }
  template<typename T>
  Vector<T>::Vector(std::vector<T> v) : _v(v), _n(_v.size())
  {

  }

  template<typename T>
  unsigned int Vector<T>::get_dim() const
  {
    return _n;
  }

  template<typename T>
  T& Vector<T>::operator()(const unsigned int& i)
  {
    return this->_v[i];
  }
  template<typename T>
  const T& Vector<T>::operator()(const unsigned int& i) const
  {
    return this->_v[i];
  }

  template<typename T>
  Vector<T> zeroes(unsigned int n)
  {
    std::vector<T> vec(n,0.0);
    Vector<T> v(vec);
    return v;
  }

}
