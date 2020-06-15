//  Vector class for Entropic Trajectories and EDMC
#pragma once

#include <vector>
#include <iostream>
#include <stdio.h>
#include <complex>

namespace ET
{
  template<typename T>
  class Vector
  {
  public:
    Vector();
    ~Vector();
    Vector(unsigned int n);
    Vector(std::vector<T> _v);

    //  Getters
    unsigned int get_dim() const;

    //  Operator overloads
    //  Access operators
    T& operator()(const unsigned int& i);
    const T& operator()(const unsigned int& i) const;

  private:
    unsigned int _n;
    std::vector<T> _v;
  };

  //  Special zero vector
  template<typename T>
  Vector<T> zeroes(unsigned int n);

  template class Vector<double>;
}
