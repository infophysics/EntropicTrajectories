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
  Vector<T>::Vector(uint64_t dim) : _dim(dim)
  {

  }
  template<typename T>
  Vector<T>::Vector(std::string name, uint64_t dim) : _dim(dim), _name(name)
  {

  }
  template<typename T>
  Vector<T>::Vector(std::vector<T> vec) : _vec(vec), _dim(_vec.size())
  {

  }
  template<typename T>
  Vector<T>::Vector(std::string name, std::vector<T> vec)
  : _vec(vec), _dim(_vec.size()), _name(name)
  {

  }

  template<typename T>
  uint64_t Vector<T>::getDim() const
  {
    return _dim;
  }
  template<typename T>
  std::vector<T> Vector<T>::getVec() const
  {
    return _vec;
  }
  template<typename T>
  std::string Vector<T>::getName() const
  {
    return _name;
  }
  template<typename T>
  void Vector<T>::setDim(uint64_t dim)
  {
    _dim = dim;
  }
  template<typename T>
  void Vector<T>::setVec(std::vector<T> vec)
  {
    _vec = vec;
    _dim = vec.size();
  }
  template<typename T>
  void Vector<T>::setName(std::string name)
  {
    _name = name;
  }

  //  Operator overloads
  template<typename T>
  Vector<T>& Vector<T>::operator=(const Vector<T>& vector)
  {
    if (&vector == this)
      return *this;

    _dim = vector.getDim();
    _name = vector.getName();
    _vec.resize(_dim);
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] = vector(i);
    }
    return *this;
  }
  template<typename T>
  bool Vector<T>::operator==(const Vector<T>& vector) const
  {
    if (_dim != vector.getDim())
      return false;
    for (uint64_t i = 0; i < _dim; i++) {
        if (vector(i) != _vec[i])
          return false;
    }
    return true;
  }
  template<typename T>
  bool Vector<T>::operator!=(const Vector<T>& vector) const
  {
    if (_dim != vector.getDim())
      return true;
    for (uint64_t i = 0; i < _dim; i++) {
        if (vector(i) != _vec[i])
          return true;
    }
    return false;
  }
  template<typename T>
  Vector<T> Vector<T>::operator-() const
  {
    std::vector<T> vec(_dim);
    for (uint64_t i = 0; i < _dim; i++)
    {
      vec[i] = -1*_vec[i];
    }
    return Vector<T>(_name,vec);
  }
  template<typename T>
  Vector<T> Vector<T>::operator+(const Vector<T>& vector) const
  {
    if(_dim != vector.getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + vector.getName() + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] + vector(i);
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const Vector<T>& vector)
  {
    if(_dim != vector.getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + vector.getName() + ")";
    _name = name;
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] += vector(i);
    }
  }
  template<typename T>
  Vector<T> Vector<T>::operator-(const Vector<T>& vector) const
  {
    if(_dim != vector.getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + vector.getName() + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] - vector(i);
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const Vector<T>& vector)
  {
    if(_dim != vector.getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + vector.getName() + ")";
    _name = name;
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] -= vector(i);
    }
  }
  template<typename T>
  //  Scalar product
  T Vector<T>::dot(const Vector<T>* vector)
  {
    if(_dim != vector->getDim())
    {    std::vector<T> vec(_dim,0.0);

      std::cout << "Vectors incompatible!" << std::endl;
      return 0;
    }
    T result = 0;
    for (uint64_t i = 0; i < _dim; i++) {
        result += _vec[i]*(*vector)(i);
    }
    return result;
  }
  template<typename T>
  //  Scalar operators
  Vector<T> Vector<T>::operator+(const T& s) const
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] + s;
    }
    return v;
  }

  template<typename T>
  Vector<T> Vector<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] - s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] * s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator/(const T& s) const
  {
    if (s == 0)
    {
      std::cout << "Division by zero!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint64_t i = 0; i < _dim; i++) {
        v(i) = _vec[i]/s;
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + ")";
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] += s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] *= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator/=(const T& s)
  {
    if (s == 0)
    {
      std::cout << "Division by zero!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    for (uint64_t i = 0; i < _dim; i++) {
        _vec[i] /= s;
    }
    return *this;
  }

  template<typename T>
  T& Vector<T>::operator()(const uint64_t& i)
  {
    return this->_vec[i];
  }
  template<typename T>
  const T& Vector<T>::operator()(const uint64_t& i) const
  {
    return this->_vec[i];
  }

  template<typename T>
  Vector<T> zeroes(uint64_t dim)
  {
    std::vector<T> vec(dim,0.0);
    Vector<T> v(vec);
    return v;
  }

}
