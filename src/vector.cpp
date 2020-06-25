//------------------------------------------------------------------------------
//  vector.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee is hereby granted.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//------------------------------------------------------------------------------
#include <vector>
#include <iostream>
#include <stdio.h>
#include <complex>
#include "vector.h"

namespace ET
{

  template<typename T>
  Vector<T>::Vector() : _dim(0), _name(" ")
  {
  }
  template<typename T>
  Vector<T>::~Vector()
  {
  }
  template<typename T>
  Vector<T>::Vector(uint32_t dim) : _dim(dim), _name(" ")
  {
  }
  template<typename T>
  Vector<T>::Vector(std::string name, uint32_t dim)
  : _dim(dim), _name(name)
  {
    _vec.resize(_dim,0.0);
  }
  template<typename T>
  Vector<T>::Vector(std::vector<T> vec)
  : _dim(vec.size()), _name(" "), _vec(vec)
  {
  }
  template<typename T>
  Vector<T>::Vector(std::string name, std::vector<T> vec)
  : _dim(vec.size()), _name(name), _vec(std::move(vec))
  {
  }
  template<typename T>
  Vector<T>::Vector(uint32_t dim, const T& init)
  : _dim(dim), _name(" ")
  {
    std::vector<T> vec(_dim,init);
    _vec = vec;
  }
  template<typename T>
  Vector<T>::Vector(std::string name, uint32_t dim, const T& init)
  : _dim(dim), _name(name)
  {
    std::vector<T> vec(_dim,init);
    _vec = vec;
  }

  template<typename T>
  uint32_t Vector<T>::getDim() const
  {
    return _dim;
  }
  template<typename T>
  std::vector<T> Vector<T>::getVec() const
  {
    return _vec;
  }
  template<typename T>
  std::vector<T>* Vector<T>::accessVec()
  {
    return &_vec;
  }
  template<typename T>
  std::string Vector<T>::getName() const
  {
    return _name;
  }
  template<typename T>
  void Vector<T>::setDim(uint32_t dim)
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
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] = vector(i);
    }
    return *this;
  }
  template<typename T>
  bool Vector<T>::operator==(const Vector<T>& vector) const
  {
    if (_dim != vector.getDim())
      return false;
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        if (vector(i) != _vec[i])
          return true;
    }
    return false;
  }
  template<typename T>
  Vector<T> Vector<T>::operator-() const
  {
    std::vector<T> vec(_dim);
    for (uint32_t i = 0; i < _dim; i++)
    {
      vec[i] = -1*_vec[i];
    }
    std::string name = "-" + _name;
    return Vector<T>(name,vec);
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
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] += vector(i);
    }
    return *this;
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
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] -= vector(i);
    }
    return *this;
  }
  template<typename T>
  //  Scalar product
  T Vector<T>::dot(const Vector<T>* vector)
  {
    if(_dim != vector->getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return 0;
    }
    T result = 0;
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] + s;
    }
    return v;
  }

  template<typename T>
  Vector<T> Vector<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint32_t i = 0; i < _dim; i++) {
        v(i) = _vec[i] - s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        v(i) = _vec[i]/s;
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] += s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++) {
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
    for (uint32_t i = 0; i < _dim; i++) {
        _vec[i] /= s;
    }
    return *this;
  }

  template<typename T>
  T& Vector<T>::operator()(const uint32_t& i)
  {
    return this->_vec[i];
  }
  template<typename T>
  const T& Vector<T>::operator()(const uint32_t& i) const
  {
    return this->_vec[i];
  }

  template<typename T>
  Vector<T> zeroes(uint32_t dim)
  {
    std::vector<T> vec(dim,0.0);
    Vector<T> v(vec);
    return v;
  }

  template<typename T>
  std::string Vector<T>::summary()
  {
    std::stringstream s;
    s.str("");
    s.clear();
    std::string sum = "dim: " + std::to_string(_dim)
                    +  ", type: "
                    + type_name<decltype(_vec[0])>();
    if (_name != " ")
    {
      sum +=  ", name: '" + _name + "'";
    }
    if (_vec.size() == 0)
    {
      sum += "\n[  empty  ]";
      return sum;
    }
    sum += "\n[ ";
    if (_dim < 10)
    {
      if (_vec[0] >= 0.0)
          sum += " ";
      for (uint32_t i = 0; i < _dim; i++)
      {
        sum += scientific_not(this->_vec[i],3);
        if (i < _dim-1)
        {
          if (_vec[i+1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
        }
      }
    }
    else
    {
      if (_vec[0] >= 0.0)
          sum += " ";
      sum += scientific_not(this->_vec[0],3);
      if (_vec[1] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->_vec[1],3);
      if (_vec[2] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->_vec[2],3);
      sum += "   ";
      sum += "...   ";
      if (_vec[_dim-3] >= 0.0)
        sum += " ";
      sum += scientific_not(this->_vec[_dim-3],3);
      if (_vec[_dim-2] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->_vec[_dim-2],3);
      if (_vec[_dim-1] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->_vec[_dim-1],3);
    }
    sum += "  ]";
    return sum;
  }

}
