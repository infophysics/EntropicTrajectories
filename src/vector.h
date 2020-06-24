//  Vector class for Entropic Trajectories and EDMC
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <complex>
#include <lapacke.h>
#include <cblas.h>
#include <stdint.h>

#include "utils.h"

namespace ET
{
  template<typename T>
  class Vector
  {
  public:
    Vector();
    ~Vector();
    Vector(uint32_t dim);
    Vector(std::string name, uint32_t dim);
    Vector(std::vector<T> vec);
    Vector(std::string name, std::vector<T> vec);
    Vector(uint32_t dim, const T& init);
    Vector(std::string name, uint32_t dim, const T& init);

    //  Getters
    uint32_t getDim() const;
    std::vector<T> getVec() const;
    std::string getName() const;
    //  Setters
    void setDim(uint32_t dim);
    void setVec(std::vector<T> vec);
    void setName(std::string name);

    //  Operator overloads
    Vector<T>& operator=(const Vector<T>& vector);
    bool operator==(const Vector<T>& vector) const;
    bool operator!=(const Vector<T>& vector) const;
    Vector<T> operator-() const;
    Vector<T> operator+(const Vector<T>& vector) const;
    Vector<T>& operator+=(const Vector<T>& vector);
    Vector<T> operator-(const Vector<T>& vector) const;
    Vector<T>& operator-=(const Vector<T>& vector);
    //  Scalar product
    T dot(const Vector<T>* vector);
    //  Scalar operators
    Vector<T> operator+(const T& s) const;
    Vector<T> operator-(const T& s) const;
    Vector<T> operator*(const T& s) const;
    Vector<T> operator/(const T& s) const;
    Vector<T>& operator+=(const T& s);
    Vector<T>& operator-=(const T& s);
    Vector<T>& operator*=(const T& s);
    Vector<T>& operator/=(const T& s);
    //  Overloads of scalar operations from the left.
    //  since we are trying to friend a template argument,
    //  the friend method must be defined within the class block.
    friend Vector<T> operator+(T s, const Vector<T>& vector)
    {
      uint32_t dim = vector.getDim();
      std::string name = "(" + std::to_string(s) + " + "  + vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (uint32_t i = 0; i < dim; i++) {
          v(i) = vector(i) + s;
      }
      return v;
    }
    friend Vector<T> operator-(T s, const Vector<T>& vector)
    {
      uint32_t dim = vector.getDim();
      std::string name = "(" + std::to_string(s) + " - "  + vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (uint32_t i = 0; i < dim; i++) {
          v(i) = s - vector(i);
      }
      return v;
    }
    friend Vector<T> operator*(T s, const Vector<T>& vector)
    {
      uint32_t dim = vector.getDim();
      std::string name = "(" + std::to_string(s) + " * "  + vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (uint32_t i = 0; i < dim; i++) {
          v(i) = vector(i) * s;
      }
      return v;
    }
    friend Vector<T> operator/(T s, const Vector<T>& vector)
    {
      uint32_t dim = vector.getDim();
      std::string name = "(" + std::to_string(s) + " / "  + vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      std::vector<T> vec(dim);
      for (uint32_t i = 0; i < dim; i++)
      {
        if (vector(i) == 0)
        {
          return v;
        }
        else
        {
          vec[i] = s / vector(i);
        }
      }
      v.setVec(vec);
      return v;
    }

    //  Access operators
    T& operator()(const uint32_t& i);
    const T& operator()(const uint32_t& i) const;

    std::string summary();
  private:
    //  dimension
    uint32_t _dim;
    //  container for the coefficients in R^n
    std::vector<T> _vec;
    //  name
    std::string _name;
  };

  //  Special zero vector
  template<typename T>
  Vector<T> zeroes(uint32_t dim);

  template class Vector<double>;
}
