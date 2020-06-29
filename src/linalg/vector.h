//------------------------------------------------------------------------------
//  vector.h
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

//------------------------------------------------------------------------------
//  Vector class
//
//------------------------------------------------------------------------------
namespace ET
{
  //--------------------------------------------------------------------------
  //  This vector class acts as a general container for (n)-dimensional
  //  vectors.  It wraps several methods from BLAS level one.
  //--------------------------------------------------------------------------
  template<typename T>
  class Vector
  {
  public:
    Vector();
    ~Vector();
    Vector(const Vector<T>& vector);
    Vector(uint32_t dim);
    Vector(std::string name, uint32_t dim);
    Vector(std::vector<T> vec);
    Vector(std::string name, std::vector<T> vec);
    Vector(uint32_t dim, const T& init);
    Vector(std::string name, uint32_t dim, const T& init);

    //  Getters
    uint32_t getDim() const;
    //  get const reference to vec
    std::vector<T> getVec() const;
    //  get access to vec
    std::vector<T>* accessVec();
    //  get access to beginning of vec
    T* data();
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
    T operator*(const Vector<T>& vector) const;
    T dot(const Vector<T>& vector) const;
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
    //  conatiner for message status
    int _flag;
    //  container for messages
    std::string _info;
  };

  //----------------------------------------------------------------------------
  //  Special vector initializers
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> zeroes(uint32_t dim);
  template<typename T>
  Vector<T> ones(uint32_t dim);
  //----------------------------------------------------------------------------

  template class Vector<double>;
  template class Vector<float>;

  //----------------------------------------------------------------------------
  //  Level 1 BLAS methods
  //----------------------------------------------------------------------------
  void DSWAP(Vector<double>& v, Vector<double>& u);
  void DSCAL(Vector<double>& v, const double& scale);
  Vector<double> DCOPY(Vector<double>& v);
  void DCOPY(Vector<double>& v, Vector<double>& u);
  void DAXPY(Vector<double>& v, const double& scale, Vector<double>& u);
  double DDOT(Vector<double>& v, Vector<double>& u);
  double DNRM2(Vector<double>& v);
  double DASUM(Vector<double>& v);
  uint32_t IDAMAX(Vector<double>& v);
  uint32_t IDAMIN(Vector<double>& v);
  //----------------------------------------------------------------------------
}
