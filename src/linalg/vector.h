//------------------------------------------------------------------------------
//  t_vector.h
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
  //  This t_vector class acts as a general container for (n)-dimensional
  //  t_vectors.  It wraps several methods from BLAS level one.
  //--------------------------------------------------------------------------
  template<typename T>
  class Vector
  {
  public:
    Vector();
    ~Vector();
    Vector(const Vector<T>& t_vector);
    Vector(size_t t_dim);
    Vector(std::string t_name, size_t t_dim);
    Vector(std::vector<T> t_vec);
    Vector(std::string t_name, std::vector<T> t_vec);
    Vector(size_t t_dim, const T& t_init);
    Vector(std::string t_name, size_t t_dim, const T& t_init);

    //  Getters
    size_t getDim() const;
    //  There are several ways we can interact with the m_vec attribute.
    //  Some methods, such as Level one BLAS, will require access to
    //  the pointer for the first entry in m_vec.  The data() method does
    //  just this.
    std::vector<T> getVec() const;  //  get const reference to m_vec
    std::vector<T>* accessVec();    //  get access to m_vec
    T* data();                      //  get access to beginning of m_vec
    std::string getName() const;
    //  Setters
    void setDim(size_t t_dim);
    void setVec(std::vector<T> t_vec);
    void setName(std::string t_name);
    void setFlag(int t_flag);
    void setInfo(std::string t_info);

    //  Operator overloads
    Vector<T>& operator=(const Vector<T>& t_vector);
    bool operator==(const Vector<T>& t_vector) const;
    bool operator!=(const Vector<T>& t_vector) const;
    Vector<T> operator-() const;
    Vector<T> operator+(const Vector<T>& t_vector) const;
    Vector<T>& operator+=(const Vector<T>& t_vector);
    Vector<T> operator-(const Vector<T>& t_vector) const;
    Vector<T>& operator-=(const Vector<T>& t_vector);
    //  Scalar product
    T operator*(const Vector<T>& t_vector) const;
    T dot(const Vector<T>& t_vector) const;
    //  Scalar operators
    Vector<T> operator+(const T& t_s) const;
    Vector<T> operator-(const T& t_s) const;
    Vector<T> operator*(const T& t_s) const;
    Vector<T> operator/(const T& t_s) const;
    Vector<T>& operator+=(const T& t_s);
    Vector<T>& operator-=(const T& t_s);
    Vector<T>& operator*=(const T& t_s);
    Vector<T>& operator/=(const T& t_s);
    //  Overloads of scalar operations from the left.
    //  since we are trying to friend a template argument,
    //  the friend method must be defined within the class block.
    friend Vector<T> operator+(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " + "
                         + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_vector(i) + t_s;
      }
      return v;
    }
    friend Vector<T> operator-(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " - "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_s - t_vector(i);
      }
      return v;
    }
    friend Vector<T> operator*(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " * "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_vector(i) * t_s;
      }
      return v;
    }
    friend Vector<T> operator/(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " / "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      std::vector<T> vec(dim);
      for (auto i = 0; i < dim; i++) {
        if (t_vector(i) == 0) {
          return v;
        }
        else {
         vec[i] = t_s / t_vector(i);
        }
      }
      v.setVec(vec);
      return v;
    }

    //  Access operators
    T& operator()(const size_t& i);
    const T& operator()(const size_t& i) const;

    const std::string summary();

  private:
    size_t m_dim{0};            //  dimension
    std::vector<T> m_vec{{0}};  //  container for the coefficients in R^n
    std::string m_name{""};     //  t_name
    int m_flag{0};              //  conatiner for message status
    std::string m_info{""};     //  container for messages
  };

  //----------------------------------------------------------------------------
  //  Special t_vector initializers
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> zeroes(size_t t_dim);
  template<typename T>
  Vector<T> ones(size_t t_dim);
  //----------------------------------------------------------------------------

  template class Vector<double>;
  template class Vector<float>;

  //----------------------------------------------------------------------------
  //  Level 1 BLAS methods
  //----------------------------------------------------------------------------
  void DSWAP(Vector<double>& t_v, Vector<double>& t_u);
  void DSCAL(Vector<double>& t_v, const double& t_scale);
  Vector<double> DCOPY(Vector<double>& t_v);
  void DCOPY(Vector<double>& t_v, Vector<double>& t_u);
  void DAXPY(Vector<double>& t_v, const double& scale, Vector<double>& t_u);
  double DDOT(Vector<double>& t_v, Vector<double>& t_u);
  double DNRM2(Vector<double>& t_v);
  double DASUM(Vector<double>& t_v);
  size_t IDAMAX(Vector<double>& t_v);
  size_t IDAMIN(Vector<double>& t_v);
  //----------------------------------------------------------------------------
}
