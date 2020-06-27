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
//------------------------------------------------------------------------------
//  vector.cpp
//
//  This file contains the instantiations of the definitions given in vector.h.
//  The various methods are listed below,
//  ----------------------------------------------------------------------------
//  Line no.  |   Method
//  ----------------------------------------------------------------------------
//  ()        |   Vector()
//  ()        |   ~Vector()
//  ()        |   Vector(uin32_t)
//  ()        |   Vector(std::string, uint32_t)
//  ()        |   Vector(std::vector<T>)
//  ()        |   Vector(std::string, std::vector<T>)
//  ()        |   Vector(uint32_t, const T&)
//  ()        |   Vector(std:string, uint32_t, const T&)
//  ()        |   getDim()
//  ()        |   getVec()
//  ()        |   accessVec()
//  ()        |   getName()
//  ()        |   setDim(uint32_t)
//  ()        |   setVec(std::vector<T>)
//  ()        |   setName(std::string)
//  ()        |   operator=(const Vector<T>&)
//  ()        |   operator==(const Vector<T>&)
//  ()        |   operator!=(const Vector<T>&)
//  ()        |   operator-()
//  ()        |   operator+(const Vector<T>&)
//  ()        |   operator+=(const Vector<T>&)
//  ()        |   operator-(const Vector<T>&)
//  ()        |   operator-=(const Vector<T>&)
//  ()        |   operator*(const Vector<T>*)
//  ()        |   dot(const Vector<T>*)
//  ()        |   operator+(const T&)
//  ()        |   operator-(const T&)
//  ()        |   operator*(const T&)
//  ()        |   operator/(const T&)
//  ()        |   operator+=(const T&)
//  ()        |   operator-=(const T&)
//  ()        |   operator*=(const T&)
//  ()        |   operator/=(const T&)
//  ()        |   operator()(uint32_t&)
//  ()        |   operator()(uint32_t&)
//  ()        |   summary()
//------------------------------------------------------------------------------

namespace ET
{
  //----------------------------------------------------------------------------
  //  Vector constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = " ", and _dim = 0
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector() : _dim(0), _name(" ")
  {
  }
  //----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::~Vector()
  {
  }
  //----------------------------------------------------------------------------
  //  Copy constructor
  //    Does not delete the copied object.
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector(const Vector<T>& vector)
  {
    _vec = vector.getVec();
    _dim = vector.getDim();
    _name = vector.getName();
  }
  //----------------------------------------------------------------------------
  //  Constructors with various sets of arguments, such as,
  //    std::string                 name,
  //    uint32_t                    dim,
  //    std::vector<T>              vector,
  //    const T&                    initial value for all elements.
  //  If the dimension is specified, the internal array element '_vec'
  //  is resized accordingly, otherwise it is left uninstantiated.
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector(uint32_t dim) : _dim(dim), _name(" ")
  {
    _vec.resize(_dim);
  }
  template<typename T>
  Vector<T>::Vector(std::string name, uint32_t dim)
  : _dim(dim), _name(name)
  {
    _vec.resize(_dim,0.0);
  }
  //  Notice that the following methods do not MOVE the vectors so that they
  //  change ownership.  Instead they are copied into _vec.
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
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //  Each attribute comes with its own setters and getters.
  //----------------------------------------------------------------------------
  template<typename T>
  uint32_t Vector<T>::getDim() const
  {
    return _dim;
  }
  //  When 'getVec()' is called it will usually create a copy of _vec.
  template<typename T>
  std::vector<T> Vector<T>::getVec() const
  {
    return _vec;
  }
  //  In order to return the '_vec' attribute so that it can be manipulated
  //  by other methods, such as those which utilize BLAS and LAPACK functions,
  //  we use 'accessVec()' to return a pointer to '_vec'.
  template<typename T>
  std::vector<T>* Vector<T>::accessVec()
  {
    return &_vec;
  }
  //  get access to beginning of _vec
  template<typename T>
  T* Vector<T>::data()
  {
    return _vec.data();
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
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>& Vector<T>::operator=(const Vector<T>& vector)
  {
    if (&vector == this)
    {
      return *this;
    }
    _dim = vector.getDim();
    _name = vector.getName();
    _vec.resize(_dim);
    for (uint32_t i = 0; i < _dim; i++)
    {
        _vec[i] = vector(i);
    }
    return *this;
  }
  template<typename T>
  bool Vector<T>::operator==(const Vector<T>& vector) const
  {
    if (_dim != vector.getDim())
    {
        return false;
    }
    for (uint32_t i = 0; i < _dim; i++)
    {
        if (vector(i) != _vec[i])
        {
          return false;
        }
    }
    return true;
  }
  template<typename T>
  bool Vector<T>::operator!=(const Vector<T>& vector) const
  {
    if (_dim != vector.getDim())
    {
      return true;
    }
    for (uint32_t i = 0; i < _dim; i++)
    {
        if (vector(i) != _vec[i])
        {
          return true;
        }
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
    for (uint32_t i = 0; i < _dim; i++)
    {
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
    for (uint32_t i = 0; i < _dim; i++)
    {
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
    for (uint32_t i = 0; i < _dim; i++)
    {
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
    for (uint32_t i = 0; i < _dim; i++)
    {
        _vec[i] -= vector(i);
    }
    return *this;
  }
  template<typename T>
  T Vector<T>::operator*(const Vector<T>& vector) const
  {
    return dot(vector);
  }
  template<typename T>
  T Vector<T>::dot(const Vector<T>& vector) const
  {
    if(_dim != vector.getDim())
    {
      std::cout << "Vectors incompatible!" << std::endl;
      return 0;
    }
    T result = 0;
    for (uint32_t i = 0; i < _dim; i++)
    {
        result += _vec[i]*vector(i);
    }
    return result;
  }
  template<typename T>
  //  Scalar operators
  Vector<T> Vector<T>::operator+(const T& s) const
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint32_t i = 0; i < _dim; i++)
    {
        v(i) = _vec[i] + s;
    }
    return v;
  }

  template<typename T>
  Vector<T> Vector<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint32_t i = 0; i < _dim; i++)
    {
        v(i) = _vec[i] - s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Vector<T> v(name, _dim, 0.0);
    for (uint32_t i = 0; i < _dim; i++)
    {
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
    for (uint32_t i = 0; i < _dim; i++)
    {
        v(i) = _vec[i]/s;
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++)
    {
        _vec[i] += s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++)
    {
        _vec[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    for (uint32_t i = 0; i < _dim; i++)
    {
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
    for (uint32_t i = 0; i < _dim; i++)
    {
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
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various methods
  //----------------------------------------------------------------------------
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
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various instantiators
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> zeroes(uint32_t dim)
  {
    return Vector<T>(dim,0.0);
  }
  template<typename T>
  Vector<T> ones(uint32_t dim)
  {
    return Vector<T>(dim,1.0);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 1 BLAS methods
  //  vector-vector operations
  //----------------------------------------------------------------------------
  //  DSWAP - swap the contents of two vectors
  //  Arguments:  v     - (n)-dim vector
  //              u     - (n)-dim vector
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DSWAP(Vector<double>& v, Vector<double>& u)
  {
    if (v.getDim() != u.getDim())
    {
      std::cout << "Vectors are incompatible!" << std::endl;
      return;
    }
    cblas_dswap(v.getDim(),//  dimension of the vectors
                v.data(),  //  pointer to elements of v
                1,         //  increment for elements of v
                u.data(),  //  pointer to elements of u
                1);        //  increment for elements of u
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DSCAL - scales a vector by a number
  //  Arguments:  v     - (n)-dim vector
  //              scale - double
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DSCAL(Vector<double>& v, const double& scale)
  {
    cblas_dscal(v.getDim(),//  dimension of the vector
                scale,     //  value to scale the vector by
                v.data(),  //  pointer to the elements of v
                1);        //  increment for the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DCOPY - returns a copy of a vector
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    Vector<double>
  //----------------------------------------------------------------------------
  Vector<double> DCOPY(Vector<double>& v)
  {
    std::vector<double> v_copy(v.getDim());
    cblas_dcopy(v.getDim(),   //  dimension of the vector
                v.data(),     //  pointer to the elements of v
                1,            //  increment for the elements of v
                v_copy.data(),//  pointer to the elements of v_copy
                1);           //  increment for the elements of v_copy
    std::string name;
    if (v.getName() != " ")
    {
      name = v.getName() + "_copy";
    }
    else
    {
      name = " ";
    }
    return Vector<double>(name,v_copy);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DCOPY - copies the contents of v into u
  //  Arguments:  v     - (n)-dim vector
  //              u     - (n)-dim vector
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DCOPY(Vector<double>& v, Vector<double>& u)
  {
    if (v.getDim() != u.getDim())
    {
      std::cout << "Vectors are incompatible!" << std::endl;
      return;
    }
    cblas_dcopy(v.getDim(),//  dimension of the vector
                v.data(),  //  pointer to the elements of v
                1,         //  increment for the elements of v
                u.data(),  //  pointer to the elements of u
                1);        //  increment for the elements of u
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DAXPY - returns a scalar multiplied by a vector added to another
  //  Arguments:  v     - (n)-dim vector
  //              scale - double
  //              u     - (n)-dim vector
  //
  //  Returns:    Vector<double>
  //----------------------------------------------------------------------------
  void DAXPY(Vector<double>& v, const double& scale, Vector<double>& u)
  {
    if (v.getDim() != u.getDim())
    {
      std::cout << "Vectors are imcompatible!" << std::endl;
      return;
    }
    cblas_daxpy(v.getDim(),//  dimension of the vectors
                scale,     //  scalar to multiply v
                v.data(),  //  pointer to the elements of v
                1,         //  increment for the elements of v
                u.data(),  //  pointer to the elements of u
                1);        //  increment for the elements of u
    std::string name;
    if (u.getName() != " " && v.getName() != " ")
    {
      name = u.getName() + " + " + std::to_string(scale)
             + " * " + v.getName();;
    }
    else if (v.getName() != " " && u.getName() == " ")
    {
      name = std::to_string(scale) + " * " + v.getName();
    }
    u.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DDOT - returns the dot product between two vectors
  //  Arguments:  v     - (n)-dim vector
  //              u     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DDOT(Vector<double>& v, Vector<double>& u)
  {
    if (v.getDim() != u.getDim())
    {
      std::cout << "Vectors are imcompatible!" << std::endl;
      return 0;
    }
    return cblas_ddot(v.getDim(),//  dimension of the vectors
                      v.data(),  //  pointer to the elements of v
                      1,         //  increment of the elements of v
                      u.data(),  //  pointer to the elements of u
                      1);        //  increment of the elements of u
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DNRM2 - returns the Euclidean norm of a vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DNRM2(Vector<double>& v)
  {
    return cblas_dnrm2(v.getDim(),//  dimension of the vectors
                       v.data(),  //  pointer to the elements of v
                       1);        //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DASUM - returns the sum of the magnitues of the coefficients of a vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DASUM(Vector<double>& v)
  {
    return cblas_dasum(v.getDim(),//  dimension of the vectors
                       v.data(),  //  pointer to the elements of v
                       1);        //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  IDAMAX - returns index to the largest element of the vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    uint32_t
  //----------------------------------------------------------------------------
  uint32_t IDAMAX(Vector<double>& v)
  {
    return cblas_idamax(v.getDim(),//  dimension of the vectors
                        v.data(),  //  pointer to the elements of v
                        1);        //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  IDAMIN - returns index to the smallest element of the vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    uint32_t
  //----------------------------------------------------------------------------
  uint32_t IDAMIN(Vector<double>& v)
  {
    return cblas_idamin(v.getDim(),//  dimension of the vectors
                        v.data(),  //  pointer to the elements of v
                        1);        //  increment of the elements of v
  }
  //----------------------------------------------------------------------------
}
