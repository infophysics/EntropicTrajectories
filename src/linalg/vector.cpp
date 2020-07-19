//------------------------------------------------------------------------------
//  vector.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]

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
  //----------------------------------------------------------------------------
  //  Vector constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = " ", and m_dim = 0
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector() : m_dim(0), m_name(" ")
  {
    //std::cout << "\nVector created at location " << this;
  }
  //----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::~Vector()
  {
    //std::cout << "\nVector at location " << this << " destroyed.";
  }
  //----------------------------------------------------------------------------
  //  Copy constructor
  //    Does not delete the copied object.
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector(const Vector<T>& t_vector)
  {
    m_vec = t_vector.getVec();
    m_dim = t_vector.getDim();
    m_name = t_vector.getName();
  }
  //----------------------------------------------------------------------------
  //  Constructors with various sets of arguments, such as,
  //    std::string                 name,
  //    size_t                      dim,
  //    std::vector<T>              t_vector,
  //    const T&                    initial value for all elements.
  //  If the dimension is specified, the internal array element 'm_vec'
  //  is resized accordingly, otherwise it is left uninstantiated.
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>::Vector(size_t t_dim) : m_dim(t_dim), m_name(" ")
  {
    //std::cout << "\nVector created at location " << this;
    m_vec.resize(m_dim);
  }
  template<typename T>
  Vector<T>::Vector(std::string t_name, size_t t_dim)
  : m_dim(t_dim), m_name(t_name)
  {
    //std::cout << "\nVector created at location " << this;
    m_vec.resize(m_dim,0.0);
  }
  //  Notice that the following methods do not MOVE the t_vectors so that they
  //  change ownership.  Instead they are copied into m_vec.
  template<typename T>
  Vector<T>::Vector(std::vector<T> t_vec)
  : m_dim(t_vec.size()), m_name(" "), m_vec(t_vec)
  {
    //std::cout << "\nVector created at location " << this;
  }
  template<typename T>
  Vector<T>::Vector(std::string t_name, std::vector<T> t_vec)
  : m_dim(t_vec.size()), m_name(t_name), m_vec(t_vec)
  {
    //std::cout << "\nVector created at location " << this;
  }
  template<typename T>
  Vector<T>::Vector(size_t t_dim, const T& t_init)
  : m_dim(t_dim), m_name(" ")
  {
    //std::cout << "\nVector created at location " << this;
    std::vector<T> vec(m_dim,t_init);
    m_vec = vec;
  }
  template<typename T>
  Vector<T>::Vector(std::string t_name, size_t t_dim, const T& t_init)
  : m_dim(t_dim), m_name(t_name)
  {
    //std::cout << "\nVector created at location " << this;
    std::vector<T> vec(m_dim,t_init);
    m_vec = vec;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //  Each attribute comes with its own setters and getters.
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Vector<T>::getDim() const
  {
    return m_dim;
  }
  //  When 'getVec()' is called it will usually create a copy of m_vec.
  template<typename T>
  std::vector<T> Vector<T>::getVec() const
  {
    return m_vec;
  }
  //  In order to return the 'm_vec' attribute so that it can be manipulated
  //  by other methods, such as those which utilize BLAS and LAPACK functions,
  //  we use 'accessVec()' to return a pointer to 'm_vec'.
  template<typename T>
  std::vector<T>* Vector<T>::accessVec()
  {
    return &m_vec;
  }
  //  get access to beginning of m_vec
  template<typename T>
  T* Vector<T>::data()
  {
    return m_vec.data();
  }
  template<typename T>
  std::string Vector<T>::getName() const
  {
    return m_name;
  }
  template<typename T>
  void Vector<T>::setDim(size_t t_dim)
  {
    m_dim = t_dim;
  }
  template<typename T>
  void Vector<T>::setVec(std::vector<T> t_vec)
  {
    m_vec = t_vec;
    m_dim = t_vec.size();
  }
  template<typename T>
  void Vector<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }
  template<typename T>
  void Vector<T>::setFlag(int t_flag)
  {
    m_flag = t_flag;
  }
  template<typename T>
  void Vector<T>::setInfo(std::string t_info)
  {
    m_info = t_info;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>& Vector<T>::operator=(const Vector<T>& t_vector)
  {
    if (&t_vector == this) {
      return *this;
    }
    m_dim = t_vector.getDim();
    m_name = t_vector.getName();
    m_vec.resize(m_dim);
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] = t_vector(i);
    }
    return *this;
  }
  template<typename T>
  bool Vector<T>::operator==(const Vector<T>& t_vector) const
  {
    if (m_dim != t_vector.getDim()) {
        return false;
    }
    for (auto i = 0; i < m_dim; i++) {
        if (t_vector(i) != m_vec[i]) {
          return false;
        }
    }
    return true;
  }
  template<typename T>
  bool Vector<T>::operator!=(const Vector<T>& t_vector) const
  {
    if (m_dim != t_vector.getDim()) {
      return true;
    }
    for (auto i = 0; i < m_dim; i++) {
        if (t_vector(i) != m_vec[i]) {
          return true;
        }
    }
    return false;
  }
  template<typename T>
  Vector<T> Vector<T>::operator-() const
  {
    std::vector<T> vec(m_dim);
    for (auto i = 0; i < m_dim; i++) {
      vec[i] = -1*m_vec[i];
    }
    std::string name = "-" + m_name;
    return Vector<T>(name,vec);
  }
  template<typename T>
  Vector<T> Vector<T>::operator+(const Vector<T>& t_vector) const
  {
    if(m_dim != t_vector.getDim()) {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " + " + t_vector.getName() + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i] + t_vector(i);
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const Vector<T>& t_vector)
  {
    if(m_dim != t_vector.getDim()) {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " + " + t_vector.getName() + ")";
    m_name = name;
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] += t_vector(i);
    }
    return *this;
  }
  template<typename T>
  Vector<T> Vector<T>::operator-(const Vector<T>& t_vector) const
  {
    if(m_dim != t_vector.getDim()) {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " - " + t_vector.getName() + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i] - t_vector(i);
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const Vector<T>& t_vector)
  {
    if(m_dim != t_vector.getDim()) {
      std::cout << "Vectors incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " - " + t_vector.getName() + ")";
    m_name = name;
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] -= t_vector(i);
    }
    return *this;
  }
  template<typename T>
  T Vector<T>::operator*(const Vector<T>& t_vector) const
  {
    return dot(t_vector);
  }
  template<typename T>
  T Vector<T>::dot(const Vector<T>& t_vector) const
  {
    if(m_dim != t_vector.getDim()) {
      std::cout << "Vectors incompatible!" << std::endl;
      return 0;
    }
    T result = 0;
    for (auto i = 0; i < m_dim; i++) {
        result += m_vec[i]*t_vector(i);
    }
    return result;
  }
  template<typename T>
  //  Scalar operators
  Vector<T> Vector<T>::operator+(const T& s) const
  {
    std::string name = "(" + m_name + " + " + std::to_string(s) + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i] + s;
    }
    return v;
  }

  template<typename T>
  Vector<T> Vector<T>::operator-(const T& s) const
  {
    std::string name = "(" + m_name + " - " + std::to_string(s) + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i] - s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator*(const T& s) const
  {
    std::string name = "(" + m_name + " * " + std::to_string(s) + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i] * s;
    }
    return v;
  }
  template<typename T>
  Vector<T> Vector<T>::operator/(const T& s) const
  {
    if (s == 0) {
      std::cout << "Division by zero!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " / " + std::to_string(s) + ")";
    Vector<T> v(name, m_dim, 0.0);
    for (auto i = 0; i < m_dim; i++) {
        v(i) = m_vec[i]/s;
    }
    return v;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator+=(const T& s)
  {
    std::string name = "(" + m_name + " + " + std::to_string(s) + ")";
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] += s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator-=(const T& s)
  {
    std::string name = "(" + m_name + " - " + std::to_string(s) + ")";
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator*=(const T& s)
  {
    std::string name = "(" + m_name + " * " + std::to_string(s) + ")";
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] *= s;
    }
    return *this;
  }
  template<typename T>
  Vector<T>& Vector<T>::operator/=(const T& s)
  {
    if (s == 0) {
      std::cout << "Division by zero!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " / " + std::to_string(s) + ")";
    for (auto i = 0; i < m_dim; i++) {
        m_vec[i] /= s;
    }
    return *this;
  }
  template<typename T>
  T& Vector<T>::operator()(const size_t& i)
  {
    return this->m_vec[i];
  }
  template<typename T>
  const T& Vector<T>::operator()(const size_t& i) const
  {
    return this->m_vec[i];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various methods
  //----------------------------------------------------------------------------
  template<typename T>
  const std::string Vector<T>::summary()
  {
    std::stringstream s;
    s.str("");
    s.clear();
    std::string sum = "dim: " + std::to_string(m_dim)
                    +  ", type: "
                    + type_name<decltype(m_vec[0])>();
    if (m_name != " ") {
      sum +=  ", name: '" + m_name + "'";
    }
    if (m_vec.size() == 0) {
      sum += "\n[  empty  ]";
      return sum;
    }
    sum += "\n[ ";
    if (m_dim < 10) {
      if (m_vec[0] >= 0.0)
          sum += " ";
      for (auto i = 0; i < m_dim; i++) {
        sum += scientific_not(this->m_vec[i],3);
        if (i < m_dim-1) {
          if (m_vec[i+1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
        }
      }
    }
    else {
      if (m_vec[0] >= 0.0)
          sum += " ";
      sum += scientific_not(this->m_vec[0],3);
      if (m_vec[1] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->m_vec[1],3);
      if (m_vec[2] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->m_vec[2],3);
      sum += "   ";
      sum += "...   ";
      if (m_vec[m_dim-3] >= 0.0)
        sum += " ";
      sum += scientific_not(this->m_vec[m_dim-3],3);
      if (m_vec[m_dim-2] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->m_vec[m_dim-2],3);
      if (m_vec[m_dim-1] >= 0.0)
        sum += "   ";
      else
        sum += "  ";
      sum += scientific_not(this->m_vec[m_dim-1],3);
    }
    sum += "  ]";
    return sum;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various instantiators
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> zeroes(size_t t_dim)
  {
    return Vector<T>(t_dim,0.0);
  }
  template<typename T>
  Vector<T> ones(size_t t_dim)
  {
    return Vector<T>(t_dim,1.0);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 1 BLAS methods
  //  Vector-Vector operations
  //----------------------------------------------------------------------------
  //  DSWAP - swap the contents of two Vectors
  //  Arguments:  v     - (n)-dim Vector
  //              u     - (n)-dim Vector
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DSWAP(Vector<double>& t_v, Vector<double>& t_u)
  {
    if (t_v.getDim() != t_u.getDim()) {
      std::cout << "Vectors are incompatible!" << std::endl;
      return;
    }
    cblas_dswap(t_v.getDim(),//  dimension of the vectors
                t_v.data(),  //  pointer to elements of v
                1,           //  increment for elements of v
                t_u.data(),  //  pointer to elements of u
                1);          //  increment for elements of u
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DSCAL - scales a vector by a number
  //  Arguments:  v     - (n)-dim vector
  //              scale - double
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DSCAL(Vector<double>& t_v, const double& t_scale)
  {
    cblas_dscal(t_v.getDim(),//  dimension of the t_vector
                t_scale,     //  value to scale the t_vector by
                t_v.data(),  //  pointer to the elements of v
                1);          //  increment for the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DCOPY - returns a copy of a vector
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    Vector<double>
  //----------------------------------------------------------------------------
  Vector<double> DCOPY(Vector<double>& t_v)
  {
    std::vector<double> v_copy(t_v.getDim());
    cblas_dcopy(t_v.getDim(), //  dimension of the t_vector
                t_v.data(),   //  pointer to the elements of v
                1,            //  increment for the elements of v
                v_copy.data(),//  pointer to the elements of v_copy
                1);           //  increment for the elements of v_copy
    std::string name;
    if (t_v.getName() != " ") {
      name = t_v.getName() + "_copy";
    }
    else {
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
  void DCOPY(Vector<double>& t_v, Vector<double>& t_u)
  {
    if (t_v.getDim() != t_u.getDim()) {
      std::cout << "Vectors are incompatible!" << std::endl;
      return;
    }
    cblas_dcopy(t_v.getDim(),//  dimension of the t_vector
                t_v.data(),  //  pointer to the elements of v
                1,           //  increment for the elements of v
                t_u.data(),  //  pointer to the elements of u
                1);          //  increment for the elements of u
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
  void DAXPY(Vector<double>& t_v, const double& t_scale, Vector<double>& t_u)
  {
    if (t_v.getDim() != t_u.getDim()) {
      std::cout << "Vectors are imcompatible!" << std::endl;
      return;
    }
    cblas_daxpy(t_v.getDim(),//  dimension of the t_vectors
                t_scale,     //  scalar to multiply v
                t_v.data(),  //  pointer to the elements of v
                1,           //  increment for the elements of v
                t_u.data(),  //  pointer to the elements of u
                1);          //  increment for the elements of u
    std::string name;
    if (t_u.getName() != " " && t_v.getName() != " ") {
      name = t_u.getName() + " + " + std::to_string(t_scale)
             + " * " + t_v.getName();;
    }
    else if (t_v.getName() != " " && t_u.getName() == " ") {
      name = std::to_string(t_scale) + " * " + t_v.getName();
    }
    t_u.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DDOT - returns the dot product between two vectors
  //  Arguments:  v     - (n)-dim vector
  //              u     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DDOT(Vector<double>& t_v, Vector<double>& t_u)
  {
    if (t_v.getDim() != t_u.getDim()) {
      std::cout << "Vectors are imcompatible!" << std::endl;
      return 0;
    }
    return cblas_ddot(t_v.getDim(),//  dimension of the t_vectors
                      t_v.data(),  //  pointer to the elements of v
                      1,           //  increment of the elements of v
                      t_u.data(),  //  pointer to the elements of u
                      1);          //  increment of the elements of u
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DNRM2 - returns the Euclidean norm of a vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DNRM2(Vector<double>& t_v)
  {
    return cblas_dnrm2(t_v.getDim(),//  dimension of the t_vectors
                       t_v.data(),  //  pointer to the elements of v
                       1);          //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DASUM - returns the sum of the magnitues of the coefficients of a vector
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    double
  //----------------------------------------------------------------------------
  double DASUM(Vector<double>& t_v)
  {
    return cblas_dasum(t_v.getDim(),//  dimension of the t_vectors
                       t_v.data(),  //  pointer to the elements of v
                       1);          //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  IDAMAX - returns index to the largest element of the vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    size_t
  //----------------------------------------------------------------------------
  size_t IDAMAX(Vector<double>& t_v)
  {
    return cblas_idamax(t_v.getDim(),//  dimension of the t_vectors
                        t_v.data(),  //  pointer to the elements of v
                        1);          //  increment of the elements of v
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  IDAMIN - returns index to the smallest element of the vector v
  //  Arguments:  v     - (n)-dim vector
  //
  //  Returns:    size_t
  //----------------------------------------------------------------------------
  size_t IDAMIN(Vector<double>& t_v)
  {
    // this is broken for now
    return 0;
    // return cblas_idamin(t_v.getDim(),//  dimension of the t_vectors
    //                     t_v.data(),  //  pointer to the elements of v
    //                     1);          //  increment of the elements of v
  }
  //----------------------------------------------------------------------------
}
